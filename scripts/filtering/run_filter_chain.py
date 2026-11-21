#!/usr/bin/env python3
"""Run the four-step chain over the final cohorts and write one filtered table each.

The chain itself is filter_stats/filter_chain.py; this is the registry of what to run it on and
where the results go. Two modes, because the two halves of the work live in different places:

  judge    reads a cohort's PDB tree and writes the topology verdicts for the representatives only.
           This is the part that needs the models, so it runs where they are (Hive). It judges one
           model per trajectory rather than all eleven recycles, which is the whole reason the
           topology screen is the last step and not the first.

  chain    reads the distance table, the convergence table and those verdicts, and writes
           <base>_chain_qc.csv plus a <base>_chain_qc.json holding the per-step counts.

Paths come from VGIC_DATA_ROOT, as everywhere else in the repository.

  python scripts/filtering/run_filter_chain.py judge --cohort kv21_f412l_c99 --model-dir <tree>
  python scripts/filtering/run_filter_chain.py chain --cohort kv21_f412l_c99
  python scripts/filtering/run_filter_chain.py chain --all
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "filter_stats"))

import filter_chain as fc  # noqa: E402

DATA_ROOT = Path(os.environ.get("VGIC_DATA_ROOT", Path.home() / ".cache" / "vgic-mutants")).expanduser()

# One entry per cohort the paper reports: Kv2.1 at 99 columns, Nav1.5 with both corrected masks,
# Cav1.2 at 50%. Every mask is 50%; Nav1.5 is the only channel with two mask designs.
COHORTS = {
    "kv21_wt_c99": ("Kv2.1", "kv21", "26-11-20_Kv2.1_wt_masked_common99_50pct"),
    "kv21_l403a_c99": ("Kv2.1", "kv21", "26-11-19_Kv2.1_l403a_masked_common99_50pct"),
    "kv21_f412l_c99": ("Kv2.1", "kv21", "26-11-20_Kv2.1_f412l_masked_common99_50pct"),
    "nav15_wt_v1": ("Nav1.5", "nav15", "26-11-10_Nav15_wt_masked_v1_50pct"),
    "nav15_qqq_v1": ("Nav1.5", "nav15", "26-11-10_Nav15_qqq_masked_v1_50pct"),
    "nav15_wt_v2": ("Nav1.5", "nav15", "26-11-10_Nav15_wt_masked_v2_50pct"),
    "nav15_qqq_v2": ("Nav1.5", "nav15", "26-11-10_Nav15_qqq_masked_v2_50pct"),
    "cav12_wt": ("Cav1.2", "cav12", "26-11-10_Cav12_wt_masked_v3_50pct"),
    "cav12_g402s": ("Cav1.2", "cav12", "26-11-10_Cav12_g402s_masked_v3_50pct"),
    "cav12_g406r": ("Cav1.2", "cav12", "26-11-10_Cav12_g406r_masked_v3_50pct"),
}


def paths(key: str) -> dict:
    channel, folder, base = COHORTS[key]
    here = DATA_ROOT / folder / "dataDistances"
    return {"channel": channel, "base": base,
            "distances": here / f"{base}_AF2_distances.csv",
            "convergence": here / f"{base}_all_ok_models.csv",
            "topology": here / f"{base}_topology_verdicts.csv",
            "output": here / f"{base}_chain_qc.csv",
            "summary": here / f"{base}_chain_qc.json"}


def judge(key: str, model_dir: Path) -> None:
    """Topology verdicts for one cohort's convergence representatives. Run where the models are."""
    p = paths(key)
    convergence = pd.read_csv(p["convergence"], low_memory=False)
    names = fc.basenames(convergence["pdb_file"])
    wanted = set(names[convergence["earliest_converged_selected"].astype(str).str.lower() == "true"])
    print(f"{key}: {len(wanted)} representatives to judge (one per surviving trajectory)")
    verdicts = fc.judge_topology(fc.collect_models(model_dir, wanted), chan=fc.get_channel(p["channel"]))
    verdicts.to_csv(p["topology"], index=False)
    failed = int((verdicts.error.astype(str) != "").sum())
    print(f"{key}: {int(verdicts.topology_ok.sum())}/{len(verdicts)} pass, {failed} could not be judged; "
          f"wrote {p['topology']}")


def chain(key: str, *, allow_missing_topology: bool = False) -> dict:
    p = paths(key)
    missing = [name for name in ("distances", "convergence") if not p[name].exists()]
    if missing:
        raise FileNotFoundError(f"{key}: {', '.join(str(p[m]) for m in missing)} not found under {DATA_ROOT}")
    verdicts = None
    if p["topology"].exists():
        verdicts = pd.read_csv(p["topology"])
    elif not allow_missing_topology:
        raise FileNotFoundError(
            f"{key}: no topology verdicts at {p['topology']}. Run the judge step where the models "
            f"are, or pass --allow-missing-topology to run the first three steps only.")

    if verdicts is None:
        print(f"{key}: WARNING running without the topology step; the chain is incomplete")
        distances = pd.read_csv(p["distances"], low_memory=False)
        convergence = pd.read_csv(p["convergence"], low_memory=False)
        chan = fc.get_channel(p["channel"])
        kept, counts_conv = fc.apply_convergence(distances, convergence)
        kept, counts_ring = fc.apply_filter_ring(kept, chan)
        counts = {"channel": chan.name, "convergence": counts_conv, "filter_ring": counts_ring,
                  "topology": "not run"}
        if chan.interface_qc:
            before = kept["trajectory"].nunique()
            kept = fc.apply_kv21_interface_qc(kept, threshold_A=fc.INTERFACE_MAX_A, trajectory_level=True)
            counts["interface_qc"] = {"trajectories_rejected": before - kept["trajectory"].nunique()}
        representatives = kept.loc[kept["is_representative"]]
        counts["surviving_trajectories"] = representatives["trajectory"].nunique()
        counts["complete"] = False
        representatives.to_csv(p["output"], index=False)
    else:
        result = fc.run_chain(pd.read_csv(p["distances"], low_memory=False),
                              pd.read_csv(p["convergence"], low_memory=False),
                              channel=p["channel"], topology_verdicts=verdicts, label=key)
        counts = {**result.counts, "complete": True}
        result.representatives.to_csv(p["output"], index=False)
    p["summary"].write_text(json.dumps({"cohort": key, **counts}, indent=2, default=str) + "\n")
    print(f"{key}: wrote {p['output']}")
    return counts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mode", choices=("judge", "chain", "list"))
    parser.add_argument("--cohort", action="append", help="a key from the registry; repeatable")
    parser.add_argument("--all", action="store_true", help="every cohort in the registry")
    parser.add_argument("--model-dir", type=Path, help="judge: the cohort's PDB tree")
    parser.add_argument("--allow-missing-topology", action="store_true",
                        help="chain: run the first three steps where no verdicts exist yet")
    args = parser.parse_args(argv)

    if args.mode == "list":
        for key, (channel, _, base) in COHORTS.items():
            p = paths(key)
            state = "ready" if p["distances"].exists() else "distances not here"
            if p["distances"].exists() and not p["topology"].exists():
                state = "needs the judge step"
            print(f"  {key:16s} {channel:7s} {base:46s} {state}")
        return 0

    keys = list(COHORTS) if args.all else (args.cohort or [])
    if not keys:
        parser.error("name at least one --cohort, or pass --all")
    if args.mode == "judge":
        if not args.model_dir:
            parser.error("the judge step needs --model-dir")
        for key in keys:
            judge(key, args.model_dir)
        return 0

    failures = {}
    for key in keys:
        try:
            chain(key, allow_missing_topology=args.allow_missing_topology)
        except FileNotFoundError as exc:
            failures[key] = str(exc)
    for key, why in failures.items():
        print(f"not run: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
