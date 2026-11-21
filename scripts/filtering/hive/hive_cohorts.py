#!/usr/bin/env python3
"""The ten final cohorts as they sit on Hive, and the two steps that run there.

The chain is filter_stats/filter_chain.py; this only says where each cohort's three inputs live on
Quobyte and where its results go. Two steps, because only one of them needs the models:

  judge   one array task per cohort: read that cohort's convergence representatives out of its PDB
          tree and write their topology verdicts. One model per trajectory, not all eleven
          recycles, which is why the topology screen is the last step of the chain and not the
          first.
  chain   every cohort in one job: convergence, filter ring, interface QC (Kv2.1 only) and then
          those verdicts, writing <cohort>_chain_qc.csv and a .json of the per-step counts.

  python hive_cohorts.py judge --index $SLURM_ARRAY_TASK_ID
  python hive_cohorts.py chain --all
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                      # the copied filter_stats/ and shared/ live here
sys.path.insert(0, str(HERE / "filter_stats"))

import filter_chain as fc  # noqa: E402

ROOT = Path("/quobyte/yarovoygrp/ahgz/vgic_mutants")
OUT = ROOT / "analysis" / "filter_chain"


def _cohort(key, channel, distances, convergence, models) -> dict:
    return {"key": key, "channel": channel, "distances": ROOT / distances,
            "convergence": ROOT / convergence, "models": ROOT / models,
            "topology": OUT / f"{key}_topology_verdicts.csv",
            "output": OUT / f"{key}_chain_qc.csv", "summary": OUT / f"{key}_chain_qc.json"}


# Kv2.1 at 99 columns, Nav1.5 with both corrected masks, Cav1.2 at 50%. Every mask is 50%; Nav1.5
# is the only channel with two mask designs. The order here is the array index order.
COHORTS = [
    _cohort("kv21_wt_c99", "Kv2.1",
            "Kv2.1/dataDistances/26-11-20_Kv2.1_wt_masked_common99_50pct_AF2_distances.csv",
            "Kv2.1/rmsd_convergence_filtering/Kv21_rmsd_convergence_wt_masked_common99_50pct/all_ok_models.csv",
            "Kv2.1/wt/masked_common99_50pct/models"),
    _cohort("kv21_l403a_c99", "Kv2.1",
            "Kv2.1/dataDistances/26-11-19_Kv2.1_l403a_masked_common99_50pct_AF2_distances.csv",
            "Kv2.1/rmsd_convergence_filtering/Kv21_rmsd_convergence_l403a_masked_common99_50pct/all_ok_models.csv",
            "Kv2.1/l403a/masked_common99_50pct/models"),
    _cohort("kv21_f412l_c99", "Kv2.1",
            "Kv2.1/dataDistances/26-11-20_Kv2.1_f412l_masked_common99_50pct_AF2_distances.csv",
            "Kv2.1/rmsd_convergence_filtering/Kv21_rmsd_convergence_f412l_masked_common99_50pct/all_ok_models.csv",
            "Kv2.1/f412l/masked_common99_50pct/models"),
    _cohort("nav15_wt_v1", "Nav1.5",
            "Nav1.5/dataDistances/26-11-10_Nav15_wt_masked_v1_50pct_AF2_distances.csv",
            "Nav1.5/Nav15_rmsd_convergence_wt_masked_v1_50pct/all_ok_models.csv",
            "Nav1.5/wt/masked_v1_50pct/models"),
    _cohort("nav15_qqq_v1", "Nav1.5",
            "Nav1.5/dataDistances/26-11-10_Nav15_qqq_masked_v1_50pct_AF2_distances.csv",
            "Nav1.5/Nav15_rmsd_convergence_qqq_masked_v1_50pct/all_ok_models.csv",
            "Nav1.5/qqq/masked_v1_50pct/models"),
    _cohort("nav15_wt_v2", "Nav1.5",
            "Nav1.5/dataDistances/26-11-10_Nav15_wt_masked_v2_50pct_AF2_distances.csv",
            "Nav1.5/Nav15_rmsd_convergence_wt_masked_v2_50pct/all_ok_models.csv",
            "Nav1.5/wt/masked_v2_50pct/models"),
    _cohort("nav15_qqq_v2", "Nav1.5",
            "Nav1.5/dataDistances/26-11-10_Nav15_qqq_masked_v2_50pct_AF2_distances.csv",
            "Nav1.5/Nav15_rmsd_convergence_qqq_masked_v2_50pct/all_ok_models.csv",
            "Nav1.5/qqq/masked_v2_50pct/models"),
    _cohort("cav12_wt", "Cav1.2",
            "Cav1.2/dataDistances/26-11-10_Cav12_wt_masked_v3_50pct_AF2_distances.csv",
            "Cav1.2/Cav12_rmsd_convergence_wt_masked_v3_50pct/all_ok_models.csv",
            "Cav1.2/wt/masked_v3_50pct/models"),
    _cohort("cav12_g402s", "Cav1.2",
            "Cav1.2/dataDistances/26-11-10_Cav12_g402s_masked_v3_50pct_AF2_distances.csv",
            "Cav1.2/Cav12_rmsd_convergence_g402s_masked_v3_50pct/all_ok_models.csv",
            "Cav1.2/g402s/masked_v3_50pct/models"),
    _cohort("cav12_g406r", "Cav1.2",
            "Cav1.2/dataDistances/26-11-10_Cav12_g406r_masked_v3_50pct_AF2_distances.csv",
            "Cav1.2/Cav12_rmsd_convergence_g406r_masked_v3_50pct/all_ok_models.csv",
            "Cav1.2/g406r/masked_v3_50pct/models"),
]
BY_KEY = {c["key"]: c for c in COHORTS}


def representatives(cohort: dict) -> set:
    convergence = pd.read_csv(cohort["convergence"], low_memory=False)
    names = fc.basenames(convergence["pdb_file"])
    chosen = convergence["earliest_converged_selected"].astype(str).str.strip().str.lower() == "true"
    return set(names[chosen])


def judge(cohort: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wanted = representatives(cohort)
    print(f"{cohort['key']}: {len(wanted)} representatives to judge from {cohort['models']}", flush=True)
    verdicts = fc.judge_topology(fc.collect_models(cohort["models"], wanted),
                                 chan=fc.get_channel(cohort["channel"]))
    verdicts.to_csv(cohort["topology"], index=False)
    failed = int((verdicts.error.astype(str) != "").sum())
    print(f"{cohort['key']}: {int(verdicts.topology_ok.sum())}/{len(verdicts)} pass the topology screen, "
          f"{failed} could not be judged; wrote {cohort['topology']}", flush=True)


def chain(cohort: dict) -> dict:
    result = fc.run_chain(pd.read_csv(cohort["distances"], low_memory=False),
                          pd.read_csv(cohort["convergence"], low_memory=False),
                          channel=cohort["channel"],
                          topology_verdicts=pd.read_csv(cohort["topology"]),
                          label=cohort["key"])
    result.representatives.to_csv(cohort["output"], index=False)
    counts = {"cohort": cohort["key"], **result.counts}
    cohort["summary"].write_text(json.dumps(counts, indent=2, default=str) + "\n")
    print(f"{cohort['key']}: wrote {cohort['output']}", flush=True)
    return counts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mode", choices=("judge", "chain", "list"))
    parser.add_argument("--index", type=int, help="array index into the cohort list")
    parser.add_argument("--cohort", action="append", help="a cohort key; repeatable")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args(argv)

    if args.mode == "list":
        for i, c in enumerate(COHORTS):
            print(f"  {i:2d}  {c['key']:16s} {c['channel']:7s} "
                  f"{'models ok' if c['models'].is_dir() else 'MODELS MISSING':14s} "
                  f"{'distances ok' if c['distances'].exists() else 'DISTANCES MISSING'}")
        return 0

    if args.index is not None:
        chosen = [COHORTS[args.index]]
    elif args.all:
        chosen = COHORTS
    elif args.cohort:
        chosen = [BY_KEY[k] for k in args.cohort]
    else:
        parser.error("give --index, --cohort or --all")

    if args.mode == "judge":
        for cohort in chosen:
            judge(cohort)
        return 0

    missing = [c["key"] for c in chosen if not c["topology"].exists()]
    if missing:
        print(f"no topology verdicts yet for {missing}; run the judge step first", flush=True)
    for cohort in chosen:
        if cohort["topology"].exists():
            chain(cohort)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
