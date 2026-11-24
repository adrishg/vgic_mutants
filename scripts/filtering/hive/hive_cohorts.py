#!/usr/bin/env python3
"""The ten final cohorts as they sit on Hive, and the two steps that run there.

The chain is filter_stats/filter_chain.py; this only says where each cohort's three inputs live on
Quobyte and where its results go. Two steps, because only one of them needs the models:

  measure one array task per cohort: read every model the distance table names and write its
          pore-VSD interface distance - the lower S6 against its own S4-S5 linker, per subunit.
          This one walks all eleven recycles, because the check rejects a trajectory on any of
          them, and it is written out once so that trying a different cutoff costs nothing.
  judge   one array task per cohort: read that cohort's convergence representatives out of its PDB
          tree and write their topology verdicts. One model per trajectory, not all eleven
          recycles, which is why the topology screen is the last step of the chain and not the
          first.
  omit    every cohort in one job: the trajectories our chain rejects, written per dataset as the
          condition,model,seed,reason CSV that Brandon's [screen] omit key reads. That is how the
          order we settled on reaches his statistics: his screen would otherwise be the only
          filter, and it runs before anything that needs a whole trajectory.
  sweep   every cohort in one job: the same inputs run through several filtering setups, so the
          cost of a setup can be read off one table instead of argued about. Nothing is recomputed
          between setups; they are all questions asked of the four tables the steps above wrote.
  verdicts every cohort in one job: the per-prediction pass/fail table, one CSV per dataset,
          holding every prediction with what each step measured and which step removed it. This is
          the deliverable someone outside the project reads; `omit` is the same walk written for
          Brandon's [screen] omit key.
  chain   every cohort in one job: convergence, filter ring, interface QC (Kv2.1 only) and then
          those verdicts, writing <cohort>_chain_qc.csv and a .json of the per-step counts.

  python hive_cohorts.py measure --index $SLURM_ARRAY_TASK_ID
  python hive_cohorts.py judge --index $SLURM_ARRAY_TASK_ID
  python hive_cohorts.py chain --all
  python hive_cohorts.py sweep --all
  python hive_cohorts.py omit --all
  python hive_cohorts.py verdicts --all
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
            "interface": OUT / f"{key}_interface.csv",
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
    # The two Nav1.5 vanilla cohorts are appended rather than filed with the other Nav1.5 entries
    # so that the array indices of the first ten, which earlier jobs and their logs refer to, do
    # not move. Their models sit in the ColabFold output folder itself, not in a models/ subfolder.
    _cohort("nav15_wt_vanilla", "Nav1.5",
            "Nav1.5/dataDistances/26-07-27_Nav15_wt_vanillaAF2_distances_all_all.csv",
            "Nav1.5/rmsd_convergence_filtering/Nav15_rmsd_convergence_wt_vanilla/all_ok_models.csv",
            "Nav1.5/wt/vanilla"),
    _cohort("nav15_qqq_vanilla", "Nav1.5",
            "Nav1.5/dataDistances/26-07-27_Nav15_qqq_vanilla_AF2_distances_all_all.csv",
            "Nav1.5/rmsd_convergence_filtering/Nav15_rmsd_convergence_qqq_vanilla/all_ok_models.csv",
            "Nav1.5/qqq/vanilla"),
]
BY_KEY = {c["key"]: c for c in COHORTS}


def representatives(cohort: dict) -> set:
    convergence = pd.read_csv(cohort["convergence"], low_memory=False)
    names = fc.basenames(convergence["pdb_file"])
    chosen = convergence["earliest_converged_selected"].astype(str).str.strip().str.lower() == "true"
    return set(names[chosen])


def measure(cohort: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    distances = pd.read_csv(cohort["distances"], usecols=["pdb_file"], low_memory=False)
    wanted = set(fc.basenames(distances["pdb_file"]))
    print(f"{cohort['key']}: measuring {len(wanted)} models under {cohort['models']}", flush=True)
    frame = fc.measure_interface(fc.collect_models(cohort["models"], wanted),
                                 fc.get_channel(cohort["channel"]))
    frame.to_csv(cohort["interface"], index=False)
    widest = frame["s6_linker_max_A"]
    print(f"{cohort['key']}: S6-to-linker median={widest.median():.2f} A, "
          f"95th={widest.quantile(0.95):.2f} A, max={widest.max():.2f} A, "
          f"{int(frame.error.astype(str).ne('').sum())} unreadable; wrote {cohort['interface']}", flush=True)


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


def chain(cohort: dict, **kwargs) -> dict:
    measures = pd.read_csv(cohort["interface"]) if cohort["interface"].exists() else None
    result = fc.run_chain(pd.read_csv(cohort["distances"], low_memory=False),
                          pd.read_csv(cohort["convergence"], low_memory=False),
                          channel=cohort["channel"],
                          topology_verdicts=pd.read_csv(cohort["topology"]),
                          interface_measures=measures,
                          label=cohort["key"], **kwargs)
    result.representatives.to_csv(cohort["output"], index=False)
    counts = {"cohort": cohort["key"], **result.counts}
    cohort["summary"].write_text(json.dumps(counts, indent=2, default=str) + "\n")
    print(f"{cohort['key']}: wrote {cohort['output']}", flush=True)
    return counts


# Each setup is (name, interface_mode, s6-linker cutoff in A, whether the topology screen runs).
# "published" is the 27 A Kv2.1 check and simply does not run for the other two channels, which is
# the asymmetry the sweep is meant to expose.
SETUPS = (
    ("convergence_only", "none", None, False),
    ("ring", "none", None, False),
    ("ring_topology", "none", None, True),
    ("ring_published", "published", None, False),
    ("ring_published_topology", "published", None, True),
    ("ring_s6linker10", "s6_linker", 10.0, False),
    ("ring_s6linker12", "s6_linker", 12.0, False),
    ("ring_s6linker16", "s6_linker", 16.0, False),
    ("ring_s6linker20", "s6_linker", 20.0, False),
    ("ring_s6linker12_topology", "s6_linker", 12.0, True),
    ("ring_published_s6linker12_topology", "both", 12.0, True),
)


# cohort key -> (dataset folder, the condition name its config uses)
FOR_CONFIG = {
    "kv21_wt_c99": ("kv21_masked99", "WT"),
    "kv21_l403a_c99": ("kv21_masked99", "L403A"),
    "kv21_f412l_c99": ("kv21_masked99", "F412L"),
    "nav15_wt_vanilla": ("nav15_vanilla", "WT"),
    "nav15_qqq_vanilla": ("nav15_vanilla", "QQQ"),
    "nav15_wt_v1": ("nav15_masked_v1", "WT"),
    "nav15_qqq_v1": ("nav15_masked_v1", "QQQ"),
    "nav15_wt_v2": ("nav15_masked_v2", "WT"),
    "nav15_qqq_v2": ("nav15_masked_v2", "QQQ"),
    "cav12_wt": ("cav12_masked_v3", "WT"),
    "cav12_g402s": ("cav12_masked_v3", "G402S"),
    "cav12_g406r": ("cav12_masked_v3", "G406R"),
}


def omit_rows(cohort: dict) -> list:
    """Every trajectory the chain rejects, with the step that rejected it.

    The steps run in the settled order and each one is asked what it took, so a trajectory is named
    for the first thing wrong with it rather than for whichever check happened to be last."""
    dataset, condition = FOR_CONFIG[cohort["key"]]
    chan = fc.get_channel(cohort["channel"])
    distances = pd.read_csv(cohort["distances"], low_memory=False)
    convergence = pd.read_csv(cohort["convergence"], low_memory=False)
    measures = pd.read_csv(cohort["interface"]) if cohort["interface"].exists() else None
    verdicts = pd.read_csv(cohort["topology"])

    before = set(fc.trajectories(distances["pdb_file"]))
    kept, _ = fc.apply_convergence(distances, convergence)
    reasons = {t: "did not converge" for t in before - set(kept["trajectory"])}

    def take(new, label):
        nonlocal kept
        for t in set(kept["trajectory"]) - set(new["trajectory"]):
            reasons[t] = label
        kept = new

    take(fc.apply_filter_ring(kept, chan)[0], "selectivity filter ring")
    if measures is not None:
        take(fc.apply_s6_linker_qc(kept, measures, threshold_A=chan.s6_linker_max_A)[0],
             "S6 off its own linker")
    if chan.interface_qc:
        from shared.dataset_selection import apply_kv21_interface_qc           # noqa: PLC0415
        take(apply_kv21_interface_qc(kept, threshold_A=fc.INTERFACE_MAX_A, trajectory_level=True),
             "pore-VSD interface 27A")
    failed = set(verdicts.loc[~fc._truth(verdicts["topology_ok"]), "model_file"].map(
        lambda n: fc.RECYCLE_SUFFIX.sub(".pdb", str(n))))
    for t in set(kept["trajectory"]) & failed:
        reasons[t] = "topology screen"

    rows = []
    for trajectory, reason in sorted(reasons.items()):
        m = fc.MODEL_META.search(trajectory.replace(".pdb", ".r0.pdb"))
        if m is None:
            raise ValueError(f"{cohort['key']}: cannot read model and seed from {trajectory!r}")
        rows.append({"dataset": dataset, "condition": condition,
                     "model": int(m.group(1)), "seed": int(m.group(2)), "reason": reason})
    return rows


def omit(chosen: list) -> None:
    """One CSV per dataset, holding every one of its conditions' rejected trajectories."""
    rows = []
    for cohort in chosen:
        if cohort["key"] not in FOR_CONFIG or not cohort["topology"].exists():
            print(f"{cohort['key']}: no chain results, skipped", flush=True)
            continue
        got = omit_rows(cohort)
        print(f"{cohort['key']}: {len(got)} trajectories omitted", flush=True)
        rows += got
    frame = pd.DataFrame(rows, columns=["dataset", "condition", "model", "seed", "reason"])
    # Every dataset gets a file, including the ones the chain rejects nothing in. Writing only the
    # non-empty ones leaves the config pointing at a path that does not exist, and "no rejects" is
    # a result worth recording rather than an absence worth skipping.
    for dataset in sorted({d for d, _c in FOR_CONFIG.values()}):
        part = frame.loc[frame["dataset"] == dataset, ["condition", "model", "seed", "reason"]]
        path = OUT / f"{dataset}_omit.csv"
        part.to_csv(path, index=False)
        print(f"wrote {path} ({len(part)} rows)", flush=True)


def sweep(chosen: list) -> pd.DataFrame:
    """Every cohort through every setup, one row each, counting what each step costs.

    The convergence-only row is there as the denominator: a step that removes nothing is not a
    clean ensemble, it is a step with no discrimination on this data, and the only way to tell the
    two apart is to see the same step bite somewhere else."""
    rows = []
    for cohort in chosen:
        distances = pd.read_csv(cohort["distances"], low_memory=False)
        convergence = pd.read_csv(cohort["convergence"], low_memory=False)
        verdicts = pd.read_csv(cohort["topology"]) if cohort["topology"].exists() else None
        measures = pd.read_csv(cohort["interface"]) if cohort["interface"].exists() else None
        for name, mode, cutoff, topology in SETUPS:
            if mode in ("s6_linker", "both") and measures is None:
                continue
            if topology and verdicts is None:
                continue
            ring = 1e9 if name == "convergence_only" else None
            result = fc.run_chain(
                distances, convergence, channel=cohort["channel"],
                topology_verdicts=verdicts if topology else _everything_passes(verdicts),
                interface_measures=measures, interface_mode=mode, s6_linker_threshold_A=cutoff,
                ring_threshold_A=ring, label=f"{cohort['key']}/{name}")
            c = result.counts
            rows.append({"cohort": cohort["key"], "channel": cohort["channel"], "setup": name,
                         "trajectories_in": c["convergence"]["trajectories_in"],
                         "converged": c["convergence"]["trajectories_out"],
                         "lost_ring": c["filter_ring"]["trajectories_rejected"],
                         "lost_published": c["interface_qc"]["trajectories_rejected"],
                         "lost_s6linker": c["s6_linker_qc"]["trajectories_rejected"],
                         "lost_topology": (c["topology"]["representatives_in"]
                                           - c["topology"]["representatives_out"]),
                         "surviving": c["surviving_trajectories"]})
    return pd.DataFrame(rows)


def verdict_rows(cohort: dict) -> pd.DataFrame:
    """The cohort's per-prediction verdict table: what each step measured and what took each loss."""
    dataset, condition = FOR_CONFIG[cohort["key"]]
    measures = pd.read_csv(cohort["interface"]) if cohort["interface"].exists() else None
    return fc.chain_verdicts(pd.read_csv(cohort["distances"], low_memory=False),
                             pd.read_csv(cohort["convergence"], low_memory=False),
                             channel=cohort["channel"],
                             topology_verdicts=pd.read_csv(cohort["topology"]),
                             interface_measures=measures,
                             dataset=dataset, condition=condition)


def verdicts(chosen: list, out_dir: Path = None) -> None:
    """One CSV per dataset, every prediction of every condition in it, passed and removed alike.

    This is the file someone outside the project reads, so it is written per dataset rather than
    per cohort: a dataset is the unit a protocol is compared in, and its conditions belong in one
    table."""
    out_dir = OUT if out_dir is None else out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    per_dataset = {}
    for cohort in chosen:
        if cohort["key"] not in FOR_CONFIG or not cohort["topology"].exists():
            print(f"{cohort['key']}: no topology verdicts, skipped", flush=True)
            continue
        frame = verdict_rows(cohort)
        passed = int((frame["status"] == "passed").sum())
        print(f"{cohort['key']}: {passed}/{len(frame)} predictions pass "
              f"({frame['removed_by'].replace('', 'passed').value_counts().to_dict()})", flush=True)
        per_dataset.setdefault(FOR_CONFIG[cohort["key"]][0], []).append(frame)
    for dataset, frames in sorted(per_dataset.items()):
        table = pd.concat(frames, ignore_index=True)
        path = out_dir / f"{dataset}_filter_verdicts.csv"
        table.to_csv(path, index=False)
        print(f"wrote {path} ({len(table)} predictions, "
              f"{int((table['status'] == 'passed').sum())} passed)", flush=True)


def _everything_passes(verdicts):
    """The verdict table with every model passing, which is how a setup leaves the screen out."""
    if verdicts is None:
        return None
    passing = verdicts.copy()
    passing["topology_ok"] = True
    return passing


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mode", choices=("measure", "judge", "chain", "sweep", "omit", "verdicts", "list"))
    parser.add_argument("--index", type=int, help="array index into the cohort list")
    parser.add_argument("--cohort", action="append", help="a cohort key; repeatable")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--out", help="verdicts only: write the per-dataset CSVs here instead")
    args = parser.parse_args(argv)

    if args.mode == "list":
        for i, c in enumerate(COHORTS):
            print(f"  {i:2d}  {c['key']:16s} {c['channel']:7s} "
                  f"{'models ok' if c['models'].is_dir() else 'MODELS MISSING':14s} "
                  f"{'interface ok' if c['interface'].exists() else 'no interface':13s} "
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

    if args.mode == "measure":
        for cohort in chosen:
            measure(cohort)
        return 0

    if args.mode == "judge":
        for cohort in chosen:
            judge(cohort)
        return 0

    if args.mode == "omit":
        omit(chosen)
        return 0

    if args.mode == "verdicts":
        verdicts(chosen, Path(args.out) if args.out else None)
        return 0

    if args.mode == "sweep":
        table = sweep(chosen)
        path = OUT / "filtering_setup_sweep.csv"
        table.to_csv(path, index=False)
        print(table.to_string(index=False), flush=True)
        print(f"wrote {path}", flush=True)
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
