#!/usr/bin/env python3
"""What changed between the Nav1.5 set sent out in November and the one the full chain keeps.

The set that went to the collaborator is not reconstructed here. It is read from the package that
was actually built and shipped,

    Nav1.5/share/Nav15_50pct_filtered_package_2026-11-11/

whose builder (Nav1.5/package_nav15_50pct.slurm) copies, per condition, exactly the `pdb_file`
column of

    Nav15_rmsd_convergence_<condition>/earliest_converged_one_per_seed_model.csv

into models/<condition>/. So the old criterion is **RMSD convergence alone, one representative per
trajectory**: the earliest recycle from which every later successive-recycle transition stays under
3 A. 500 models per condition.

That is the same unit the chain works in, which makes the comparison a set comparison rather than
an argument: the chain starts from the same convergence step and then applies the selectivity-filter
ring, the pore-VSD interface check and the topology screen. Two things are therefore checked rather
than assumed:

  - that the package's convergence tables are still byte-identical to the ones on Hive, so the
    convergence step itself has not moved under us;
  - that `earliest_converged_one_per_seed_model.csv` holds exactly the rows `all_ok_models.csv`
    marks `earliest_converged_selected`, which is what the chain reads.

If both hold, every difference is the chain's doing and nothing else's.

Writes one per-model CSV per dataset plus a summary, and prints the summary.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/quobyte/yarovoygrp/ahgz/vgic_mutants")
PACKAGE = ROOT / "Nav1.5" / "share" / "Nav15_50pct_filtered_package_2026-11-11"
CHAIN = ROOT / "analysis" / "filter_chain"
RECYCLE = __import__("re").compile(r"\.r(\d+)\.pdb$")

# package condition folder -> (our dataset, our condition, the Hive convergence folder)
CONDITIONS = {
    "wt_masked_v1_50pct": ("nav15_masked_v1", "WT", "Nav15_rmsd_convergence_wt_masked_v1_50pct"),
    "qqq_masked_v1_50pct": ("nav15_masked_v1", "QQQ", "Nav15_rmsd_convergence_qqq_masked_v1_50pct"),
    "wt_masked_v2_50pct": ("nav15_masked_v2", "WT", "Nav15_rmsd_convergence_wt_masked_v2_50pct"),
    "qqq_masked_v2_50pct": ("nav15_masked_v2", "QQQ", "Nav15_rmsd_convergence_qqq_masked_v2_50pct"),
}
# Datasets the chain now covers that the November package never held.
NOT_IN_PACKAGE = ("nav15_vanilla",)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def recycle_of(name: str):
    m = RECYCLE.search(str(name))
    return int(m.group(1)) if m else None


def old_set(folder: str) -> tuple:
    """The models the package shipped for one condition, and the two provenance checks."""
    here = PACKAGE / "convergence" / folder
    selection = pd.read_csv(here / "earliest_converged_one_per_seed_model.csv", low_memory=False)
    shipped = sorted(p.name for p in (PACKAGE / "models" / folder).glob("*.pdb"))
    listed = sorted(Path(str(p)).name for p in selection["pdb_file"])

    live = ROOT / "Nav1.5" / CONDITIONS[folder][2]
    checks = {
        "package_selection_matches_shipped_models": listed == shipped,
        "package_all_ok_matches_hive": (
            digest(here / "all_ok_models.csv") == digest(live / "all_ok_models.csv")),
        "package_selection_matches_hive": (
            digest(here / "earliest_converged_one_per_seed_model.csv")
            == digest(live / "earliest_converged_one_per_seed_model.csv")),
    }
    # the chain reads earliest_converged_selected out of all_ok_models.csv; it must be this same set
    all_ok = pd.read_csv(live / "all_ok_models.csv", low_memory=False)
    chosen = all_ok["earliest_converged_selected"].astype(str).str.strip().str.lower() == "true"
    checks["selection_equals_earliest_converged_selected"] = (
        sorted(Path(str(p)).name for p in all_ok.loc[chosen, "pdb_file"]) == listed)
    return set(shipped), checks


def compare(dataset: str, verdicts: pd.DataFrame) -> tuple:
    rows, summary, checks = [], [], {}
    for folder, (ds, condition, _live) in CONDITIONS.items():
        if ds != dataset:
            continue
        old, folder_checks = old_set(folder)
        checks[folder] = folder_checks

        here = verdicts[verdicts["condition"] == condition]
        passed = here[here["status"] == "passed"]
        new = set(passed["representative_model_file"])
        by_representative = here.set_index("representative_model_file")

        for name in sorted(old | new):
            in_old, in_new = name in old, name in new
            row = by_representative.loc[name] if name in by_representative.index else None
            rows.append({
                "dataset": dataset, "condition": condition, "model_file": name,
                "model": None if row is None else row["model"],
                "seed": None if row is None else row["seed"],
                "recycle": recycle_of(name),
                "in_november_package": in_old, "in_new_filter_chain": in_new,
                "change": ("unchanged_kept" if in_old and in_new else
                           "removed_by_new_chain" if in_old else "added_by_new_chain"),
                "removed_by": "" if row is None else row["removed_by"],
                "removed_because": "" if row is None else row["removed_because"],
            })
        summary.append({
            "dataset": dataset, "condition": condition,
            "november_package_models": len(old), "new_chain_models": len(new),
            "in_both": len(old & new), "removed_by_new_chain": len(old - new),
            "added_by_new_chain": len(new - old),
            "november_recycle0": sum(recycle_of(n) == 0 for n in old),
            "new_recycle0": sum(recycle_of(n) == 0 for n in new),
        })
    return pd.DataFrame(rows), pd.DataFrame(summary), checks


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, default=CHAIN)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    summaries, failed = [], []
    for dataset in sorted({d for d, _c, _l in CONDITIONS.values()}):
        verdicts = pd.read_csv(CHAIN / f"{dataset}_filter_verdicts.csv", low_memory=False)
        table, summary, checks = compare(dataset, verdicts)
        for folder, result in checks.items():
            for name, ok in result.items():
                print(f"  provenance {folder} {name}: {'OK' if ok else 'MISMATCH'}", flush=True)
                if not ok:
                    failed.append(f"{folder}/{name}")
        path = args.out / f"{dataset}_old_vs_new.csv"
        table.to_csv(path, index=False)
        print(f"wrote {path} ({len(table)} models)", flush=True)
        summaries.append(summary)

    for dataset in NOT_IN_PACKAGE:
        print(f"{dataset}: no November package - that protocol was never sent, so there is nothing "
              f"to compare it against", flush=True)

    summary = pd.concat(summaries, ignore_index=True)
    path = args.out / "nav15_old_vs_new_summary.csv"
    summary.to_csv(path, index=False)
    print(summary.to_string(index=False), flush=True)
    print(f"wrote {path}", flush=True)
    if failed:
        print(f"PROVENANCE CHECKS FAILED: {failed}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
