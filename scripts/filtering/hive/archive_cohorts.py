#!/usr/bin/env python3
"""The eighteen condition directories that step 1 of Brandon's pipeline needs as ColabFold archives.

His build_dv.py reads each condition with `tarfile.open(path, "r|gz")` - a real .tar.gz, streamed,
one ColabFold job per archive. Our predictions are loose PDB trees on Quobyte, so each condition is
packed once here and the configs point at the result.

Only the saved-recycle files go in (`*.r<N>.pdb`). The unsuffixed final model of each trajectory is
left out: it duplicates the last recycle, and step 1 reads recycles.

  python archive_cohorts.py list
  python archive_cohorts.py pack --index $SLURM_ARRAY_TASK_ID
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path("/quobyte/yarovoygrp/ahgz/vgic_mutants")
OUT = ROOT / "analysis" / "archives"

# (channel folder, condition, protocol folder, dataset the config will call it)
COHORTS = [
    ("Kv2.1", "wt", "vanilla", "kv21_vanilla"),
    ("Kv2.1", "l403a", "vanilla", "kv21_vanilla"),
    ("Kv2.1", "f412l", "vanilla", "kv21_vanilla"),
    ("Kv2.1", "wt", "masked_common99_50pct", "kv21_masked99"),
    ("Kv2.1", "l403a", "masked_common99_50pct", "kv21_masked99"),
    ("Kv2.1", "f412l", "masked_common99_50pct", "kv21_masked99"),
    ("Nav1.5", "wt", "vanilla", "nav15_vanilla"),
    ("Nav1.5", "qqq", "vanilla", "nav15_vanilla"),
    ("Nav1.5", "wt", "masked_v1_50pct", "nav15_masked_v1"),
    ("Nav1.5", "qqq", "masked_v1_50pct", "nav15_masked_v1"),
    ("Nav1.5", "wt", "masked_v2_50pct", "nav15_masked_v2"),
    ("Nav1.5", "qqq", "masked_v2_50pct", "nav15_masked_v2"),
    ("Cav1.2", "wt", "vanilla", "cav12_vanilla"),
    ("Cav1.2", "g402s", "vanilla", "cav12_vanilla"),
    ("Cav1.2", "g406r", "vanilla", "cav12_vanilla"),
    ("Cav1.2", "wt", "masked_v3_50pct", "cav12_masked_v3"),
    ("Cav1.2", "g402s", "masked_v3_50pct", "cav12_masked_v3"),
    ("Cav1.2", "g406r", "masked_v3_50pct", "cav12_masked_v3"),
]


def entry(i: int) -> dict:
    channel, condition, protocol, dataset = COHORTS[i]
    source = ROOT / channel / condition / protocol
    models = source / "models" if (source / "models").is_dir() else source
    return {"channel": channel, "condition": condition, "protocol": protocol, "dataset": dataset,
            "models": models, "archive": OUT / f"{dataset}_{condition}.tar.gz"}


def pack(i: int) -> None:
    e = entry(i)
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted(p.relative_to(e["models"]) for p in e["models"].rglob("*.pdb")
                   if ".r" in p.name.rsplit("_", 1)[-1])
    if not files:
        raise SystemExit(f"{e['archive'].name}: no recycle files under {e['models']}")
    jobs = {p.name.split("_unrelaxed")[0] for p in files}
    if len(jobs) != 1:
        raise SystemExit(f"{e['archive'].name}: {len(jobs)} ColabFold jobs {sorted(jobs)[:3]}; "
                         f"step 1 wants one per archive (D-019)")
    print(f"{e['archive'].name}: {len(files)} recycle files, job {jobs.pop()}", flush=True)
    listing = e["archive"].with_suffix(".files")        # one writer per file, and tar reads it
    listing.write_text("\n".join(str(p) for p in files) + "\n")
    # Written to a .part first so a killed job never leaves a half archive that looks finished.
    part = e["archive"].with_suffix(".part")
    subprocess.run(["tar", "-czf", str(part), "-C", str(e["models"]), "-T", str(listing)], check=True)
    part.rename(e["archive"])
    listing.unlink()
    print(f"{e['archive'].name}: {e['archive'].stat().st_size / 2**30:.2f} GiB", flush=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mode", choices=("pack", "list"))
    parser.add_argument("--index", type=int)
    args = parser.parse_args(argv)
    if args.mode == "list":
        for i in range(len(COHORTS)):
            e = entry(i)
            print(f"  {i:2d}  {e['dataset']:16s} {e['condition']:7s} "
                  f"{'models ok' if e['models'].is_dir() else 'MISSING':10s} "
                  f"{'packed' if e['archive'].exists() else '-'}")
        return 0
    if args.index is None:
        parser.error("pack needs --index")
    pack(args.index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
