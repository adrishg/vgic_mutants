#!/usr/bin/env python3
"""Stage the six current Kv2.1 cohorts for the conformation pipeline (v2.0.1).

Writes one manifest (three vanilla cohorts, three masked_common99_50pct cohorts)
and links every per-structure part that already exists with status ``ok`` into
the new output directory, so the array job only computes what is missing.
Parts are keyed by sha1 of the resolved model path, exactly as the pipeline does.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

ROOT = Path("/quobyte/yarovoygrp/ahgz/vgic_mutants")
KV = ROOT / "Kv2.1"
OUT = Path(sys.argv[1])
PARTS = OUT / "parts"
PARTS.mkdir(parents=True, exist_ok=True)

# Existing part stores, in order of preference.
STORES = [
    KV / "analysis/kv21_conformation_all_v201/parts",           # vanilla, August run
    KV / "conformation_c99/l403a_masked_common99_50pct/parts",  # L403A 99 positions
    KV / "conformation_c99/wt_masked_common99_50pct/parts",     # WT 99 positions (failed run)
]

rows = []
for condition in ("wt", "l403a", "f412l"):
    vanilla = KV / condition / "vanilla"
    found = sorted(p for p in vanilla.glob("*.pdb"))
    rows += [(condition, "vanilla", str(p), str(vanilla)) for p in found]
    masked = KV / condition / "masked_common99_50pct" / "models"
    found = sorted(masked.rglob("*.pdb"))
    rows += [(condition, "masked", str(p), str(masked)) for p in found]

with (OUT / "manifest_six_cohorts.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["condition", "protocol", "source_path", "discovery_root"])
    writer.writerows(rows)

counts = Counter((r[0], r[1]) for r in rows)
linked, missing, bad = Counter(), Counter(), Counter()
targets = [(r[0], r[1], hashlib.sha1(str(Path(r[2]).resolve()).encode()).hexdigest()[:16]) for r in rows]
targets += [("experimental", "8SD3", "8SD3"), ("experimental", "8SDA", "8SDA")]
for condition, protocol, sid in targets:
    key = (condition, protocol)
    dest = PARTS / f"{sid}.json"
    if dest.exists() or dest.is_symlink():
        linked[key] += 1
        continue
    for store in STORES:
        source = store / f"{sid}.json"
        if not source.is_file():
            continue
        try:
            status = json.load(source.open()).get("status")
        except Exception:  # unreadable part: recompute
            status = "unreadable"
        if status == "ok":
            os.symlink(source, dest)
            linked[key] += 1
            break
        bad[key] += 1
    else:
        missing[key] += 1

report = {
    "manifest_rows": {f"{k[0]}_{k[1]}": v for k, v in sorted(counts.items())},
    "parts_linked_ok": {f"{k[0]}_{k[1]}": v for k, v in sorted(linked.items())},
    "parts_to_compute": {f"{k[0]}_{k[1]}": v for k, v in sorted(missing.items())},
    "existing_parts_rejected_not_ok": {f"{k[0]}_{k[1]}": v for k, v in sorted(bad.items())},
}
(OUT / "prepare_report.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
