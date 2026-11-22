#!/usr/bin/env python3
"""Build the Kv2.1 interface-QC tables for whichever masked cohort is registered.

``build_kv21_structural_interface_qc_50pct.py`` names its six input bases in a
literal table, which was fine while one masked cohort was ever current. The
99-position replacement run changed those filenames, so this driver reads the
masked bases from ``analysis/active_masked_cohorts.tsv`` instead and reuses the
other script unchanged for the filtering itself. The vanilla bases are not
registered anywhere and stay as that script defines them.

Usage mirrors the script it wraps::

    VGIC_DATA_ROOT=/path/to/data python3 \\
        scripts/filtering/build_kv21_structural_interface_qc_from_registry.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_kv21_structural_interface_qc_50pct as builder  # noqa: E402
from shared.dataset_selection import kv21_active_masked_distance_path  # noqa: E402


def registered_masked_bases() -> dict[tuple[str, str], str]:
    """The masked base stem per condition, as the active-cohort registry gives it."""
    bases = {}
    for condition in ("WT", "L403A", "F412L"):
        path = Path(kv21_active_masked_distance_path(REPO_ROOT, condition))
        stem = path.stem
        while stem.endswith("_all"):
            stem = stem[:-4]
        bases[(condition, "masked")] = stem
    return bases


def main() -> None:
    for key, base in registered_masked_bases().items():
        builder.BASES[key] = base
    for (condition, protocol), base in sorted(builder.BASES.items()):
        print(f"base {condition:6s} {protocol:8s} {base}")
    builder.main()


if __name__ == "__main__":
    main()
