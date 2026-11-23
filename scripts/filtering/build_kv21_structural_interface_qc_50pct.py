#!/usr/bin/env python3
"""Build Kv2.1 structural/interface-QC distance tables for the 50% release.

The same three filters are applied to every vanilla and masked_v3_50pct
cohort so that both protocols share one definition:

1. All-OK-3 convergence (the ``*_all_ok_rmsd_3A.csv`` input tables).
2. G377 selectivity-filter ring integrity: a trajectory is rejected when any
   retained recycle has a G377 interchain C-alpha distance above
   ``--g377-max`` (default 20 Å; no vanilla or 50% trajectory exceeds 16.6 Å).
3. Pore–VSD interface integrity: ``shared.dataset_selection.apply_kv21_interface_qc``
   with the published 27 Å trajectory-level cutoff.

The stable-core alignment QC used for the earlier masked cohort requires the
Hive RMSD v2 table and is not applied here, to either protocol.

Outputs are written beside the inputs in ``$VGIC_DATA_ROOT/kv21/dataDistances``:
``<base>_all_ok_rmsd_3A_structural_qc.csv``,
``<base>_all_ok_rmsd_3A_structural_interface_qc.csv`` and the matching
``<base>_extra_all_ok_rmsd_3A_structural_interface_qc.csv`` companion tables.
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

from shared import data_access  # noqa: E402
from shared.dataset_selection import (  # noqa: E402
    apply_kv21_interface_qc,
    kv21_active_masked_distance_path,
)

BASES = {
    ("WT", "vanilla"): "26-02-11_Kv2.1_wt_vanillaAF2_distances",
    ("L403A", "vanilla"): "26-02-11_Kv2.1_l403a_vanillaAF2_distances",
    ("F412L", "vanilla"): "26-02-11_Kv2.1_f412l_vanillaAF2_distances",
}
# The masked bases come from the registry, so the builder follows the active
# cohort instead of naming one.
for _condition in ("WT", "L403A", "F412L"):
    BASES[(_condition, "masked")] = kv21_active_masked_distance_path(
        REPO_ROOT, _condition
    ).stem


def trajectory(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.r\d+\.pdb$", ".pdb", regex=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--g377-max", type=float, default=20.0)
    parser.add_argument("--interface-max", type=float, default=27.0)
    args = parser.parse_args()

    root = os.environ.get("VGIC_DATA_ROOT")
    if not root:
        raise SystemExit("Set VGIC_DATA_ROOT to the external data root.")
    out_dir = Path(root) / "kv21" / "dataDistances"
    out_dir.mkdir(parents=True, exist_ok=True)

    audit = []
    for (condition, protocol), base in BASES.items():
        main_rel = f"kv21/dataDistances/{base}_all_ok_rmsd_3A.csv"
        frame = pd.read_csv(data_access.resolve(main_rel), low_memory=False)
        # Masked releases ship a filtered companion; the vanilla companion is
        # the complete 26-02-10 export and is restricted below by pdb_file.
        extra_candidates = [
            f"kv21/dataDistances/{base}_extra_all_ok_rmsd_3A.csv",
            f"kv21/dataDistances/{base.replace('26-02-11', '26-02-10')}_extra.csv",
        ]
        for extra_rel in extra_candidates:
            try:
                extra = pd.read_csv(data_access.resolve(extra_rel), low_memory=False)
                break
            except (data_access.DataUnavailable, FileNotFoundError):
                continue
        else:
            raise FileNotFoundError(f"{base}: no companion extra table found")
        frame.attrs["dataset_label"] = f"Kv2.1 {condition} {protocol}"

        g377 = [c for c in frame.columns if c.startswith("CA_CA_") and c.count("GLY377") == 2]
        if len(g377) != 6:
            raise KeyError(f"{base}: expected six G377 ring columns, found {len(g377)}")
        traj = trajectory(frame["pdb_file"])
        ring = frame[g377].apply(pd.to_numeric, errors="coerce").max(axis=1)
        bad_ring = set(traj[ring.gt(args.g377_max)])
        structural = frame.loc[~traj.isin(bad_ring)].copy()
        structural.attrs.update(frame.attrs)
        structural.to_csv(out_dir / f"{base}_all_ok_rmsd_3A_structural_qc.csv", index=False)

        interface = apply_kv21_interface_qc(
            structural, threshold_A=args.interface_max, trajectory_level=True
        )
        interface.to_csv(
            out_dir / f"{base}_all_ok_rmsd_3A_structural_interface_qc.csv", index=False
        )
        keep = set(interface["pdb_file"])
        extra_qc = extra.loc[extra["pdb_file"].isin(keep)]
        if len(extra_qc) != len(interface):
            raise ValueError(f"{base}: extra table does not cover every retained model")
        extra_qc.to_csv(
            out_dir / f"{base}_extra_all_ok_rmsd_3A_structural_interface_qc.csv", index=False
        )
        audit.append({
            "condition": condition,
            "protocol": protocol,
            "all_ok_3_rows": len(frame),
            "all_ok_3_trajectories": int(traj.nunique()),
            "g377_rejected_trajectories": len(bad_ring),
            "interface_rejected_trajectories": int(
                interface.attrs["kv21_interface_qc_rejected_trajectories"]
            ),
            "final_rows": len(interface),
            "final_trajectories": int(trajectory(interface["pdb_file"]).nunique()),
        })

    summary = pd.DataFrame(audit)
    summary.to_csv(out_dir / "kv21_50pct_structural_interface_qc_summary.csv", index=False)
    print(summary.to_string(index=False))
    (out_dir / "kv21_50pct_structural_interface_qc_settings.json").write_text(
        json.dumps({"g377_max_A": args.g377_max, "interface_max_A": args.interface_max}, indent=2)
    )


if __name__ == "__main__":
    main()
