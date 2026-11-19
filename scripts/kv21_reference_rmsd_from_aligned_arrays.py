#!/usr/bin/env python3
"""Kv2.1 model-to-reference Cα RMSD (8SD3, 8SDA) from the aligned 50% arrays.

Implements the Supplementary Table S3 definitions on the merged aligned-Cα
arrays of the 50% release (all six cohorts, chain correspondence already
resolved by the cyclic, adjacency-preserving mapping):

* each tetramer is fitted once to the reference on the 184–422 Cα core of all
  four chains (``selected_core_postfit_rmsd_A``);
* ``whole_matched_tetramer`` is the Cα RMSD over every residue resolved in both
  model and reference, without refitting;
* regional ``__ca__core_aligned`` values reuse that same transformation;
  ``__ca__local_aligned`` values refit on the region itself.

Backbone (N/C/O) variants are not available from Cα arrays and are omitted.
Residue numbers are model/construct numbering (displayed + 2).

Output: ``$VGIC_DATA_ROOT/kv21/dataRMSD/Kv21_50pct_all_models_vs_8SD3_8SDA_RMSD.csv.gz``
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

CORE = [(184, 422)]
REGIONS = {
    "s6_bundle_working": [(389, 422)],
    "distal_s6_working": [(398, 414)],
    "l403_region": [(398, 409)],
    "f412_region": [(407, 417)],
    "hydrophobic_nexus": [(314, 318), (327, 331), (401, 405), (410, 414)],
}


def residue_mask(raw: np.ndarray, ranges) -> np.ndarray:
    mask = np.zeros(raw.shape, bool)
    for start, end in ranges:
        mask |= (raw >= start) & (raw <= end)
    return mask


def kabsch(moving: np.ndarray, fixed: np.ndarray):
    mc, fc = moving.mean(0), fixed.mean(0)
    h = (moving - mc).T @ (fixed - fc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, fc - mc @ rot.T


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1)))) if len(a) else np.nan


def main() -> None:
    root = Path(os.environ["VGIC_DATA_ROOT"]) / "kv21"
    merged = root / "dataRMSF" / "merged"
    meta = pd.read_csv(merged / "kv21_alignment_metadata.csv", low_memory=False)
    coords = np.load(merged / "kv21_aligned_ca_coordinates.npy", mmap_mode="r")
    present = np.load(merged / "kv21_aligned_ca_present.npy", mmap_mode="r")
    refs = np.load(root / "dataRMSF" / "references" / "kv21_aligned_references.npz", allow_pickle=False)
    raw = refs["raw_residue_numbers"].astype(int)
    chains = refs["canonical_chains"].astype(str)
    core_res = residue_mask(raw, CORE)
    region_res = {name: residue_mask(raw, r) for name, r in REGIONS.items()}

    rows = []
    for i, m in meta.iterrows():
        base = {
            "dataset": m.dataset, "sequence_condition": m.sequence_condition,
            "protocol": "masked" if "mask" in str(m.protocol) else "vanilla",
            "pdb_file": Path(str(m.pdb_file)).name, "model_path": m.model_path,
            "mapping_core_ca_rmsd_diagnostic_A_value": m.best_mapping_core_ca_rmsd_A,
        }
        if not bool(m.alignment_success):
            for ref_id in refs["reference_ids"].astype(str):
                rows.append({**base, "reference_id": ref_id, "analysis_status": "alignment_failed",
                             "analysis_error": m.alignment_error})
            continue
        xyz = np.asarray(coords[i], dtype=np.float64)          # (4, n, 3)
        ok = np.asarray(present[i], bool) & np.isfinite(xyz).all(-1)
        for r, ref_id in enumerate(refs["reference_ids"].astype(str)):
            ref = refs["coords"][r].astype(np.float64)
            both = ok & refs["present"][r].astype(bool)
            core = both & core_res[None, :]
            rot, shift = kabsch(xyz[core], ref[core])
            fitted = xyz @ rot.T + shift
            row = {**base, "reference_id": ref_id, "analysis_status": "ok", "analysis_error": "",
                   "selected_core_postfit_rmsd_A": rmsd(fitted[core], ref[core]),
                   "selected_core_matched_atoms": int(core.sum()),
                   "whole_matched_tetramer__ca__core_aligned_rmsd_A": rmsd(fitted[both], ref[both]),
                   "whole_matched_tetramer__ca__core_aligned_matched_atoms": int(both.sum())}
            for name, res in region_res.items():
                sel = both & res[None, :]
                row[f"{name}__ca__core_aligned_rmsd_A"] = rmsd(fitted[sel], ref[sel])
                row[f"{name}__ca__core_aligned_matched_atoms"] = int(sel.sum())
                if sel.sum() >= 3:
                    lr, ls = kabsch(xyz[sel], ref[sel])
                    row[f"{name}__ca__local_aligned_rmsd_A"] = rmsd(xyz[sel] @ lr.T + ls, ref[sel])
                for c, chain in enumerate(chains):
                    sel_c = sel[c]
                    row[f"{name}__chain_{chain}__ca__core_aligned_rmsd_A"] = rmsd(fitted[c][sel_c], ref[c][sel_c])
            rows.append(row)
        if i % 5000 == 0:
            print(f"{i}/{len(meta)}", flush=True)

    out = root / "dataRMSD"
    out.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows)
    path = out / "Kv21_50pct_all_models_vs_8SD3_8SDA_RMSD.csv.gz"
    table.to_csv(path, index=False)
    print(path, table.shape)


if __name__ == "__main__":
    main()
