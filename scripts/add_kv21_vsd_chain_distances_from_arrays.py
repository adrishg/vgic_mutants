#!/usr/bin/env python3
"""Add chain-resolved Kv2.1 F236–R308 Cα distances to the 50% QC tables.

The distance exports contain only the chain-A F236–R308 (model 238–310)
measurement. Intrachain Cα–Cα distances are invariant to the common
superposition, so the missing B–D values are taken from the merged aligned-Cα
arrays of the 50% release. Canonical chains are translated back to model chain
labels with the recorded ``chain_mapping`` (``model->canonical``), and the
recomputed chain-A value is checked against the exported column before any
column is written.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

FIRST, SECOND = 238, 310
SUFFIXES = ("_all_ok_rmsd_3A_structural_interface_qc.csv",)


def column(chain: str) -> str:
    return f"CA_CA_{chain}_ARG310_CA-{chain}_PHE238_CA"


def main() -> None:
    root = Path(os.environ["VGIC_DATA_ROOT"]) / "kv21"
    merged = root / "dataRMSF" / "merged"
    meta = pd.read_csv(merged / "kv21_alignment_metadata.csv", usecols=["pdb_file", "chain_mapping", "alignment_success"], low_memory=False)
    coords = np.load(merged / "kv21_aligned_ca_coordinates.npy", mmap_mode="r")
    present = np.load(merged / "kv21_aligned_ca_present.npy", mmap_mode="r")
    raw = np.load(merged / "kv21_raw_residue_numbers.npy")
    canonical = list(np.load(merged / "kv21_canonical_chains.npy", allow_pickle=True).astype(str))
    i1, i2 = int(np.where(raw == FIRST)[0][0]), int(np.where(raw == SECOND)[0][0])

    values = {}
    for index, row in meta.iterrows():
        if not bool(row.alignment_success):
            continue
        xyz = np.asarray(coords[index][:, [i1, i2]], dtype=float)
        ok = np.asarray(present[index][:, [i1, i2]], bool).all(axis=1)
        dist = np.where(ok, np.linalg.norm(xyz[:, 0] - xyz[:, 1], axis=1), np.nan)
        model_of = {canon: model for model, canon in (pair.split("->") for pair in str(row.chain_mapping).split(";"))}
        values[Path(str(row.pdb_file)).name] = {model_of[c]: dist[k] for k, c in enumerate(canonical)}

    for path in sorted((root / "dataDistances").glob("26-*_Kv2.1_wt_*" + SUFFIXES[0])):  # Fig. 3 uses WT only
        if "_extra_" in path.name:
            continue
        frame = pd.read_csv(path, low_memory=False)
        names = frame["pdb_file"].map(lambda p: Path(str(p)).name)
        computed = pd.DataFrame([values.get(n, {}) for n in names])
        check = (computed["A"] - pd.to_numeric(frame[column("A")], errors="coerce")).abs()
        if check.max() > 0.02:  # exported distances are rounded to 0.01 Å
            raise ValueError(f"{path.name}: chain-A check failed (max |Δ|={check.max():.3f} Å)")
        for chain in "BCD":
            frame[column(chain)] = computed[chain].round(2).to_numpy()
        frame.to_csv(path, index=False)
        print(f"{path.name}: chain-A max |Δ| {check.max():.3f} Å; added B–D for {len(frame)} rows")


if __name__ == "__main__":
    main()
