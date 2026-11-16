#!/usr/bin/env python3
"""Generate NaV1.5 IFM/receptor distance tables from local PDB files.

The model numbering is the sequence-verified AlphaFold numbering used by the
publication tables: IFM F/Q1170 to N1343 and N1449. Shortest distances use
non-hydrogen atoms and the companion C-alpha columns use the named C-alpha
atoms. The script has no dependency on a cluster-only analysis package.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


def parse_pdb(path: Path, chain: str) -> tuple[dict[tuple[str, int], np.ndarray], dict[tuple[str, int], np.ndarray]]:
    heavy: dict[tuple[str, int], list[np.ndarray]] = {}
    ca: dict[tuple[str, int], np.ndarray] = {}
    with path.open(errors="ignore") as handle:
        for line in handle:
            if line[:6].strip() not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            if line[16] not in {" ", "A"} or (chain and line[21].strip() != chain):
                continue
            try:
                key = (line[17:20].strip(), int(line[22:26]))
                atom = line[12:16].strip()
                element = line[76:78].strip() or atom[:1]
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except (IndexError, ValueError):
                continue
            if element.upper() != "H":
                heavy.setdefault(key, []).append(xyz)
            if atom == "CA":
                ca[key] = xyz
    return {key: np.vstack(value) for key, value in heavy.items()}, ca


def minimum_distance(first: np.ndarray, second: np.ndarray) -> float:
    delta = first[:, None, :] - second[None, :, :]
    return float(np.sqrt(np.einsum("ijk,ijk->ij", delta, delta).min()))


def calculate(path: Path, chain: str, ifm_residue: str) -> dict[str, object]:
    heavy, ca = parse_pdb(path, chain)
    row: dict[str, object] = {"pdb_file": path.name}
    for receptor in (1343, 1449):
        first = heavy.get((ifm_residue, 1170))
        second = heavy.get(("ASN", receptor))
        first_ca = ca.get((ifm_residue, 1170))
        second_ca = ca.get(("ASN", receptor))
        short_label = f"shortest_{ifm_residue}1170-ASN{receptor}"
        ca_label = f"CA_{ifm_residue}1170_CA-ASN{receptor}_CA"
        row[short_label] = minimum_distance(first, second) if first is not None and second is not None else np.nan
        row[ca_label] = math.dist(first_ca, second_ca) if first_ca is not None and second_ca is not None else np.nan
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--ifm-residue", choices=("PHE", "GLN"), required=True)
    parser.add_argument("--chain", default="A")
    args = parser.parse_args()
    paths = sorted(args.input_dir.rglob("*.pdb"))
    if not paths:
        raise SystemExit(f"No PDB files found below {args.input_dir}")
    rows = [calculate(path, args.chain, args.ifm_residue) for path in paths]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
