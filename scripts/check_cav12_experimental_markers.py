#!/usr/bin/env python3
"""Recompute the Cav1.2 experimental C-alpha markers by sequence alignment and compare them with
the literals the cross-channel figure reads from Cav12_distanceDistribution_vsExperimental.ipynb.

Why this exists. The modeled construct is UniProt Q13936 residues 1-1685, which carries both members
of two mutually exclusive exon pairs in tandem (model 932-951 / 952-971 and 1297-1324 / 1325-1352).
8FD7 chain K holds one member of each pair and is numbered continuously, so its numbering against
the model is piecewise: +0 for model 112-931, -20 for model 952-1296, -48 for model 1325-1658. The
blocks follow the exons, not the repeats: DIII S1 is still +0 and DIV S1-S2 is -20. 8WE6 and 8HLP
keep Q13936 numbering. A model residue inside an exon copy the structure lacks (F940 in 8FD7 and
8WE6) has no sequence-identical counterpart and is reported as unmapped rather than guessed.

Usage: python scripts/check_cav12_experimental_markers.py [--output table.csv]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from shared import experimental_overlays  # noqa: E402

CHAINS = {"8HLP": "A", "8WE6": "A", "8FD7": "K"}
LITERAL_INDEX = {"8HLP": 0, "8WE6": 7, "8FD7": 8}   # position in each exp_distances_list
MODEL_LENGTH = 1685
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
         "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def model_sequence() -> str:
    lines = (ROOT / "datasets/uniprot/Q13936.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))[:MODEL_LENGTH]


def ca_atoms(path: Path, chain: str):
    numbers, letters, coords = [], [], {}
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA" \
                and line[16] in " A":
            number = int(line[22:26])
            if number in coords:
                continue
            numbers.append(number)
            letters.append(THREE.get(line[17:20], "X"))
            coords[number] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return numbers, "".join(letters), coords


def model_to_structure(model: str, numbers, sequence: str) -> dict:
    """Model position -> structure residue number, identical residues only."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score, aligner.end_gap_score = -12, -1, 0
    alignment = aligner.align(model, sequence)[0]
    mapping = {}
    for (m0, m1), (s0, s1) in zip(*alignment.aligned):
        for step in range(m1 - m0):
            if model[m0 + step] == sequence[s0 + step]:
                mapping[m0 + step + 1] = numbers[s0 + step]
    return mapping


def offset_blocks(mapping: dict):
    blocks = []
    for position in sorted(mapping):
        offset = mapping[position] - position
        if blocks and blocks[-1][2] == offset:
            blocks[-1][1] = position
        else:
            blocks.append([position, position, offset])
    return [block for block in blocks if block[1] - block[0] >= 5]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    model = model_sequence()
    notebook = ROOT / "cav12/Cav12_distanceDistribution_vsExperimental.ipynb"
    literal_lists = experimental_overlays._literal_assignments(notebook, "exp_distances_list")
    rows = []
    for pdb_id, chain in CHAINS.items():
        numbers, sequence, coords = ca_atoms(ROOT / f"cav12/experimental/{pdb_id}.pdb", chain)
        mapping = model_to_structure(model, numbers, sequence)
        print(pdb_id, "chain", chain, "offset blocks (model first, model last, offset):",
              offset_blocks(mapping))
        for literals in literal_lists:
            for alias, values in literals[LITERAL_INDEX[pdb_id]].items():
                first, second = (int(re.sub(r"\D", "", part)) for part in alias.split("-"))
                a, b = mapping.get(first), mapping.get(second)
                recomputed = (float(np.linalg.norm(coords[a] - coords[b]))
                              if a is not None and b is not None else np.nan)
                rows.append({"structure": pdb_id, "alias_model_numbering": alias,
                             "structure_residue_1": a, "structure_residue_2": b,
                             "notebook_literal_A": float(values[0]),
                             "recomputed_A": round(recomputed, 2) if np.isfinite(recomputed) else np.nan})
    table = pd.DataFrame(rows)
    table["difference_A"] = (table.recomputed_A - table.notebook_literal_A).round(2)
    table["status"] = np.where(table.recomputed_A.isna(), "model residue has no counterpart",
                               np.where(table.difference_A.abs() <= 0.06, "agrees", "DIFFERS"))
    print(table.to_string(index=False))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
