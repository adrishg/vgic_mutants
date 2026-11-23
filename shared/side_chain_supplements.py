"""One atom definition for the Nav1.5 IFM and Cav1.2 R406 contact distances.

The vanilla ``shortest_*`` columns (2026-02 and 2026-07 exports) were written by
``measure_shortest_distance`` of the cluster analysis package, which uses only
the terminal atoms of each side chain (Arg NH1/NH2, Asp OD1/OD2, Phe CZ, ...).
The November 2026 masked exports used every non-hydrogen atom. Vanilla-versus-
masked contact comparisons therefore mixed two definitions.

Both protocols were recomputed from the predicted structures with one rule,
side-chain non-hydrogen atoms (every non-hydrogen atom except N, CA, C, O and
OXT). The slim tables live in ``<channel>/dataDistances/side_chain_20261123/``
of the data root and hold every model of a cohort; the retained models are
still the rows of the final-QC table they are applied to.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

from . import data_access

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

_PATTERNS = (
    (
        re.compile(r"cav1\.?2_g406r_(vanilla|masked_v3_50pct)", re.I),
        "cav12/dataDistances/side_chain_20261123/Cav12_g406r_{arm}_R406_side_chain.csv",
        None,
    ),
    (
        re.compile(r"nav15_(wt|qqq)_(vanilla|masked_v1_50pct|masked_v2_50pct)", re.I),
        "nav15/dataDistances/side_chain_20261123/Nav15_{sequence}_{arm}_IFM_side_chain.csv",
        "sequence",
    ),
)


def supplement_path(table_name: str) -> str | None:
    """Return the data-root path of the side-chain table for a distance table."""
    name = Path(str(table_name)).name
    for pattern, template, leading in _PATTERNS:
        match = pattern.search(name)
        if not match:
            continue
        groups = [group.lower() for group in match.groups()]
        if leading:
            return template.format(sequence=groups[0], arm=groups[1])
        return template.format(arm=groups[0])
    return None


def apply_side_chain_supplement(
    frame: pd.DataFrame, table_name: str, *, add_missing: bool = True
) -> pd.DataFrame:
    """Replace the contact columns of ``frame`` by their side-chain values.

    ``table_name`` is the file the frame was read from. Frames of tables that
    have no supplement, or that carry no ``pdb_file`` column, are returned
    unchanged. Every row of ``frame`` must be present in the supplement.
    """
    relative = supplement_path(table_name)
    if relative is None or "pdb_file" not in frame.columns:
        return frame
    table = pd.read_csv(data_access.resolve(relative))
    table.index = table.pop("pdb_file").map(lambda value: Path(str(value)).name)
    models = frame["pdb_file"].map(lambda value: Path(str(value)).name)
    absent = set(models) - set(table.index)
    if absent:
        raise KeyError(
            f"{len(absent)} models of {Path(str(table_name)).name} are not in {relative}"
        )
    columns = [
        column for column in table.columns if add_missing or column in frame.columns
    ]
    values = table.loc[models, columns].set_axis(frame.index)
    out = frame.drop(columns=[column for column in columns if column in frame.columns])
    return pd.concat([out, values], axis=1)


def side_chain_atoms(pdb_path: str | Path, chain: str = "A") -> dict:
    """Side-chain non-hydrogen coordinates keyed by ``(residue name, number)``."""
    residues: dict[tuple[str, int], list[tuple[float, float, float]]] = {}
    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[21] != chain or line[16] not in " A":
                continue
            name = line[12:16].strip()
            element = line[76:78].strip().upper() or name.lstrip("0123456789")[:1].upper()
            if element == "H" or name in BACKBONE_ATOMS:
                continue
            key = (line[17:20].strip(), int(line[22:26]))
            residues.setdefault(key, []).append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            )
    return residues


def shortest_side_chain_distance(residues: dict, first, second) -> float:
    """Minimum distance between the side-chain atoms of two residue keys."""
    if first not in residues or second not in residues:
        return float("nan")
    return min(math.dist(a, b) for a in residues[first] for b in residues[second])
