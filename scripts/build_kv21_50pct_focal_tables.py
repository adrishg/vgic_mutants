#!/usr/bin/env python3
"""Kv2.1 focal tables for the 50% release, built from the filtered distance tables.

Both protocols use the ``*_all_ok_rmsd_3A_structural_interface_qc`` tables
(All-OK-3 convergence, G377 ring and 27 Å pore–VSD interface checks).

Outputs (in ``$VGIC_DATA_ROOT``):

* ``kv21/dataDistances/analysis/L403A_E423_N179_all_structure_distances.csv``:
  per-structure E423–N179 Cα distance in every subunit (model labels
  GLU425/ASN181), with the maximum and the number of subunits above the
  experimentally derived 12.841 Å cutoff.
* ``kv21/dataDistances/analysis/F412L_nexus_shortest_contacts_long.csv``:
  per-structure shortest heavy-atom distance from L412 (model 414) to the
  same-subunit L316 (318) and to the closest neighboring-subunit L329 (331)
  and L403 (405), using the column selection of ``shared.kv21_nexus_analysis``.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from shared import data_access  # noqa: E402
from shared.dataset_selection import filtered_masked_distance_path  # noqa: E402
from shared.kv21_nexus_analysis import _nexus_distance_columns  # noqa: E402

SHIFT_CUTOFF_A = 12.84133701269344
L403_COLUMNS = [f"CA_CA_{c}_GLU425_CA-{c}_ASN181_CA" for c in "ABCD"]
VANILLA = "kv21/dataDistances/26-02-11_Kv2.1_{}_vanillaAF2_distances_all_ok_rmsd_3A_structural_interface_qc.csv"
CONTACT_LABELS = {
    318: "L412–L316\n(same subunit)",
    331: "L412–L329\n(neighbor)",
    405: "L412–L403\n(neighbor)",
}


def source(condition: str, protocol: str) -> Path:
    if protocol == "vanilla":
        return data_access.resolve(VANILLA.format(condition.lower()))
    return data_access.resolve(filtered_masked_distance_path(REPO, "kv21", condition))


def main() -> None:
    out = Path(os.environ["VGIC_DATA_ROOT"]) / "kv21" / "dataDistances" / "analysis"
    out.mkdir(parents=True, exist_ok=True)

    parts = []
    for protocol in ("vanilla", "masked"):
        frame = pd.read_csv(source("L403A", protocol), usecols=["pdb_file", *L403_COLUMNS])
        values = frame[L403_COLUMNS].apply(pd.to_numeric, errors="coerce")
        table = pd.DataFrame({"condition": protocol, "pdb_file": frame.pdb_file.map(lambda p: Path(str(p)).name)})
        table[L403_COLUMNS] = values.to_numpy()
        for c, column in zip("ABCD", L403_COLUMNS):
            table[f"chain_{c}"] = values[column].to_numpy()
        table["maximum_distance_A"] = values.max(axis=1).to_numpy()
        table["minimum_distance_A"] = values.min(axis=1).to_numpy()
        table["mean_distance_A"] = values.mean(axis=1).to_numpy()
        table["chain_with_maximum"] = values.idxmax(axis=1).str.extract(r"CA_CA_([A-D])_")[0].to_numpy()
        table["number_shifted_subunits"] = values.gt(SHIFT_CUTOFF_A).sum(axis=1).to_numpy()
        parts.append(table)
    pd.concat(parts, ignore_index=True).to_csv(out / "L403A_E423_N179_all_structure_distances.csv", index=False)

    records = []
    for protocol in ("vanilla", "masked"):
        path = source("F412L", protocol)
        header = pd.read_csv(path, nrows=0).columns
        columns = _nexus_distance_columns(header)
        frame = pd.read_csv(path, usecols=["pdb_file", *columns], low_memory=False)
        for partner, label in CONTACT_LABELS.items():
            chosen = [c for c in columns if re.search(fr"(?:^|[A-Z]){partner}[A-D](?:-|$)", c)]
            nearest = frame[chosen].apply(pd.to_numeric, errors="coerce").min(axis=1)
            records.append(pd.DataFrame({
                "Protocol": protocol.capitalize(),
                "pdb_file": frame.pdb_file.map(lambda p: Path(str(p)).name),
                "model_path": frame.pdb_file,
                "Contact": label,
                "Shortest heavy-atom distance (Å)": nearest,
                "number_of_available_pairs": len(chosen),
            }))
    contacts = pd.concat(records, ignore_index=True)
    contacts.to_csv(out / "F412L_nexus_shortest_contacts_long.csv", index=False)
    print(contacts.groupby(["Protocol", "Contact"])["Shortest heavy-atom distance (Å)"].describe().round(2))


if __name__ == "__main__":
    main()
