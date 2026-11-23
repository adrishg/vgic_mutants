#!/usr/bin/env python3
"""Which atom set does a ``shortest_*`` column use? Direct check on single models.

For each model and residue pair, print the minimum distance over (a) every
non-hydrogen atom and (b) side-chain non-hydrogen atoms only (no N, CA, C, O,
OXT), next to the C-alpha distance. Compare by eye with the table values.
"""
import sys
from itertools import product
from pathlib import Path

import numpy as np

BACKBONE = {"N", "CA", "C", "O", "OXT"}
PAIRS = [("A", 405, "B", 414), ("A", 414, "A", 315), ("A", 405, "B", 418),
         ("A", 414, "D", 331), ("A", 414, "A", 318)]


def atoms(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        name, chain, number = line[12:16].strip(), line[21], int(line[22:26])
        element = (line[76:78].strip() or name[0]).upper()
        if element in {"H", "D"}:
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        out.setdefault((chain, number), []).append((name, line[17:20], xyz))
    return out


for path in sys.argv[1:]:
    table = atoms(path)
    print(Path(path).name)
    for c1, n1, c2, n2 in PAIRS:
        a, b = table[(c1, n1)], table[(c2, n2)]
        every = min(np.linalg.norm(x[2] - y[2]) for x, y in product(a, b))
        sa = [x for x in a if x[0] not in BACKBONE]
        sb = [y for y in b if y[0] not in BACKBONE]
        side = min(np.linalg.norm(x[2] - y[2]) for x, y in product(sa, sb)) if sa and sb else float("nan")
        ca = np.linalg.norm([x[2] for x in a if x[0] == "CA"][0] - [y[2] for y in b if y[0] == "CA"][0])
        print(f"  {a[0][1]}{n1}{c1}-{b[0][1]}{n2}{c2}: all_heavy={every:.2f} side_chain={side:.2f} CA={ca:.2f}")
