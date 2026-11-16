#!/usr/bin/env python3
"""Generate Cav1.2 distance tables using an existing table schema.

This generator reuses the current distance-table headers so
the new masked ensembles have identical main and extra measurement columns.
Distances are calculated on chain A from all ``*.pdb`` files below the input
directory.  Shortest distances use non-hydrogen atoms; C-alpha distances use
the named residue C-alpha atoms.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import tempfile
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import Iterable

import numpy as np


SHORT_PATTERN = re.compile(r"^shortest_([A-Z]{3})(\d+)-([A-Z]{3})(\d+)$")
CA_PATTERN = re.compile(r"^CA_([A-Z]{3})(\d+)_CA-([A-Z]{3})(\d+)_CA$")

# Keep measurements in template order.  The source templates interleave
# shortest-atom and C-alpha columns; separating these lists would shift values
# under the wrong headers when a row is written.
_MEASUREMENTS: list[tuple[str, tuple[str, int, str, int]]] = []


def parse_pdb(path: Path, chain: str) -> tuple[dict[tuple[str, int], np.ndarray], dict[tuple[str, int], np.ndarray]]:
    heavy: dict[tuple[str, int], list[np.ndarray]] = {}
    ca: dict[tuple[str, int], np.ndarray] = {}
    with path.open(errors="ignore") as handle:
        for line in handle:
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            if line[21].strip() != chain or line[16] not in {" ", "A"}:
                continue
            try:
                residue = (line[17:20].strip(), int(line[22:26]))
                atom = line[12:16].strip()
                element = line[76:78].strip() or atom[:1]
                xyz = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                )
            except (IndexError, ValueError):
                continue
            if element.upper() != "H":
                heavy.setdefault(residue, []).append(xyz)
            if atom == "CA":
                ca[residue] = xyz
    return {key: np.vstack(value) for key, value in heavy.items()}, ca


def minimum_distance(first: np.ndarray, second: np.ndarray) -> float:
    delta = first[:, None, :] - second[None, :, :]
    return float(np.sqrt(np.einsum("ijk,ijk->ij", delta, delta).min()))


def calculate_row(path: Path, chain: str) -> tuple[list[str], int]:
    heavy, ca = parse_pdb(path, chain)
    row = [path.name]
    missing = 0
    for kind, (residue_a, number_a, residue_b, number_b) in _MEASUREMENTS:
        coordinates = ca if kind == "ca" else heavy
        first = coordinates.get((residue_a, number_a))
        second = coordinates.get((residue_b, number_b))
        if first is None or second is None:
            row.append("")
            missing += 1
        elif kind == "ca":
            row.append(f"{math.dist(first, second):.2f}".rstrip("0").rstrip("."))
        else:
            row.append(f"{minimum_distance(first, second):.1f}")
    return row, missing


def read_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return next(csv.reader(handle))


def configure_columns(headers: Iterable[str]) -> None:
    global _MEASUREMENTS
    _MEASUREMENTS = []
    for label in headers:
        if label == "pdb_file":
            continue
        match = SHORT_PATTERN.fullmatch(label)
        if match:
            _MEASUREMENTS.append(
                ("shortest", (match.group(1), int(match.group(2)), match.group(3), int(match.group(4))))
            )
            continue
        match = CA_PATTERN.fullmatch(label)
        if match:
            _MEASUREMENTS.append(
                ("ca", (match.group(1), int(match.group(2)), match.group(3), int(match.group(4))))
            )
            continue
        raise ValueError(f"Unrecognized distance column: {label}")


def write_table(
    paths: list[Path],
    headers: list[str],
    output: Path,
    chain: str,
    workers: int,
) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".part", dir=output.parent)
    os.close(fd)
    missing = 0
    rows = 0
    try:
        with open(temporary_name, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(headers)
            with ProcessPoolExecutor(max_workers=workers) as pool:
                calculate = partial(calculate_row, chain=chain)
                for row, row_missing in pool.map(
                    calculate, paths, chunksize=16
                ):
                    writer.writerow(row)
                    rows += 1
                    missing += row_missing
        if rows != len(paths):
            raise RuntimeError(f"Wrote {rows} rows for {len(paths)} input PDB files")
        if missing:
            raise RuntimeError(f"Found {missing} missing distance measurements in {output.name}")
        os.replace(temporary_name, output)
    except Exception:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--main-template", type=Path, required=True)
    parser.add_argument("--extra-template", type=Path, required=True)
    parser.add_argument("--main-output", type=Path, required=True)
    parser.add_argument("--extra-output", type=Path, required=True)
    parser.add_argument("--chain", default="A")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    paths = sorted(args.input_dir.rglob("*.pdb"))
    if not paths:
        raise SystemExit(f"No PDB files found below {args.input_dir}")
    if args.main_output.exists() or args.extra_output.exists():
        raise SystemExit("Refusing to overwrite an existing output")

    main_headers = read_header(args.main_template)
    extra_headers = read_header(args.extra_template)
    for headers in (main_headers, extra_headers):
        if headers[0] != "pdb_file":
            raise SystemExit("Distance templates must begin with pdb_file")
        configure_columns(headers[1:])

    # The main and extra tables use different schemas, so write them in two
    # passes while retaining the exact source-template column order.
    configure_columns(main_headers[1:])
    main_rows = write_table(paths, main_headers, args.main_output, args.chain, args.workers)
    configure_columns(extra_headers[1:])
    extra_rows = write_table(paths, extra_headers, args.extra_output, args.chain, args.workers)
    if main_rows != extra_rows:
        raise SystemExit("Main and extra tables have different row counts")
    print(f"Wrote {main_rows} rows to {args.main_output} and {args.extra_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
