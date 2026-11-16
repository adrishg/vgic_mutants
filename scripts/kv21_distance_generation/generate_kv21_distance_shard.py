#!/usr/bin/env python3
"""Generate a shard of full-tetramer Kv2.1 distance tables.

The existing Kv2.1 CSV headers are the schema.  Both shortest-heavy-atom and
C-alpha columns are parsed directly from those headers, including chain labels,
so the generated table retains the historical column order and all A-D and
inter-subunit measurements.  Shards are independent and can be concatenated
after successful completion.
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

import numpy as np


SHORT = re.compile(r"^shortest_([A-Z]{3})(\d+)([A-Z])-([A-Z]{3})(\d+)([A-Z])$")
CA = re.compile(
    r"^CA_CA_([A-Z])_([A-Z]{3})(\d+)_CA-([A-Z])_([A-Z]{3})(\d+)_CA$"
)
MEASUREMENTS: list[tuple[str, tuple[str, int, str, str, int, str]]] = []


def parse_pdb(path: Path):
    heavy: dict[tuple[str, str, int], list[np.ndarray]] = {}
    ca: dict[tuple[str, str, int], np.ndarray] = {}
    with path.open(errors="ignore") as handle:
        for line in handle:
            if line[:6].strip() not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            if line[16] not in {" ", "A"}:
                continue
            try:
                chain = line[21].strip()
                residue = line[17:20].strip()
                number = int(line[22:26])
                atom = line[12:16].strip()
                element = line[76:78].strip() or atom[:1]
                xyz = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                )
            except (IndexError, ValueError):
                continue
            key = (chain, residue, number)
            if element.upper() != "H":
                heavy.setdefault(key, []).append(xyz)
            if atom == "CA":
                ca[key] = xyz
    return {key: np.vstack(value) for key, value in heavy.items()}, ca


def minimum_distance(first: np.ndarray, second: np.ndarray) -> float:
    delta = first[:, None, :] - second[None, :, :]
    return float(np.sqrt(np.einsum("ijk,ijk->ij", delta, delta).min()))


def calculate_row(path: Path) -> tuple[list[str], int]:
    heavy, ca = parse_pdb(path)
    row = [path.name]
    missing = 0
    for kind, (res_a, num_a, chain_a, res_b, num_b, chain_b) in MEASUREMENTS:
        coordinates = ca if kind == "ca" else heavy
        first = coordinates.get((chain_a, res_a, num_a))
        second = coordinates.get((chain_b, res_b, num_b))
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


def configure(headers: list[str]) -> None:
    global MEASUREMENTS
    MEASUREMENTS = []
    for label in headers:
        if label == "pdb_file":
            continue
        match = SHORT.fullmatch(label)
        if match:
            MEASUREMENTS.append(("shortest", (match[1], int(match[2]), match[3], match[4], int(match[5]), match[6])))
            continue
        match = CA.fullmatch(label)
        if match:
            MEASUREMENTS.append(("ca", (match[2], int(match[3]), match[1], match[5], int(match[6]), match[4])))
            continue
        raise ValueError(f"Unrecognized Kv2.1 distance column: {label}")


def write_rows(paths: list[Path], headers: list[str], output: Path, workers: int) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".part", dir=output.parent)
    os.close(fd)
    rows = missing = 0
    try:
        with open(temporary, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(headers)
            with ProcessPoolExecutor(max_workers=workers) as pool:
                for row, row_missing in pool.map(calculate_row, paths, chunksize=8):
                    writer.writerow(row)
                    rows += 1
                    missing += row_missing
        if rows != len(paths):
            raise RuntimeError(f"Wrote {rows} rows for {len(paths)} inputs")
        if missing:
            raise RuntimeError(f"Missing measurements: {missing} cell values in {output}")
        os.replace(temporary, output)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--main-template", type=Path, required=True)
    parser.add_argument("--extra-template", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise SystemExit("Invalid shard index")
    paths = sorted(args.input_dir.rglob("*.pdb"))
    if not paths:
        raise SystemExit(f"No PDB files below {args.input_dir}")
    paths = paths[args.shard_index :: args.shard_count]
    main_headers = read_header(args.main_template)
    extra_headers = read_header(args.extra_template)
    for headers in (main_headers, extra_headers):
        if headers[0] != "pdb_file":
            raise SystemExit("Distance templates must begin with pdb_file")
    configure(main_headers[1:])
    main_output = args.output_dir / f"distances_main_part_{args.shard_index:02d}.csv"
    extra_output = args.output_dir / f"distances_extra_part_{args.shard_index:02d}.csv"
    if main_output.exists() or extra_output.exists():
        raise SystemExit("Refusing to overwrite existing shard")
    n_main = write_rows(paths, main_headers, main_output, args.workers)
    configure(extra_headers[1:])
    n_extra = write_rows(paths, extra_headers, extra_output, args.workers)
    if n_main != n_extra:
        raise SystemExit("Main/extra row counts differ")
    print(f"shard={args.shard_index} inputs={len(paths)} rows={n_main}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
