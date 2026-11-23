#!/usr/bin/env python3
"""Write the seven dataset configs Brandon's steps 1 and 2 read, in UniProt numbering.

One dataset per channel and protocol, with that channel's conditions inside it, because the
comparison of interest is mutant against wild type within one protocol. Vanilla and masked are
separate datasets rather than a second factor: they are different predictions of the same thing,
not two levels of one treatment.

Every residue range is converted from the models' numbering with filter_stats/numbering.py. Nothing
here is written by hand in UniProt numbers; see that module for why hand-writing them is the one
mistake that runs to completion and reports the wrong residues.

  python scripts/filtering/configs/make_configs.py --out datasets
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "filter_stats"))

import channels as ch                                                  # noqa: E402
import numbering as nb                                                 # noqa: E402

QUOBYTE = "/quobyte/yarovoygrp/ahgz/vgic_mutants"
ARCHIVES = f"{QUOBYTE}/analysis/archives"
RESULTS = f"{QUOBYTE}/analysis/stats"

FIT_SEGMENTS = ("S5", "pore_helix", "filter", "S6_upper")
MEASURE_SEGMENTS = ("S1", "S2", "S3", "S4", "linker", "S5", "pore_helix", "filter",
                    "S6_upper", "S6_lower")
# Region names take letters, digits and '-' only: they become column names like spread_<region>.
REGIONS = {"vsd": ("S1", "S2", "S3", "S4"), "s4s5-linker": ("linker",),
           "pore": ("S5", "pore_helix", "filter", "S6-upper".replace("-", "_")),
           "filter": ("filter",), "gate": ("S6_lower",)}

REFERENCES = {"Kv2.1": ("kv21/experimental/8SD3.pdb", ["A", "B", "C", "D"], []),
              "Nav1.5": ("nav15/experimental/7DTC.pdb", ["A"], []),
              "Cav1.2": ("cav12/experimental/8WE6.pdb", ["A"], ["C", "D"])}

# dataset -> (channel, protocol folder, control, [(condition, markers)])
DATASETS = {
    "kv21_vanilla": ("Kv2.1", "vanilla", "WT",
                     [("WT", {407: "L", 416: "F"}), ("L403A", {407: "A", 416: "F"}),
                      ("F412L", {407: "L", 416: "L"})]),
    "kv21_masked99": ("Kv2.1", "masked_common99_50pct", "WT",
                      [("WT", {407: "L", 416: "F"}), ("L403A", {407: "A", 416: "F"}),
                       ("F412L", {407: "L", 416: "L"})]),
    "nav15_vanilla": ("Nav1.5", "vanilla", "WT",
                      [("WT", {1485: "I", 1486: "F", 1487: "M"}),
                       ("QQQ", {1485: "Q", 1486: "Q", 1487: "Q"})]),
    "nav15_masked_v1": ("Nav1.5", "masked_v1_50pct", "WT",
                        [("WT", {1485: "I", 1486: "F", 1487: "M"}),
                         ("QQQ", {1485: "Q", 1486: "Q", 1487: "Q"})]),
    "nav15_masked_v2": ("Nav1.5", "masked_v2_50pct", "WT",
                        [("WT", {1485: "I", 1486: "F", 1487: "M"}),
                         ("QQQ", {1485: "Q", 1486: "Q", 1487: "Q"})]),
    "cav12_vanilla": ("Cav1.2", "vanilla", "WT",
                      [("WT", {402: "G", 406: "G"}), ("G402S", {402: "S", 406: "G"}),
                       ("G406R", {402: "G", 406: "R"})]),
    "cav12_masked_v3": ("Cav1.2", "masked_v3_50pct", "WT",
                        [("WT", {402: "G", 406: "G"}), ("G402S", {402: "S", 406: "G"}),
                         ("G406R", {402: "G", 406: "R"})]),
}
ARCHIVE_CONDITION = {"WT": "wt", "L403A": "l403a", "F412L": "f412l", "QQQ": "qqq",
                     "G402S": "g402s", "G406R": "g406r"}


def _ranges(channel: str, names) -> list:
    """The UniProt ranges those segments cover, over every repeat, merged and in order."""
    out = []
    for table in ch.channel(channel).segments:
        for name in names:
            if name in table:
                out += nb.to_uniprot_ranges(channel, table[name])
    out.sort()
    merged = []
    for pair in out:
        if merged and pair[0] <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], pair[1])
        else:
            merged.append(list(pair))
    return merged


def _without(ranges, positions) -> list:
    """The ranges with those residues taken out, splitting a range where one sits inside it.

    The models are superposed on the fit set, so a residue the treatment changes cannot be in it
    (D-030): fitting on it would let the mutation move every other residue. Kv2.1 is where this
    bites - the L403A site is UniProt 407, the last residue of S6_upper, which is otherwise exactly
    the kind of stable helix a fit set wants."""
    drop = set(positions)
    out = []
    for first, last in ranges:
        run = []
        for n in range(first, last + 1):
            if n in drop:
                if run:
                    out.append([run[0], run[-1]])
                run = []
            else:
                run.append(n)
        if run:
            out.append([run[0], run[-1]])
    return out


def _pairs(ranges) -> str:
    return "[" + ", ".join(f"[{a}, {b}]" for a, b in ranges) + "]"


def _segment_tables(channel: str) -> str:
    """[segments] for a homotetramer, or four [[segments]] for a pseudo-tetramer."""
    order = ("S1", "S2", "S3", "S4", "linker", "S5", "pore_helix", "filter", "S6_upper", "S6_lower")
    tables = ch.channel(channel).segments
    header = "[segments]" if len(tables) == 1 else "[[segments]]"
    out = []
    for table in tables:
        lines = [header]
        for name in order:
            # A segment has to be one unbroken span, so where a truncation splits it the longest
            # mapped run is used. Spanning the gap instead would claim residues the construct does
            # not have, and the screen checks every segment against the measure set.
            spans = nb.to_uniprot_ranges(channel, table[name])
            first, last = max(spans, key=lambda r: r[1] - r[0])
            lines.append(f"{name} = [{first}, {last}]")
        out.append("\n".join(lines))
    return "\n\n".join(out)


def write(dataset: str, out_dir: Path) -> Path:
    channel, protocol, control, conditions = DATASETS[dataset]
    folder = out_dir / dataset
    folder.mkdir(parents=True, exist_ok=True)
    ref_path, ref_chains, ignore = REFERENCES[channel]
    treated = sorted({p for _n, markers in conditions for p in markers})
    chains = list(ch.channel(channel).chains)

    body = [f'# {dataset}: {channel}, {protocol}. Written by make_configs.py - edit that, not this.',
            f'# Ranges are UniProt ({nb.ACCESSIONS[channel]}), converted from the models\' own',
            f'# numbering by filter_stats/numbering.py. The mutation names are labels: L403A is',
            f'# UniProt 407, F412L is 416, QQQ is the IFM motif at 1485-1487.',
            '',
            f'name = "{dataset}"',
            f'uniprot_fasta = "../uniprot/{nb.ACCESSIONS[channel]}.fasta"',
            f'output_dir = "{RESULTS}/{dataset}"',
            f'model_chains = {chains!r}'.replace("'", '"'),
            '',
            '[reference]',
            f'path = "../../{ref_path}"',
            f'chains = {ref_chains!r}'.replace("'", '"')]
    if ignore:
        body.append(f'ignore_chains = {ignore!r}'.replace("'", '"'))
    body += ['',
             '[ranges]',
             f'fit = {_pairs(_without(_ranges(channel, FIT_SEGMENTS), treated))}',
             f'measure = {_pairs(_ranges(channel, MEASURE_SEGMENTS))}',
             f'treated = {_pairs([[p, p] for p in treated])}',
             '',
             '[regions]']
    for region, names in REGIONS.items():
        body.append(f'{region} = {_pairs(_ranges(channel, names))}')
    body += ['',
             '[design]',
             'factors = ["genotype"]',
             f'control = "{control}"']
    for name, markers in conditions:
        archive = f"{ARCHIVES}/{dataset}_{ARCHIVE_CONDITION[name]}.tar.gz"
        marker_text = ", ".join(f'"{p}" = "{letter}"' for p, letter in sorted(markers.items()))
        body += ['',
                 '[[conditions]]',
                 f'name = "{name}"',
                 f'archive = "{archive}"',
                 f'levels = {{ "genotype" = "{name}" }}',
                 f'markers = {{ {marker_text} }}']
    body += ['',
             '[screen]',
             '# The checks that need no reference. own_place and filter_quarter are left off here:',
             '# they compare against the reference structure, and the chain we settled on judges',
             '# topology reference-free. first_recycle = 1 because recycle 0 is not analysed.',
             'domain_swapped = true',
             'wiring = true',
             'box_corners = true',
             's6_on_own_s5 = true',
             'own_place = false',
             'filter_quarter = false',
             'first_recycle = 1',
             '',
             _segment_tables(channel),
             '']
    path = folder / f"{dataset}.toml"
    path.write_text("\n".join(body))
    return path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, default=REPO / "datasets")
    args = parser.parse_args(argv)
    for dataset in DATASETS:
        print(f"wrote {write(dataset, args.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
