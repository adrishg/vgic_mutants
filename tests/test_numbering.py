"""Model residue numbers to UniProt: the maps the dataset configs are written with."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "filter_stats"))

import channels  # noqa: E402
import numbering as nb  # noqa: E402


def test_every_anchor_lands_where_the_model_files_say_it_does():
    """The anchors were read off the model PDBs; this is what stops the tables drifting from them."""
    for channel, anchors in nb.ANCHORS.items():
        for model, _residue, uniprot in anchors:
            assert nb.to_uniprot(channel, model) == uniprot, f"{channel} model {model}"


def test_the_mutation_names_are_labels_and_not_residue_numbers():
    """L403A is model 405 / UniProt 407, and F412L is model 414 / UniProt 416.

    Writing either name into a config as a residue number points at the wrong residue and still
    runs, which is the whole reason this module exists."""
    assert nb.to_uniprot("Kv2.1", 405) == 407
    assert nb.to_uniprot("Kv2.1", 414) == 416
    assert nb.to_uniprot("Kv2.1", 403) != 403


def test_a_span_across_a_truncation_becomes_several_ranges():
    """Repeat II's S5 runs into the gap the construct's truncation leaves, so it cannot be one pair.

    Returning a single pair here would silently take in the residues the construct does not have."""
    s5 = channels.channel("Nav1.5").segments[1]["S5"]
    ranges = nb.to_uniprot_ranges("Nav1.5", s5)
    assert len(ranges) >= 1
    assert ranges[0][0] == nb.to_uniprot("Nav1.5", s5[0])
    covered = sum(b - a + 1 for a, b in ranges)
    assert covered <= s5[1] - s5[0] + 1, "a conversion may lose residues, never invent them"


def test_unmapped_residues_have_no_number():
    assert nb.to_uniprot("Kv2.1", 1) is None               # the construct's cloning tag
    assert nb.to_uniprot("Kv2.1", 601) is None             # past the construct
    assert nb.to_uniprot("Cav1.2", 1) == 1                 # an exact prefix, nothing lost


def test_nav15_maps_every_model_residue():
    """The rat construct is loop-truncated but contiguous: the truncations are gaps in rat
    numbering, not gaps in the model's. Against human Q14524 it looked gapped, which is what a
    wrong reference sequence does to an alignment."""
    unmapped = [n for n in range(1, 1573) if nb.to_uniprot("Nav1.5", n) is None]
    assert unmapped == []
    assert nb.to_uniprot("Nav1.5", 461) == 461             # the last residue before the first cut
    assert nb.to_uniprot("Nav1.5", 462) == 659             # and the first after it


@pytest.mark.parametrize("channel", sorted(nb.BLOCKS))
def test_blocks_are_ordered_and_do_not_overlap(channel):
    blocks = nb.BLOCKS[channel]
    for (a0, a1, _), (b0, _b1, _o) in zip(blocks, blocks[1:]):
        assert a0 <= a1 < b0, f"{channel}: {a0}-{a1} then {b0}"


@pytest.mark.parametrize("channel", sorted(nb.BLOCKS))
def test_every_segment_the_chain_uses_has_a_uniprot_range(channel):
    """A segment that converts to nothing would become an empty region in the config."""
    for table in channels.channel(channel).segments:
        for name, span in table.items():
            assert nb.to_uniprot_ranges(channel, span), f"{channel} {name} {span}"


@pytest.mark.parametrize("channel", sorted(nb.BLOCKS))
def test_segments_run_in_order_without_sharing_a_residue(channel):
    """The screen rejects a table where one residue belongs to two segments.

    Deposited HELIX records share boundaries - Cav1.2 repeat IV annotates S3 as starting where S2
    ends, and several repeats run the filter into the pore helix - so channels.py trims them."""
    order = ("S1", "S2", "S3", "S4", "linker", "S5", "pore_helix", "filter", "S6_upper", "S6_lower")
    for index, table in enumerate(channels.channel(channel).segments, start=1):
        previous = None
        for name in order:
            first, last = table[name]
            assert first <= last, f"{channel} repeat {index} {name} is empty"
            if previous is not None:
                assert first > previous, f"{channel} repeat {index} {name} overlaps the one before"
            previous = last
