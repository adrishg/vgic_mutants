#!/usr/bin/env python3
"""Model residue numbers to UniProt, for the three constructs these cohorts were predicted from.

Brandon's pipeline is written entirely in UniProt numbering: the fit, measure and treated ranges,
the regions, the mutation markers. Our cohorts are not. Every channel here uses a different
construct, and two of the three are numbered unlike the reference sequence, so a config written in
the models' own numbers would point at the wrong residues and still run.

The maps below were derived by aligning each construct's FASTA to its UniProt sequence and then
checked against the model files themselves, which is the only check that settles it:

  Kv2.1    rat P15387, one block, UniProt = model + 2. The construct is rat, not human: rat covers
           598 of 600 residues in a single block where human Q14721 breaks into four. The filter
           glycine the distance tables call G377 is UniProt G379, and the TVGYG motif sits at
           UniProt 375-379.
  Nav1.5   human Q14524, piecewise. The construct is loop-truncated, so the offset is +0 at the very
           N-terminus, -1 through repeat I, +196 through repeat II's voltage sensor, +194 from
           repeat II's pore onward, and +316 through repeats III and IV. DEKA lands on the canonical
           D372/E898/K1419/A1711 and the repeat II gating charges on R808/R811/R814, which is what
           says the map is right.
  Cav1.2   human Q13936, one block, UniProt = model. The "short" construct is an exact prefix.

A third numbering exists and is not used here: the mutation names. L403A is model residue 405 and
UniProt 407; F412L is model 414 and UniProt 416. Those names follow the classic rat drk1 numbering,
four lower than current UniProt. Nav1.5's QQQ is the IFM inactivation motif, model 1169-1171 and
UniProt I1485/F1486/M1487. Cav1.2's G402S and G406R are the only ones that mean what they say in
every numbering. Read a mutation name as a label, never as a residue number.

Where a construct residue falls in a gap between blocks it has no UniProt number at all, so a range
converts to a list of ranges rather than one pair. That is why the regions in the configs are
written as lists.
"""
from __future__ import annotations

ACCESSIONS = {"Kv2.1": "P15387", "Nav1.5": "Q14524", "Cav1.2": "Q13936"}

# (first model residue, last model residue, what to add to get the UniProt number)
BLOCKS = {
    "Kv2.1": ((3, 600, +2),),
    "Nav1.5": ((5, 35, 0), (49, 73, -1), (75, 120, -1), (122, 134, -1), (136, 202, -1),
               (204, 286, -1), (288, 302, -1), (323, 461, -1), (462, 483, +196), (498, 512, +196),
               (524, 546, +196), (548, 605, +196), (607, 673, +196), (677, 759, +194),
               (761, 789, +194), (849, 862, +194), (872, 883, +316), (885, 1066, +316),
               (1074, 1125, +316), (1137, 1240, +316), (1242, 1262, +316), (1264, 1414, +316),
               (1426, 1495, +316), (1498, 1572, +316)),
    "Cav1.2": ((1, 1685, 0),),
}

# Checked against the model PDBs: model number, the residue there, and the UniProt number it is.
ANCHORS = {
    "Kv2.1": ((377, "G", 379), (405, "L", 407), (414, "F", 416)),
    "Nav1.5": ((373, "D", 372), (704, "E", 898), (1103, "K", 1419), (1395, "A", 1711),
               (1169, "I", 1485), (1170, "F", 1486), (1171, "M", 1487), (612, "R", 808)),
    "Cav1.2": ((363, "E", 363), (706, "E", 706), (1135, "E", 1135), (1464, "E", 1464),
               (402, "G", 402), (406, "G", 406)),
}


def to_uniprot(channel: str, number: int):
    """The UniProt number of one model residue, or None where the construct has no counterpart."""
    for first, last, offset in BLOCKS[channel]:
        if first <= number <= last:
            return number + offset
    return None


def to_uniprot_ranges(channel: str, span) -> list:
    """A model range as the UniProt ranges it covers, in order, skipping what does not map.

    One range can become several: a truncated construct runs two UniProt stretches together, and
    writing the span as a single pair would quietly take in the residues between them."""
    first, last = int(span[0]), int(span[1])
    out = []
    for low, high, offset in BLOCKS[channel]:
        a, b = max(first, low), min(last, high)
        if a <= b:
            out.append([a + offset, b + offset])
    out.sort()
    merged = []
    for pair in out:
        if merged and pair[0] == merged[-1][1] + 1:
            merged[-1][1] = pair[1]
        else:
            merged.append(list(pair))
    return merged
