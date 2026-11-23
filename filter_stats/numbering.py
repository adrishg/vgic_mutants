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
  Nav1.5   rat P15389, three blocks: +0 through model 461, +197 from 462 to 871, +318 from 872 to
           the end. The construct is rat, not human, which is what the model files say and what the
           prediction job names say too: they carry 6uz0, the rat cardiac channel, and 6UZ3's own
           SEQADV records the same loop deletions from residue 462 that our block boundary falls on.
           Against human Q14524 the same sequence breaks into twenty-four blocks and still leaves 53
           residues mismatched; against rat it is three blocks covering all 1572 with no mismatch at
           all. DEKA lands on D373/E901/K1421/A1713 and the IFM motif on 1487-1489, where 7FBS, the
           rat QQQ structure, records its three engineered glutamines.
  Cav1.2   human Q13936, one block, UniProt = model. The "short" construct is an exact prefix.

A third numbering exists and is not used here: the mutation names. L403A is model residue 405 and
UniProt 407; F412L is model 414 and UniProt 416. Those names follow the classic rat drk1 numbering,
four lower than current UniProt. Nav1.5's QQQ is the IFM inactivation motif, model 1169-1171 and
UniProt I1487/F1488/M1489. Cav1.2's G402S and G406R are the only ones that mean what they say in
every numbering. Read a mutation name as a label, never as a residue number.

Where a construct residue falls in a gap between blocks it has no UniProt number at all, so a range
converts to a list of ranges rather than one pair. That is why the regions in the configs are
written as lists.
"""
from __future__ import annotations

ACCESSIONS = {"Kv2.1": "P15387", "Nav1.5": "P15389", "Cav1.2": "Q13936"}

# (first model residue, last model residue, what to add to get the UniProt number)
BLOCKS = {
    "Kv2.1": ((3, 600, +2),),
    "Nav1.5": ((1, 461, 0), (462, 871, +197), (872, 1572, +318)),
    "Cav1.2": ((1, 1685, 0),),
}

# Checked against the model PDBs: model number, the residue there, and the UniProt number it is.
ANCHORS = {
    "Kv2.1": ((377, "G", 379), (405, "L", 407), (414, "F", 416)),
    "Nav1.5": ((373, "D", 373), (704, "E", 901), (1103, "K", 1421), (1395, "A", 1713),
               (1169, "I", 1487), (1170, "F", 1488), (1171, "M", 1489), (612, "R", 809)),
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
