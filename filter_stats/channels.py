#!/usr/bin/env python3
"""What each channel needs for the quality-control chain: its filter ring, and its repeats.

Kv2.1 is a homotetramer, four chains of one sequence. Nav1.5 and Cav1.2 are pseudo-tetramers, one
chain with four non-identical repeats, so everything the chain does per subunit has to be said four
times over, in the cohort's own numbering.

Two things measured on 2026-11-21, both in 26-11/filterComparison, decide how much of the Kv2.1
chain transfers:

1. The selectivity-filter ring generalizes as a definition and is nearly inert as a filter. EEEE
   (E363/E706/E1135/E1464) in Cav1.2 and DEKA (D373/E704/K1103/A1395) in Nav1.5 are the same six
   inter-subunit Ca-Ca distances the G377 ring is, and the columns are already in the cohort tables.
   But the widest of the six never exceeded 10.8 A in any Cav1.2 cohort or 12.0 A in any Nav1.5
   cohort, and the G377 step rejected no trajectory at all in either 99-column Kv2.1 cohort. The
   ring is a tripwire worth keeping because it costs nothing, not the step that finds broken models.

2. The topology screen is what generalizes as a filter, with one caveat. `wiring` and
   `s6_on_own_s5` passed 4/4 on every structure tried - both Kv2.1 experimental structures, Cav1.2
   8WE6 and 8HLP, 15 Cav1.2 models and 19 Nav1.5 models. `box_corners` needs care: its 45-135 deg
   window assumes a square, which a homotetramer's lower S6 ring is (8SD3 gives 90/90/90/90) and a
   pseudo-tetramer's is not. With the lower-S6 window below each gate residue the Cav1.2 ring reads
   86/90/92/93 deg and passes, and Nav1.5 reads about 70/75/90/125 deg and mostly passes - 17 of 19
   models, the two failures sitting at 137-138 deg, at the edge rather than misfolded. So the
   window is a per-channel number, not a constant, and S6_LOWER_LENGTH matters: the angles are very
   sensitive to a lower-S6 range that runs past the modeled helix.

A third finding, measured on 2026-11-21 across fourteen experimental structures, replaces the
Kv2.1 interface check for the other two channels. The published 27 A cutoff is the widest of the
K429/E425/K422-to-N181/V184 contacts, which is the end of S6 against the residues just before S1 of
the same subunit. That pairing is a contact only in Kv2.1 (8.6 A in 8SD3); in Cav1.2 and Nav1.5 the
same two elements sit 25 to 50 A apart, so a cutoff on them would be a threshold on a non-contact.
What does hold everywhere is the lower S6 against its own S4-S5 linker, the joint the gate actually
hangs from: 5.1-6.9 A over six Kv2.1 structures, 4.8-7.6 A over two Cav1.2, 5.6-9.2 A over six
Nav1.5, and in every one of them the subunit's own linker is nearer than any other subunit's. That
is the measure `s6_linker_max_A` thresholds, and 12 A sits well clear of the 9.2 A widest contact
any experimental structure shows.

A note on Nav1.5 numbering. The cohorts use a loop-truncated construct, and its offset against the
experimental structures is piecewise - DI -1, DII +196 through the voltage sensor but +194 from the
filter on, DIII and DIV +316. So a Nav1.5 segment table has to be built in the models' own numbering,
which is what SEGMENTS holds; an experimental structure plus one constant will not do.
"""
from __future__ import annotations

from dataclasses import dataclass

S6_LOWER_LENGTH = 12              # residues below the gate: a fixed window, not the rest of the helix


VSD_SEGMENTS = ("S1", "S2", "S3", "S4", "linker")
PORE_SEGMENTS = ("S5", "pore_helix", "filter", "S6_upper", "S6_lower")


@dataclass(frozen=True)
class Channel:
    """One channel's chain layout, filter ring and repeat segments."""
    name: str
    chains: tuple                 # the chain ids a model holds
    ring_position: int = 0        # homotetramer: one residue number, the same in all four chains
    ring_residues: tuple = ()     # pseudo-tetramer: the four filter residues, in model numbering
    ring_max_A: float = 20.0
    segments: tuple = ()          # one table per repeat; a homotetramer gives one table for all four
    interface_qc: bool = False    # the published Kv2.1 pore-VSD cutoff applies to this channel
    s6_linker_max_A: float = 12.0  # the gate-to-S4-S5-linker contact may open this far and no further
    experimental_offsets: tuple = ()  # per repeat, (voltage sensor, pore): add to a model number

    def experimental_segments(self) -> tuple:
        """The segment tables written in the experimental structures' numbering instead of the models'.

        Calibrating a cutoff means measuring the same thing on a solved structure that the check
        measures on a model, and for Nav1.5 that is not one subtraction: the construct's offset
        against 7DTC is piecewise, and within repeat II it differs either side of the filter. Taking
        the offsets from here is what keeps the calibration and the check the same measurement."""
        if not self.experimental_offsets:
            return self.segments
        tables = self.segments if len(self.segments) > 1 else self.segments * len(self.experimental_offsets)
        out = []
        for table, (vsd, pore) in zip(tables, self.experimental_offsets):
            shifted = {}
            for key, (first, last) in table.items():
                shift = vsd if key in VSD_SEGMENTS else pore
                shifted[key] = (first + shift, last + shift)
            out.append(shifted)
        return tuple(out)

    def subunits(self) -> tuple:
        """(chain, segment table) once per subunit, whichever way the channel is built.

        A homotetramer has four chains sharing one table; a pseudo-tetramer has one chain with four.
        Everything measured per subunit goes through here so it does not have to care which."""
        if len(self.segments) == 1:
            return tuple((c, self.segments[0]) for c in self.chains)
        return tuple((self.chains[0], s) for s in self.segments)


def _repeat(s1, s2, s3, s4, linker, s5, pore_helix, filter_, s6_start, gate) -> dict:
    """One repeat, with the lower S6 taken as a fixed window below the gate residue.

    Deposited HELIX records sometimes share a boundary residue between neighbors - the filter is
    annotated as the end of the pore helix, S3 as starting where S2 stops - so each segment is
    trimmed to start after the one before it. A residue in two segments at once is not a
    segmentation, and the screen rejects the table outright."""
    table = {"S1": s1, "S2": s2, "S3": s3, "S4": s4, "linker": linker, "S5": s5,
             "pore_helix": pore_helix, "filter": filter_,
             "S6_upper": (s6_start, gate), "S6_lower": (gate + 1, gate + S6_LOWER_LENGTH)}
    out, previous_end = {}, None
    for name, (first, last) in table.items():
        if previous_end is not None and first <= previous_end:
            first = previous_end + 1
        if first > last:
            raise ValueError(f"segment {name} {table[name]} is swallowed by the one before it")
        out[name] = (first, last)
        previous_end = last
    return out


# Kv2.1: 8SD3's own HELIX records, plus 2 for the model numbering. S6 breaks at L403/P404, which is
# both the mutation site and the published kink, so the lower S6 is the published range rather than a
# window below a gate residue - this is the segmentation the filter-order comparison was run with.
KV21_SEGMENTS = ({"S1": (183, 207), "S2": (226, 248), "S3": (253, 279), "S4": (291, 311),
                  "linker": (312, 325), "S5": (326, 352), "pore_helix": (361, 374),
                  "filter": (375, 380), "S6_upper": (385, 405), "S6_lower": (406, 430)},)

# Cav1.2: 8WE6's HELIX records, which are already in the cohorts' numbering. Each element is the
# helix holding its anchor - the S2 charge-transfer phenylalanine, the S4 gating charges, the EEEE
# glutamate, the S6 gate residue.
CAV12_SEGMENTS = (
    _repeat((125, 143), (153, 182), (193, 217), (232, 252), (253, 267), (268, 290),
            (350, 361), (362, 367), 381, 401),
    _repeat((513, 542), (551, 577), (586, 608), (613, 633), (634, 646), (647, 674),
            (691, 705), (705, 710), 730, 749),
    _repeat((902, 919), (928, 950), (955, 975), (1011, 1033), (1036, 1049), (1051, 1072),
            (1120, 1133), (1134, 1139), 1162, 1182),
    _repeat((1227, 1256), (1264, 1301), (1301, 1322), (1372, 1392), (1393, 1408), (1410, 1431),
            (1449, 1462), (1463, 1468), 1499, 1516),
)

# Nav1.5: 7DTC's HELIX records mapped into the cohorts' construct numbering by the piecewise offsets
# in the module docstring, each one checked against the anchor residue's identity.
NAV15_SEGMENTS = (
    _repeat((131, 150), (158, 178), (192, 208), (219, 232), (236, 249), (251, 274),
            (358, 371), (372, 377), 390, 415),
    _repeat((520, 540), (548, 573), (585, 602), (612, 624), (628, 641), (647, 676),
            (690, 703), (703, 708), 717, 742),
    _repeat((888, 905), (918, 946), (955, 977), (983, 1002), (1003, 1014), (1015, 1041),
            (1088, 1102), (1102, 1107), 1130, 1154),
    _repeat((1213, 1231), (1237, 1266), (1275, 1293), (1305, 1322), (1324, 1336), (1337, 1363),
            (1380, 1393), (1394, 1399), 1428, 1455),
)

CHANNELS = {
    "Kv2.1": Channel(name="Kv2.1", chains=("A", "B", "C", "D"), ring_position=377,
                     ring_max_A=20.0, segments=KV21_SEGMENTS, interface_qc=True,
                     s6_linker_max_A=12.0, experimental_offsets=((-2, -2),) * 4),
    "Cav1.2": Channel(name="Cav1.2", chains=("A",),
                      ring_residues=("GLU363", "GLU706", "GLU1135", "GLU1464"),
                      ring_max_A=20.0, segments=CAV12_SEGMENTS, interface_qc=False,
                      s6_linker_max_A=12.0, experimental_offsets=((0, 0),) * 4),
    "Nav1.5": Channel(name="Nav1.5", chains=("A",),
                      ring_residues=("ASP373", "GLU704", "LYS1103", "ALA1395"),
                      ring_max_A=20.0, segments=NAV15_SEGMENTS, interface_qc=False,
                      s6_linker_max_A=12.0,
                      experimental_offsets=((1, 1), (196, 194), (316, 316), (316, 316))),
}

ALIASES = {"kv21": "Kv2.1", "kv2.1": "Kv2.1", "cav12": "Cav1.2", "cav1.2": "Cav1.2",
           "nav15": "Nav1.5", "nav1.5": "Nav1.5"}


def channel(name: str) -> Channel:
    """One channel by any of the spellings the repository uses for it."""
    key = ALIASES.get(str(name).strip().lower(), str(name).strip())
    if key not in CHANNELS:
        raise KeyError(f"unknown channel {name!r}; known: {sorted(CHANNELS)}")
    return CHANNELS[key]
