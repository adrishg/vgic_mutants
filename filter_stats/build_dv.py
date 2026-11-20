#!/usr/bin/env python3
"""Step 1: structural-diversity dependent variables (DVs) from ColabFold predictions.

Run:  python -B pipeline/build_dv.py --config datasets/<dataset>/<dataset>.toml

Order of operations:
 1. Reads the dataset config (`config.py`): conditions and their archives, chains, the fit,
    measure and treated residue ranges in UniProt numbering, an optional reference structure.
 2. Numbers the reference structure, if any: each chain's sequence is aligned to the UniProt
    sequence, and each run of consecutively numbered residues takes one offset.
    Structure letters that differ from UniProt stop step 1 only inside the compared set
    (fit ∪ measure); every difference is reported.
 3. Reads each condition's ColabFold archive (one job; files `...model_<M>_seed_<S>.r<N>.pdb` at
    any depth) into `parse/`. This step is slow.
    Predictions missing a recycle are omitted; this removes early termination runs.
 4. Numbers the models by aligning their sequence to UniProt, and checks the
    markers (expected residues at mutation sites) and every other residue against UniProt.
 5. Superposes every model on the frame. The frame is either a reference structure, or without one, a model chosen
    over the fit set. Each model's chain labelling is chosen so that the models of a
    group match a common frame.
 6. Computes the DVs per group (condition × recycle) over the measure set: spread (each model's RMS
    distance to the group's average model); with a reference also spread_ref, deviation and offset
    over the measured residues the structure resolves. Checks two identities: bias-variance decomposition and.
 7. Writes `aligned/` (coords.npy, sites.csv, dv.csv, manifest.json) and `report/`
    (data_report.json, alignment.aln, step1.log) under the config's output folder.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import re
import subprocess
import sys
import tarfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

from config import Config, ConfigError, load_config, resolved

HERE = Path(__file__).resolve().parent                  # pipeline/
PROJECT = HERE.parent


class ParseError(Exception):
    """A coordinate file cannot be read as the config describes (chains, residue numbering)."""


class ArchiveError(Exception):
    """An archive is not what the config says it is (duplicates, jobs, sequences, markers, recycles)."""


class NumberingError(Exception):
    """The sequences cannot be numbered in UniProt numbering without changing a computed number."""


class IdentityError(Exception):
    """The DVs do not satisfy the identities they must satisfy; the outputs are wrong."""


STEP1_ERRORS = (ConfigError, ParseError, ArchiveError, NumberingError, IdentityError)
# step 2's files in <output>/regions/ (region_dv.py); step 1 removes them when it rebuilds aligned/ (D-061 item 1)
STEP2_FILES = ("manifest.json", "dv.csv", "region_report.json", "step2.log")


# =========================================================================== reading Cα coordinates
#
# A residue counts when gemmi's residue table calls its name an amino acid (ATOM or HETATM records
# alike) and its Cα atom is a carbon. Its letter is its parent's (MSE -> M). Ions, waters, sugars
# and lipids never count, even when an atom is named CA (D-039). A chain counts when it has at least
# one such residue; chains without any are skipped and reported (D-045).

@dataclass
class ChainCA:
    resnums: np.ndarray              # residue numbers, in file order
    seq: str                         # one letter per residue
    xyz: np.ndarray                  # (n, 3) Cα coordinates
    modified: list = field(default_factory=list)   # [(resnum, name)]: modified amino acids (D-039 item 4)


class Chains(dict):
    """{chain ID: ChainCA}, in the order the chains were asked for, plus what was skipped."""

    def __init__(self, items, ignored, absent_ignored):
        super().__init__(items)
        self.ignored = ignored                    # [(chain, reason)], reason "no amino acids" | "ignore_chains"
        self.absent_ignored = absent_ignored      # ignore_chains entries the file does not contain


@dataclass
class _Atom:
    chain: str
    resnum: int
    icode: str
    resname: str
    altloc: str
    occupancy: float
    element: str
    xyz: tuple


_RESIDUE_TABLE: dict = {}


def _residue_info(name: str):
    """(is an amino acid, one-letter code, is modified) from gemmi's residue table, cached."""
    if name not in _RESIDUE_TABLE:
        info = gemmi.find_tabulated_residue(name)
        if info is None or not info.is_amino_acid():
            _RESIDUE_TABLE[name] = (False, "", False)
        else:
            code = info.one_letter_code
            letter = code.upper() if code.strip() and code.isalpha() else "X"
            _RESIDUE_TABLE[name] = (True, letter, not info.is_standard())
    return _RESIDUE_TABLE[name]


_CA_LINE = re.compile(r"^(?:ATOM  |HETATM).{6}(?: CA |CA  |  CA)[^\n]*", re.M)
_CHAIN_OF_LINE = re.compile(r"^(?:ATOM  |HETATM).{15}(.)", re.M)


def _pdb_atoms(text: str):
    """Atoms named CA, and every chain ID that has coordinates, from PDB text (first model only)."""
    end = text.find("\nENDMDL")
    if end >= 0:
        text = text[:end]
    atoms = []
    for match in _CA_LINE.finditer(text):
        line = match.group(0).ljust(80)
        name = line[12:16]
        element = line[76:78].strip().upper() or (name[1] if name[0] == " " else name[:2].strip())
        occupancy = line[54:60].strip()
        atoms.append(_Atom(chain=line[21], resnum=int(line[22:26]), icode=line[26], resname=line[17:20].strip(),
                           altloc=line[16], occupancy=float(occupancy) if occupancy else 1.0,
                           element=element.upper(),
                           xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return atoms, set(_CHAIN_OF_LINE.findall(text))


def _gemmi_atoms(path: Path):
    """The same from any file gemmi reads (used for mmCIF): author chain IDs and numbers (D-039)."""
    st = gemmi.read_structure(str(path))
    atoms, present = [], set()
    for chain in st[0]:
        present.add(chain.name)
        for res in chain:
            for atom in res:
                if atom.name != "CA":
                    continue
                alt = atom.altloc if atom.altloc not in ("\0", "") else " "
                atoms.append(_Atom(chain=chain.name, resnum=res.seqid.num, icode=res.seqid.icode or " ",
                                   resname=res.name, altloc=alt, occupancy=atom.occ,
                                   element=atom.element.name.upper(), xyz=(atom.pos.x, atom.pos.y, atom.pos.z)))
    return atoms, present


def _chain_residues(chain: str, atoms: list, source: str) -> ChainCA:
    """One chain's residues from its Cα atoms: one Cα per residue (highest-occupancy altloc)."""
    groups: dict = {}
    for a in atoms:
        groups.setdefault((a.resnum, a.icode), []).append(a)
    resnums, letters, xyz, modified = [], [], [], []
    for (num, icode), group in groups.items():
        if icode.strip():
            raise ParseError(f"{source}: chain {chain} residue {num}{icode} has an insertion code; "
                             f"renumber the file so every residue has its own number")
        if len(group) > 1:
            alts = [a.altloc for a in group]
            if " " in alts or len(set(alts)) < len(alts):
                raise ParseError(f"{source}: chain {chain} has residue number {num} more than once "
                                 f"(a duplicate residue number)")
        best = max(group, key=lambda a: a.occupancy)          # altlocs: the highest occupancy wins (S-8)
        _, letter, is_modified = _residue_info(best.resname)
        resnums.append(num)
        letters.append(letter)
        xyz.append(best.xyz)
        if is_modified:
            modified.append((num, best.resname))
    return ChainCA(np.array(resnums, int), "".join(letters), np.array(xyz, float).reshape(-1, 3), modified)


def _select_chains(atoms, present, chains, ignore_chains, source) -> Chains:
    by_chain: dict = {}
    for a in atoms:
        is_amino, _, _ = _residue_info(a.resname)
        if is_amino and a.element == "C":                    # a calcium named CA is not a Cα (D-039 item 3)
            by_chain.setdefault(a.chain, []).append(a)
    protein = set(by_chain)
    missing = [c for c in chains if c not in protein]
    if missing:
        raise ParseError(f"{source}: chain(s) {missing} are listed in the config but the file has no "
                         f"amino acids in them (chains with protein: {sorted(protein)})")
    unlisted = sorted(protein - set(chains) - set(ignore_chains))
    if unlisted:
        raise ParseError(f"{source}: protein chain(s) {unlisted} are not listed in the config "
                         f"(listed: {list(chains)}); list them under ignore_chains to skip them")
    ignored = []
    for c in sorted(present - set(chains)):
        if c in ignore_chains:
            ignored.append((c, "ignore_chains"))
        elif c not in protein:
            ignored.append((c, "no amino acids"))
    absent = [c for c in ignore_chains if c not in present]
    return Chains({c: _chain_residues(c, by_chain[c], source) for c in chains}, ignored, absent)


def read_text_ca(text: str, chains) -> Chains:
    """Cα coordinates of the listed chains from PDB text (a model file)."""
    atoms, present = _pdb_atoms(text)
    return _select_chains(atoms, present, list(chains), [], "model file")


def read_structure_ca(path, chains, ignore_chains=()) -> Chains:
    """Cα coordinates of the listed chains from a PDB or mmCIF file, optionally gzipped (`.gz`), given
    as a str or Path. A file that cannot be read stops step 1 with a ParseError naming it (D-057)."""
    path = Path(path)
    try:
        raw = path.read_bytes()
        text = (gzip.decompress(raw) if path.suffix.lower() == ".gz" else raw).decode()
        inner = Path(path.stem) if path.suffix.lower() == ".gz" else path
        if inner.suffix.lower() in (".cif", ".mmcif") or text.lstrip().startswith("data_"):
            atoms, present = _gemmi_atoms(path)                   # gemmi reads .gz itself
        else:
            atoms, present = _pdb_atoms(text)
    except (OSError, EOFError, UnicodeDecodeError, ValueError, RuntimeError) as exc:
        raise ParseError(f"{path.name}: cannot be read as a PDB or mmCIF file ({type(exc).__name__}: {exc})") from None
    return _select_chains(atoms, present, list(chains), list(ignore_chains), path.name)


# =========================================================================== archives

MODEL_FILE = re.compile(r"model_(\d+)_seed_(\d+)\.r(\d+)\.pdb$")     # one saved recycle (D-019 item 2)
JOB_PREFIX = re.compile(r"^(.+?)_(?:un)?relaxed_rank_\d+")


@dataclass
class Archive:
    job: str
    labels: pd.DataFrame             # model, seed, recycle, member; row i <-> coords[i]
    coords: np.ndarray               # (n models, n chains, sequence length, 3), file chain order
    sequence: str                    # the models' sequence (every chain, every model)
    omitted: list                    # [(model, seed)]: incomplete predictions, left out (D-011)
    max_recycle: int
    ignored_chains: list             # [{"chain", "reason"}], merged over the model files (D-054 item 3)
    num_recycles: int | None = None  # from the job's config.json, when the archive has one
    n_files: int = 0


def _job_name(base: str, hit) -> str:
    m = JOB_PREFIX.match(base)
    return m.group(1) if m else base[:hit.start()].rstrip("_")


def scan_archive(path, chains) -> Archive:
    """Read every saved-recycle model file of a ColabFold archive, streaming (never seeking)."""
    path = Path(path)
    names: dict = {}                 # (model, seed, recycle) -> [member names]
    jobs: dict = {}
    models: dict = {}                # (model, seed, recycle) -> (n chains, L, 3)
    ignored: set = set()
    first = None                     # (key, sequence) of the first model read
    num_recycles = None
    with tarfile.open(path, "r|gz") as tf:
        for member in tf:
            if not member.isfile():
                continue
            base = member.name.rsplit("/", 1)[-1]
            if base.startswith("._"):
                continue                                          # macOS AppleDouble companion (D-057)
            if base == "config.json":
                try:
                    num_recycles = json.loads(tf.extractfile(member).read()).get("num_recycles")
                except (ValueError, AttributeError):
                    pass
                continue
            hit = MODEL_FILE.search(base)
            if not hit:
                continue                                          # scores, the bare final copy, ...
            key = tuple(int(x) for x in hit.groups())
            names.setdefault(key, []).append(member.name)
            if len(names[key]) > 1:
                continue                                          # reported below with every duplicate
            jobs.setdefault(_job_name(base, hit), []).append(member.name)
            try:
                ca = read_text_ca(tf.extractfile(member).read().decode(), chains)
            except ParseError as exc:
                raise ParseError(f"{path.name}: {member.name}: {exc}") from None
            except (UnicodeDecodeError, ValueError) as exc:       # not a readable PDB file (D-057)
                raise ParseError(f"{path.name}: {member.name}: cannot be read as a PDB file "
                                 f"({type(exc).__name__}: {exc})") from None
            c0 = chains[0]
            for c in chains[1:]:                                  # D-041 check 1
                if ca[c].seq != ca[c0].seq or not np.array_equal(ca[c].resnums, ca[c0].resnums):
                    raise ArchiveError(f"{path.name}: model {key[0]} seed {key[1]} recycle {key[2]}: chain {c}'s "
                                       f"sequence differs from chain {c0}'s{_first_difference(ca[c0].seq, ca[c].seq)}")
            if first is None:
                first = (key, ca[c0].seq, ca[c0].resnums)
            elif ca[c0].seq != first[1] or not np.array_equal(ca[c0].resnums, first[2]):    # D-041 check 2
                raise ArchiveError(f"{path.name}: the models' sequences differ: model {key[0]} seed {key[1]} "
                                   f"recycle {key[2]} vs model {first[0][0]} seed {first[0][1]} recycle "
                                   f"{first[0][2]}{_first_difference(first[1], ca[c0].seq)}")
            models[key] = np.stack([ca[c].xyz for c in chains]).astype(np.float32)
            ignored |= set(ca.ignored)
    duplicates = {k: v for k, v in names.items() if len(v) > 1}
    if duplicates:                                                # D-026: list every duplicate
        lines = [f"  model {k[0]} seed {k[1]} recycle {k[2]}: {v}" for k, v in sorted(duplicates.items())]
        raise ArchiveError(f"{path.name}: duplicate model files:\n" + "\n".join(lines))
    if len(jobs) > 1:
        raise ArchiveError(f"{path.name}: holds more than one ColabFold job: {sorted(jobs)} (D-019: one job per archive)")
    if not models:
        raise ArchiveError(f"{path.name}: no saved-recycle model files (…model_<M>_seed_<S>.r<N>.pdb)")

    max_recycle = max(k[2] for k in models)
    by_prediction: dict = {}
    for m, s, r in models:
        by_prediction.setdefault((m, s), set()).add(r)
    full = set(range(max_recycle + 1))
    omitted = sorted(p for p, recs in by_prediction.items() if recs != full)      # D-011
    keep = sorted(k for k in models if (k[0], k[1]) not in omitted)
    if not keep:                                                  # D-057 item 1: never vanish silently
        found = sorted({k[2] for k in models})
        raise ArchiveError(f"{path.name}: no complete prediction: every one of its {len(by_prediction)} "
                           f"predictions lacks one of recycles 0-{max_recycle} (recycles found: {found})")
    labels = pd.DataFrame([{"model": m, "seed": s, "recycle": r, "member": names[(m, s, r)][0]}
                           for m, s, r in keep], columns=["model", "seed", "recycle", "member"])
    coords = np.stack([models[k] for k in keep])
    return Archive(job=next(iter(jobs)), labels=labels, coords=coords, sequence=first[1], omitted=omitted,
                   max_recycle=max_recycle,
                   ignored_chains=[{"chain": c, "reason": r} for c, r in sorted(ignored)],
                   num_recycles=num_recycles, n_files=sum(len(v) for v in names.values()))


def _first_difference(a: str, b: str) -> str:
    for i, (x, y) in enumerate(zip(a, b), 1):
        if x != y:
            return f" (first at residue {i}: {x} vs {y})"
    return f" (lengths {len(a)} vs {len(b)})" if len(a) != len(b) else ""


# =========================================================================== sequence alignment
#
# Local alignment with identity scoring (D-038): +1 identical, -1 different; a gap costs -5 for its
# first residue and -1 for each further one. Residues outside the best-scoring stretch get no
# number. Ambiguity is defined by the score (D-054 item 6): a position is ambiguous when two
# alignments with the best score treat it differently — paired with different residues, or paired
# in one and unpaired in the other — including alignments that begin or end with a stretch scoring
# exactly 0. It is computed here from forward and backward scores, not from a library's list.

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 1.0, -1.0, -5.0, -1.0
_NEG = -np.inf


@dataclass
class Alignment:
    pairs: dict                      # {1-based position in a: position in b}, one best alignment
    ambiguous: set                   # positions in a that two best alignments treat differently
    partners: dict                   # {position in a: set of positions in b, None = unpaired}, over all best alignments
    score: float


def _forward(S):
    """Best score of an alignment ending at each cell, in three states (Gotoh):
    M ends with a[i] paired to b[j], X with a[i] against a gap, Y with b[j] against a gap."""
    n, m = S.shape
    M = np.full((n + 1, m + 1), _NEG)
    X = np.full((n + 1, m + 1), _NEG)
    Y = np.full((n + 1, m + 1), _NEG)
    j = np.arange(m + 1)
    for i in range(1, n + 1):
        before = np.maximum(np.maximum(M[i - 1, :-1], X[i - 1, :-1]), Y[i - 1, :-1])
        M[i, 1:] = S[i - 1] + np.maximum(before, 0.0)             # 0: a local alignment may start here
        X[i, 1:] = np.maximum(np.maximum(M[i - 1, 1:] + GAP_OPEN, X[i - 1, 1:] + GAP_EXTEND), Y[i - 1, 1:] + GAP_OPEN)
        # Y[i, j] = max(M[i, j-1] + open, X[i, j-1] + open, Y[i, j-1] + extend), unrolled along the row:
        opened = np.maximum(M[i], X[i]) + GAP_OPEN
        best = np.maximum.accumulate(opened - GAP_EXTEND * j)
        Y[i, 1:] = GAP_EXTEND * (j[1:] - 1) + best[:-1]
    return M, X, Y


def _backward(S):
    """Best score still to come after each cell in each state (the mirror image of _forward)."""
    n, m = S.shape
    BM = np.full((n + 2, m + 2), _NEG)
    BX = np.full((n + 2, m + 2), _NEG)
    BY = np.full((n + 2, m + 2), _NEG)
    j = np.arange(m + 1)
    for i in range(n, -1, -1):
        diag = np.full(m + 1, _NEG)
        if i < n:
            diag[:m] = S[i] + BM[i + 1, 1:m + 1]                  # pair (i+1, j+1) next
        down = BX[i + 1, :m + 1]                                  # a[i+1] against a gap next
        # BY[i, j] = max(diag, down + open, BY[i, j+1] + extend), unrolled right to left:
        start = np.maximum(diag, down + GAP_OPEN)
        best = np.maximum.accumulate((start + GAP_EXTEND * j)[::-1])[::-1]
        BY[i, :m + 1] = best - GAP_EXTEND * j
        right = BY[i, 1:m + 2]                                    # b[j+1] against a gap next
        BM[i, :m + 1] = np.maximum(np.maximum(0.0, diag), np.maximum(down + GAP_OPEN, right + GAP_OPEN))
        BX[i, :m + 1] = np.maximum(diag, np.maximum(down + GAP_EXTEND, right + GAP_OPEN))
    return BM[:n + 1, :m + 1], BX[:n + 1, :m + 1], BY[:n + 1, :m + 1]


def align_sequences(a: str, b: str) -> Alignment:
    """Align a to b (D-038); positions are 1-based."""
    n, m = len(a), len(b)
    if not n or not m:
        return Alignment({}, set(), {i: {None} for i in range(1, n + 1)}, 0.0)
    A = np.frombuffer(a.encode(), np.uint8)
    B = np.frombuffer(b.encode(), np.uint8)
    S = np.where(A[:, None] == B[None, :], MATCH, MISMATCH)
    M, X, Y = _forward(S)
    best = M[1:, 1:].max()
    if best <= 0:
        return Alignment({}, set(), {i: {None} for i in range(1, n + 1)}, 0.0)
    BM, BX, _ = _backward(S)
    on_m = M + BM == best                                         # a[i] paired with b[j] in some best alignment
    on_x = X + BX == best                                         # a[i] against a gap in some best alignment
    starts = np.zeros((n + 1, m + 1), bool)
    starts[1:, 1:] = S + BM[1:, 1:] == best                       # some best alignment begins with pair (i, j)
    first_rows = np.nonzero(starts.any(axis=1))[0]
    last_rows = np.nonzero((M == best).any(axis=1))[0]
    latest_start, earliest_end = first_rows.max(), last_rows.min()
    partners, ambiguous = {}, set()
    for i in range(1, n + 1):
        p = set((np.nonzero(on_m[i, 1:])[0] + 1).tolist())
        if on_x[i].any() or latest_start > i or earliest_end < i:
            p.add(None)                                           # unpaired in some best alignment
        partners[i] = p
        if len(p) > 1:
            ambiguous.add(i)
    return Alignment(_traceback(S, M, X, Y, best), ambiguous, partners, float(best))


def _traceback(S, M, X, Y, best) -> dict:
    """One best alignment: the first best end cell, preferring to stop as soon as the score allows."""
    i, j = (int(v) for v in np.argwhere(M == best)[0])
    state, pairs = "M", {}
    while True:
        if state == "M":
            pairs[i] = j
            before = M[i, j] - S[i - 1, j - 1]
            if before == 0:
                break
            state = next(s for s, T in (("M", M), ("X", X), ("Y", Y)) if T[i - 1, j - 1] == before)
            i, j = i - 1, j - 1
        elif state == "X":
            v = X[i, j]
            state = next(s for s, prev in (("M", M[i - 1, j] + GAP_OPEN), ("X", X[i - 1, j] + GAP_EXTEND),
                                           ("Y", Y[i - 1, j] + GAP_OPEN)) if prev == v)
            i -= 1
        else:
            v = Y[i, j]
            state = next(s for s, prev in (("M", M[i, j - 1] + GAP_OPEN), ("Y", Y[i, j - 1] + GAP_EXTEND),
                                           ("X", X[i, j - 1] + GAP_OPEN)) if prev == v)
            j -= 1
    return dict(sorted(pairs.items()))


def _ranges_of(numbers) -> list:
    """Sorted integers -> [[start, end], ...] of consecutive stretches."""
    out = []
    for u in sorted(numbers):
        if out and u == out[-1][1] + 1:
            out[-1][1] = u
        else:
            out.append([u, u])
    return out


def _gaps(numbered) -> list:
    """UniProt ranges inside the numbered span that are not numbered."""
    numbered = set(numbered)
    if not numbered:
        return []
    return _ranges_of(set(range(min(numbered), max(numbered) + 1)) - numbered)


# =========================================================================== numbering the models

@dataclass
class ModelNumbering:
    to_uniprot: dict                 # model residue number -> UniProt number
    report: dict                     # identity, mismatches, unnumbered, gaps (D-043 item 6)


def number_condition(cond, sequence: str, uniprot: str, compared: set) -> ModelNumbering:
    """Number one condition's models (they share one sequence) and check them against UniProt."""
    aln = align_sequences(sequence, uniprot)
    pairs = aln.pairs
    model_of = {u: i for i, u in pairs.items()}
    for u, want in sorted(cond.markers.items()):                  # the archive holds what the config says
        i = model_of.get(u)
        if i is None and any(u in p for p in aln.partners.values()):
            raise NumberingError(f"condition {cond.name!r}: the marker at UniProt {u} is left unpaired: an equally "
                                 f"good alignment ends before it (it is within two residues of an end of the "
                                 f"construct, where its mismatch costs what the next match gains). A marker this "
                                 f"close to an end cannot be checked (a known limitation, D-057)")
        if i is None:
            raise ArchiveError(f"condition {cond.name!r}: marker at UniProt {u}: the models have no residue there "
                               f"(is the marker's position right, and is this the right archive?)")
        got = sequence[i - 1]
        if got != want:
            raise ArchiveError(f"condition {cond.name!r}: marker mismatch at UniProt {u}: the config expects "
                               f"{want}, the models show {got}. Is this the right archive for {cond.name!r}?")
    mismatches, unlisted = [], []
    for i, u in pairs.items():
        got, expected = sequence[i - 1], uniprot[u - 1]
        if got != expected:
            listed = cond.markers.get(u) == got
            mismatches.append({"uniprot": u, "expected": expected, "found": got, "listed": listed})
            if not listed:
                unlisted.append(f"UniProt {u}: {expected} -> {got} (model residue {i})")
    if unlisted:                                                  # D-041 check 4: models are checked everywhere
        raise NumberingError(f"condition {cond.name!r}: model residues differ from UniProt (mismatch) away "
                             f"from the markers: " + "; ".join(unlisted) + ". A change made in every condition "
                             f"is listed as a marker with the same residue in every condition.")
    tied = sorted(i for i in aln.ambiguous if any(p is not None and p in compared for p in aln.partners[i]))
    if tied:                                                      # D-038 rule 2
        detail = "; ".join(f"model residue {i} -> UniProt {sorted(p for p in aln.partners[i] if p)}" for i in tied)
        raise NumberingError(f"condition {cond.name!r}: equally good alignments number fit or measure residues "
                             f"differently (ambiguous): {detail}")
    lacking = sorted(u for u in compared if u not in model_of)
    if lacking:                                                   # D-038 item 3
        raise NumberingError(f"condition {cond.name!r}: fit or measure residues {_ranges_of(lacking)} have no "
                             f"model counterpart (the models lack them)")
    identical = sum(sequence[i - 1] == uniprot[u - 1] for i, u in pairs.items())
    return ModelNumbering(pairs, {"identity": identical / len(pairs), "mismatches": mismatches,
                                  "unnumbered": [i for i in range(1, len(sequence) + 1) if i not in pairs],
                                  "gaps": _gaps(pairs.values())})


def check_conditions_agree(conf: Config, sequences: dict, numbering: dict, compared: set) -> None:
    """D-041 check 3: conditions differ only at listed marker sites, and are numbered alike."""
    ref = conf.control
    for cond in conf.conditions:
        if cond.name == ref:
            continue
        a, b = sequences[ref], sequences[cond.name]
        if len(a) != len(b):
            raise NumberingError(f"conditions {ref!r} and {cond.name!r} have models of different length "
                                 f"({len(a)} vs {len(b)} residues): the constructs differ, not only at markers")
        marker_sites = set(cond.markers) | set(conf.condition(ref).markers)
        outside = [i for i in range(1, len(a) + 1) if a[i - 1] != b[i - 1]
                   and numbering[ref].to_uniprot.get(i) not in marker_sites]
        if outside:
            raise NumberingError(f"conditions {ref!r} and {cond.name!r} differ at model residue(s) {outside}, "
                                 f"which are not listed marker sites")
        ua = {u: i for i, u in numbering[ref].to_uniprot.items()}
        ub = {u: i for i, u in numbering[cond.name].to_uniprot.items()}
        differ = sorted(u for u in compared if ua.get(u) != ub.get(u))
        if differ:
            raise NumberingError(f"conditions {ref!r} and {cond.name!r} are numbered differently at UniProt {differ}")


# =========================================================================== numbering the structure
#
# Each structure chain is aligned to the full UniProt sequence (D-050). A run of consecutively
# numbered residues takes one offset (structure number -> UniProt number) from its aligned part,
# so every residue of the run is numbered, its ends included (D-044). A run with no aligned part
# takes the offset of its nearest aligned neighbour if exactly one such offset makes every letter
# match (D-047). Problems stop step 1 only inside the compared set (D-051, D-052).

@dataclass
class StructureChain:
    numbered: dict                   # UniProt number -> index into the chain's ChainCA arrays
    report: dict
    problems: list                   # messages for problems inside the compared set


def _runs(resnums) -> list:
    """Indices of each stretch of consecutively numbered residues."""
    runs = []
    for k in range(len(resnums)):
        if runs and resnums[k] == resnums[runs[-1][-1]] + 1:
            runs[-1].append(k)
        else:
            runs.append([k])
    return runs


def number_structure_chain(chain: str, ca: ChainCA, uniprot: str, mutations: dict, compared: set) -> StructureChain:
    n_uni = len(uniprot)
    aln = align_sequences(ca.seq, uniprot)
    runs = _runs(ca.resnums)
    offset, state = [None] * len(runs), [None] * len(runs)
    for r, ks in enumerate(runs):
        candidates, always = set(), []
        for k in ks:
            p = aln.partners.get(k + 1, {None})
            candidates |= {u - ca.resnums[k] for u in p if u is not None}
            if None not in p:
                always.append((k, p))                             # paired in every best alignment
        valid = {d for d in candidates if all(ca.resnums[k] + d in p for k, p in always)}
        if len(valid) == 1:
            offset[r], state[r] = valid.pop(), "aligned"
        elif not candidates or len(valid) > 1:
            state[r] = "unplaced"                                 # no aligned part, or a tie: D-047
        else:
            state[r] = "contradiction"                            # its own numbers disagree with the sequence

    def letters_fit(ks, d):
        for k in ks:
            u = ca.resnums[k] + d
            if not 1 <= u <= n_uni or ca.seq[k] not in (uniprot[u - 1], mutations.get(u)):
                return False
        return True

    aligned = [r for r in range(len(runs)) if state[r] == "aligned"]
    problems, unnumbered = [], []
    for r, ks in enumerate(runs):
        if state[r] == "aligned":
            continue
        before = max((x for x in aligned if x < r), default=None)
        after = min((x for x in aligned if x > r), default=None)
        if state[r] == "unplaced":
            fits = {offset[x] for x in (before, after) if x is not None and letters_fit(ks, offset[x])}
            if len(fits) == 1:
                offset[r], state[r] = fits.pop(), "placed"
                continue
        lo = offset[before] + ca.resnums[runs[before][-1]] if before is not None else 0
        hi = offset[after] + ca.resnums[runs[after][0]] if after is not None else n_uni + 1
        span = f"structure residues {ca.resnums[ks[0]]}-{ca.resnums[ks[-1]]} of chain {chain}"
        if any(lo < u < hi for u in compared):
            why = ("number consecutively across residues the sequence needs a gap for (a numbering "
                   "contradiction)" if state[r] == "contradiction" else
                   "form a run that cannot be placed: neither neighbouring run's offset makes its letters "
                   "match UniProt (was it renumbered?)")
            problems.append(f"{span} {why}; it must lie between UniProt {lo} and {hi}, which includes fit or "
                            f"measure residues")
        unnumbered += [int(ca.resnums[k]) for k in ks]

    numbered, where = {}, {}
    for r, ks in enumerate(runs):
        if offset[r] is None:
            continue
        for k in ks:
            u = int(ca.resnums[k] + offset[r])
            if not 1 <= u <= n_uni:
                unnumbered.append(int(ca.resnums[k]))             # e.g. a tag numbered below 1
            else:
                where.setdefault(u, []).append(k)
    for u, ks in where.items():
        if len(ks) == 1:
            numbered[u] = ks[0]
        else:
            if u in compared:
                problems.append(f"chain {chain}: structure residues {[int(ca.resnums[k]) for k in ks]} are all "
                                f"numbered as UniProt {u}")
            unnumbered += [int(ca.resnums[k]) for k in ks]

    mismatches = {"compared": [], "outside": []}
    for u in sorted(numbered):
        found, expected = ca.seq[numbered[u]], uniprot[u - 1]
        if found != expected:
            listed = mutations.get(u) == found
            mismatches["compared" if u in compared else "outside"].append(
                {"uniprot": u, "expected": expected, "found": found, "listed": listed})
            if u in compared and not listed:
                problems.append(f"chain {chain}: structure residue {int(ca.resnums[numbered[u]])} (UniProt {u}) is "
                                f"{found} where UniProt has {expected} (an unlisted mismatch in the fit or measure "
                                f"set); if the structure really differs there, list it under [reference] mutations")
    identical = sum(ca.seq[k] == uniprot[u - 1] for u, k in numbered.items())
    report = {"identity": identical / len(numbered) if numbered else 0.0,
              "mismatches": mismatches,
              "listed": [{"uniprot": u, "found": ca.seq[numbered[u]] if u in numbered else None}
                         for u in sorted(mutations)],
              "unnumbered": sorted(unnumbered),
              "gaps": _gaps(numbered),
              "modified": [{"resnum": int(num), "name": name} for num, name in ca.modified]}
    return StructureChain(numbered, report, problems)


@dataclass
class StructureNumbering:
    chains: Chains
    per_chain: dict                  # chain -> StructureChain
    ignored_chains: list             # [{"chain", "reason"}]
    no_effect: list                  # config entries that changed nothing (D-052)


def number_structure(conf: Config, compared: set, fit: set) -> StructureNumbering:
    ref = conf.reference
    chains = read_structure_ca(ref.path, ref.chains, ref.ignore_chains)
    per_chain = {c: number_structure_chain(c, chains[c], conf.uniprot, ref.mutations, compared) for c in ref.chains}
    problems = [p for c in ref.chains for p in per_chain[c].problems]
    if problems:
        raise NumberingError(f"{ref.path.name}: the structure cannot be numbered where it matters:\n  "
                             + "\n  ".join(problems))
    lacking = [f"chain {c}: {_ranges_of(fit - set(per_chain[c].numbered))}" for c in ref.chains
               if fit - set(per_chain[c].numbered)]
    if lacking:
        raise NumberingError(f"{ref.path.name}: fit-set residues are missing from the structure (every fit "
                             f"residue must be resolved in every chain): " + "; ".join(lacking))
    no_effect = []
    for u, letter in sorted(ref.mutations.items()):
        shown = [per_chain[c].report["listed"] for c in ref.chains]
        found = [next(x["found"] for x in s if x["uniprot"] == u) for s in shown]
        if letter == conf.uniprot[u - 1]:
            no_effect.append({"key": "reference.mutations", "item": str(u), "reason": "same as UniProt"})
        elif all(f is None for f in found):
            no_effect.append({"key": "reference.mutations", "item": str(u), "reason": "not resolved"})
        elif letter not in found:
            no_effect.append({"key": "reference.mutations", "item": str(u), "reason": "not shown"})
    for c in chains.absent_ignored:
        no_effect.append({"key": "reference.ignore_chains", "item": c, "reason": "not in the structure"})
    return StructureNumbering(chains, per_chain, [{"chain": c, "reason": r} for c, r in chains.ignored], no_effect)


# =========================================================================== geometry
#
# A model's chain labels carry no meaning, so each model is matched to the frame by trying the
# chain mappings that keep ring neighbours together: the n cyclic shifts in both directions around
# the ring (one direction only for two chains). The fit set picks the direction (the other one is
# upside down); a progressive start gives each model a first shift (D-055); D-016 then lets each
# model of a group take the shift that best matches the group's average model over the measure
# set, and D-035 turns the whole group to best match the frame.

def ring_order(centroids: np.ndarray) -> np.ndarray:
    """Order of the chains around the ring's axis, from their centroids."""
    n = len(centroids)
    if n <= 2:
        return np.arange(n)
    c = centroids - centroids.mean(0)
    axis = np.linalg.svd(c)[2][-1]                                # normal to the plane of the centroids
    e1 = c[0] - np.dot(c[0], axis) * axis
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    return np.argsort(np.arctan2(c @ e2, c @ e1))


def fit_transform(P: np.ndarray, Q: np.ndarray):
    """Rotation R and centres (p0, q0) with (P - p0) @ R + q0 closest to Q (Kabsch), and the RMSD."""
    p0, q0 = P.mean(0), Q.mean(0)
    Pc, Qc = P - p0, Q - q0
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = 1.0 if np.linalg.det(U @ Vt) >= 0 else -1.0               # a rotation, never a mirror image
    R = U @ np.diag([1.0, 1.0, d]) @ Vt
    msd = max((Pc ** 2).sum() + (Qc ** 2).sum() - 2.0 * (S[0] + S[1] + d * S[2]), 0.0) / len(P)
    return R, p0, q0, float(np.sqrt(msd))


def _mappings(frame_ring, model_ring) -> list:
    """Candidate mappings as (direction, shift, perm); perm[slot] = model chain in that frame slot."""
    n = len(frame_ring)
    directions = [model_ring] if n <= 2 else [model_ring, model_ring[::-1]]
    out = []
    for d, direction in enumerate(directions):
        for k in range(n):
            perm = np.empty(n, int)
            perm[frame_ring] = np.roll(direction, -k)             # ring position p holds direction[(p + k) % n]
            out.append((d, k, perm))
    return out


@dataclass
class Labelled:
    perms: np.ndarray                # (N, n): model chain in each frame slot
    coords: np.ndarray               # (N, n, sites, 3): fitted coordinates
    core_fit_rmsd: np.ndarray        # (N,)
    turn: int | None
    margin: float | None


def label_group(X: np.ndarray, frame: np.ndarray, fit_t: np.ndarray, meas_t: np.ndarray,
                compare: np.ndarray) -> Labelled:
    """Fit one group's models (X: N × n chains × sites × 3) onto the frame (n × sites × 3).

    fit_t, meas_t: which sites are in the fit and measure sets; compare: (n × sites) the sites
    used to compare the group's average model with the frame (D-035 item 3). The models must come
    in (model, seed) order: the labelling's start is built in that order (D-055).
    """
    N, n = X.shape[:2]
    frame_fit = frame[:, fit_t].reshape(-1, 3)

    def fitted(m, perm, sites=slice(None)):
        R, p0, q0, rmsd = fit_transform(X[m][perm][:, fit_t].reshape(-1, 3), frame_fit)
        return (X[m][perm][:, sites] - p0) @ R + q0, rmsd

    if n == 1:                                                    # a single chain: no labelling (D-016 item 5)
        out = [fitted(m, np.zeros(1, int)) for m in range(N)]
        return Labelled(np.zeros((N, 1), int), np.array([y for y, _ in out]), np.array([r for _, r in out]), None, None)

    frame_ring = ring_order(frame[:, fit_t].mean(1))
    shifts, choice = [], np.zeros(N, int)
    for m in range(N):                                            # the fit set picks the direction (which side is up)
        cands = _mappings(frame_ring, ring_order(X[m][:, fit_t].mean(1)))
        rmsds = [fit_transform(X[m][p][:, fit_t].reshape(-1, 3), frame_fit)[3] for _, _, p in cands]
        d_best, k_best, _ = cands[int(np.argmin(rmsds))]
        shifts.append([p for d, _, p in cands if d == d_best])    # index = shift
        choice[m] = k_best
    measured = np.array([[fitted(m, shifts[m][k], meas_t)[0] for k in range(n)] for m in range(N)])   # N × n × n × s × 3

    # D-055: a progressive start. With a symmetric core the core-best shifts are effectively random,
    # and a group that starts split evenly between labellings never moves; instead each model, in
    # the order given ((model, seed)), takes the shift closest to the average of those placed so far.
    placed = measured[0, choice[0]].copy()
    for m in range(1, N):
        choice[m] = int(((measured[m] - placed / m) ** 2).sum(axis=(1, 2, 3)).argmin())
        placed += measured[m, choice[m]]

    for _ in range(1000):                                         # D-016: agree with the group's average
        average = measured[np.arange(N), choice].mean(0)
        d2 = ((measured - average) ** 2).sum(axis=(2, 3, 4))      # N × n
        current = d2[np.arange(N), choice]
        better = d2.argmin(axis=1)
        moves = d2[np.arange(N), better] < current - 1e-9 * (1.0 + current)
        if not moves.any():
            break
        choice[moves] = better[moves]
    else:
        raise RuntimeError("the labelling did not settle")        # cannot happen: each round lowers the total

    cmp_in_meas = compare[:, meas_t]
    target = frame[:, meas_t][cmp_in_meas]
    dists = []
    for t in range(n):                                            # D-035: turn the whole group
        average = measured[np.arange(N), (choice + t) % n].mean(0)
        dists.append(float(np.sqrt(((average[cmp_in_meas] - target) ** 2).sum(-1).mean())))
    turn = 0 if dists[0] <= min(dists) + 1e-9 else int(np.argmin(dists))    # an exact tie keeps D-016's result
    margin = min(d for t, d in enumerate(dists) if t != turn) - dists[turn]
    final = (choice + turn) % n
    perms = np.array([shifts[m][final[m]] for m in range(N)])
    out = [fitted(m, perms[m]) for m in range(N)]
    return Labelled(perms, np.array([y for y, _ in out]), np.array([r for _, r in out]), turn, float(margin))


# =========================================================================== DVs and their identities

def group_dvs(Y: np.ndarray, meas_cols: np.ndarray, ref_cols: np.ndarray | None, frame_flat: np.ndarray | None) -> dict:
    """DVs of one group from its fitted coordinates Y (N × sites × 3), all chains' sites flattened."""
    x = Y[:, meas_cols]
    out = {"spread": np.sqrt(((x - x.mean(0)) ** 2).sum(-1).mean(-1))}
    if ref_cols is not None:
        xr, fr = Y[:, ref_cols], frame_flat[ref_cols]
        out["spread_ref"] = np.sqrt(((xr - xr.mean(0)) ** 2).sum(-1).mean(-1))
        out["deviation"] = np.sqrt(((xr - fr) ** 2).sum(-1).mean(-1))
        out["offset"] = np.full(len(Y), np.sqrt(((xr.mean(0) - fr) ** 2).sum(-1).mean()))
    return out


def check_identities(dv: pd.DataFrame, coords: np.ndarray, sites: pd.DataFrame, rtol: float = 1e-6,
                     by=("condition", "recycle")) -> None:
    """Raise IdentityError unless, in every group (condition × recycle, or the columns `by`; step 2 groups
    by condition × network × recycle, D-070):
    - mean(deviation²) = mean(spread_ref²) + offset² (with a reference), and
    - the mean pairwise squared distance between models equals 2·mean(spread²) on the in_measure
      sites' coordinates (and, with a reference, the same for spread_ref on the in_measure and
      resolved sites). The second check reads the coordinates, so it also catches wrong rows."""
    from scipy.spatial.distance import pdist
    x = np.asarray(coords, np.float64).reshape(len(dv), len(sites), 3)
    meas = sites.in_measure.to_numpy(bool)
    checks = [("spread", meas)]
    with_ref = "spread_ref" in dv.columns
    if with_ref:
        checks.append(("spread_ref", meas & sites.resolved.to_numpy(bool)))
    failures = []

    def close(a, b):
        return abs(a - b) <= rtol * max(abs(a), abs(b)) + 1e-12

    for key, idx in dv.groupby(list(by)).indices.items():
        g = dv.iloc[idx]
        name = " ".join(f"{k} {v}" for k, v in zip(by, key if isinstance(key, tuple) else (key,)))
        if with_ref:
            lhs = (g.deviation ** 2).mean()
            rhs = (g.spread_ref ** 2).mean() + g.offset.iloc[0] ** 2
            if not close(lhs, rhs) or g.offset.nunique() != 1:
                failures.append(f"{name}: mean(deviation²) {lhs:.6g} != mean(spread_ref²) + offset² {rhs:.6g}")
        for column, mask in checks:
            v = x[idx][:, mask].reshape(len(idx), -1)
            pairwise = pdist(v, "sqeuclidean").sum() / mask.sum() / len(idx) ** 2 if len(idx) > 1 else 0.0
            spread2 = (g[column] ** 2).mean()
            if not close(spread2, pairwise):
                failures.append(f"{name}: mean({column}²) {spread2:.6g} != pairwise identity {pairwise:.6g}")
    if failures:
        raise IdentityError("identity checks failed:\n  " + "\n  ".join(failures))


# =========================================================================== provenance

def file_digest(path: Path) -> dict:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return {"sha256": h.hexdigest(), "size": Path(path).stat().st_size}


def code_fingerprint() -> str:
    """SHA-256 of build_dv.py's bytes followed by config.py's bytes (D-036 item 3)."""
    return hashlib.sha256((HERE / "build_dv.py").read_bytes() + (HERE / "config.py").read_bytes()).hexdigest()


def _git(*args) -> str:
    try:
        return subprocess.run(["git", "-C", str(PROJECT), *args], capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return ""


def _versions() -> dict:
    import Bio
    import scipy
    return {"python": platform.python_version(), "numpy": np.__version__, "gemmi": gemmi.__version__,
            "biopython": Bio.__version__, "pandas": pd.__version__, "scipy": scipy.__version__}


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _write_json(path: Path, data) -> None:
    path.write_text(json.dumps(_jsonable(data), indent=2) + "\n")


class Log:
    """Prints a step's messages and keeps them for its log file (report/step1.log, regions/step2.log)."""

    def __init__(self):
        self.lines = []

    def __call__(self, *parts):
        text = " ".join(str(p) for p in parts)
        print(text, flush=True)
        self.lines.append(text)


# =========================================================================== the parse layer (D-024)

def _parse_layer(conf: Config, cond, digest: dict, code_sha: str, out: Path, log: Log):
    """Read an archive, or reuse its parse/ files when archive, chains and code are unchanged."""
    npy, csv, rec_path = (out / "parse" / f"{cond.name}{ext}" for ext in (".npy", ".csv", ".json"))
    if rec_path.exists() and npy.exists() and csv.exists():
        rec = json.loads(rec_path.read_text())
        if (rec.get("archive") == digest and rec.get("model_chains") == conf.model_chains
                and rec.get("code_sha256") == code_sha):
            log(f"  {cond.name}: parse/ reused ({rec['n_models']} models)")
            archive = Archive(job=rec["job"], labels=pd.read_csv(csv), coords=np.load(npy), sequence=rec["sequence"],
                              omitted=[tuple(p) for p in rec["omitted"]], max_recycle=rec["max_recycle"],
                              ignored_chains=rec["ignored_chains"], num_recycles=rec.get("num_recycles"),
                              n_files=rec.get("n_files", 0))
            return archive, "reused"
    rec_path.unlink(missing_ok=True)                              # the record is written last
    log(f"  {cond.name}: reading {cond.archive} ...")
    try:
        archive = scan_archive(cond.archive, conf.model_chains)
    except (ArchiveError, ParseError) as exc:
        raise type(exc)(f"condition {cond.name!r}: {exc}") from None
    np.save(npy, archive.coords)
    archive.labels.to_csv(csv, index=False)
    _write_json(rec_path, {"archive": digest, "model_chains": conf.model_chains, "code_sha256": code_sha,
                           "job": archive.job, "sequence": archive.sequence, "max_recycle": archive.max_recycle,
                           "omitted": archive.omitted, "ignored_chains": archive.ignored_chains,
                           "num_recycles": archive.num_recycles, "n_files": archive.n_files,
                           "n_models": len(archive.labels)})
    log(f"  {cond.name}: {len(archive.labels)} models in {len(archive.labels.groupby(['model', 'seed']))} "
        f"complete predictions; {len(archive.omitted)} omitted")
    return archive, "built"


# =========================================================================== the frame without a reference (D-040)

def _frame_from_config(conf: Config, archives: dict) -> tuple:
    f = conf.frame
    a = archives[f.condition]
    recycle = a.max_recycle if f.recycle is None else f.recycle
    available = sorted({(int(m), int(s)) for m, s in zip(a.labels.model, a.labels.seed)})
    if (f.model, f.seed) in a.omitted:
        raise ConfigError(f"[frame] model {f.model} seed {f.seed} of {f.condition!r} is an omitted (incomplete) "
                          f"prediction; complete predictions: {available}")
    if (f.model, f.seed) not in available:
        raise ConfigError(f"[frame] model {f.model} seed {f.seed} is not in {f.condition!r}; complete "
                          f"predictions (model, seed): {available}")
    if not 0 <= recycle <= a.max_recycle:
        raise ConfigError(f"[frame] recycle {recycle} does not exist; {f.condition!r} has recycles 0-{a.max_recycle}")
    row = int(np.nonzero(((a.labels.model == f.model) & (a.labels.seed == f.seed) & (a.labels.recycle == recycle)).to_numpy())[0][0])
    return f.condition, row, recycle


def choose_frame_model(conf: Config, archives: dict, site_X: dict, masks: dict, log: Log) -> tuple:
    """The control's lowest-spread complete model at its last recycle, confirmed by refitting."""
    a = archives[conf.control]
    recycle = a.max_recycle
    rows = np.nonzero((a.labels.recycle == recycle).to_numpy())[0]
    keys = [(int(a.labels.model.iloc[r]), int(a.labels.seed.iloc[r])) for r in rows]
    target, info = lowest_spread_pick(site_X[conf.control][rows], keys, masks, log)
    return conf.control, int(rows[target]), recycle, info


def lowest_spread_pick(X: np.ndarray, keys: list, masks: dict, log: Log, detail: bool = False) -> tuple:
    """D-040's representative model of one group (X: models × chains × sites × 3, keys: (model, seed)):
    fit every model onto a target (first the lowest model, then seed), pick the lowest spread, and repeat
    with the pick as the target until the pick no longer changes. Returns (index into X, info); with
    `detail`, info also holds the pick's spread, the runner-up and the gap (D-040 item 6). The screen
    (screen.py) uses it again on the models that pass, when there is no structure (D-072)."""
    order = sorted(range(len(keys)), key=lambda i: keys[i])
    target = order[0]                                             # provisional: lowest model, then seed
    provisional, targets, rounds = keys[target], [target], 0
    meas_cols = np.tile(masks["meas_t"], X.shape[1])
    while True:
        rounds += 1
        lab = label_group(X, X[target], masks["fit_t"], masks["meas_t"], masks["frame_compare"])
        stored = lab.coords.reshape(len(keys), -1, 3).astype(np.float32).astype(np.float64)   # as in dv.csv
        spread = group_dvs(stored, meas_cols, None, None)["spread"]
        ranked = sorted(range(len(keys)), key=lambda i: (round(float(spread[i]), 12), keys[i]))
        pick = ranked[0]
        log(f"  frame pick, round {rounds}: target model {keys[target][0]} seed {keys[target][1]} -> lowest spread "
            f"model {keys[pick][0]} seed {keys[pick][1]} ({spread[pick]:.4f} Å)")
        if pick == target:
            break
        if pick in targets:                                       # D-046 item 2
            path = " -> ".join(f"m{keys[t][0]}s{keys[t][1]}" for t in targets + [pick])
            raise ConfigError(f"the frame pick never settles ({path}); name a frame model under [frame]")
        targets.append(pick)
        target = pick
    info = {"provisional": {"model": provisional[0], "seed": provisional[1]}, "rounds": rounds}
    if detail:
        info["spread"] = float(spread[target])
        if len(ranked) > 1:
            r = ranked[1]
            info["runner_up"] = {"model": keys[r][0], "seed": keys[r][1], "spread": float(spread[r])}
            info["gap"] = float(spread[r] - spread[target])
        else:
            info["runner_up"], info["gap"] = None, None
    return target, info


# =========================================================================== the alignment file (D-038 item 4)

def write_alignment(path: Path, uniprot: str, rows: list) -> None:
    """Clustal file in UniProt columns; rows = [(name, [(UniProt number or None, letter), ...])].
    Residues without a UniProt number sit in gapped UniProt columns next to their neighbours. Written
    here rather than with Biopython, whose writer cuts names to 30 characters (D-057 item 8)."""
    n = len(uniprot)
    placed = []
    for name, residues in rows:
        blocks, letters, last = {}, {}, 0
        numbers = [u for u, _ in residues if u is not None]
        for u, letter in residues:
            if u is not None and u > last:
                letters[u], last = letter, u
            else:                                                  # before the first numbered residue: just before it
                block = last if last or not numbers else numbers[0] - 1
                blocks.setdefault(block, []).append(letter)
        placed.append((name, blocks, letters))
    width = {b: max(len(bl.get(b, [])) for _, bl, _ in placed) for b in range(n + 1)}
    records = []
    for name, blocks, letters in placed:
        text = []
        for b in range(n + 1):
            ins = blocks.get(b, [])
            text.append("".join(ins) + "-" * (width[b] - len(ins)))
            if b < n:
                text.append(letters.get(b + 1, "-"))
        records.append((name, "".join(text)))
    name_width = max(len(name) for name, _ in records) + 4
    lines = ["CLUSTAL W multiple sequence alignment (step 1, build_dv.py; columns follow UniProt)", "", ""]
    for start in range(0, len(records[0][1]), 60):
        lines += [f"{name:<{name_width}}{seq[start:start + 60]}" for name, seq in records]
        lines.append("")
    path.write_text("\n".join(lines) + "\n")


# =========================================================================== step 1

def _residues(ranges) -> set:
    return {u for a, b in ranges for u in range(a, b + 1)}


def build(conf: Config) -> dict:
    """Run step 1 for a loaded config; returns {"parse": {condition: "built"|"reused"},
    "aligned": "built"|"reused", ...}. Raises one of STEP1_ERRORS on a problem."""
    log = Log()
    out = Path(conf.output_dir)
    for sub in ("parse", "aligned", "report"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    for old in ("data_report.json", "alignment.aln"):             # D-057 item 7: no stale report
        (out / "report" / old).unlink(missing_ok=True)
    try:
        return _build(conf, out, log)
    except Exception as exc:
        log.lines.append(f"step 1 stopped ({type(exc).__name__}): {exc}")      # in the log, not printed twice
        raise
    finally:
        (out / "report" / "step1.log").write_text("\n".join(log.lines) + "\n")


def _build(conf: Config, out: Path, log: Log) -> dict:
    fit, measure = _residues(conf.fit), _residues(conf.measure)
    compared = fit | measure
    uniprot = conf.uniprot
    log(f"step 1: {conf.name} ({conf.path})")
    log(f"  UniProt {len(uniprot)} residues; fit {len(fit)}, measure {len(measure)} residues per chain; "
        f"{len(conf.conditions)} conditions, control {conf.control!r}")

    structure = None                                              # D-050: before the archives are read
    if conf.reference is not None:
        structure = number_structure(conf, compared, fit)
        for c in conf.reference.chains:
            r = structure.per_chain[c].report
            log(f"  structure chain {c}: identity {r['identity']:.3f}, {len(r['mismatches']['compared'])} "
                f"mismatches in the compared set, {len(r['mismatches']['outside'])} outside")

    digests = {c.name: file_digest(c.archive) for c in conf.conditions}
    by_digest: dict = {}
    for name, d in digests.items():
        by_digest.setdefault(d["sha256"], []).append(name)
    copies = [names for names in by_digest.values() if len(names) > 1]
    if copies:                                                    # D-037 item 2: a renamed copy
        raise ArchiveError(f"conditions {copies[0]} have identical archives (the same checksum)")

    code_sha = code_fingerprint()
    archives, parse_status = {}, {}
    log("archives:")
    for cond in conf.conditions:
        archives[cond.name], parse_status[cond.name] = _parse_layer(conf, cond, digests[cond.name], code_sha, out, log)
    highest = {name: a.max_recycle for name, a in archives.items()}
    if len(set(highest.values())) > 1:                           # D-042
        detail = ", ".join(f"{n}: r{h}" + (f" (config.json num_recycles {archives[n].num_recycles})"
                                           if archives[n].num_recycles is not None else "") for n, h in highest.items())
        raise ArchiveError(f"conditions have different highest recycles ({detail}); every condition needs the same")

    numbering = {}
    for cond in conf.conditions:                                  # D-041, D-038
        numbering[cond.name] = number_condition(cond, archives[cond.name].sequence, uniprot, compared)
    check_conditions_agree(conf, {n: a.sequence for n, a in archives.items()}, numbering, compared)

    # the site table (D-025): fit ∪ measure residues of every chain, in frame-slot order
    us = sorted(compared)
    n_chains = len(conf.model_chains)
    model_of = {u: i for i, u in numbering[conf.control].to_uniprot.items()}
    frame_chains = conf.reference.chains if structure else conf.model_chains
    rows = []
    for j in range(n_chains):
        numbered = structure.per_chain[frame_chains[j]].numbered if structure else None
        ca = structure.chains[frame_chains[j]] if structure else None
        for u in us:
            have = numbered is None or u in numbered
            rows.append({"chain_index": j, "chain": frame_chains[j], "model_resnum": model_of[u], "uniprot_resnum": u,
                         "ref_resnum": int(ca.resnums[numbered[u]]) if structure and u in numbered else pd.NA,
                         "in_fit": u in fit, "in_measure": u in measure, "resolved": have})
    sites = pd.DataFrame(rows).astype({"ref_resnum": "Int64"})
    fit_t = np.array([u in fit for u in us])
    meas_t = np.array([u in measure for u in us])
    resolved_js = sites.resolved.to_numpy(bool).reshape(n_chains, len(us))
    model_idx = np.array([model_of[u] - 1 for u in us])
    site_X = {name: a.coords[:, :, model_idx, :].astype(np.float64) for name, a in archives.items()}
    for a in archives.values():
        a.coords = None                                           # only the sites are needed from here (D-057)

    archive_digests = {name: digests[name] for name in digests}
    inputs = {"uniprot_fasta": file_digest(conf.uniprot_fasta),
              "reference": file_digest(conf.reference.path) if conf.reference else None}
    man_path = out / "aligned" / "manifest.json"
    aligned_files = [out / "aligned" / f for f in ("coords.npy", "sites.csv", "dv.csv")]
    reuse = False
    if man_path.exists() and all(p.exists() for p in aligned_files):
        old = json.loads(man_path.read_text())
        reuse = (old.get("config_text") == conf.text and old.get("archives") == archive_digests
                 and old.get("inputs") == inputs and old.get("code", {}).get("sha256") == code_sha)
    if reuse:
        log("aligned/ reused: config, archives, structure, FASTA and code unchanged")
        frame_info, groups = old["frame"], old["groups"]
    else:
        man_path.unlink(missing_ok=True)                          # D-054 item 1: before anything is written
        for name in STEP2_FILES:                                  # D-061 item 1: they belong to the old aligned/
            (out / "regions" / name).unlink(missing_ok=True)
        frame_info, groups = _aligned_layer(conf, out, archives, site_X, sites, structure, fit_t, meas_t,
                                            resolved_js, us, frame_chains, log)

    _write_report(conf, out, archives, digests, numbering, structure, frame_info, groups)
    if not reuse or "built" in parse_status.values():
        git_commit = _git("rev-parse", "HEAD").strip() or "unknown"
        _write_json(man_path, {
            "config_text": conf.text, "config": resolved(conf), "archives": archive_digests, "inputs": inputs,
            "code": {"sha256": code_sha, "git_commit": git_commit, "git_dirty": bool(_git("status", "--porcelain").strip())},
            "versions": _versions(), "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "frame": frame_info, "groups": groups})
    log(f"done: {out}")
    return {"parse": parse_status, "aligned": "reused" if reuse else "built", "output_dir": str(out)}


def structure_frame(structure: StructureNumbering, frame_chains, us) -> np.ndarray:
    """The structure's Cα coordinates at every site: frame chains × UniProt residues `us` × xyz, NaN where
    the structure does not resolve the residue. Step 2 (`region_dv.py`) rebuilds the frame with it."""
    frame = np.full((len(frame_chains), len(us), 3), np.nan)
    for j, c in enumerate(frame_chains):
        ca, numbered = structure.chains[c], structure.per_chain[c].numbered
        for t, u in enumerate(us):
            if u in numbered:
                frame[j, t] = ca.xyz[numbered[u]]
    return frame


def _aligned_layer(conf, out, archives, site_X, sites, structure, fit_t, meas_t, resolved_js, us, frame_chains, log):
    """Fit every model onto the frame, compute the DVs, write aligned/ (all but the manifest)."""
    n_chains = len(conf.model_chains)
    meas_cols = np.tile(meas_t, n_chains)
    if structure is not None:
        frame = structure_frame(structure, frame_chains, us)
        compare = resolved_js & meas_t[None, :]
        ref_cols = compare.reshape(-1)
        frame_info = {"kind": "reference"}
    else:
        masks = {"fit_t": fit_t, "meas_t": meas_t, "frame_compare": np.ones((n_chains, len(us)), bool) & meas_t[None, :]}
        if conf.frame is not None:
            cond_name, row, recycle = _frame_from_config(conf, archives)
            extra = None
        else:
            cond_name, row, recycle, extra = choose_frame_model(conf, archives, site_X, masks, log)
        frame = site_X[cond_name][row]
        a = archives[cond_name]
        frame_key = (int(a.labels.model.iloc[row]), int(a.labels.seed.iloc[row]))
        compare = masks["frame_compare"]
        ref_cols = None
        frame_info = {"kind": "model", "chosen_by": "config" if conf.frame is not None else "lowest_spread",
                      "condition": cond_name, "model": frame_key[0], "seed": frame_key[1], "recycle": int(recycle)}
        if extra:
            frame_info.update(extra)

    dv_rows, coord_rows, groups = [], [], []
    for cond in conf.conditions:
        a = archives[cond.name]
        for recycle in sorted(a.labels.recycle.unique()):
            idx = np.nonzero((a.labels.recycle == recycle).to_numpy())[0]
            idx = idx[np.lexsort((a.labels.seed.to_numpy()[idx], a.labels.model.to_numpy()[idx]))]
            lab = label_group(site_X[cond.name][idx], frame, fit_t, meas_t, compare)
            Y32 = lab.coords.reshape(len(idx), -1, 3).astype(np.float32)        # what coords.npy will hold
            dvs = group_dvs(Y32.astype(np.float64), meas_cols, ref_cols, frame.reshape(-1, 3) if ref_cols is not None else None)
            for i, r in enumerate(idx):
                m, s = int(a.labels.model.iloc[r]), int(a.labels.seed.iloc[r])
                row = {"condition": cond.name, **{f: cond.levels[f] for f in conf.factors}, "model": m, "seed": s,
                       "recycle": int(recycle), "prediction": f"{cond.name}_m{m}_s{s}",
                       "mapping": ",".join(conf.model_chains[c] for c in lab.perms[i]),
                       "core_fit_rmsd": float(lab.core_fit_rmsd[i])}
                row.update({k: float(v[i]) for k, v in dvs.items()})
                dv_rows.append(row)
            coord_rows.append(Y32.reshape(len(idx), -1))
            groups.append({"condition": cond.name, "recycle": int(recycle), "n_models": int(len(idx)),
                           "turn": lab.turn, "margin": lab.margin})
            log(f"  {cond.name} recycle {recycle}: {len(idx)} models, spread {dvs['spread'].mean():.3f} Å"
                + (f", deviation {dvs['deviation'].mean():.3f} Å" if "deviation" in dvs else "")
                + (f", turn {lab.turn} (margin {lab.margin:.3f} Å)" if lab.turn is not None else ""))
    dv = pd.DataFrame(dv_rows)
    if frame_info["kind"] == "model":
        g = dv[(dv.condition == frame_info["condition"]) & (dv.recycle == frame_info["recycle"])]
        is_frame = (g.model == frame_info["model"]) & (g.seed == frame_info["seed"])
        frame_info["spread"] = float(g[is_frame].spread.iloc[0])
        if frame_info["chosen_by"] == "lowest_spread":
            ranked = g.sort_values(["spread", "model", "seed"])
            runner = ranked[~((ranked.model == frame_info["model"]) & (ranked.seed == frame_info["seed"]))]
            if len(runner):
                r0 = runner.iloc[0]
                frame_info["runner_up"] = {"model": int(r0.model), "seed": int(r0.seed), "spread": float(r0.spread)}
                frame_info["gap"] = float(r0.spread) - frame_info["spread"]
            else:
                frame_info["runner_up"], frame_info["gap"] = None, None

    np.save(out / "aligned" / "coords.npy", np.concatenate(coord_rows).astype(np.float32))
    sites.to_csv(out / "aligned" / "sites.csv", index=False)
    dv.to_csv(out / "aligned" / "dv.csv", index=False)
    text_columns = {c: str for c in ["condition", *conf.factors]}
    check_identities(pd.read_csv(out / "aligned" / "dv.csv", dtype=text_columns, keep_default_na=False),
                     np.load(out / "aligned" / "coords.npy"), pd.read_csv(out / "aligned" / "sites.csv"))
    log("identity checks passed: mean(deviation²) = mean(spread_ref²) + offset²; pairwise = 2·mean(spread²)")
    return frame_info, groups


def _write_report(conf, out, archives, digests, numbering, structure, frame_info, groups) -> None:
    report_archives = {}
    for cond in conf.conditions:
        a = archives[cond.name]
        report_archives[cond.name] = {
            "job": a.job, "sha256": digests[cond.name]["sha256"], "size": digests[cond.name]["size"],
            "n_predictions": int(len(a.labels.groupby(["model", "seed"]))), "n_models": int(len(a.labels)),
            "recycles": sorted(int(r) for r in a.labels.recycle.unique()),
            "seeds": sorted(int(s) for s in a.labels.seed.unique()),
            "omitted": [list(p) for p in a.omitted],
            "markers": {str(u): aa for u, aa in sorted(cond.markers.items())},
            "ignored_chains": a.ignored_chains}
    alignment = {"conditions": {c.name: numbering[c.name].report for c in conf.conditions}}
    if structure is not None:
        alignment["reference"] = {c: structure.per_chain[c].report for c in conf.reference.chains}
        alignment["ignored_chains"] = structure.ignored_chains
    _write_json(out / "report" / "data_report.json", {
        "archives": report_archives, "alignment": alignment,
        "no_effect": structure.no_effect if structure is not None else [],
        "frame": frame_info, "groups": groups})

    rows = [("UniProt", [(u, aa) for u, aa in enumerate(conf.uniprot, 1)])]
    for cond in conf.conditions:
        seq, to_u = archives[cond.name].sequence, numbering[cond.name].to_uniprot
        rows.append((cond.name, [(to_u.get(i), aa) for i, aa in enumerate(seq, 1)]))
    if structure is not None:
        for c in conf.reference.chains:
            ca, numbered = structure.chains[c], structure.per_chain[c].numbered
            u_of = {k: u for u, k in numbered.items()}
            rows.append((f"ref_{c}", [(u_of.get(k), aa) for k, aa in enumerate(ca.seq)]))
    write_alignment(out / "report" / "alignment.aln", conf.uniprot, rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Step 1 of the pipeline: reads each condition's ColabFold archive, numbers models and "
                    "reference structure in UniProt numbering, superposes every model on the frame over the "
                    "fit set, and writes the structural-diversity DVs (spread; with a structure also "
                    "spread_ref, deviation, offset) and a data report to the config's output folder.")
    parser.add_argument("--config", required=True, metavar="datasets/<dataset>/<dataset>.toml",
                        help="the dataset config (TOML)")
    args = parser.parse_args(argv)
    try:
        summary = build(load_config(args.config))
    except STEP1_ERRORS as exc:
        print(f"step 1 stopped ({type(exc).__name__}): {exc}", file=sys.stderr)
        return 1
    print(f"parse/: {summary['parse']}; aligned/: {summary['aligned']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
