"""Load and check a dataset config for steps 1 and 2 (`build_dv.py`, `region_dv.py`).

A dataset config is a TOML file, usually `datasets/<dataset>/<dataset>.toml` (D-029). It names
every input file, the conditions and their archives, the chains, and the residue ranges. Paths are
relative to the config file; residue ranges are inclusive and in UniProt numbering (D-023).

    name = "<dataset>"
    uniprot_fasta = "<UniProt FASTA>"
    output_dir = "<folder>"                 # optional; default <project>/results/<name>
    model_chains = ["A", "B", "C", "D"]
    [reference]                             # optional: an experimental structure (D-006)
    path = "<reference structure file>"     # PDB or mmCIF; author chain IDs and numbers (D-039)
    chains = ["A", "B", "C", "D"]           # same count as model_chains
    mutations = { "<UniProt number>" = "<letter>" }   # optional: known differences (D-041); X = unknown (UNK)
    ignore_chains = ["<chain>"]             # optional: chains of the structure to skip (D-045)
    [ranges]
    fit = [[<start>, <end>], ...]           # the core the models are superposed on (D-014)
    measure = [[<start>, <end>], ...]       # where structural diversity is measured
    treated = [[<start>, <end>], ...]       # residues the treatment changes; may be empty (D-030)
    [regions]                               # optional: regions step 2 measures separately (D-060)
    <region> = [[<start>, <end>], ...]      # inside the measure set; letters, digits and - only
    [design]
    factors = ["<factor>", ...]
    control = "<condition>"
    [frame]                                 # optional, only without [reference] (D-040)
    model = <model>                         # also: condition (default the control),
    seed = <seed>                           #       recycle (default the last)
    [[conditions]]
    name = "<condition>"
    archive = "<archive .tar.gz>"
    levels = { "<factor>" = "<level>", ... }      # text that reads back as itself, e.g. "WT" (D-057)
    markers = { "<UniProt number>" = "<letter>" }     # expected residue at each mutation site
    [screen]                                # optional: the topology screen (D-063, screen.py)
    domain_swapped = true                   # required: each voltage sensor on the next pore domain
                                            # clockwise seen from outside (true) or on its own (false)
    wiring = true                           # optional switches, defaults as shown
    own_place = true
    box_corners = true
    s6_on_own_s5 = true
    filter_quarter = false
    first_recycle = 1                       # lowest recycle judged by the screen and carried into the
                                            # statistics; 0 includes recycle 0 (D-063)
    omit = "<file>"                         # optional CSV: condition,model,seed,reason (removed by hand)
    [segments]                              # with [screen]: one subunit's segments, for four chains
    S1 = [<start>, <end>]                   # ... S2, S3, S4, linker, S5, pore_helix, filter,
                                            # S6_upper, S6_lower, in this order, inside the measure set
    [[segments]]                            # or, for one chain with four repeats (Nav, Cav): four
                                            # tables, one per repeat, in order
    [statistics]                            # optional (row 4, fit_model.py)
    min_models_per_group = 10               # groups (condition × network × recycle) with fewer models are
                                            # left out of the fit (D-068); 0 keeps every group

`load_config(path)` returns a `Config` or raises `ConfigError` naming the offending value. It
checks everything that can be checked without reading the archives or the structure (D-037,
D-040, D-041, D-045, D-046, D-054); it reads only the config and the UniProt FASTA, and creates
nothing.
"""
from __future__ import annotations

import io
import os
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]          # this file lives in <project>/pipeline/

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")               # letters allowed in markers
STRUCTURE_LETTERS = AMINO_ACIDS | {"X"}                 # [reference] mutations may also list X (UNK), D-057
# names become file names (D-046): letters, digits, _ - . only, with at least one letter or digit (D-057)
SAFE_NAME = re.compile(r"[A-Za-z0-9_.-]*[A-Za-z0-9][A-Za-z0-9_.-]*")
# region names become column names spread_<region>, spread_ref_<region>, ... and later file names:
# letters, digits and -, with at least one letter or digit; no _, so a column name reads one way only (D-060)
REGION_NAME = re.compile(r"[A-Za-z0-9-]*[A-Za-z0-9][A-Za-z0-9-]*")
REGION_DVS = ("spread", "spread_ref", "deviation", "offset")    # step 2's columns per region: <dv>_<region>
# dv.csv's own columns; a factor with one of these names would clash with them (D-037 item 4).
RESERVED = {"condition", "model", "seed", "recycle", "prediction", "mapping", "core_fit_rmsd", *REGION_DVS}
TOP_KEYS = {"name", "uniprot_fasta", "output_dir", "model_chains", "reference", "ranges", "regions", "design",
            "frame", "conditions", "screen", "segments", "statistics"}
MIN_MODELS_DEFAULT = 10                                  # D-068
FIRST_RECYCLE_DEFAULT = 1                                # D-063; 0 judges recycle 0 as well
# the topology screen (D-063): the ten segments of a subunit, in sequence order, and the checks' switches
ELEMENTS = ("S1", "S2", "S3", "S4", "linker", "S5", "pore_helix", "filter", "S6_upper", "S6_lower")
SCREEN_SWITCHES = {"wiring": True, "own_place": True, "box_corners": True, "s6_on_own_s5": True,
                   "filter_quarter": False}
OMIT_COLUMNS = ["condition", "model", "seed", "reason"]


class ConfigError(ValueError):
    """The dataset config is invalid; the message names the offending value."""


@dataclass
class Reference:
    path: Path
    chains: list
    mutations: dict = field(default_factory=dict)        # UniProt number -> letter the structure shows
    ignore_chains: list = field(default_factory=list)


@dataclass
class Condition:
    name: str
    levels: dict                                         # factor -> level
    archive: Path
    markers: dict = field(default_factory=dict)          # UniProt number -> expected letter


@dataclass
class Frame:
    condition: str
    model: int
    seed: int
    recycle: int | None = None                           # None = the last recycle


@dataclass
class Screen:                                            # the switches' defaults are SCREEN_SWITCHES
    domain_swapped: bool
    wiring: bool
    own_place: bool
    box_corners: bool
    s6_on_own_s5: bool
    filter_quarter: bool
    first_recycle: int = FIRST_RECYCLE_DEFAULT           # lowest recycle judged (D-063)
    omit: Path | None = None                             # CSV of predictions removed by hand
    omitted: list = field(default_factory=list)          # its rows: (condition, model, seed, reason)


@dataclass
class Config:
    path: Path                                           # the config file
    text: str                                            # its TOML text, exactly as written (D-043)
    name: str
    output_dir: Path
    uniprot_fasta: Path
    uniprot: str                                         # the UniProt sequence read from the FASTA
    model_chains: list
    reference: Reference | None
    factors: list
    conditions: list
    control: str
    fit: list                                            # [(start, end)], inclusive, UniProt numbering
    measure: list
    treated: list
    frame: Frame | None
    regions: dict = field(default_factory=dict)          # name -> [(start, end)], in the order written (D-060)
    screen: Screen | None = None                         # the topology screen's settings (D-063)
    segments: list = field(default_factory=list)         # one table (four chains) or four (one chain)
    min_models_per_group: int = MIN_MODELS_DEFAULT       # the statistics' minimum group size (D-068)

    def condition(self, name: str) -> Condition:
        return next(c for c in self.conditions if c.name == name)


# --------------------------------------------------------------------------- small checks

def _need(table: dict, key: str, where: str):
    if key not in table:
        raise ConfigError(f"{where}: missing '{key}'")
    return table[key]


def _no_unknown_keys(table: dict, allowed: set, where: str) -> None:
    unknown = sorted(set(table) - allowed)
    if unknown:
        raise ConfigError(f"{where}: unknown key(s) {unknown} (allowed: {sorted(allowed)})")


def _is_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _string_list(value, where: str, allow_empty: bool = False) -> list:
    if not isinstance(value, list) or not all(isinstance(x, str) and x for x in value):
        raise ConfigError(f"{where} must be a list of non-empty strings, got {value!r}")
    if not value and not allow_empty:
        raise ConfigError(f"{where} must not be empty")
    repeated = sorted({x for x in value if value.count(x) > 1})
    if repeated:
        raise ConfigError(f"{where} lists {repeated} more than once")
    return list(value)


def _ranges(value, where: str, n_uniprot: int, allow_empty: bool) -> list:
    """[[start, end], ...] -> [(start, end)]: integer pairs with 1 <= start <= end <= UniProt length."""
    if not isinstance(value, list):
        raise ConfigError(f"{where} must be a list of [start, end] pairs, got {value!r}")
    if not value and not allow_empty:
        raise ConfigError(f"{where} must not be empty")
    out = []
    for r in value:
        if not (isinstance(r, list) and len(r) == 2 and all(_is_int(x) for x in r)):
            raise ConfigError(f"{where}: {r!r} is not a pair of integers [start, end]")
        start, end = r
        if not 1 <= start <= end <= n_uniprot:
            raise ConfigError(f"{where}: {r!r} must satisfy 1 <= start <= end <= {n_uniprot} "
                              f"(the UniProt sequence's length)")
        out.append((start, end))
    return out


def _positions(value, where: str, n_uniprot: int, letters=AMINO_ACIDS) -> dict:
    """{"<UniProt number>": "<letter>"} -> {int: letter}, checked against the UniProt sequence."""
    if not isinstance(value, dict):
        raise ConfigError(f"{where} must be a table of UniProt number = letter, got {value!r}")
    out, written = {}, {}
    for key, letter in value.items():
        if not re.fullmatch(r"[0-9]+", str(key)) or not 1 <= int(key) <= n_uniprot:
            raise ConfigError(f"{where}: position {key!r} is not a UniProt residue number (1-{n_uniprot})")
        if int(key) in written:                                          # e.g. "50" and "050" (D-057)
            raise ConfigError(f"{where}: UniProt position {int(key)} is written more than once "
                              f"({written[int(key)]!r} and {key!r})")
        if not (isinstance(letter, str) and letter in letters):
            allowed = "the 20 amino-acid letters" + (" or X (unknown)" if "X" in letters else "")
            raise ConfigError(f"{where}: {key} = {letter!r} is not one of {allowed}")
        out[int(key)] = letter
        written[int(key)] = key
    return out


def _reads_back(value) -> bool:
    """True when a CSV reader (pandas, as later steps use) reads the value back as the same text, not
    as a number, a true/false value or a missing value such as NA (D-057 item 2)."""
    import pandas as pd
    if not isinstance(value, str):
        return False
    back = pd.read_csv(io.StringIO(pd.DataFrame({"v": [value]}).to_csv(index=False)))["v"]
    return pd.api.types.is_string_dtype(back.dtype) and back.iloc[0] == value


def _residues(ranges) -> set:
    return {u for a, b in ranges for u in range(a, b + 1)}


def _as_ranges(numbers) -> list:
    """Integers -> [[start, end], ...] of consecutive stretches, for messages."""
    out = []
    for u in sorted(numbers):
        if out and u == out[-1][1] + 1:
            out[-1][1] = u
        else:
            out.append([u, u])
    return out


def region_columns(name: str) -> list:
    """The columns step 2 may add for a region: spread_<name>, spread_ref_<name>, ..."""
    return [f"{dv}_{name}" for dv in REGION_DVS]


def _regions(value, measure: set, factors: list, n_uniprot: int) -> dict:
    """[regions] -> {name: [(start, end)]}, checked against the measure set and the factors (D-060)."""
    if not isinstance(value, dict):
        raise ConfigError(f"[regions] must be a table of <region> = [[start, end], ...], got {value!r}")
    out, by_case = {}, {}
    for name, ranges in value.items():
        if not REGION_NAME.fullmatch(name):
            raise ConfigError(f"[regions]: name {name!r} may contain only letters, digits and '-', with at least "
                              f"one letter or digit (no '_': region columns are named like spread_ref_<region>)")
        if name.lower() == "ref":
            raise ConfigError(f"[regions]: the name {name!r} is not allowed: the column spread_ref is step 1's")
        if not _reads_back(name):
            raise ConfigError(f"[regions]: name {name!r} would not read back from a table as the same text (a CSV "
                              f"reader takes it for a number, true/false or a missing value); rename it")
        if name.lower() in by_case:
            raise ConfigError(f"[regions]: names {by_case[name.lower()]!r} and {name!r} differ only in case; later "
                              f"steps use region names in file names, and some file systems ignore case")
        by_case[name.lower()] = name
        where = f"[regions] {name}"
        rs = _ranges(ranges, where, n_uniprot, allow_empty=False)
        outside = _residues(rs) - measure
        if outside:
            raise ConfigError(f"{where}: residue(s) {_as_ranges(outside)} are outside the measure set; step 1 "
                              f"keeps coordinates only there, so every region must lie inside [ranges] measure")
        out[name] = rs
    clash = sorted(set(factors) & {c for name in out for c in region_columns(name)})
    if clash:
        raise ConfigError(f"[design] factor(s) {clash} are named like region columns of step 2's dv.csv; rename them")
    return out


def _existing_file(base: Path, value, where: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ConfigError(f"{where} must be a path, got {value!r}")
    path = base / value
    if not path.is_file():
        raise ConfigError(f"{where}: file not found: {path}")
    return path


def read_fasta(path: Path) -> str:
    """The one sequence of a FASTA file, upper case, without whitespace."""
    lines = path.read_text().splitlines()
    headers = [l for l in lines if l.startswith(">")]
    if len(headers) != 1 or not lines[0].startswith(">"):
        raise ConfigError(f"uniprot_fasta: {path} must hold exactly one sequence (found {len(headers)} headers)")
    seq = "".join(l.strip() for l in lines[1:]).upper()
    bad = sorted(set(seq) - AMINO_ACIDS - {"U", "X"})
    if not seq or bad:
        raise ConfigError(f"uniprot_fasta: {path} holds no sequence or unexpected letters {bad}")
    return seq


def _segments(value, measure: set, n_chains: int, n_uniprot: int) -> list:
    """[segments] (one table) or [[segments]] (four tables) -> [{element: (start, end)}], checked."""
    tables = [value] if isinstance(value, dict) else value
    if not isinstance(tables, list) or not all(isinstance(t, dict) for t in tables):
        raise ConfigError("[segments] must be one table, or four [[segments]] tables (one per repeat)")
    if not ((n_chains == 4 and len(tables) == 1) or (n_chains == 1 and len(tables) == 4)):
        raise ConfigError(f"[segments]: the screen works on four-fold channels only: four chains with one "
                          f"[segments] table, or one chain with four [[segments]] tables (one per repeat); "
                          f"this config has {n_chains} chain(s) and {len(tables)} table(s)")
    out, last = [], 0
    for i, t in enumerate(tables):
        where = "[segments]" if len(tables) == 1 else f"[[segments]] #{i + 1}"
        _no_unknown_keys(t, set(ELEMENTS), where)
        table = {}
        for name in ELEMENTS:
            r = _need(t, name, where)
            (start, end), = _ranges([r], f"{where} {name}", n_uniprot, allow_empty=False)
            outside = sorted(set(range(start, end + 1)) - measure)
            if outside:
                raise ConfigError(f"{where} {name} {list(r)}: residue(s) {_ranges_text(outside)} are outside the "
                                  f"measure set; the screen reads step 1's measured coordinates")
            if start <= last:
                raise ConfigError(f"{where} {name} {list(r)} overlaps or comes before the segment before it; "
                                  f"segments must follow the sequence order {list(ELEMENTS)} without overlap")
            table[name] = (start, end)
            last = end
        out.append(table)
    return out


def _ranges_text(residues: list) -> str:
    runs, start = [], residues[0]
    for a, b in zip(residues, residues[1:] + [None]):
        if b != a + 1:
            runs.append(f"{start}-{a}" if a != start else f"{a}")
            start = b
    return ", ".join(runs)


def _screen(value, base: Path, conditions: list) -> Screen:
    """[screen] -> Screen, with its omit file read and checked against the conditions."""
    if not isinstance(value, dict):
        raise ConfigError("[screen] must be a table")
    _no_unknown_keys(value, {"domain_swapped", "omit", "first_recycle", *SCREEN_SWITCHES}, "[screen]")
    swapped = _need(value, "domain_swapped", "[screen]")
    if not isinstance(swapped, bool):
        raise ConfigError(f"[screen] domain_swapped must be true or false, got {swapped!r}")
    switches = {}
    for key, default in SCREEN_SWITCHES.items():
        v = value.get(key, default)
        if not isinstance(v, bool):
            raise ConfigError(f"[screen] {key} must be true or false, got {v!r}")
        switches[key] = v
    first = value.get("first_recycle", FIRST_RECYCLE_DEFAULT)
    if not isinstance(first, int) or isinstance(first, bool) or first < 0:
        raise ConfigError(f"[screen] first_recycle must be a whole number, 0 or more, got {first!r}")
    omit, omitted = None, []
    if "omit" in value:
        omit = _existing_file(base, value["omit"], "[screen] omit")
        import csv
        with open(omit, newline="") as fh:
            rows = list(csv.reader(fh))
        if not rows or [c.strip() for c in rows[0]] != OMIT_COLUMNS:
            raise ConfigError(f"[screen] omit: {omit} must start with the header {','.join(OMIT_COLUMNS)} "
                              f"(each named prediction needs a reason)")
        names = [c.name for c in conditions]
        for k, row in enumerate(rows[1:], 2):
            if not any(x.strip() for x in row):
                continue
            if len(row) != 4:
                raise ConfigError(f"[screen] omit: {omit} line {k} has {len(row)} fields, not 4")
            cond, model, seed, reason = (x.strip() for x in row)
            if cond not in names:
                raise ConfigError(f"[screen] omit: {omit} line {k}: condition {cond!r} is not one of {names}")
            if not (re.fullmatch(r"[0-9]+", model) and re.fullmatch(r"[0-9]+", seed)):
                raise ConfigError(f"[screen] omit: {omit} line {k}: model and seed must be whole numbers")
            if not reason:
                raise ConfigError(f"[screen] omit: {omit} line {k}: the reason is empty")
            if any((cond, int(model), int(seed)) == r[:3] for r in omitted):
                raise ConfigError(f"[screen] omit: {omit} line {k} names condition {cond!r} model {model} seed {seed} "
                                  f"a second time")
            omitted.append((cond, int(model), int(seed), reason))
    return Screen(domain_swapped=swapped, first_recycle=first, omit=omit, omitted=omitted, **switches)


# --------------------------------------------------------------------------- loading

def load_config(path) -> Config:
    """Read and check a dataset config (TOML). Raises ConfigError naming the offending value."""
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"config file not found: {path}")
    text = path.read_text()
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        raise ConfigError(f"{path}: not valid TOML ({exc})") from exc
    base = path.parent
    _no_unknown_keys(data, TOP_KEYS, "config")

    name = _need(data, "name", "config")
    if not (isinstance(name, str) and SAFE_NAME.fullmatch(name)):
        raise ConfigError(f"name {name!r} may contain only letters, digits, '_', '-' and '.', "
                          f"and needs at least one letter or digit")

    fasta = _existing_file(base, _need(data, "uniprot_fasta", "config"), "uniprot_fasta")
    uniprot = read_fasta(fasta)
    n_uni = len(uniprot)

    if "output_dir" in data:
        if not isinstance(data["output_dir"], str) or not data["output_dir"]:
            raise ConfigError(f"output_dir must be a path, got {data['output_dir']!r}")
        output_dir = base / data["output_dir"]                       # relative to the config (D-046)
    else:
        output_dir = PROJECT / "results" / name                      # D-027
    data_dir = (PROJECT / "data").resolve()
    if output_dir.resolve() == data_dir or data_dir in output_dir.resolve().parents:
        raise ConfigError(f"output_dir {output_dir} is inside the project's data/ folder, which holds "
                          f"inputs only; choose another folder")

    model_chains = _string_list(_need(data, "model_chains", "config"), "model_chains")

    reference = None
    if "reference" in data:
        r = data["reference"]
        if not isinstance(r, dict):
            raise ConfigError("[reference] must be a table")
        _no_unknown_keys(r, {"path", "chains", "mutations", "ignore_chains"}, "[reference]")
        chains = _string_list(_need(r, "chains", "[reference]"), "[reference] chains")
        if len(chains) != len(model_chains):
            raise ConfigError(f"[reference] chains {chains} must have as many chains as model_chains "
                              f"{model_chains} (D-028)")
        ignore = _string_list(r.get("ignore_chains", []), "[reference] ignore_chains", allow_empty=True)
        both = sorted(set(chains) & set(ignore))
        if both:
            raise ConfigError(f"[reference] chain(s) {both} are listed in both chains and ignore_chains")
        reference = Reference(path=_existing_file(base, _need(r, "path", "[reference]"), "[reference] path"),
                              chains=chains,
                              mutations=_positions(r.get("mutations", {}), "[reference] mutations", n_uni,
                                                   letters=STRUCTURE_LETTERS),
                              ignore_chains=ignore)

    ranges = _need(data, "ranges", "config")
    if not isinstance(ranges, dict):
        raise ConfigError("[ranges] must be a table")
    _no_unknown_keys(ranges, {"fit", "measure", "treated"}, "[ranges]")
    fit = _ranges(_need(ranges, "fit", "[ranges]"), "[ranges] fit", n_uni, allow_empty=False)
    measure = _ranges(_need(ranges, "measure", "[ranges]"), "[ranges] measure", n_uni, allow_empty=False)
    treated = _ranges(ranges.get("treated", []), "[ranges] treated", n_uni, allow_empty=True)
    shared = sorted(_residues(fit) & _residues(treated))
    if shared:                                                       # D-030, D-037 item 7
        raise ConfigError(f"[ranges] fit and treated share residue(s) {shared}: the fit set must not "
                          f"contain treated residues (D-030)")

    design = _need(data, "design", "config")
    if not isinstance(design, dict):
        raise ConfigError("[design] must be a table")
    _no_unknown_keys(design, {"factors", "control"}, "[design]")
    factors = _string_list(_need(design, "factors", "[design]"), "[design] factors")
    reserved = sorted(set(factors) & RESERVED)
    if reserved:
        raise ConfigError(f"[design] factor(s) {reserved} are reserved: they are dv.csv column names")
    regions = _regions(data.get("regions", {}), _residues(measure), factors, n_uni)

    raw_conditions = _need(data, "conditions", "config")
    if not isinstance(raw_conditions, list) or not raw_conditions:
        raise ConfigError("[[conditions]] must list at least one condition")
    conditions = []
    for i, c in enumerate(raw_conditions):
        where = f"[[conditions]] #{i + 1}"
        if not isinstance(c, dict):
            raise ConfigError(f"{where} must be a table")
        _no_unknown_keys(c, {"name", "archive", "levels", "markers"}, where)
        cname = _need(c, "name", where)
        if not (isinstance(cname, str) and SAFE_NAME.fullmatch(cname)):
            raise ConfigError(f"{where}: name {cname!r} may contain only letters, digits, '_', '-' and '.', "
                              f"with at least one letter or digit (condition names become file names)")
        if not _reads_back(cname):
            raise ConfigError(f"{where}: name {cname!r} would not read back from dv.csv as the same text (a CSV "
                              f"reader takes it for a number, true/false or a missing value); rename it")
        where = f"condition {cname!r}"
        levels = _need(c, "levels", where)
        if not isinstance(levels, dict) or set(levels) != set(factors):
            raise ConfigError(f"{where}: levels {levels!r} must give exactly one level for each factor {factors}")
        unreadable = {f: v for f, v in levels.items() if not _reads_back(v)}
        if unreadable:
            raise ConfigError(f"{where}: level(s) {unreadable} would not read back from dv.csv as the same text "
                              f"(levels are text: a CSV reader takes these for a number, true/false or a "
                              f"missing value); write them as text with a letter, e.g. \"dose10\"")
        conditions.append(Condition(name=cname, levels=dict(levels),
                                    archive=_existing_file(base, _need(c, "archive", where), f"{where} archive"),
                                    markers=_positions(c.get("markers", {}), f"{where} markers", n_uni)))

    names = [c.name for c in conditions]
    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:                                                     # D-037 item 1
        raise ConfigError(f"condition name(s) {repeated} used more than once")
    by_case = {}
    for n in names:                                                  # D-054 item 2
        by_case.setdefault(n.lower(), []).append(n)
    clash = [v for v in by_case.values() if len(v) > 1]
    if clash:
        raise ConfigError(f"condition names {clash[0]} differ only in case; they would share files in "
                          f"parse/ on file systems that ignore case (macOS, Windows)")
    for i, a in enumerate(conditions):                               # D-037 item 2: one file, two conditions
        for b in conditions[i + 1:]:
            if a.archive.resolve() == b.archive.resolve() or os.path.samefile(a.archive, b.archive):
                raise ConfigError(f"conditions {a.name!r} and {b.name!r} name the same archive file {a.archive}")
            if all(a.levels[f] == b.levels[f] for f in factors):     # D-037 item 3
                raise ConfigError(f"conditions {a.name!r} and {b.name!r} have identical levels {a.levels}")

    control = _need(design, "control", "[design]")
    if control not in names:
        raise ConfigError(f"[design] control {control!r} is not one of the conditions {names}")

    frame = None
    if "frame" in data:
        f = data["frame"]
        if reference is not None:                                    # D-040
            raise ConfigError("[frame] names a frame model, which is used only without [reference]; "
                              "with a reference structure the structure is the frame")
        if not isinstance(f, dict):
            raise ConfigError("[frame] must be a table")
        _no_unknown_keys(f, {"condition", "model", "seed", "recycle"}, "[frame]")
        fcond = f.get("condition", control)
        if fcond not in names:
            raise ConfigError(f"[frame] condition {fcond!r} is not one of the conditions {names}")
        for key in ("model", "seed"):
            if not _is_int(_need(f, key, "[frame]")):
                raise ConfigError(f"[frame] {key} must be an integer, got {f[key]!r}")
        if "recycle" in f and not _is_int(f["recycle"]):
            raise ConfigError(f"[frame] recycle must be an integer, got {f['recycle']!r}")
        if "screen" in data:                                       # screen review 2, F3
            first = data["screen"].get("first_recycle", FIRST_RECYCLE_DEFAULT) if isinstance(data["screen"], dict) \
                else FIRST_RECYCLE_DEFAULT
            if _is_int(first) and _is_int(f.get("recycle")) and int(f["recycle"]) < int(first):
                raise ConfigError(f"[frame] recycle {f['recycle']} is below [screen] first_recycle ({first}), so it is "
                                  f"not screened and with a [screen] and no structure it cannot be the representative "
                                  f"model; name a later recycle (or leave recycle out for the last)")
        frame = Frame(condition=fcond, model=f["model"], seed=f["seed"], recycle=f.get("recycle"))

    screen = _screen(data["screen"], base, conditions) if "screen" in data else None
    segments = []
    if "segments" in data:
        segments = _segments(data["segments"], _residues(measure), len(model_chains), n_uni)
    elif screen is not None:
        raise ConfigError("[screen] needs the subunit's segments: add a [segments] table (four chains) or four "
                          "[[segments]] tables (one chain with four repeats)")

    min_models = MIN_MODELS_DEFAULT
    if "statistics" in data:
        st = data["statistics"]
        if not isinstance(st, dict):
            raise ConfigError("[statistics] must be a table")
        _no_unknown_keys(st, {"min_models_per_group"}, "[statistics]")
        min_models = st.get("min_models_per_group", MIN_MODELS_DEFAULT)
        if not (_is_int(min_models) and min_models >= 0):
            raise ConfigError(f"[statistics] min_models_per_group must be a whole number, 0 or more, got {min_models!r}")

    return Config(path=path, text=text, name=name, output_dir=output_dir, uniprot_fasta=fasta,
                  uniprot=uniprot, model_chains=model_chains, reference=reference, factors=factors,
                  conditions=conditions, control=control, fit=fit, measure=measure, treated=treated,
                  frame=frame, regions=regions, screen=screen, segments=segments, min_models_per_group=min_models)


def resolved(config: Config) -> dict:
    """The config's values with defaults filled in, in the TOML's layout (for manifest.json, D-043)."""
    ref = config.reference
    return {
        "name": config.name,
        "uniprot_fasta": str(config.uniprot_fasta),
        "output_dir": str(config.output_dir),
        "model_chains": list(config.model_chains),
        "reference": None if ref is None else {
            "path": str(ref.path), "chains": list(ref.chains),
            "mutations": {str(u): aa for u, aa in sorted(ref.mutations.items())},
            "ignore_chains": list(ref.ignore_chains)},
        "ranges": {key: [list(r) for r in getattr(config, key)] for key in ("fit", "measure", "treated")},
        "regions": {name: [list(r) for r in ranges] for name, ranges in config.regions.items()},
        "design": {"factors": list(config.factors), "control": config.control},
        "frame": None if config.frame is None else {
            "condition": config.frame.condition, "model": config.frame.model,
            "seed": config.frame.seed, "recycle": config.frame.recycle},
        "conditions": [{"name": c.name, "archive": str(c.archive), "levels": dict(c.levels),
                        "markers": {str(u): aa for u, aa in sorted(c.markers.items())}}
                       for c in config.conditions],
        "screen": screen_settings(config),
        "statistics": {"min_models_per_group": config.min_models_per_group},
    }


def screen_settings(config: Config) -> dict | None:
    """The screen's settings with defaults filled in, the predictions named in its omit file, and its
    segments (for the screen's manifest, which step 2 compares with the config); None without [screen].
    The omit file's path is left out, so how it was typed does not matter; its content is here."""
    s = config.screen
    if s is None:
        return None
    return {"domain_swapped": s.domain_swapped, "first_recycle": s.first_recycle,
            **{k: getattr(s, k) for k in SCREEN_SWITCHES},
            "omitted": [list(r) for r in s.omitted],
            "segments": [{k: list(v) for k, v in t.items()} for t in config.segments]}
