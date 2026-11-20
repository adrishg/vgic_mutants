#!/usr/bin/env python3
"""The topology screen: removes models whose fold is wrong (D-063, D-072; plan generalization/notes/screen-plan.md).

Run:  python -B pipeline/screen.py --config datasets/<dataset>/<dataset>.toml
It runs after step 1 (`build_dv.py`) and before step 2 (`region_dv.py`), which then uses the kept
models only.

Each model (one prediction at one recycle) is judged on its own, at every recycle except 0 (D-063
3h). The checks, exactly as first built and validated (D-063, "Decision"; thresholds fixed here, 3d):
1. wiring: every voltage sensor (S1 start to S4 end) touches the pore domain (S5 start to S6 end)
   it should: that pore domain has strictly more Cα pairs within 8 Å than any other (ties and no
   contact fail).
2. own_place: in slices every 2 Å along the model's own pore axis, each element of every subunit is
   nearer its own reference place than any other element's place. The filter is not tested, though
   its reference places still bound the others' (D-063, "Packing-arrangement test").
3. box_corners: the four lower S6s, projected on the model's own membrane plane, keep a box: going
   round the ring, each turns between 45° and 135° to the next.
4. s6_on_own_s5: each S6 (S6 upper start to S6 lower end) has more Cα pairs within 8 Å with its own
   subunit's S5 than with any other subunit's S5.
5. filter_quarter (off by default): each filter's centre is closer, in angle around the pore axis,
   to its own reference filter than to any other.

With a structure in the config, it is the reference (D-063 3e, D-072 item 1):
- wiring expects the structure's partners;
- the box goes round in the structure's order;
- the structure must agree with `domain_swapped` and pass the switched-on box and S6 checks itself
  (D-072 item 3);
- segment residues the structure lacks are recorded and left out of every model before it is
  checked, so model and structure are compared on the same residues. A switched-on check that needs
  a segment the structure lacks entirely stops the run.
Without a structure (D-072 item 2):
- only wiring (each voltage sensor on the pore domain 90° clockwise seen from outside, or on its
  own, from `domain_swapped`), box_corners and s6_on_own_s5 run, each judged on the model itself;
- own_place and filter_quarter are not run;
- the representative model (D-040: the control's model at its last recycle closest to the group's
  average, or the config's [frame]) is then picked among the models that pass, and step 2 fits every
  model onto it.

A subunit is a chain of a four-chain channel, or a repeat of a one-chain channel (Nav, Cav: 3f),
cut out by the config's segments. The screen works on four-fold channels only (3g).

Writes <output>/screen/: verdicts.csv (one row per model: each check's value and result, whether the
user named the prediction, the first failing reason, kept), summary.json (kept per condition, recycle
and network; removed per condition, recycle, network and reason; the named predictions; the
structure's missing residues), manifest.json (written last: the settings, the reference, checksums of
step 1's files and of verdicts.csv, and of this file's code, so step 2 can tell whether the screen is
up to date) and screen.log. Running the screen removes step 2's results, which belonged to the old
verdicts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from build_dv import STEP1_ERRORS, STEP2_FILES, Log, _git, _residues, _versions, _write_json, file_digest, lowest_spread_pick
from config import ELEMENTS, SCREEN_SWITCHES, Config, ConfigError, load_config, screen_settings
from region_dv import SCREEN_CODE_FILES, Step1ResultsError, frame_sites, read_step1

HERE = Path(__file__).resolve().parent
CHECKS = tuple(SCREEN_SWITCHES)          # the config's switches, in the order a first failing reason is named
NEED_REFERENCE = ("own_place", "filter_quarter")         # not run without a structure (D-072 item 2)
CONTACT = 8.0                     # Å: a Cα pair closer than this is a contact (wiring, s6_on_own_s5)
SLICE = 2.0                       # Å between own-place slices, on a grid through depth 0
LINKER_HALF = 3.0                 # Å: the flat linker is tested in slices this close to its mean depth
BOX = (45.0, 135.0)               # degrees: allowed turn from one lower S6 to the next
LINKER, FILTER = ELEMENTS.index("linker"), ELEMENTS.index("filter")
PARTS = (*ELEMENTS, "vsd", "pd", "s5", "s6")              # every coordinate set a subunit holds
STEP1_FILES = ("dv.csv", "sites.csv", "coords.npy")
SCREEN_FILES = ("manifest.json", "verdicts.csv", "summary.json", "screen.log")


class ScreenError(Exception):
    """The screen cannot run: not a four-fold channel, a reference that is unfit or disagrees with
    domain_swapped, or no model to pick the representative model from."""


STEP_ERRORS = (ScreenError, Step1ResultsError, *STEP1_ERRORS)


def code_sha256() -> str:
    """SHA-256 of the code the verdicts depend on (this file, and the step 1 and config code it calls):
    step 2 re-checks it, because changed code may change the verdicts."""
    return hashlib.sha256(b"".join((HERE / f).read_bytes() for f in SCREEN_CODE_FILES)).hexdigest()


# =========================================================================== subunits

@dataclass
class Subunit:
    """One subunit's Cα coordinates (NaN where unresolved) per element, and its residue numbers."""
    el: dict                  # element -> (k, 3)
    res: dict                 # element, "vsd" and "pd" -> UniProt numbers of their rows
    vsd: np.ndarray           # S1 start to S4 end
    pd: np.ndarray            # S5 start to S6 lower end
    s5: np.ndarray
    s6: np.ndarray            # S6 upper start to S6 lower end

    def part(self, name: str) -> np.ndarray:
        return self.el[name] if name in self.el else getattr(self, name)


def split(xyz, u, segments) -> list:
    """Four subunits from coordinates (chains × residues × 3) with UniProt numbers u: four chains
    with one segment table, or one chain with four tables (one per repeat)."""
    xyz, u = np.asarray(xyz, dtype=float), np.asarray(u)
    tables = [segments] if isinstance(segments, dict) else list(segments)
    if xyz.shape[0] == 4 and len(tables) == 1:
        pairs = [(xyz[c], tables[0]) for c in range(4)]
    elif xyz.shape[0] == 1 and len(tables) == 4:
        pairs = [(xyz[0], t) for t in tables]
    else:
        raise ScreenError(f"the screen works on four-fold channels only: four chains with one segment table, or "
                          f"one chain with four repeats (four tables); got {xyz.shape[0]} chain(s) and "
                          f"{len(tables)} table(s)")

    def m(a, b):
        return (u >= a) & (u <= b)

    return [Subunit(el={k: X[m(*t[k])] for k in ELEMENTS},
                    res={**{k: u[m(*t[k])] for k in ELEMENTS}, "vsd": u[m(t["S1"][0], t["S4"][1])],
                         "pd": u[m(t["S5"][0], t["S6_lower"][1])]},
                    vsd=X[m(t["S1"][0], t["S4"][1])], pd=X[m(t["S5"][0], t["S6_lower"][1])],
                    s5=X[m(*t["S5"])], s6=X[m(t["S6_upper"][0], t["S6_lower"][1])]) for X, t in pairs]


def _masked(subs: list, present: list) -> list:
    """The subunits with NaN wherever the reference lacks the residue (D-072 item 1)."""
    out = []
    for s, p in zip(subs, present):
        cut = {name: np.where(p[name][:, None], s.part(name), np.nan) for name in PARTS}
        out.append(Subunit(el={k: cut[k] for k in ELEMENTS}, res=s.res, vsd=cut["vsd"], pd=cut["pd"],
                           s5=cut["s5"], s6=cut["s6"]))
    return out


def _clean(Q: np.ndarray) -> np.ndarray:
    return Q[~np.isnan(Q).any(1)]


def frame(subs: list) -> tuple:
    """centre, unit normal (the fourfold axis, from the linkers to the filters, i.e. towards the
    outside), and an in-plane basis e1, e2 with (e1, e2, normal) right-handed."""
    pdc = np.array([np.nanmean(s.pd, axis=0) for s in subs])
    centre = pdc.mean(0)
    n = np.linalg.svd(pdc - centre)[2][-1]
    filt = np.nanmean(np.concatenate([s.el["filter"] for s in subs]), axis=0)
    link = np.nanmean(np.concatenate([s.el["linker"] for s in subs]), axis=0)
    if (filt - link) @ n < 0:
        n = -n
    e1 = np.cross(n, [1.0, 0, 0]) if abs(n[0]) < 0.9 else np.cross(n, [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    return centre, n, e1, np.cross(n, e1)


def _wrap(a):
    return (np.asarray(a) + 180.0) % 360.0 - 180.0


def _contacts(A: list, B: list) -> np.ndarray:
    """K[i, j] = Cα pairs closer than CONTACT between A[i] and B[j]; NaN never counts."""
    K = np.zeros((len(A), len(B)), int)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            K[i, j] = int((np.linalg.norm(_clean(a)[:, None] - _clean(b)[None], axis=-1) < CONTACT).sum())
    return K


def _top(K: np.ndarray) -> np.ndarray:
    """Per row, the column with strictly most contacts, or -1 (a tie, or no contact)."""
    s = np.sort(K, axis=1)
    return np.where((s[:, -1] > s[:, -2]) & (s[:, -1] > 0), K.argmax(1), -1)


def _ring(subs: list, fr: tuple) -> tuple:
    """From the pore domains' angles about the axis, seen from outside: each subunit's clockwise and
    anticlockwise neighbour, and the angles (degrees)."""
    centre, _, e1, e2 = fr
    pdc = np.array([np.nanmean(s.pd, axis=0) for s in subs]) - centre
    ang = np.degrees(np.arctan2(pdc @ e2, pdc @ e1))
    turn = _wrap(ang[None, :] - ang[:, None])                    # turn[i, j]: from i to j
    np.fill_diagonal(turn, np.nan)
    cw = np.array([int(np.nanargmin(np.abs(turn[i] + 90))) for i in range(4)])
    acw = np.array([int(np.nanargmin(np.abs(turn[i] - 90))) for i in range(4)])
    return cw, acw, ang


# =========================================================================== own place

def _elements(subs: list, fr: tuple) -> dict:
    """(subunit, element) -> (point on its axis line, unit direction, lowest, highest and mean depth)."""
    centre, n = fr[0], fr[1]
    el = {}
    for c, s in enumerate(subs):
        for k, name in enumerate(ELEMENTS):
            Q = _clean(s.el[name])
            if len(Q) < 2:
                continue
            p0 = Q.mean(0)
            d = np.linalg.svd(Q - p0)[2][0]
            z = (Q - centre) @ n
            el[c, k] = (p0, d, z.min(), z.max(), z.mean())
    return el


def _positions(points: np.ndarray, depths: np.ndarray, fr: tuple, el: dict) -> np.ndarray:
    """Where each tested (slice, subunit, element) crosses its slice, in the in-plane basis; the flat
    linker at its centroid."""
    centre, n, e1, e2 = fr
    out = np.full((len(points), 2), np.nan)
    for i, (si, c, k) in enumerate(points):
        if (c, k) not in el:
            continue
        p0, d = el[c, k][0], el[c, k][1]
        if k == LINKER or abs(d @ n) < 1e-6:
            v = p0 - centre
        else:
            t = (centre @ n + depths[si] - p0 @ n) / (d @ n)
            v = p0 + t * d - centre
        out[i] = (v @ e1, v @ e2)
    return out


def _place(pos: np.ndarray, ref: np.ndarray) -> tuple:
    """Turn positions about the axis onto the reference's: the circular median of the angle differences
    (from their circular mean, the differences are wrapped and their median taken)."""
    a_ref = np.arctan2(ref[:, 1], ref[:, 0])
    a_mod = np.arctan2(pos[:, 1], pos[:, 0])
    ok = (np.linalg.norm(ref, axis=1) > 1) & (np.linalg.norm(pos, axis=1) > 1)
    diff = (a_ref - a_mod)[ok]
    mean = np.arctan2(np.sin(diff).mean(), np.cos(diff).mean())
    theta = mean + np.median((diff - mean + np.pi) % (2 * np.pi) - np.pi)
    c, s = np.cos(theta), np.sin(theta)
    return pos @ np.array([[c, s], [-s, c]]), theta


def _arrows(subs: list, fr: tuple) -> list:
    """Each lower S6 as a 2D arrow from its first to its last Cα, projected on the membrane plane."""
    centre, _, e1, e2 = fr
    out = []
    for s in subs:
        Q = _clean(s.el["S6_lower"])
        f = lambda p: np.array([(p - centre) @ e1, (p - centre) @ e2])
        out.append((f(Q[0]), f(Q[-1])))
    return out


def _corners(arrows: list, succ) -> np.ndarray:
    """Per lower S6, how far (degrees) its turn to the next falls outside BOX (0 inside)."""
    v = []
    for c in range(4):
        (a, b), (p, q) = arrows[c], arrows[succ[c]]
        r, s = b - a, q - p
        turn = float(_wrap(np.degrees(np.arctan2(s[1], s[0]) - np.arctan2(r[1], r[0]))))
        v.append(max(0.0, BOX[0] - turn, turn - BOX[1]))
    return np.array(v)


def _filter_angles(subs: list, fr: tuple, theta: float) -> np.ndarray:
    centre, _, e1, e2 = fr
    cs, sn = np.cos(theta), np.sin(theta)
    out = []
    for s in subs:
        v = np.nanmean(s.el["filter"], axis=0) - centre
        a, b = v @ e1, v @ e2
        out.append(np.degrees(np.arctan2(a * sn + b * cs, a * cs - b * sn)))
    return np.array(out)


def _s6_on_own(subs: list) -> np.ndarray:
    K = _contacts([s.s6 for s in subs], [s.s5 for s in subs])
    return np.array([K[i, i] > np.delete(K[i], i).max() for i in range(4)])


# =========================================================================== the reference

@dataclass
class Reference:
    checks: tuple             # the checks it was prepared for (the switched-on ones)
    present: list             # per subunit: part -> bool per residue, False where the structure lacks it
    missing: dict             # subunit (1-4) -> segment, or "between_segments" -> UniProt numbers it lacks
    partner: np.ndarray | None = None   # wiring: the pore domain each voltage sensor should touch
    depths: np.ndarray | None = None    # own_place: slice depths (Å from the pore domain's centre)
    points: np.ndarray | None = None    # own_place: (P, 3) tested (slice index, subunit, element)
    places: np.ndarray | None = None    # own_place: (P, 2) the reference's place of each point
    slices: list = field(default_factory=list)
    untested: list = field(default_factory=list)   # own_place: (subunit, element) the structure lacks
    succ: np.ndarray | None = None      # box_corners: the lower S6 after each, going round the ring
    filters: np.ndarray | None = None   # filter_quarter: each filter's angle about the axis (degrees)


def prepare_reference(subs: list, domain_swapped: bool, checks=CHECKS) -> Reference:
    """Everything the switched-on `checks` compare a model with; ScreenError if the structure is unfit
    for them."""
    checks = tuple(c for c in CHECKS if c in checks)
    present = [{name: ~np.isnan(s.part(name)).any(1) for name in PARTS} for s in subs]
    missing = {}
    for c, s in enumerate(subs):
        gaps = {k: [int(x) for x in s.res[k][~present[c][k]]] for k in ELEMENTS if (~present[c][k]).any()}
        in_segments = set().union(*(set(s.res[k].tolist()) for k in ELEMENTS))
        loops = sorted({int(x) for part in ("vsd", "pd") for x in s.res[part][~present[c][part]]} - in_segments)
        if loops:                                 # loops inside the voltage sensor or pore domain (wiring)
            gaps["between_segments"] = loops
        if gaps:
            missing[c + 1] = gaps
    ref = Reference(checks=checks, present=present, missing=missing)
    uses_frame = [c for c in checks if c != "s6_on_own_s5"]      # every check but the S6 one uses the axis
    if uses_frame:                                              # review 2, F1: what the axis needs
        lacks = [f"subunit {c + 1}'s pore domain" for c, s in enumerate(subs) if len(_clean(s.pd)) < 3]
        lacks += [f"the {k} of every subunit" for k in ("linker", "filter")
                  if all(len(_clean(s.el[k])) == 0 for s in subs)]
        if lacks:
            raise ScreenError(f"the structure lacks {', '.join(lacks)}, which the screen needs to find the pore axis "
                              f"and which side is outside (used by {uses_frame}); use a structure that has them, or "
                              f"switch those checks off")
    if "wiring" in checks:
        gone = [c + 1 for c, s in enumerate(subs) if len(_clean(s.vsd)) == 0]
        if gone:
            raise ScreenError(f"the structure lacks the voltage sensor of subunit(s) {gone}, which wiring needs; switch "
                              f"wiring off, or use a structure that has it")
    fr = frame(subs)

    def lacking(names, need=2):
        return [(c + 1, k) for c, s in enumerate(subs) for k in names if len(_clean(s.part(k))) < need]

    if "wiring" in checks:
        partner = _top(_contacts([s.vsd for s in subs], [s.pd for s in subs]))
        if (partner < 0).any() or sorted(partner) != [0, 1, 2, 3]:
            raise ScreenError(f"the structure's voltage sensors do not each touch one pore domain (partners "
                              f"{(partner + 1).tolist()}, 0 = a tie or no contact); it cannot serve as the reference "
                              f"for wiring")
        cw, _, _ = _ring(subs, fr)
        want = cw if domain_swapped else np.arange(4)
        if not (partner == want).all():
            raise ScreenError(f"the structure disagrees with domain_swapped = {str(domain_swapped).lower()}: its voltage "
                              f"sensors touch pore domains {(partner + 1).tolist()} (subunits 1-4); "
                              + ("swapped means each on the next one clockwise seen from outside"
                                 if domain_swapped else "not swapped means each on its own"))
        ref.partner = partner

    if "own_place" in checks or "filter_quarter" in checks:      # the filter quarter uses own place's turn
        el = _elements(subs, fr)
        if "own_place" in checks:
            ref.untested = [(c + 1, ELEMENTS[k]) for c in range(4) for k in range(len(ELEMENTS)) if (c, k) not in el]
        lo = min(v[2] for v in el.values())
        hi = max(v[3] for v in el.values())
        depths = np.arange(np.floor(lo / SLICE) * SLICE, np.ceil(hi / SLICE) * SLICE + SLICE / 2, SLICE)
        points = []
        for si, z in enumerate(depths):
            for (c, k), (p0, d, zlo, zhi, zm) in el.items():
                if (k == LINKER and abs(z - zm) <= LINKER_HALF) or (k != LINKER and zlo <= z <= zhi):
                    points.append((si, c, k))
        points = np.array(points)
        count = np.bincount(points[:, 0], minlength=len(depths))
        points = points[count[points[:, 0]] >= 2]                # a slice needs at least two elements
        ref.depths, ref.points = depths, points
        ref.places = _positions(points, depths, fr, el)
        ref.slices = [np.nonzero(points[:, 0] == si)[0] for si in range(len(depths)) if (points[:, 0] == si).any()]

    if "box_corners" in checks:
        gone = lacking(["S6_lower"])
        if gone:
            raise ScreenError(f"the structure lacks the lower S6 of subunit(s) {[c for c, _ in gone]}, which box_corners "
                              f"needs; switch box_corners off, or use a structure that has it")
        _, acw, _ = _ring(subs, fr)
        if sorted(acw) != [0, 1, 2, 3] or (_corners(_arrows(subs, fr), acw) > 0).any():
            raise ScreenError("the structure fails its own box_corners check (its lower S6s do not form a box); "
                              "it cannot serve as the reference; switch box_corners off if it does not suit this channel")
        ref.succ = acw

    if "s6_on_own_s5" in checks:
        gone = lacking(["s5", "s6"], need=1)
        if gone:
            raise ScreenError(f"the structure lacks {gone} (subunit, part), which s6_on_own_s5 needs; switch it off, "
                              f"or use a structure that has them")
        if not _s6_on_own(subs).all():
            raise ScreenError("the structure fails its own s6_on_own_s5 check (an S6 packs on another subunit's S5); "
                              "it cannot serve as the reference")

    if "filter_quarter" in checks:
        gone = lacking(["filter"], need=1)
        if gone:
            raise ScreenError(f"the structure lacks the filter of subunit(s) {[c for c, _ in gone]}, which "
                              f"filter_quarter needs; switch filter_quarter off, or use a structure that has it")
        ref.filters = _filter_angles(subs, fr, 0.0)
    return ref


# =========================================================================== one model

def check(subs: list, ref: Reference) -> dict:
    """The value and result of each check the reference was prepared for, on the residues the structure
    has (D-072 item 1)."""
    subs = _masked(subs, ref.present)
    fr = frame(subs)
    out = {}
    if ref.partner is not None:
        wired = int((_top(_contacts([s.vsd for s in subs], [s.pd for s in subs])) == ref.partner).sum())
        out["wiring"], out["wiring_ok"] = wired, wired == 4
    theta = 0.0
    if ref.points is not None:
        pos, theta = _place(_positions(ref.points, ref.depths, fr, _elements(subs, fr)), ref.places)
    if "own_place" in ref.checks:
        v = np.full(len(ref.points), np.inf)
        for idx in ref.slices:
            D = np.linalg.norm(pos[idx][:, None] - ref.places[idx][None], axis=-1)
            own = D[np.arange(len(idx)), np.arange(len(idx))].copy()
            D[np.arange(len(idx)), np.arange(len(idx))] = np.inf
            v[idx] = own - D.min(1)
        v = np.where(np.isnan(v), np.inf, v)                      # an element the model lacks fails
        out["own_place"] = float(v[ref.points[:, 2] != FILTER].max())
        out["own_place_ok"] = bool(out["own_place"] <= 0)
    if ref.succ is not None:
        corners = _corners(_arrows(subs, fr), ref.succ)
        out["box_corners"], out["box_corners_ok"] = float(corners.max()), bool(corners.max() <= 0)
    if "s6_on_own_s5" in ref.checks:
        own = int(_s6_on_own(subs).sum())
        out["s6_on_own_s5"], out["s6_on_own_s5_ok"] = own, own == 4
    if ref.filters is not None:
        d = np.abs(_wrap(_filter_angles(subs, fr, theta)[:, None] - ref.filters[None]))
        q = max(float(d[c, c] - np.delete(d[c], c).min()) for c in range(4))
        out["filter_quarter"], out["filter_quarter_ok"] = q, q <= 0
    return out


def check_without_reference(subs: list, domain_swapped: bool, checks=CHECKS) -> dict:
    """wiring, box_corners and s6_on_own_s5, each judged on the model itself (no structure, D-072 item 2)."""
    fr = frame(subs)
    cw, acw, _ = _ring(subs, fr)
    ring = sorted(cw) == [0, 1, 2, 3] and sorted(acw) == [0, 1, 2, 3]   # each subunit once around the ring (review 2, F2)
    out = {}
    if "wiring" in checks:
        want = cw if domain_swapped else np.arange(4)
        wired = int((_top(_contacts([s.vsd for s in subs], [s.pd for s in subs])) == want).sum())
        out["wiring"], out["wiring_ok"] = wired, ring and wired == 4
    if "box_corners" in checks:
        corners = _corners(_arrows(subs, fr), acw) if ring else np.array([np.inf])
        out["box_corners"], out["box_corners_ok"] = float(corners.max()), bool(corners.max() <= 0)
    if "s6_on_own_s5" in checks:
        own = int(_s6_on_own(subs).sum())
        out["s6_on_own_s5"], out["s6_on_own_s5_ok"] = own, own == 4
    return out


# =========================================================================== the step

def build(conf: Config) -> dict:
    """Run the screen for a loaded config; returns {"models": n, "kept": n, "removed": n, "output_dir": str}.
    Raises one of STEP_ERRORS on a problem."""
    if conf.screen is None:
        raise ConfigError("the config has no [screen] table; add one to run the screen (see config.py)")
    log = Log()
    out = Path(conf.output_dir)
    folder = out / "screen"
    folder.mkdir(parents=True, exist_ok=True)
    results = [name for name in SCREEN_FILES if name != "screen.log"]      # the manifest first
    for name in results:
        (folder / name).unlink(missing_ok=True)
    for name in STEP2_FILES:                       # step 2's results belonged to the old verdicts (review B2)
        (out / "regions" / name).unlink(missing_ok=True)
    try:
        return _build(conf, folder, log)
    except BaseException as exc:
        for name in results:
            (folder / name).unlink(missing_ok=True)
        log.lines.append(f"screen stopped ({type(exc).__name__}): {exc}")
        raise
    finally:
        (folder / "screen.log").write_text("\n".join(log.lines) + "\n")


def _parsed_sites(conf: Config, step1: dict, cname: str) -> tuple:
    """Step 1's parse/ files of one condition: (labels, coordinates at step 1's sites, as float64 on use)."""
    parse = Path(conf.output_dir) / "parse"
    try:
        lab = pd.read_csv(parse / f"{cname}.csv")
        P = np.load(parse / f"{cname}.npy", mmap_mode="r")
    except (OSError, ValueError) as exc:
        raise Step1ResultsError(f"step 1's parse/ files of {cname!r} are missing or damaged ({exc}); re-run step 1") from None
    s0 = step1["sites"][step1["sites"].chain_index == 0]
    return lab, P, s0.model_resnum.to_numpy() - 1


def pick_representative(conf: Config, step1: dict, verdicts: pd.DataFrame, log: Log) -> dict:
    """Without a structure (D-072 item 2): the representative model among the models that pass. The
    config's [frame] if it names one (it must pass), otherwise D-040's lowest-spread model of the
    control's kept models at its last recycle."""
    kept = verdicts[verdicts.kept]
    if conf.frame is not None:
        f = conf.frame
        recycle = f.recycle if f.recycle is not None else int(verdicts[verdicts.condition == f.condition].recycle.max())
        ok = ((kept.condition == f.condition) & (kept.model == f.model) & (kept.seed == f.seed)
              & (kept.recycle == recycle)).any()
        if not ok:
            raise ScreenError(f"the [frame] model (condition {f.condition!r}, model {f.model}, seed {f.seed}, recycle "
                              f"{recycle}) fails the screen or is named in the omit file; name another one")
        return {"kind": "model", "chosen_by": "config", "condition": f.condition, "model": f.model, "seed": f.seed,
                "recycle": int(recycle)}
    cname = conf.control
    recycle = int(verdicts[verdicts.condition == cname].recycle.max())
    group = kept[(kept.condition == cname) & (kept.recycle == recycle)].sort_values(["model", "seed"])
    if group.empty:
        raise ScreenError(f"no model of the control {cname!r} at recycle {recycle} passes the screen, so there is no "
                          f"representative model; name one under [frame]")
    lab, P, model_idx = _parsed_sites(conf, step1, cname)
    pos = {(int(m), int(s), int(r)): i for i, (m, s, r) in enumerate(zip(lab.model, lab.seed, lab.recycle))}
    rows = [pos[(int(m), int(s), recycle)] for m, s in zip(group.model, group.seed)]
    us = step1["us"]
    meas_t = np.array([u in _residues(conf.measure) for u in us])
    masks = {"fit_t": np.array([u in _residues(conf.fit) for u in us]), "meas_t": meas_t,
             "frame_compare": np.ones((len(conf.model_chains), len(us)), bool) & meas_t[None, :]}
    keys = list(zip(group.model.astype(int), group.seed.astype(int)))
    target, info = lowest_spread_pick(np.asarray(P[rows], dtype=np.float64)[:, :, model_idx, :], keys, masks, log,
                                      detail=True)
    return {"kind": "model", "chosen_by": "lowest_spread_after_screen", "condition": cname, "model": keys[target][0],
            "seed": keys[target][1], "recycle": recycle, **info}


def _build(conf: Config, folder: Path, log: Log) -> dict:
    log(f"screen: {conf.name} ({conf.path})")
    settings = conf.screen
    step1 = read_step1(conf)
    dv, us = step1["dv"], np.array(step1["us"])
    n_chains = len(conf.model_chains)
    X = step1.pop("coords").astype(np.float64).reshape(len(dv), n_chains, len(us), 3)
    screened = np.nonzero((dv.recycle >= 1).to_numpy())[0]
    if len(screened) == 0:
        raise ScreenError("step 1's results hold no model after recycle 0; there is nothing to screen")
    present = set(zip(dv.condition, dv.model.astype(int), dv.seed.astype(int)))
    named = {(c, m, s): r for c, m, s, r in settings.omitted}
    absent = sorted(set(named) - present)
    if absent:                                                  # before the slow part (review E2)
        raise ConfigError(f"[screen] omit names prediction(s) {absent} that step 1's results do not hold")

    on = [c for c in CHECKS if getattr(settings, c)]
    if conf.reference is not None:
        run = on
        ref_xyz = frame_sites(conf, list(us)).reshape(n_chains, len(us), 3)
        ref = prepare_reference(split(ref_xyz, us, conf.segments), settings.domain_swapped, run)
        judge = lambda subs: check(subs, ref)
        log(f"  reference: the structure {conf.reference.path}; it passes its own switched-on checks")
        if ref.missing:
            log(f"  segment residues the structure lacks, left out of every model (D-072): {ref.missing}")
        if ref.untested:
            log(f"  own place cannot test (the structure lacks them): {ref.untested}")
    else:
        run = [c for c in on if c not in NEED_REFERENCE]
        judge = lambda subs: check_without_reference(subs, settings.domain_swapped, run)
        ref = None
        log(f"  no structure in the config: running {run} on each model itself; not run: "
            f"{[c for c in on if c in NEED_REFERENCE]} (D-072)")

    rows = []
    for i in screened:
        key = (dv.condition.iloc[i], int(dv.model.iloc[i]), int(dv.seed.iloc[i]))
        v = judge(split(X[i], us, conf.segments))
        reason = "named" if key in named else next((c for c in run if not v[f"{c}_ok"]), "")
        row = {"condition": key[0], "model": key[1], "seed": key[2], "recycle": int(dv.recycle.iloc[i])}
        for c in CHECKS:
            row[c], row[f"{c}_ok"] = v.get(c, np.nan), v.get(f"{c}_ok", np.nan)
        rows.append({**row, "named": key in named, "reason": reason, "kept": reason == ""})
    verdicts = pd.DataFrame(rows)

    if conf.reference is not None:
        reference = {"kind": "reference", "missing_residues": {str(k): v for k, v in ref.missing.items()},
                     "own_place_untested": [list(t) for t in ref.untested]}
    else:
        reference = pick_representative(conf, step1, verdicts, log)
        log(f"  representative model, picked among the models that pass (D-072): {reference}")

    summary = {"checks_run": run, "checks_not_run": [c for c in CHECKS if c not in run], "reference": reference,
               "models": len(verdicts), "models_kept": int(verdicts.kept.sum()),
               "kept": {}, "removed": {}, "removed_by_reason": {},
               "named": [{"condition": c, "model": m, "seed": s, "reason": r} for c, m, s, r in settings.omitted]}
    for cond in conf.conditions:
        g = verdicts[verdicts.condition == cond.name]
        nets = sorted(g.model.unique())
        summary["kept"][cond.name], summary["removed"][cond.name] = {}, {}
        for r in sorted(g.recycle.unique()):
            here = g[g.recycle == r]
            summary["kept"][cond.name][str(r)] = {str(m): int((here.kept & (here.model == m)).sum()) for m in nets}
            summary["removed"][cond.name][str(r)] = {
                str(m): {k: int(n) for k, n in here[~here.kept & (here.model == m)].reason.value_counts().items()}
                for m in nets}
        summary["removed_by_reason"][cond.name] = {k: int(n) for k, n in g[~g.kept].reason.value_counts().items()}
        log(f"  {cond.name}: {int(g.kept.sum())} of {len(g)} models kept"
            + (f"; removed: {summary['removed_by_reason'][cond.name]}" if summary["removed_by_reason"][cond.name] else ""))

    verdicts.to_csv(folder / "verdicts.csv", index=False)
    _write_json(folder / "summary.json", summary)
    _write_json(folder / "manifest.json", {
        "config_text": conf.text, "settings": screen_settings(conf), "reference": reference,
        "omit_file": None if settings.omit is None else file_digest(settings.omit),
        "step1_manifest": {"sha256": step1["manifest_sha256"], "created": step1["manifest"]["created"]},
        "step1_files": {f"aligned/{name}": file_digest(step1["aligned"] / name) for name in STEP1_FILES},
        "verdicts": file_digest(folder / "verdicts.csv"),
        "screen_code_sha256": code_sha256(),
        "code": {"git_commit": _git("rev-parse", "HEAD").strip() or "unknown",
                 "git_dirty": bool(_git("status", "--porcelain").strip())},
        "versions": _versions(), "created": datetime.now(timezone.utc).isoformat(timespec="seconds")})
    log(f"done: {folder}")
    return {"models": len(verdicts), "kept": int(verdicts.kept.sum()), "removed": int((~verdicts.kept).sum()),
            "output_dir": str(conf.output_dir)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="The topology screen: after step 1, checks each model's fold and records which models step 2 "
                    "keeps, in <output>/screen/.")
    parser.add_argument("--config", required=True, metavar="datasets/<dataset>/<dataset>.toml",
                        help="the dataset config (TOML)")
    args = parser.parse_args(argv)
    try:
        summary = build(load_config(args.config))
    except STEP_ERRORS as exc:
        print(f"screen stopped ({type(exc).__name__}): {exc}", file=sys.stderr)
        return 1
    print(f"screen: {summary['kept']} of {summary['models']} models kept")
    return 0


if __name__ == "__main__":
    sys.exit(main())
