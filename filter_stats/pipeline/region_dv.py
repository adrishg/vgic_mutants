#!/usr/bin/env python3
"""Step 2: structural-diversity DVs per region, from step 1's results.

Run:  python -B pipeline/region_dv.py --config datasets/<dataset>/<dataset>.toml

What step 2 does, in order (decision cards in generalization/decisions/):
 1. Reads the dataset config (`config.py`), including its optional [regions] table: each region is
    a list of UniProt ranges inside the measure set (D-060 items 1-2).
 2. Checks that step 1's results (<output>/aligned/) were built from this config apart from
    [regions], and from the same UniProt FASTA and structure; otherwise it stops (D-060 item 7).
 3. Recomputes step 1's own DVs from aligned/coords.npy and checks that they agree with
    aligned/dv.csv, so the coordinates belong to the table.
 4. For every region and model, on step 1's fitted coordinates with step 1's chain labelling (D-059,
    D-016 item 3): spread; and when the structure resolves at least half of the region's sites, also
    spread_ref, deviation and offset over the resolved sites (D-060 items 3-4). Every DV that uses a
    group's average model (spread, spread_ref, offset) uses the average of the model's own condition ×
    network × recycle group: the five networks are different predictors, and each prediction's spread
    measures how much one network's repeats vary (D-070). Step 1's own columns in dv.csv are recomputed
    the same way; deviation, mapping and core_fit_rmsd are step 1's. The DVs are step 1's
    own function applied to the region's sites. Step 1's two identities are checked per region, and
    a region equal to the measure set must reproduce step 1's columns (D-060 item 5).
 5. Writes <output>/regions/: dv.csv (step 1's columns and rows, then the region columns),
    region_report.json, manifest.json (written last) and step2.log (D-060 item 6). Step 1's files
    are never changed.
With a [screen] table in the config (the topology screen, D-063, screen.py), step 2 first checks
that the screen's results in <output>/screen/ are up to date (same settings, same step 1 files),
and then uses the kept models only: at every recycle except 0, each condition × recycle group's
chain labelling and fit are redone on its kept models (step 1's labelling of a model depends on the
other models in its group, D-055), and step 1's own columns are recomputed for them. dv.csv then
holds recycle 0's models and the kept models.

The contract the tests check is in generalization/notes/regions-plan.md ("Interface tested").
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tomllib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from build_dv import (STEP1_ERRORS, STEP2_FILES, IdentityError, Log, _git, _ranges_of, _residues, _versions,
                      _write_json, check_identities, file_digest, group_dvs, label_group, number_structure,
                      structure_frame)
from config import REGION_DVS, Config, load_config, resolved, screen_settings

HERE = Path(__file__).resolve().parent                  # pipeline/
RESOLVED_SHARE = 0.5          # D-060 item 4: DVs against the structure need at least half of the region's sites
RTOL = 1e-9                 # recomputed values are the same arithmetic on the same numbers
STEP1_GROUP = ("condition", "recycle")                       # step 1's groups
GROUP = ("condition", "model", "recycle")                    # step 2's groups: one network each (D-070)


class Step1ResultsError(Exception):
    """Step 1's results are missing or damaged, or were built from another config, FASTA, structure or
    archive."""


SCREEN_CODE_FILES = ("screen.py", "build_dv.py", "config.py", "region_dv.py")   # what the screen's verdicts depend on
PARSE_ATOL = 1e-4             # Å: refitting from parse/ must give step 1's coordinates (allows another machine's rounding)


class ScreenResultsError(Exception):
    """The config has a [screen], but the screen's results are missing, damaged or out of date."""


# step 1's errors (ConfigError, ParseError, ArchiveError, NumberingError, IdentityError) can also stop step 2
STEP2_ERRORS = (Step1ResultsError, ScreenResultsError, *STEP1_ERRORS)


def percent(k: int, n: int) -> str:
    """k of n as a percentage rounded down to one decimal, so a share under half never reads 50.0% (D-062)."""
    return f"{(1000 * k // n) / 10:.1f}%"


def code_fingerprint() -> str:
    """SHA-256 of region_dv.py, build_dv.py and config.py, in that order."""
    return hashlib.sha256(b"".join((HERE / f).read_bytes() for f in ("region_dv.py", "build_dv.py", "config.py"))
                          ).hexdigest()


# =========================================================================== step 1's results

def _settings(text: str) -> dict:
    """A config's settings that step 1 uses: without [regions], [screen], [segments] and [statistics], as parsed TOML
    (so how a path was typed does not matter)."""
    data = tomllib.loads(text)
    for key in ("regions", "screen", "segments", "statistics"):
        data.pop(key, None)
    return data


def read_step1(conf: Config) -> dict:
    """Step 1's aligned/ files, after checking they belong to this config (D-060 item 7)."""
    aligned = Path(conf.output_dir) / "aligned"
    files = {name: aligned / name for name in ("manifest.json", "dv.csv", "sites.csv", "coords.npy")}
    missing = [name for name, p in files.items() if not p.exists()]
    if missing:
        raise Step1ResultsError(f"no complete step 1 results in {aligned} (missing {missing}); run step 1 first")
    manifest_bytes = files["manifest.json"].read_bytes()
    try:                                                          # D-061 item 5: a damaged record is a message
        manifest = json.loads(manifest_bytes)
        old = _settings(manifest["config_text"])
        fasta, reference = manifest["inputs"]["uniprot_fasta"], manifest["inputs"]["reference"]
        sizes = {name: d["size"] for name, d in manifest["archives"].items()}
        manifest["created"]
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        raise Step1ResultsError(f"step 1's record {files['manifest.json']} is damaged or incomplete "
                                f"({type(exc).__name__}: {exc}); re-run step 1") from None
    new = _settings(conf.text)
    if old != new:
        differ = sorted(k for k in old.keys() | new.keys() if old.get(k) != new.get(k))
        raise Step1ResultsError(f"step 1's results in {aligned} were built from a different config (it differs in "
                                f"{differ}, not counting [regions], [screen], [segments] and [statistics]); re-run step 1 first")
    if fasta != file_digest(conf.uniprot_fasta):
        raise Step1ResultsError(f"the UniProt FASTA {conf.uniprot_fasta} changed after step 1 ran; re-run step 1 first")
    if conf.reference is not None and reference != file_digest(conf.reference.path):
        raise Step1ResultsError(f"the structure {conf.reference.path} changed after step 1 ran; re-run step 1 first")
    for cond in conf.conditions:                                  # D-061 item 2: by size; checksums are step 1's job
        size = cond.archive.stat().st_size
        if size != sizes.get(cond.name):
            raise Step1ResultsError(f"the archive {cond.archive} of condition {cond.name!r} changed after step 1 ran "
                                    f"(now {size} bytes; step 1 read {sizes.get(cond.name)}); re-run step 1 first")

    text_columns = {c: str for c in ["condition", *conf.factors]}
    dv = pd.read_csv(files["dv.csv"], dtype=text_columns, keep_default_na=False,
                     float_precision="round_trip")                # D-061 item 4: step 1's numbers exactly as written
    sites = pd.read_csv(files["sites.csv"])
    coords = np.load(files["coords.npy"])
    us = sorted(_residues(conf.fit) | _residues(conf.measure))
    n_chains = len(conf.model_chains)
    expected_sites = [(j, u) for j in range(n_chains) for u in us]
    if list(zip(sites.chain_index, sites.uniprot_resnum)) != expected_sites:
        raise Step1ResultsError(f"{files['sites.csv']} does not list the config's fit and measure residues for "
                                f"{n_chains} chain(s); re-run step 1")
    if coords.shape != (len(dv), 3 * len(sites)):
        raise Step1ResultsError(f"{files['coords.npy']} has shape {coords.shape}, but dv.csv has {len(dv)} rows and "
                                f"sites.csv {len(sites)} sites; re-run step 1")
    return {"manifest": manifest, "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(), "dv": dv,
            "sites": sites, "coords": coords, "us": us, "aligned": aligned}


def frame_sites(conf: Config, us: list) -> np.ndarray | None:
    """The structure's coordinates per site (sites × 3, sites.csv order), rebuilt with step 1's own
    numbering; None without a structure."""
    if conf.reference is None:
        return None
    structure = number_structure(conf, _residues(conf.fit) | _residues(conf.measure), _residues(conf.fit))
    return structure_frame(structure, conf.reference.chains, us).reshape(-1, 3)


# =========================================================================== the DVs

def dvs_over(dv: pd.DataFrame, X: np.ndarray, cols: np.ndarray, ref_cols: np.ndarray | None,
             frame: np.ndarray | None, by=GROUP) -> dict:
    """step 1's DVs (group_dvs) per group (condition × network × recycle, or `by`) over the sites `cols`;
    with `ref_cols` also the DVs against the structure over those sites. X: rows × sites × 3 (float64)."""
    names = ["spread"] + (list(REGION_DVS[1:]) if ref_cols is not None else [])
    out = {name: np.empty(len(dv)) for name in names}
    for _, idx in dv.groupby(list(by)).indices.items():
        values = group_dvs(X[idx], cols, ref_cols, frame)
        for name in names:
            out[name][idx] = values[name]
    return out


def _disagree(a, b) -> bool:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return not np.allclose(a, b, rtol=RTOL, atol=1e-12)


def check_step1_columns(dv: pd.DataFrame, X: np.ndarray, sites: pd.DataFrame, frame) -> None:
    """Recompute step 1's DVs from coords.npy; raise IdentityError unless they agree with dv.csv."""
    meas = sites.in_measure.to_numpy(bool)
    ref_cols = meas & sites.resolved.to_numpy(bool) if frame is not None else None
    again = dvs_over(dv, X, meas, ref_cols, frame, by=STEP1_GROUP)
    wrong = [c for c, v in again.items() if _disagree(v, dv[c])]
    if wrong:
        raise IdentityError(f"step 1's columns {wrong} do not follow from aligned/coords.npy: either the coordinates "
                            f"do not belong to aligned/dv.csv, or this code numbers the structure differently from "
                            f"the code that ran step 1; re-run step 1")


# =========================================================================== the screen's results

def read_screen(conf: Config, step1: dict) -> dict:
    """The screen's verdicts, after checking they are up to date: made by the current screen.py, with this
    config's [screen] and [segments] and the same omit file, from these step 1 files, and covering exactly
    step 1's models at [screen] first_recycle and above. Returns {"kept": set of (condition, model, seed,
    recycle), "first_recycle": the lowest recycle judged, "removed": n,
    "screened": n, "reference": the screen's reference, "manifest_sha256"}."""
    folder = Path(conf.output_dir) / "screen"
    man_path, verdicts_path = folder / "manifest.json", folder / "verdicts.csv"
    if not (man_path.exists() and verdicts_path.exists()):
        raise ScreenResultsError(f"the config has a [screen], but there are no screen results in {folder}; "
                                 f"run the screen (pipeline/screen.py) after step 1")
    raw = man_path.read_bytes()
    try:
        manifest = json.loads(raw)
        old_settings, old_files, old_omit = manifest["settings"], manifest["step1_files"], manifest["omit_file"]
        old_verdicts, old_code, reference = manifest["verdicts"], manifest["screen_code_sha256"], manifest["reference"]
    except (ValueError, KeyError, TypeError) as exc:
        raise ScreenResultsError(f"the screen's record {man_path} is damaged ({type(exc).__name__}: {exc}); "
                                 f"re-run the screen") from None
    if old_settings != json.loads(json.dumps(screen_settings(conf))):
        raise ScreenResultsError(f"the screen's results in {folder} were made with other [screen] or [segments] "
                                 f"settings, or another omit list, than the config's; re-run the screen")
    omit = conf.screen.omit
    if old_omit != (None if omit is None else file_digest(omit)):
        raise ScreenResultsError(f"the screen's omit file {omit} changed after the screen ran; re-run the screen")
    now = {f"aligned/{name}": file_digest(step1["aligned"] / name) for name in ("dv.csv", "sites.csv", "coords.npy")}
    if old_files != now:
        raise ScreenResultsError(f"the screen's results in {folder} were made from other step 1 results "
                                 f"(step 1 was rebuilt since); re-run the screen")
    if old_code != hashlib.sha256(b"".join((HERE / f).read_bytes() for f in SCREEN_CODE_FILES)).hexdigest():
        raise ScreenResultsError(f"the screen's code ({', '.join(SCREEN_CODE_FILES)}) changed after the screen ran; "
                                 f"re-run the screen")
    if old_verdicts != file_digest(verdicts_path):
        raise ScreenResultsError(f"{verdicts_path} changed after the screen wrote it; re-run the screen")
    v = pd.read_csv(verdicts_path, dtype={"condition": str}, keep_default_na=False)
    keys = list(zip(v.condition, v.model, v.seed, v.recycle))
    dv = step1["dv"]
    first = conf.screen.first_recycle
    screened_rows = dv.recycle >= first
    if sorted(keys) != sorted(zip(dv.condition[screened_rows], dv.model[screened_rows],
                                  dv.seed[screened_rows], dv.recycle[screened_rows])):
        raise ScreenResultsError(f"{verdicts_path} does not list exactly step 1's models at recycle {first} and "
                                 f"above; re-run the screen")
    is_kept = v.kept.astype(str).str.lower() == "true"
    return {"kept": {k for k, ok in zip(keys, is_kept) if ok}, "first_recycle": int(first),
            "removed": int((~is_kept).sum()), "screened": len(v),
            "reference": reference, "manifest_sha256": hashlib.sha256(raw).hexdigest()}


def refit_kept(conf: Config, step1: dict, kept: set, log: Log, reference: dict | None = None) -> tuple:
    """dv rows and coordinates with the screen applied. Each condition × recycle group's kept models are
    labelled and fitted again on their own (label_group, as step 1 does), and step 1's own columns are
    recomputed; a group with nothing removed keeps step 1's result, and recycles below [screen]
    first_recycle are not screened. Without a structure, the frame is the representative model the screen
    picked among the kept models (`reference`, D-072); if it is not step 1's frame model, every group, the
    unscreened recycles included, is fitted onto it. Before that, one group per condition is refitted from
    parse/ onto step 1's own frame and must give step 1's result exactly. Returns (dv, X: rows × sites × 3)."""
    dv, sites, us = step1["dv"], step1["sites"], step1["us"]
    X1 = step1["coords"].astype(np.float64).reshape(len(dv), len(sites), 3)
    out = Path(conf.output_dir)
    n_chains = len(conf.model_chains)
    fit, measure = _residues(conf.fit), _residues(conf.measure)
    fit_t = np.array([u in fit for u in us])
    meas_t = np.array([u in measure for u in us])
    model_idx = sites[sites.chain_index == 0].model_resnum.to_numpy() - 1
    meas_cols = np.tile(meas_t, n_chains)

    parsed = {}
    for cond in conf.conditions:
        try:
            lab = pd.read_csv(out / "parse" / f"{cond.name}.csv")
            P = np.load(out / "parse" / f"{cond.name}.npy", mmap_mode="r")
        except (OSError, ValueError) as exc:
            raise Step1ResultsError(f"step 1's parse/ files of {cond.name!r} are missing or damaged ({exc}); "
                                    f"re-run step 1") from None
        pos = {(int(m), int(s), int(r)): i for i, (m, s, r) in enumerate(zip(lab.model, lab.seed, lab.recycle))}
        parsed[cond.name] = (pos, P)

    def raw(cname, keys):
        pos, P = parsed[cname]
        return np.asarray(P[[pos[k] for k in keys]], dtype=np.float64)[:, :, model_idx, :]

    if conf.reference is not None:
        frame1 = frame_sites(conf, us).reshape(n_chains, len(us), 3)
        compare = sites.resolved.to_numpy(bool).reshape(n_chains, len(us)) & meas_t[None, :]
        ref_cols, frame_flat = compare.reshape(-1), frame1.reshape(-1, 3)
        frame2, new_frame = frame1, False
    else:
        compare, ref_cols, frame_flat = np.ones((n_chains, len(us)), bool) & meas_t[None, :], None, None
        f = step1["manifest"]["frame"]
        frame1 = raw(f["condition"], [(f["model"], f["seed"], f["recycle"])])[0]
        r = reference or f
        new_frame = (r["condition"], r["model"], r["seed"], r["recycle"]) != (f["condition"], f["model"], f["seed"], f["recycle"])
        frame2 = raw(r["condition"], [(r["model"], r["seed"], r["recycle"])])[0] if new_frame else frame1
        if new_frame:
            log(f"  every model is fitted onto the representative model the screen picked: {r}")

    def refit(cname, rows_in_dv, frame):
        keys = [(int(dv.model.iloc[i]), int(dv.seed.iloc[i]), int(dv.recycle.iloc[i])) for i in rows_in_dv]
        L = label_group(raw(cname, keys), frame, fit_t, meas_t, compare)
        return L, L.coords.reshape(len(keys), -1, 3).astype(np.float32).astype(np.float64)

    first_recycle = conf.screen.first_recycle if conf.screen is not None else 1
    keep_rows, new_X, verified = [], [], set()
    columns = {c: dv[c].to_numpy(dtype=object).copy() for c in ("mapping", "core_fit_rmsd", *REGION_DVS) if c in dv}
    for (cname, recycle), idx in dv.groupby(["condition", "recycle"], sort=False).indices.items():
        idx = np.sort(idx)
        if cname not in verified:            # once per condition: parse/ must give step 1's own result (review B4)
            _, Y = refit(cname, idx, frame1)
            if not np.allclose(Y, X1[idx], rtol=0, atol=PARSE_ATOL):
                raise Step1ResultsError(f"refitting {cname} recycle {recycle} from parse/ does not reproduce "
                                        f"aligned/coords.npy: the parse cache does not belong to step 1's results; "
                                        f"re-run step 1")
            verified.add(cname)
        if recycle < first_recycle:                                       # below [screen] first_recycle: unjudged
            mine = idx
        else:
            mine = idx[[(dv.condition.iloc[i], dv.model.iloc[i], dv.seed.iloc[i], dv.recycle.iloc[i]) in kept
                        for i in idx]]
        if len(mine) == len(idx) and not new_frame:                       # nothing removed: step 1's own result
            keep_rows.extend(idx)
            new_X.append(X1[idx])
            continue
        if len(mine) == 0:
            log(f"  {cname} recycle {recycle}: no model kept")
            continue
        L, Y = refit(cname, mine, frame2)
        values = group_dvs(Y, meas_cols, ref_cols, frame_flat)
        for k, i in enumerate(mine):
            columns["mapping"][i] = ",".join(conf.model_chains[c] for c in L.perms[k])
            columns["core_fit_rmsd"][i] = float(L.core_fit_rmsd[k])
            for name, v in values.items():
                columns[name][i] = float(v[k])
        keep_rows.extend(mine)
        new_X.append(Y)
        log(f"  {cname} recycle {recycle}: {len(mine)} of {len(idx)} models kept; labelling and fit redone")
    new_dv = dv.copy()
    for c, v in columns.items():
        new_dv[c] = v if c == "mapping" else v.astype(float)
    order = np.argsort(np.array(keep_rows), kind="stable")
    rows = np.array(keep_rows)[order]
    X = np.concatenate(new_X)[order]
    return new_dv.iloc[rows].reset_index(drop=True), X


# =========================================================================== step 2

def build(conf: Config) -> dict:
    """Run step 2 for a loaded config; returns {"regions": [names], "with_reference": {name: bool},
    "output_dir": str}. Raises one of STEP2_ERRORS on a problem."""
    log = Log()
    folder = Path(conf.output_dir) / "regions"
    folder.mkdir(parents=True, exist_ok=True)
    results = [name for name in STEP2_FILES if name != "step2.log"]      # the manifest first: never beside new files
    for name in results:
        (folder / name).unlink(missing_ok=True)
    try:
        return _build(conf, folder, log)
    except BaseException as exc:                                 # errors and interrupts alike (D-062 item 1)
        for name in results:                                      # D-061 item 6: nothing from a failed run
            (folder / name).unlink(missing_ok=True)
        log.lines.append(f"step 2 stopped ({type(exc).__name__}): {exc}")
        raise
    finally:
        (folder / "step2.log").write_text("\n".join(log.lines) + "\n")


def _build(conf: Config, folder: Path, log: Log) -> dict:
    log(f"step 2: {conf.name} ({conf.path})")
    step1 = read_step1(conf)
    log(f"  step 1 results: {step1['aligned']} (built {step1['manifest']['created']})")
    dv, sites, us = step1["dv"], step1["sites"], step1["us"]
    X = step1["coords"].astype(np.float64).reshape(len(dv), len(sites), 3)
    frame = frame_sites(conf, us)
    check_step1_columns(dv, X, sites, frame)
    log("  step 1's columns agree with aligned/coords.npy")
    screen = None
    if conf.screen is not None:
        screen = read_screen(conf, step1)
        log(f"  the screen's results are up to date: {screen['removed']} of {screen['screened']} screened models "
            f"removed (recycles below {screen['first_recycle']} are not screened)")
        dv, X = refit_kept(conf, step1, screen["kept"], log, screen["reference"])
    step1.pop("coords")
    meas_site = sites.in_measure.to_numpy(bool)                  # step 1's own columns, by network (D-070)
    dv = dv.copy()
    for name, v in dvs_over(dv, X, meas_site, meas_site & sites.resolved.to_numpy(bool) if frame is not None else None,
                            frame).items():
        if name != "deviation":                                  # a model's own distance: no group involved
            dv[name] = v
    log("  groups: condition × network × recycle (D-070); step 1's spread columns recomputed per network")

    fit, measure, treated = _residues(conf.fit), _residues(conf.measure), _residues(conf.treated)
    uniprot_of_site = sites.uniprot_resnum.to_numpy()
    resolved_site = sites.resolved.to_numpy(bool)
    n_chains = len(conf.model_chains)
    columns, report_regions, with_reference = {}, [], {}
    if not conf.regions:
        log("  no [regions] in the config: nothing to add; regions/dv.csv holds step 1's columns only")
    for name, ranges in conf.regions.items():
        residues = _residues(ranges)
        cols = np.isin(uniprot_of_site, sorted(residues))
        n_sites = int(cols.sum())
        n_resolved = int((cols & resolved_site).sum()) if frame is not None else None
        if frame is None:
            reason = "no structure in the config: spread only"
        elif n_resolved < RESOLVED_SHARE * n_sites:
            reason = (f"the structure resolves {n_resolved} of its {n_sites} sites ({percent(n_resolved, n_sites)}), "
                      f"under half: spread only")
        else:
            reason = None
        with_reference[name] = reason is None
        values = dvs_over(dv, X, cols, cols & resolved_site if reason is None else None, frame)
        for dv_name, v in values.items():
            columns[f"{dv_name}_{name}"] = v
        report_regions.append({
            "name": name, "ranges": [list(r) for r in ranges], "n_residues": len(residues), "n_sites": n_sites,
            "n_resolved_sites": n_resolved, "resolved_fraction": None if n_resolved is None else n_resolved / n_sites,
            "with_reference": reason is None, "reason": reason,
            "n_treated": len(residues & treated), "treated_fraction": len(residues & treated) / len(residues),
            "n_fit": len(residues & fit), "equals_measure": residues == measure})
        log(f"  {name}: UniProt {[list(r) for r in ranges]}, {len(residues)} residues per chain × {n_chains}; "
            + (f"{percent(n_resolved, n_sites)} of its sites resolved in the structure; " if n_resolved is not None else "")
            + f"{percent(len(residues & treated), len(residues))} treated; "
            + (", ".join(values) if reason is None else reason))

    out_dv = pd.concat([dv, pd.DataFrame(columns, index=dv.index)], axis=1)
    _check_regions(conf, out_dv, X, sites)
    if conf.regions:
        log("identity checks passed per region: mean(deviation²) = mean(spread_ref²) + offset²; "
            "pairwise = 2·mean(spread²)")
        last = out_dv[out_dv.recycle == out_dv.recycle.max()]
        means = last.groupby("condition", sort=False)[[f"spread_{n}" for n in conf.regions]].mean()
        log(f"mean spread at recycle {int(out_dv.recycle.max())} (Å):\n" + means.round(3).to_string())

    names = list(conf.regions)
    overlaps = [{"regions": [a, b], "n_residues": len(_residues(conf.regions[a]) & _residues(conf.regions[b]))}
                for i, a in enumerate(names) for b in names[i + 1:]]
    in_regions = set().union(*(_residues(r) for r in conf.regions.values())) if names else set()
    report = {"regions": report_regions, "measured_in_no_region": _ranges_of(measure - in_regions),
              "overlaps": [o for o in overlaps if o["n_residues"]],
              "step1": {"manifest_sha256": step1["manifest_sha256"], "created": step1["manifest"]["created"]},
              "groups": "condition × network × recycle (D-070)",
              # the frame every model was fitted onto (review 2, S3): the structure, step 1's frame model, or
              # without a structure the screen's representative model (D-072)
              "frame": ({"kind": "reference"} if conf.reference is not None
                        else (screen["reference"] if screen is not None else step1["manifest"]["frame"]))}
    if report["measured_in_no_region"] and names:
        log(f"  measured residues in no region: {report['measured_in_no_region']}")

    out_dv.to_csv(folder / "dv.csv", index=False)
    _write_json(folder / "region_report.json", report)
    _write_json(folder / "manifest.json", {
        "config_text": conf.text, "regions": resolved(conf)["regions"],
        "step1_manifest": {"sha256": step1["manifest_sha256"], "created": step1["manifest"]["created"]},
        # what later steps check regions/ against: unchanged while step 1 reuses aligned/ (D-062 item 4)
        "step1_files": {name: file_digest(step1["aligned"] / name) for name in ("dv.csv", "sites.csv", "coords.npy")},
        "screen_manifest": None if screen is None else {"sha256": screen["manifest_sha256"]},
        "code": {"sha256": code_fingerprint(), "git_commit": _git("rev-parse", "HEAD").strip() or "unknown",
                 "git_dirty": bool(_git("status", "--porcelain").strip())},
        "versions": _versions(), "created": datetime.now(timezone.utc).isoformat(timespec="seconds")})
    log(f"done: {folder}")
    return {"regions": names, "with_reference": with_reference, "output_dir": str(conf.output_dir)}


def _check_regions(conf: Config, out_dv: pd.DataFrame, X: np.ndarray, sites: pd.DataFrame) -> None:
    """Step 1's two identities per region (check_identities on the region's sites), and a region equal
    to the measure set must reproduce step 1's columns (D-060 item 5)."""
    measure = _residues(conf.measure)
    keys = out_dv[list(GROUP)]
    for name, ranges in conf.regions.items():
        present = [c for c in REGION_DVS if f"{c}_{name}" in out_dv.columns]
        region_dv = pd.concat([keys, out_dv[[f"{c}_{name}" for c in present]].set_axis(present, axis=1)], axis=1)
        region_sites = sites.assign(in_measure=sites.uniprot_resnum.isin(_residues(ranges)))
        try:
            check_identities(region_dv, X, region_sites, by=GROUP)
        except IdentityError as exc:
            raise IdentityError(f"region {name}: {exc}") from None
        if _residues(ranges) == measure:
            wrong = [c for c in present if c in out_dv.columns and _disagree(out_dv[f"{c}_{name}"], out_dv[c])]
            if wrong:
                raise IdentityError(f"region {name} equals the measure set, but its {wrong} differ from step 1's")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Step 2 of the pipeline: from step 1's results, computes the structural-diversity DVs of "
                    "each region in the config's [regions] table (spread; with a structure also spread_ref, "
                    "deviation, offset) and writes them, with a region report, to <output>/regions/.")
    parser.add_argument("--config", required=True, metavar="datasets/<dataset>/<dataset>.toml",
                        help="the dataset config (TOML)")
    args = parser.parse_args(argv)
    try:
        summary = build(load_config(args.config))
    except STEP2_ERRORS as exc:
        print(f"step 2 stopped ({type(exc).__name__}): {exc}", file=sys.stderr)
        return 1
    print(f"regions: {summary['regions']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
