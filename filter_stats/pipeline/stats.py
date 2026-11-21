#!/usr/bin/env python3
"""The statistics step: each condition against the control, for every measure and recycle.

Run:  python -B pipeline/stats.py --config datasets/<dataset>/<dataset>.toml
It reads step 2's results (<output>/regions/) and writes <output>/stats/. Decisions: D-007 (any number of
conditions against one control, two-sided), D-066 (one family), D-069 (REML), D-070 (spread against each
network's own average), D-074 (spread and offset, not deviation), D-075 (no small-group correction),
D-076 (the 5% split between measures, Dunnett within each), D-077 (one mixed model, one value per group),
D-078 (no minimum group size).

The measures: spread and offset of each region in the config (spread only where step 2 has no offset,
for example without a structure); without [regions], step 1's spread and offset. Recycle 0 is not analysed.

1. One value per group (one condition, one network, one recycle). Spread: the root-mean-square of the
   group's model spreads, that is how scattered the group's models are around their average model; a group
   of one model has no spread. Offset: the distance from the group's average model to the structure.
   Every group is used wherever its value exists.
2. One mixed model per measure, fitted by REML:
       log(value) ~ condition × recycle + network + (1 | condition:network)
   - condition × recycle: each condition's own mean at each recycle (no shape assumed across recycles);
   - network: every condition was run through the same networks, so each comparison is made within
     network, as in a paired design;
   - (1 | condition:network): a pair's groups at different recycles come from the same seeds, so they are
     not independent; the random intercept counts each pair once. The pairs are the units, and the
     networks are each condition's repeats.
3. Each condition other than the control is compared with the control at each recycle: the difference of
   their log means, reported as a ratio of geometric means. Its standard error comes from the fixed
   effects' covariance (X' V^-1 X)^-1 with the fitted variances.
4. Degrees of freedom: pairs − (conditions + networks − 1), i.e. (conditions − 1) × (networks − 1) for
   complete data. With fewer than 2 conditions or 2 networks there is nothing to test against.
5. Multiple testing: the 5% is split evenly between the measures (Bonferroni); within a measure, its
   comparisons get one two-sided Dunnett cutoff from the multivariate t with their correlations. Adjusted
   p-values and simultaneous 95% intervals use the same cutoff.

Caveats, reported with the results: the log values are about normal with one noise level per measure; all
recycles of a pair are equally alike; groups count equally whatever their size; each average model is
estimated from its own group's models, so spread reads slightly low and offset slightly high in small
groups (not corrected).

Writes <output>/stats/: comparisons.csv (one row per measure × condition × recycle), models.json (per
measure: groups, pairs, degrees of freedom, variances, cutoff; and the counts a methods text needs),
manifest.json (written last) and stats.log.
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
from scipy.optimize import brentq
from scipy.stats import multivariate_t
from scipy.stats import t as t_dist
from statsmodels.regression.mixed_linear_model import MixedLM

from build_dv import STEP1_ERRORS, Log, _git, _versions, _write_json, file_digest
from config import Config, load_config

HERE = Path(__file__).resolve().parent
ALPHA = 0.05
SEED = 0                          # the multivariate t's integration is randomised; a fixed seed repeats it
STATS_FILES = ("manifest.json", "comparisons.csv", "models.json", "stats.log")
MODEL = "log(value) ~ condition * recycle + network + (1 | condition:network)"


class StatsInputError(Exception):
    """Step 2's results are missing, or out of date with the config, step 1 or the screen."""


class StatsDesignError(Exception):
    """The design leaves nothing to test against (fewer than 2 conditions or networks, no degrees of freedom)."""


STEP_ERRORS = (StatsInputError, StatsDesignError, *STEP1_ERRORS)


# =========================================================================== one value per group

OFFSET_TOLERANCE = 1e-6           # an offset is the group's one distance to the structure: the same for every model


def _group_offset(x: np.ndarray, measure: str, cond: str, net: str, rec: int) -> float:
    """The group's offset, which is one value shared by all its models (D-074).

    The group's mean is returned rather than an arbitrary member's, and a group whose members disagree by
    more than OFFSET_TOLERANCE is an error rather than a silent pick of whichever model groupby ordered
    first: that disagreement means the grouping key no longer matches the region the offset was measured
    over, and the whole measure would be wrong in a way nothing downstream could see."""
    finite = x[np.isfinite(x)]
    if len(finite) == 0:
        return float("nan")
    spread = float(finite.max() - finite.min())
    if spread > OFFSET_TOLERANCE:
        raise StatsInputError(
            f"{measure}: the group (condition {cond!r}, network {net!r}, recycle {rec}) holds {len(finite)} models "
            f"whose offsets differ by {spread:.3g} A, more than the {OFFSET_TOLERANCE:g} A an offset is allowed to "
            f"vary within a group. An offset is one distance per group, so this means the grouping key and the "
            f"region the offset was measured over no longer agree.")
    return float(finite.mean())


def group_values(dv: pd.DataFrame, measures: list, first_recycle: int = 1) -> pd.DataFrame:
    """One row per group (condition, network, recycle >= first_recycle): pair, number of models and each
    measure's value. first_recycle must match the screen's, or the statistics analyse a different set
    from the one the screen passed (D-063)."""
    d = dv[dv.recycle >= first_recycle]
    rows = []
    for (cond, net, rec), g in d.groupby(["condition", "model", "recycle"], sort=True):
        row = {"condition": cond, "network": net, "recycle": int(rec), "pair": f"{cond}:{net}", "n": len(g)}
        for m in measures:
            x = g[m].to_numpy(dtype=float)
            if m.startswith("offset"):
                row[m] = _group_offset(x, m, cond, net, rec)                # one value for the whole group
            else:
                row[m] = float(np.sqrt(np.mean(x ** 2))) if len(g) >= 2 else np.nan
        rows.append(row)
    return pd.DataFrame(rows, columns=["condition", "network", "recycle", "pair", "n", *measures])


# =========================================================================== the mixed model

def fe_covariance(X: np.ndarray, pairs: np.ndarray, var_pair: float, var_noise: float) -> np.ndarray:
    """(X' V^-1 X)^-1, V = var_noise·I + var_pair within each pair: the fixed effects' covariance given the
    fitted variances. For a pair of k rows, V^-1 = (I − var_pair / (var_noise + k·var_pair) · 11') / var_noise."""
    A = np.zeros((X.shape[1], X.shape[1]))
    for p in pd.unique(pairs):
        Xi = X[pairs == p]
        s = Xi.sum(axis=0)
        A += (Xi.T @ Xi - var_pair / (var_noise + len(Xi) * var_pair) * np.outer(s, s)) / var_noise
    return np.linalg.inv(A)


def fit_measure(groups: pd.DataFrame, measure: str, control: str, conditions=None) -> dict:
    """Fit one measure (module docstring, items 2-4) and compare each condition with the control."""
    usable = np.isfinite(groups[measure]) & (groups[measure] > 0)
    # Why these are counted apart: a group of one model has no spread (D-078 sets no minimum group size, so
    # nothing upstream stops this), and dropping it here removes a whole condition x network x recycle cell
    # from the design rather than just changing a value. With a screen in front of the statistics that is
    # how a filter quietly edits the experiment, so the count is logged and written into models.json.
    too_small = int((~usable & (groups.n.to_numpy() < 2)).sum())
    d = groups[usable].reset_index(drop=True)
    present = list(dict.fromkeys(d.condition))
    order = [c for c in (conditions or present) if c in present]
    nets = sorted(d.network.unique())
    n_pairs, n_conds = d.pair.nunique(), len(present)
    df = n_pairs - (n_conds + len(nets) - 1)
    if control not in present or n_conds < 2 or len(nets) < 2 or df < 1:
        raise StatsDesignError(
            f"{measure}: nothing to test against: the comparisons need the control and at least one other condition, "
            f"each run through at least 2 networks ({n_conds} condition(s), {len(nets)} network(s), {df} degrees of "
            f"freedom; control {control!r} {'present' if control in present else 'missing'})")

    cells = sorted(set(zip(d.condition, d.recycle)))                  # one mean per condition × recycle ...
    col = {cell: j for j, cell in enumerate(cells)}
    X = np.zeros((len(d), len(cells) + len(nets) - 1))
    X[np.arange(len(d)), [col[k] for k in zip(d.condition, d.recycle)]] = 1.0
    for k, net in enumerate(nets[1:]):                                # ... plus each network's level
        X[:, len(cells) + k] = (d.network.to_numpy() == net)
    y = np.log(d[measure].to_numpy(dtype=float))
    pairs = d.pair.to_numpy()
    fit = MixedLM(y, X, groups=pairs).fit(reml=True)
    beta = np.asarray(fit.fe_params)
    var_pair, var_noise = float(np.asarray(fit.cov_re)[0, 0]), float(fit.scale)
    V = fe_covariance(X, pairs, var_pair, var_noise)

    rows, L = [], []
    for cond in [c for c in order if c != control]:
        for rec in sorted(d.recycle.unique()):
            if (cond, rec) in col and (control, rec) in col:
                ell = np.zeros(X.shape[1])
                ell[col[(cond, rec)]], ell[col[(control, rec)]] = 1.0, -1.0
                est, se = float(ell @ beta), float(np.sqrt(ell @ V @ ell))
                rows.append({"condition": cond, "recycle": int(rec), "log_ratio": est, "ratio": np.exp(est),
                             "SE": se, "t": est / se})
                L.append(ell)
            else:                                                     # a condition with no group here: not computable
                rows.append({"condition": cond, "recycle": int(rec), "log_ratio": np.nan, "ratio": np.nan,
                             "SE": np.nan, "t": np.nan})
    comparisons = pd.DataFrame(rows, columns=["condition", "recycle", "log_ratio", "ratio", "SE", "t"])
    if L:
        Lk = np.array(L)
        C = Lk @ V @ Lk.T
        sd = np.sqrt(np.diag(C))
        R = C / np.outer(sd, sd)
    else:
        R = np.zeros((0, 0))
    info = {"groups": int(len(d)), "groups_without_value": int(len(groups) - len(d)),
            "groups_dropped_too_small": too_small,
            "groups_dropped_other": int(len(groups) - len(d) - too_small), "pairs": int(n_pairs),
            "conditions": int(n_conds), "networks": int(len(nets)), "var_pair": var_pair, "var_noise": var_noise,
            "converged": bool(fit.converged)}
    return {"comparisons": comparisons, "R": R, "df": int(df), "info": info}


# =========================================================================== the cutoff

def block_prob(c: float, R: np.ndarray, df: float) -> float:
    """P(every |T_i| <= c) for T multivariate t with correlations R and df degrees of freedom."""
    k = len(R)
    if k == 0:
        return 1.0
    if k == 1:
        return float(1 - 2 * t_dist.sf(c, df))
    mvt = multivariate_t(loc=np.zeros(k), shape=R, df=df, allow_singular=True)
    return float(mvt.cdf(np.full(k, c), lower_limit=np.full(k, -c), random_state=SEED))


def cutoff(R: np.ndarray, df: float, alpha: float) -> float:
    """Dunnett's two-sided cutoff: the c at which P(every |T| <= c) = 1 − alpha."""
    if len(R) == 0:
        return np.nan
    target, hi = 1 - alpha, 8.0
    while block_prob(hi, R, df) < target:                             # few degrees of freedom need a high cutoff
        hi *= 2
        if hi > 1e6:
            raise StatsDesignError(f"no cutoff below {hi:g} reaches {target}")
    return float(brentq(lambda c: block_prob(c, R, df) - target, 0.1, hi, xtol=1e-6))


def adjusted_p(t_value: float, R: np.ndarray, df: float, n_measures: int) -> float:
    """The comparison's p-value adjusted within its measure (Dunnett) and between measures (× n_measures)."""
    if not np.isfinite(t_value):
        return np.nan
    return float(min(1.0, max(0.0, n_measures * (1.0 - block_prob(abs(t_value), R, df)))))


# =========================================================================== step 2's results

def _settings(text: str) -> dict:
    """A config's settings without [statistics], as parsed TOML."""
    data = tomllib.loads(text)
    data.pop("statistics", None)
    return data


def read_regions(conf: Config) -> pd.DataFrame:
    """Step 2's dv.csv, after checking it is up to date: made from this config (apart from [statistics]),
    from step 1's current files and, with a screen, from the screen's current results."""
    out = Path(conf.output_dir)
    man_path, dv_path = out / "regions" / "manifest.json", out / "regions" / "dv.csv"
    if not (man_path.exists() and dv_path.exists()):
        raise StatsInputError(f"no step 2 results in {out / 'regions'}; run step 2 (region_dv.py) first")
    try:
        m = json.loads(man_path.read_text())
        old_text, old_files, old_screen = m["config_text"], m["step1_files"], m.get("screen_manifest")
    except (ValueError, KeyError, TypeError) as exc:
        raise StatsInputError(f"step 2's record {man_path} is damaged ({type(exc).__name__}: {exc}); re-run step 2") from None
    if _settings(old_text) != _settings(conf.text):
        raise StatsInputError(f"step 2's results in {out / 'regions'} were made from a different config (apart from "
                              f"[statistics]); re-run step 2 (and the screen or step 1 if you changed theirs)")
    now = {name: file_digest(out / "aligned" / name) for name in ("dv.csv", "sites.csv", "coords.npy")}
    if old_files != now:
        raise StatsInputError("step 2's results were made from other step 1 results (step 1 was rebuilt since); "
                              "re-run step 2")
    if conf.screen is not None:
        sm = out / "screen" / "manifest.json"
        if not sm.exists() or old_screen != {"sha256": hashlib.sha256(sm.read_bytes()).hexdigest()}:
            raise StatsInputError("step 2's results were made from other screen results; re-run step 2")
    return pd.read_csv(dv_path, dtype={"condition": str}, keep_default_na=False, float_precision="round_trip")


def measures(conf: Config, dv: pd.DataFrame) -> list:
    """Spread and offset of each region (offset where step 2 has it); without [regions], step 1's."""
    if conf.regions:
        names = []
        for r in conf.regions:
            names += [m for m in (f"spread_{r}", f"offset_{r}") if m in dv.columns]
        return names
    return [m for m in ("spread", "offset") if m in dv.columns]


def counts(conf: Config, dv: pd.DataFrame, groups: pd.DataFrame) -> dict:
    """The totals a methods text needs (D-076 item 3)."""
    d = dv[dv.recycle >= 1]
    screen = None
    verdicts = Path(conf.output_dir) / "screen" / "verdicts.csv"
    if conf.screen is not None and verdicts.exists():
        v = pd.read_csv(verdicts, usecols=["kept"])
        kept = v.kept.astype(str).str.lower().isin(["true", "1"])
        screen = {"screened": int(len(v)), "kept": int(kept.sum()), "removed": int((~kept).sum())}
    return {"models": int(len(d)), "predictions": int(d.prediction.nunique()), "groups": int(len(groups)),
            "conditions": int(groups.condition.nunique()), "networks": int(groups.network.nunique()),
            "pairs": int(groups.pair.nunique()), "recycles": sorted(int(r) for r in groups.recycle.unique()),
            "smallest_group": int(groups.n.min()) if len(groups) else 0, "screen": screen}


# =========================================================================== the step

def build(conf: Config) -> dict:
    """Run the statistics step for a loaded config; returns a summary dict. Raises one of STEP_ERRORS."""
    log = Log()
    folder = Path(conf.output_dir) / "stats"
    folder.mkdir(parents=True, exist_ok=True)
    results = [name for name in STATS_FILES if name != "stats.log"]
    for name in results:
        (folder / name).unlink(missing_ok=True)
    try:
        return _build(conf, folder, log)
    except BaseException as exc:
        for name in results:
            (folder / name).unlink(missing_ok=True)
        log.lines.append(f"statistics stopped ({type(exc).__name__}): {exc}")
        raise
    finally:
        (folder / "stats.log").write_text("\n".join(log.lines) + "\n")


def _build(conf: Config, folder: Path, log: Log) -> dict:
    log(f"statistics: {conf.name} ({conf.path}); model {MODEL}")
    dv = read_regions(conf)
    names = measures(conf, dv)
    if not names:
        raise StatsInputError("step 2's results hold no spread or offset to test")
    first = conf.screen.first_recycle if conf.screen is not None else 1
    groups = group_values(dv, names, first)
    alpha_m = ALPHA / len(names)
    order = [c.name for c in conf.conditions]
    log(f"  {len(names)} measures: the 5% is split evenly, {100 * alpha_m:.3g}% each (Bonferroni), "
        f"then Dunnett within each measure")
    tables, per_measure = [], {}
    for m in names:
        f = fit_measure(groups, m, conf.control, order)
        crit = cutoff(f["R"], f["df"], alpha_m)
        c = f["comparisons"].copy()
        c.insert(0, "measure", m)
        c["lo95"] = np.exp(c.log_ratio - crit * c.SE)
        c["hi95"] = np.exp(c.log_ratio + crit * c.SE)
        c["p_adj"] = [adjusted_p(t, f["R"], f["df"], len(names)) for t in c.t]
        c["significant"] = c.p_adj < ALPHA
        c["cutoff"], c["df"] = crit, f["df"]
        tables.append(c)
        i = f["info"]
        per_measure[m] = {**i, "df": f["df"], "cutoff": crit, "comparisons": int(c.t.notna().sum()),
                          "significant": int(c.significant.sum())}
        log(f"  {m}: {i['groups']} groups ({i['groups_without_value']} without a value, of which "
            f"{i['groups_dropped_too_small']} held fewer than 2 models), {i['pairs']} pairs, "
            f"df {f['df']}; variances: pair {i['var_pair']:.4f}, noise {i['var_noise']:.4f}; cutoff |t| > {crit:.3f}; "
            f"{per_measure[m]['significant']} of {per_measure[m]['comparisons']} significant"
            + ("" if i["converged"] else "; WARNING: the fit did not converge"))
    out = pd.concat(tables, ignore_index=True)
    first = ["measure", "condition", "recycle", "ratio", "lo95", "hi95", "p_adj", "significant", "log_ratio", "SE", "t",
             "cutoff", "df"]
    out = out[first]
    out.to_csv(folder / "comparisons.csv", index=False)
    totals = counts(conf, dv, groups)
    _write_json(folder / "models.json", {"model": MODEL, "control": conf.control, "alpha": ALPHA,
                                         "n_measures": len(names), "alpha_per_measure": alpha_m,
                                         "counts": totals, "measures": per_measure})
    regions = Path(conf.output_dir) / "regions"
    _write_json(folder / "manifest.json", {
        "config_text": conf.text,
        "regions_files": {f"regions/{n}": file_digest(regions / n) for n in ("dv.csv", "manifest.json")},
        "code": {"sha256": hashlib.sha256((HERE / "stats.py").read_bytes()).hexdigest(),
                 "git_commit": _git("rev-parse", "HEAD").strip() or "unknown",
                 "git_dirty": bool(_git("status", "--porcelain").strip())},
        "versions": _versions(), "created": datetime.now(timezone.utc).isoformat(timespec="seconds")})
    n_sig = int(out.significant.sum())
    log(f"  counts: {totals['models']} models from {totals['predictions']} predictions in {totals['groups']} groups "
        f"({totals['pairs']} condition × network pairs)" + (f"; screen: {totals['screen']['removed']} of "
                                                              f"{totals['screen']['screened']} models removed"
                                                              if totals["screen"] else ""))
    log(f"  {n_sig} of {int(out.t.notna().sum())} comparisons significant at the 5% family-wise level")
    log(f"done: {folder}")
    return {"measures": names, "comparisons": len(out), "significant": n_sig, "output_dir": str(conf.output_dir)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="The statistics step: from step 2's results, each condition against the control for every "
                    "measure and recycle (one mixed model per measure; Bonferroni between measures, Dunnett within), "
                    "in <output>/stats/.")
    parser.add_argument("--config", required=True, metavar="datasets/<dataset>/<dataset>.toml",
                        help="the dataset config (TOML)")
    args = parser.parse_args(argv)
    try:
        summary = build(load_config(args.config))
    except STEP_ERRORS as exc:
        print(f"statistics stopped ({type(exc).__name__}): {exc}", file=sys.stderr)
        return 1
    print(f"statistics: {summary['significant']} of {summary['comparisons']} comparisons significant")
    return 0


if __name__ == "__main__":
    sys.exit(main())
