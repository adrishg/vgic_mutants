#!/usr/bin/env python3
"""The quality-control chain, in the one order that is defensible.

    1. all-OK-3 convergence          earliest recycle from which every later successive-recycle
                                     transition stays under 3 A; that recycle is the trajectory's
                                     representative, and it and everything after it are retained
    2. selectivity-filter ring       trajectory rejected if any retained recycle exceeds 20 A
                                     (Kv2.1 G377, Cav1.2 EEEE, Nav1.5 DEKA)
    3. pore-VSD interface            trajectory level: one bad recycle rejects the whole trajectory.
                                     Two forms. The published Kv2.1 one thresholds S6's end against
                                     the residues before S1 at 27 A, read off the distance table.
                                     The generalized one thresholds the lower S6 against its own
                                     S4-S5 linker at 12 A, measured from the models, and works for
                                     all three channels because that contact exists in all three.
    4. topology screen               Brandon Harris's reference-free checks (filter_stats/pipeline/
                                     screen.py), judging each surviving trajectory's representative

Channel differences - the ring residues, the repeat segments, whether an interface cutoff exists -
live in channels.py, which also records what was measured about how far each step transfers.

Why this order and not the cheap-screen-first one that looks more attractive. Two reasons, both
measured in 26-11/filterComparison/FILTER_ORDER_FINDINGS.md section 9, neither of them about cost:

  - Convergence is defined on a contiguous recycle series. It compares each recycle with the previous
    available one and takes the earliest recycle whose every later transition stays under 3 A. A
    topology screen that removes recycles out of the middle of a trajectory changes what "the previous
    available recycle" means, so running it first does not reorder two filters, it silently redefines
    one of them, and "converged" stops meaning what the published number means.
  - The pore-VSD interface QC is trajectory level and has to see the whole trajectory. Screening first
    deletes exactly the recycles the interface QC exists to catch, and it then finds almost nothing to
    reject: 8 L403A trajectories instead of 107, and 11 WT instead of 297. The ensemble is not better,
    the QC has been made blind. The headline metric duly rises in the WT control as well as in the
    mutant, which is permissiveness, not signal.

The topology screen judges only the representative, not all eleven recycles. That is about eleven
times cheaper, but the reason is that the representative is the unit the rest of the analysis reports.
It judges recycle 0 too: the convergence step picks r0 as the representative for over half of some
cohorts, so a screen that skipped it would let exactly those models through unjudged.

Steps 2 and 3 are absolute thresholds applied per model and then taken to the trajectory, so none of
them is a statistic over the set another has thinned and their order among themselves changes no
survivor and no count - only which name a doomed trajectory is filed under in chain_verdicts(). The
convention is the generalized S6-linker check before the Kv2.1-only published one, so that one
structural failure carries one name in all three channels. docs/FILTER_CHAIN_ORDER.md has the whole
argument and the numbers behind it.

Only the published form is a cutoff anyone has published. The generalized form was calibrated on
fourteen experimental structures (channels.py, item 3) and is the only one Nav1.5 and Cav1.2 can
take, so run_chain() offers both and the sweep in scripts/filtering/hive compares them.

A note on the distance columns: their names carry the residue name, so position 405 is LEU405 in the
wild type and ALA405 in L403A. Match them by position with ring_columns(), never by residue name.

Run:  python filter_stats/filter_chain.py --distances <cohort>_AF2_distances.csv \
          --convergence <cohort>_all_ok_models.csv --model-dir <pdb tree> --output <chain>.csv
"""
from __future__ import annotations

import argparse
import itertools
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from shared.dataset_selection import apply_kv21_interface_qc           # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from channels import Channel, channel as get_channel                   # noqa: E402

G377_RING_MAX_A = 20.0            # the published selectivity-filter ring cutoff, trajectory level
INTERFACE_MAX_A = 27.0            # the published pore-VSD interface cutoff, trajectory level
S6_LINKER_MAX_A = 12.0            # the generalized pore-VSD cutoff, trajectory level
TOPOLOGY_CHECKS = ("wiring", "box_corners", "s6_on_own_s5")            # the checks needing no reference
RECYCLE_SUFFIX = re.compile(r"\.r\d+\.pdb$")
MODEL_META = re.compile(r"model_(\d+)_seed_(\d+)\.r(\d+)\.pdb$")

# The helix boundaries are 8SD3's own HELIX records; the models are numbered two higher. S6 breaks at
# L403/P404, which is both the mutation site and the published S6 kink.
SEGMENTS_EXPERIMENTAL = {
    "S1": (181, 205), "S2": (224, 246), "S3": (251, 277), "S4": (289, 309),
    "linker": (310, 323), "S5": (324, 350), "pore_helix": (359, 372),
    "filter": (373, 378), "S6_upper": (383, 403), "S6_lower": (404, 428),
}
MODEL_NUMBERING_OFFSET = 2
SEGMENTS_MODEL = {k: (a + MODEL_NUMBERING_OFFSET, b + MODEL_NUMBERING_OFFSET)
                  for k, (a, b) in SEGMENTS_EXPERIMENTAL.items()}


# =========================================================================== small shared helpers

def basenames(column: pd.Series) -> pd.Series:
    """The file name alone, since the cohort CSVs and the PDB trees disagree about the leading path."""
    return column.astype(str).str.rsplit("/", n=1).str[-1]


def trajectories(column: pd.Series) -> pd.Series:
    """The trajectory key: one (model, seed) pair, with the recycle suffix dropped."""
    return basenames(column).str.replace(RECYCLE_SUFFIX, ".pdb", regex=True)


def ring_columns(frame: pd.DataFrame, position: int) -> list:
    """The six inter-chain CA-CA columns of a homotetramer's ring at `position`, matched by position.

    The column names carry the residue name, which differs between cohorts at any mutated position
    (LEU405 in the wild type, ALA405 in L403A), so a hardcoded name silently returns nothing for one
    of the two conditions."""
    pattern = re.compile(rf"^CA_CA_([A-D])_[A-Z]{{3}}{position}_CA-([A-D])_[A-Z]{{3}}{position}_CA$")
    found = [c for c in frame.columns if (m := pattern.match(c)) and m.group(1) != m.group(2)]
    if len(found) != 6:
        raise KeyError(f"expected the 6 inter-chain columns of the position-{position} ring, found "
                       f"{len(found)}: {found}")
    return sorted(found)


def repeat_ring_columns(frame: pd.DataFrame, residues) -> list:
    """The six columns of a pseudo-tetramer's filter ring (Nav1.5 DEKA, Cav1.2 EEEE).

    One chain, four non-identical filter residues, so the columns name a pair of residues and the
    pair can be written either way round."""
    found = []
    for first, second in itertools.combinations(residues, 2):
        options = [f"CA_{first}_CA-{second}_CA", f"CA_{second}_CA-{first}_CA"]
        here = [c for c in options if c in frame.columns]
        if not here:
            raise KeyError(f"the filter ring needs a column for {first}-{second}; tried {options}")
        found.append(here[0])
    return found


def filter_ring_columns(frame: pd.DataFrame, chan: Channel) -> list:
    """The channel's six selectivity-filter ring columns, whichever shape its models have."""
    if chan.ring_residues:
        return repeat_ring_columns(frame, chan.ring_residues)
    return ring_columns(frame, chan.ring_position)


def _truth(column: pd.Series) -> pd.Series:
    """The convergence CSV writes its flags as the strings 'True'/'False'."""
    return column.astype(str).str.strip().str.lower().eq("true")


# =========================================================================== step 1: convergence

def apply_convergence(distances: pd.DataFrame, convergence: pd.DataFrame) -> tuple:
    """Keep the retained models of each trajectory and mark its representative (all-OK-3).

    The convergence itself is computed by the published script, which walks the contiguous recycle
    series; this reads its verdict. `all_ok` marks the retained models (the representative and every
    recycle after it) and `earliest_converged_selected` marks the representative."""
    for column in ("pdb_file", "all_ok", "earliest_converged_selected"):
        if column not in convergence.columns:
            raise KeyError(f"the convergence table needs the {column!r} column")
    names = basenames(convergence["pdb_file"])
    retained = set(names[_truth(convergence["all_ok"])])
    representatives = set(names[_truth(convergence["earliest_converged_selected"])])

    out = distances.copy()
    out.attrs.update(distances.attrs)
    out["model_file"] = basenames(out["pdb_file"])
    out["trajectory"] = trajectories(out["pdb_file"])
    out["is_representative"] = out["model_file"].isin(representatives)
    kept = out.loc[out["model_file"].isin(retained)].copy()
    kept.attrs.update(out.attrs)
    return kept, {"trajectories_in": out["trajectory"].nunique(),
                  "models_in": len(out),
                  "trajectories_out": kept["trajectory"].nunique(),
                  "models_out": len(kept)}


# =========================================================================== step 2: the G377 ring

def apply_filter_ring(frame: pd.DataFrame, chan: Channel, *, threshold_A: float = None) -> tuple:
    """Reject a whole trajectory if any retained recycle opens the selectivity-filter ring too far.

    Kv2.1's G377 ring, Cav1.2's EEEE and Nav1.5's DEKA are the same six inter-subunit Ca-Ca
    distances. The step is kept for all three because it costs nothing, but it is a tripwire: it has
    never rejected anything in the cohorts measured so far (channels.py, item 1)."""
    threshold_A = chan.ring_max_A if threshold_A is None else threshold_A
    columns = filter_ring_columns(frame, chan)
    widest = frame[columns].apply(pd.to_numeric, errors="coerce").max(axis=1)
    failed = set(frame.loc[widest.gt(threshold_A), "trajectory"])
    kept = frame.loc[~frame["trajectory"].isin(failed)].copy()
    kept.attrs.update(frame.attrs)
    return kept, {"trajectories_rejected": len(failed),
                  "models_rejected": int(len(frame) - len(kept)),
                  "threshold_A": float(threshold_A),
                  "widest_retained_A": float(widest.loc[kept.index].max()) if len(kept) else float("nan")}


# =========================================================================== step 3b: S6 against its own linker

def s6_linker_distances(path: Path, chan: Channel, *, segments=None) -> list:
    """How far each subunit's lower S6 sits from its own S4-S5 linker, closest Ca-Ca pair.

    This is the pore-VSD interface check written so that all three channels can take it. The
    published Kv2.1 version thresholds the end of S6 against the residues before S1, and that
    pairing is a contact in a domain-swapped homotetramer and nowhere else - in Cav1.2 and Nav1.5
    those elements are 25 to 50 A apart (channels.py, item 3). The gate's own S4-S5 linker is the
    joint it hangs from in every voltage-gated channel, and across fourteen experimental structures
    the closest approach runs 4.8 to 9.2 A with the subunit's own linker always nearest. A model
    whose S6 has let go of that linker has come apart at the one place the pore is held."""
    coords, numbers = _load_ca(path, chains=chan.chains)
    chain_index = {c: i for i, c in enumerate(chan.chains)}
    number_index = {int(n): i for i, n in enumerate(numbers)}
    out = []
    subunits = chan.subunits() if segments is None else tuple(
        zip([c for c, _ in chan.subunits()], segments))
    for chain, segments in subunits:
        rows = coords[chain_index[chain]]
        a = np.array([rows[number_index[n]] for n in range(*_inclusive(segments["S6_lower"]))
                      if n in number_index])
        b = np.array([rows[number_index[n]] for n in range(*_inclusive(segments["linker"]))
                      if n in number_index])
        a = a[np.isfinite(a).all(axis=1)] if len(a) else a
        b = b[np.isfinite(b).all(axis=1)] if len(b) else b
        if not len(a) or not len(b):
            out.append(float("nan"))
            continue
        out.append(float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min()))
    return out


def _inclusive(span) -> tuple:
    """A (first, last) residue span as the half-open range that walks it."""
    return int(span[0]), int(span[1]) + 1


def measure_interface(paths: dict, chan: Channel) -> pd.DataFrame:
    """One row per model: each subunit's S6-to-linker distance and the widest of them.

    Walking a whole cohort's PDB tree is the expensive part, so the measure is written out once and
    every threshold is then a question asked of the table rather than of the models again."""
    rows = []
    for name, path in sorted(paths.items()):
        try:
            per = s6_linker_distances(Path(path), chan)
            error = ""
        except Exception as exc:                                       # noqa: BLE001
            per, error = [float("nan")] * len(chan.subunits()), f"{type(exc).__name__}: {exc}"
        row = {"pdb_file": name, "error": error}
        for i, value in enumerate(per, start=1):
            row[f"s6_linker_{i}_A"] = value
        finite = [v for v in per if np.isfinite(v)]
        row["s6_linker_max_A"] = max(finite) if finite else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def apply_s6_linker_qc(frame: pd.DataFrame, measures: pd.DataFrame, *, threshold_A: float = S6_LINKER_MAX_A,
                       trajectory_level: bool = True) -> tuple:
    """Reject trajectories whose S6 has let go of its linker in any retained recycle."""
    widest = frame["model_file"].map(
        measures.set_index("pdb_file")["s6_linker_max_A"].to_dict())
    unmeasured = int(widest.isna().sum())
    over = widest.gt(threshold_A)
    if trajectory_level:
        rejected = frame["trajectory"].isin(set(frame.loc[over, "trajectory"]))
    else:
        rejected = over
    kept = frame.loc[~rejected].copy()
    kept.attrs.update(frame.attrs)
    return kept, {"trajectories_rejected": frame["trajectory"].nunique() - kept["trajectory"].nunique(),
                  "models_rejected": int(rejected.sum()),
                  "models_unmeasured": unmeasured,
                  "threshold_A": float(threshold_A),
                  "widest_retained_A": float(widest.loc[kept.index].max()) if len(kept) else float("nan")}


# =========================================================================== step 4: the topology screen

def _load_ca(path: Path, chains=("A", "B", "C", "D")) -> tuple:  # noqa: D401
    """CA coordinates as (chain, residue, 3) over the union of residue numbers, plus that numbering."""
    per = {c: {} for c in chains}
    with open(path) as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            if line[21] not in per or line[16] not in (" ", "A"):      # one altloc only
                continue
            per[line[21]][int(line[22:26])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    numbers = sorted(set().union(*(set(d) for d in per.values())))
    index = {n: i for i, n in enumerate(numbers)}
    coords = np.full((len(chains), len(numbers), 3), np.nan)
    for ci, chain in enumerate(chains):
        for number, xyz in per[chain].items():
            coords[ci, index[number]] = xyz
    return coords, np.array(numbers)


def _screen_module():
    """Brandon's screen.py, imported for its check functions only.

    His pipeline wants one ColabFold .tar.gz per condition and ours are loose PDB trees on Quobyte,
    so the integration borrows the geometry and leaves the plumbing: the checks take plain coordinate
    arrays and need none of his config or archive handling."""
    here = Path(__file__).resolve().parent
    pipeline = here / "pipeline"
    if str(pipeline) not in sys.path:
        sys.path.insert(0, str(pipeline))
    try:
        import screen                                                  # noqa: PLC0415
    except ImportError:
        # His file imports three modules at module scope that the checks never use, and one of them
        # pulls in gemmi, which the ColabFold environment on Hive does not have. Stand-ins for those
        # names, and only then, so where the real pipeline can be imported it is the one used.
        sys.path.insert(0, str(here / "screen_shim"))
        import screen                                                  # noqa: PLC0415
    return screen


def judge_topology(paths: dict, *, chan: Channel = None, segments=None, chains=None,
                   domain_swapped: bool = True, checks=TOPOLOGY_CHECKS) -> pd.DataFrame:
    """Run the reference-free checks over the given {model file: path} models, one row each.

    A model the geometry cannot judge is recorded with an error and `topology_ok` False rather than
    ending the run; on a cohort of several thousand that difference is hours of work."""
    screen = _screen_module()
    if segments is None:
        segments = (list(chan.segments) if chan is not None and len(chan.segments) > 1
                    else (chan.segments[0] if chan is not None else SEGMENTS_MODEL))
    if chains is None:
        chains = chan.chains if chan is not None else ("A", "B", "C", "D")
    rows = []
    for name in sorted(paths):
        meta = MODEL_META.search(name)
        row = {"model_file": name,
               "model": int(meta.group(1)) if meta else np.nan,
               "seed": int(meta.group(2)) if meta else np.nan,
               "recycle": int(meta.group(3)) if meta else np.nan}
        try:
            coords, numbering = _load_ca(Path(paths[name]), chains)
            verdict = screen.check_without_reference(
                screen.split(coords, numbering, segments), domain_swapped, checks)
        except Exception as exc:                                       # one bad model, not a dead job
            rows.append({**row, **{c: np.nan for c in checks}, **{f"{c}_ok": False for c in checks},
                         "topology_ok": False, "error": f"{type(exc).__name__}: {exc}"})
            continue
        oks = {f"{c}_ok": bool(verdict[f"{c}_ok"]) for c in checks}
        rows.append({**row, **{c: verdict[c] for c in checks}, **oks,
                     "topology_ok": all(oks.values()), "error": ""})
    return pd.DataFrame(rows)


def collect_models(model_dir: Path, wanted: set) -> dict:
    """{model file: path} for the wanted models, found anywhere under the cohort's shard directories."""
    found = {}
    for path in Path(model_dir).rglob("*.pdb"):
        if path.name in wanted and path.name not in found:
            found[path.name] = path
    missing = sorted(wanted - set(found))
    if missing:
        raise FileNotFoundError(f"{len(missing)} model(s) not found under {model_dir}, "
                                f"for example {missing[:3]}")
    return found


def apply_topology(representatives: pd.DataFrame, *, chan: Channel = None, model_dir: Path = None,
                   verdicts: pd.DataFrame = None, **kwargs) -> tuple:
    """Screen the representatives, from a PDB tree or from a verdict table already computed."""
    if (model_dir is None) == (verdicts is None):
        raise ValueError("give either model_dir, to judge the representatives now, or verdicts, to "
                         "read a screen that has already been run")
    wanted = set(representatives["model_file"])
    if verdicts is None:
        verdicts = judge_topology(collect_models(model_dir, wanted), chan=chan, **kwargs)
    else:
        verdicts = verdicts.copy()
        key = "model_file" if "model_file" in verdicts.columns else "pdb_file"
        verdicts["model_file"] = basenames(verdicts[key])
        verdicts = verdicts.loc[verdicts["model_file"].isin(wanted)]
    passing = set(verdicts.loc[_truth(verdicts["topology_ok"]), "model_file"])
    unjudged = sorted(wanted - set(verdicts["model_file"]))
    if unjudged:
        raise KeyError(f"{len(unjudged)} representative(s) have no topology verdict, for example "
                       f"{unjudged[:3]}")
    kept = representatives.loc[representatives["model_file"].isin(passing)].copy()
    kept.attrs.update(representatives.attrs)
    return kept, verdicts, {"representatives_in": len(representatives),
                            "representatives_out": len(kept),
                            "could_not_judge": int((verdicts["error"].astype(str) != "").sum())}


# =========================================================================== the per-prediction verdict table

# The chain's steps in the order they run, each with the sentence a reader who has never seen this
# code should be able to act on. The keys are the `removed_by` values in the verdict table.
FILTER_STEPS = (
    ("convergence",
     "No recycle of this prediction reached the all-OK-3 convergence criterion, so the trajectory "
     "never settled."),
    ("selectivity_filter_ring",
     "The selectivity-filter ring opened past the cutoff in at least one retained recycle."),
    ("pore_vsd_interface_s6_linker",
     "The lower S6 let go of its own S4-S5 linker in at least one retained recycle, so the gate is "
     "no longer held by the joint it hangs from."),
    ("pore_vsd_interface_published",
     "The published Kv2.1 pore-VSD interface opened past the cutoff in at least one retained "
     "recycle. This check exists for Kv2.1 only."),
    ("topology",
     "The representative model failed at least one reference-free topology check."),
)
STEP_REASON = dict(FILTER_STEPS)

PASS, FAIL, NOT_APPLICABLE = "pass", "fail", "not_applicable"


def _per_trajectory_max(frame: pd.DataFrame, columns: list) -> dict:
    """The widest of `columns` over each trajectory's rows, as {trajectory: value}."""
    if not len(frame):
        return {}
    widest = frame[columns].apply(pd.to_numeric, errors="coerce").max(axis=1)
    return widest.groupby(frame["trajectory"]).max().to_dict()


def _published_interface_columns() -> list:
    """The 24 columns the published Kv2.1 pore-VSD check reads, in shared/dataset_selection.py."""
    return [f"CA_CA_{chain}_{first}_CA-{chain}_{second}_CA"
            for chain in "ABCD"
            for first in ("LYS429", "GLU425", "LYS422")
            for second in ("ASN181", "VAL184")]


def chain_verdicts(distances: pd.DataFrame, convergence: pd.DataFrame, *, channel: str,
                   topology_verdicts: pd.DataFrame, interface_measures=None,
                   dataset: str = "", condition: str = "",
                   ring_threshold_A: float = None, interface_threshold_A: float = INTERFACE_MAX_A,
                   s6_linker_threshold_A: float = None) -> pd.DataFrame:
    """One row per prediction: what every step measured, whether it passed, and if not, which step
    took it and why.

    This is the table a reader outside the project gets. `run_chain` answers "what survives"; this
    answers "what happened to each prediction", which is the same walk through the same steps in
    the same order, with nothing discarded along the way.

    A prediction here is one (model, seed) trajectory - the eleven recycles of one ColabFold
    prediction - because every step but the topology screen judges the trajectory as a whole.

    Measurements are taken over each trajectory's retained recycles, which is the set the filters
    themselves see. A trajectory that never converged has no retained recycles, so its later
    measurements are blank rather than zero."""
    chan = get_channel(channel)
    distances = distances.copy()
    distances.attrs["dataset_label"] = f"{dataset} {condition}".strip() or chan.name

    all_trajectories = sorted(set(trajectories(distances["pdb_file"])))
    recycles_in = trajectories(distances["pdb_file"]).value_counts().to_dict()

    kept, _ = apply_convergence(distances, convergence)
    converged = set(kept["trajectory"])
    retained_counts = kept["trajectory"].value_counts().to_dict()
    representative = (kept.loc[kept["is_representative"]]
                      .set_index("trajectory")["model_file"].to_dict())

    reasons, status_step = {}, {}
    for trajectory in set(all_trajectories) - converged:
        status_step[trajectory] = "convergence"

    ring_max = _per_trajectory_max(kept, filter_ring_columns(kept, chan))
    published_columns = _published_interface_columns()
    has_published = chan.interface_qc and all(c in kept.columns for c in published_columns)
    published_max = _per_trajectory_max(kept, published_columns) if has_published else {}
    if interface_measures is not None:
        by_model = interface_measures.set_index("pdb_file")["s6_linker_max_A"].to_dict()
        s6_values = kept["model_file"].map(by_model)
        s6_max = s6_values.groupby(kept["trajectory"]).max().to_dict()
    else:
        s6_max = {}

    def take(new, step):
        nonlocal kept
        for trajectory in set(kept["trajectory"]) - set(new["trajectory"]):
            status_step[trajectory] = step
        kept = new

    take(apply_filter_ring(kept, chan, threshold_A=ring_threshold_A)[0], "selectivity_filter_ring")
    if interface_measures is not None:
        threshold = chan.s6_linker_max_A if s6_linker_threshold_A is None else s6_linker_threshold_A
        take(apply_s6_linker_qc(kept, interface_measures, threshold_A=threshold)[0],
             "pore_vsd_interface_s6_linker")
    if chan.interface_qc:
        take(apply_kv21_interface_qc(kept, threshold_A=interface_threshold_A, trajectory_level=True),
             "pore_vsd_interface_published")

    verdicts = topology_verdicts.copy()
    key = "model_file" if "model_file" in verdicts.columns else "pdb_file"
    verdicts["model_file"] = basenames(verdicts[key])
    topology_ok = dict(zip(verdicts["model_file"], _truth(verdicts["topology_ok"])))
    check_ok = {c: dict(zip(verdicts["model_file"], _truth(verdicts[f"{c}_ok"])))
                for c in TOPOLOGY_CHECKS if f"{c}_ok" in verdicts.columns}
    unjudged = sorted(n for t, n in representative.items() if t in set(kept["trajectory"])
                      and n not in topology_ok)
    if unjudged:
        raise KeyError(f"{len(unjudged)} representative(s) reach the topology step with no verdict, for "
                       f"example {unjudged[:3]}; re-run the judge step over this cohort")
    for trajectory in sorted(set(kept["trajectory"])):
        name = representative.get(trajectory)
        if name is not None and not topology_ok.get(name, False):
            status_step[trajectory] = "topology"

    ring_threshold = chan.ring_max_A if ring_threshold_A is None else ring_threshold_A
    s6_threshold = chan.s6_linker_max_A if s6_linker_threshold_A is None else s6_linker_threshold_A

    def verdict(value, threshold, ran):
        if not ran or value is None or not np.isfinite(value):
            return NOT_APPLICABLE
        return PASS if value <= threshold else FAIL

    rows = []
    for trajectory in all_trajectories:
        meta = MODEL_META.search(trajectory.replace(".pdb", ".r0.pdb"))
        name = representative.get(trajectory)
        rep_meta = MODEL_META.search(name) if name else None
        step = status_step.get(trajectory)
        ring = ring_max.get(trajectory, float("nan"))
        s6 = s6_max.get(trajectory, float("nan"))
        pub = published_max.get(trajectory, float("nan"))
        rows.append({
            "dataset": dataset,
            "condition": condition,
            "channel": chan.name,
            "prediction": trajectory,
            "model": int(meta.group(1)) if meta else np.nan,
            "seed": int(meta.group(2)) if meta else np.nan,
            "recycles_in_distance_table": int(recycles_in.get(trajectory, 0)),
            "recycles_retained_after_convergence": int(retained_counts.get(trajectory, 0)),
            "representative_model_file": name or "",
            "representative_recycle": int(rep_meta.group(3)) if rep_meta else np.nan,
            "convergence_result": PASS if trajectory in converged else FAIL,
            "selectivity_filter_ring_max_A": ring,
            "selectivity_filter_ring_cutoff_A": float(ring_threshold),
            "selectivity_filter_ring_result": verdict(ring, ring_threshold, trajectory in converged),
            "pore_vsd_s6_linker_max_A": s6,
            "pore_vsd_s6_linker_cutoff_A": float(s6_threshold) if s6_max else np.nan,
            "pore_vsd_s6_linker_result": verdict(s6, s6_threshold, bool(s6_max)),
            "pore_vsd_published_max_A": pub,
            "pore_vsd_published_cutoff_A": float(interface_threshold_A) if has_published else np.nan,
            "pore_vsd_published_result": verdict(pub, interface_threshold_A, has_published),
            **{f"topology_{c}_result": (PASS if check_ok[c].get(name, False) else FAIL)
               if (name and name in check_ok[c]) else NOT_APPLICABLE
               for c in check_ok},
            "topology_result": (PASS if topology_ok.get(name, False) else FAIL)
                               if (name and name in topology_ok) else NOT_APPLICABLE,
            "status": "removed" if step else "passed",
            "removed_by": step or "",
            "removed_because": STEP_REASON[step] if step else "",
        })
    out = pd.DataFrame(rows)
    return out.sort_values(["model", "seed"], kind="stable").reset_index(drop=True)


# =========================================================================== the chain

@dataclass
class ChainResult:
    """What a cohort looks like after the four steps."""
    representatives: pd.DataFrame      # one row per surviving trajectory: the filtered table to analyse
    retained: pd.DataFrame             # every retained model of those trajectories
    topology: pd.DataFrame             # the per-representative topology verdicts
    counts: dict = field(default_factory=dict)


def run_chain(distances: pd.DataFrame, convergence: pd.DataFrame, *, channel: str = "Kv2.1",
              model_dir=None, topology_verdicts=None, ring_threshold_A: float = None,
              interface_threshold_A: float = INTERFACE_MAX_A, interface_measures=None,
              interface_mode: str = "auto", s6_linker_threshold_A: float = None,
              label: str = None, **topology_kwargs) -> ChainResult:
    """Pass one cohort through convergence, the filter ring, the interface QC and the topology
    screen, in that order and no other (see the module docstring for why the order is not free).

    The third step has two forms. `published` is the 27 A Kv2.1 cutoff read off the distance table,
    which exists for Kv2.1 alone. `s6_linker` is the generalized version, the lower S6 against its
    own S4-S5 linker, which needs `interface_measures` from measure_interface() and works for all
    three channels. `auto` takes the generalized one where the measures are given and falls back to
    the published one for Kv2.1, which keeps the Kv2.1 numbers we have already published intact
    when no measures are passed."""
    chan = get_channel(channel)
    label = label or f"{chan.name} cohort"
    distances = distances.copy()
    distances.attrs["dataset_label"] = label
    counts = {"channel": chan.name}

    kept, counts["convergence"] = apply_convergence(distances, convergence)
    kept, counts["filter_ring"] = apply_filter_ring(kept, chan, threshold_A=ring_threshold_A)

    mode = interface_mode
    if mode == "auto":
        mode = "s6_linker" if interface_measures is not None else ("published" if chan.interface_qc else "none")
    if mode in ("s6_linker", "both"):
        if interface_measures is None:
            raise ValueError(f"interface_mode={interface_mode!r} needs interface_measures from measure_interface()")
        threshold = chan.s6_linker_max_A if s6_linker_threshold_A is None else s6_linker_threshold_A
        kept, counts["s6_linker_qc"] = apply_s6_linker_qc(
            kept, interface_measures, threshold_A=threshold, trajectory_level=True)
    if mode in ("published", "both") and chan.interface_qc:
        before = kept["trajectory"].nunique()
        kept = apply_kv21_interface_qc(kept, threshold_A=interface_threshold_A, trajectory_level=True)
        counts["interface_qc"] = {
            "trajectories_rejected": before - kept["trajectory"].nunique(),
            "models_rejected": int(kept.attrs.get("kv21_interface_qc_rejected_rows", 0)),
            "threshold_A": float(interface_threshold_A)}
    else:
        counts["interface_qc"] = {"trajectories_rejected": 0, "models_rejected": 0,
                                  "not_run": f"interface_mode={mode} for {chan.name}"}
    counts.setdefault("s6_linker_qc", {"trajectories_rejected": 0, "models_rejected": 0,
                                       "not_run": f"interface_mode={mode}"})
    counts["interface_mode"] = mode

    representatives = kept.loc[kept["is_representative"]].copy()
    representatives.attrs.update(kept.attrs)
    representatives, verdicts, counts["topology"] = apply_topology(
        representatives, chan=chan, model_dir=model_dir, verdicts=topology_verdicts, **topology_kwargs)

    surviving = set(representatives["trajectory"])
    retained = kept.loc[kept["trajectory"].isin(surviving)].copy()
    retained.attrs.update(kept.attrs)
    counts["surviving_trajectories"] = len(surviving)
    print(f"{label}: {counts['convergence']['trajectories_in']} trajectories in, "
          f"{counts['convergence']['trajectories_out']} converged, "
          f"-{counts['filter_ring']['trajectories_rejected']} filter ring, "
          f"-{counts['interface_qc']['trajectories_rejected']} interface, "
          f"-{counts['s6_linker_qc']['trajectories_rejected']} S6-linker, "
          f"-{counts['topology']['representatives_in'] - counts['topology']['representatives_out']} topology; "
          f"{len(surviving)} surviving")
    return ChainResult(representatives=representatives, retained=retained, topology=verdicts,
                       counts=counts)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--distances", required=True, type=Path, help="the cohort's AF2 distance CSV")
    parser.add_argument("--convergence", required=True, type=Path,
                        help="the cohort's all_ok_models CSV from the published convergence script")
    parser.add_argument("--model-dir", type=Path, help="the cohort's PDB tree, to screen the representatives")
    parser.add_argument("--topology", type=Path, help="instead of --model-dir: topology verdicts already computed")
    parser.add_argument("--channel", default="Kv2.1", help="Kv2.1, Nav1.5 or Cav1.2")
    parser.add_argument("--label", default=None)
    parser.add_argument("--output", required=True, type=Path, help="where to write the filtered table")
    parser.add_argument("--topology-output", type=Path, help="where to write the topology verdicts")
    args = parser.parse_args(argv)

    result = run_chain(
        pd.read_csv(args.distances, low_memory=False),
        pd.read_csv(args.convergence, low_memory=False),
        channel=args.channel,
        model_dir=args.model_dir,
        topology_verdicts=pd.read_csv(args.topology) if args.topology else None,
        label=args.label)
    result.representatives.to_csv(args.output, index=False)
    print(f"wrote {args.output}")
    if args.topology_output:
        result.topology.to_csv(args.topology_output, index=False)
        print(f"wrote {args.topology_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
