#!/usr/bin/env python3
"""The quality-control chain, in the one order that is defensible.

    1. all-OK-3 convergence          earliest recycle from which every later successive-recycle
                                     transition stays under 3 A; that recycle is the trajectory's
                                     representative, and it and everything after it are retained
    2. selectivity-filter ring       trajectory rejected if any retained recycle exceeds 20 A
                                     (Kv2.1 G377, Cav1.2 EEEE, Nav1.5 DEKA)
    3. pore-VSD interface, 27 A      trajectory level: one bad recycle rejects the whole trajectory
                                     (Kv2.1 only - no cutoff is published for the other two)
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
    pipeline = Path(__file__).resolve().parent / "pipeline"
    if str(pipeline) not in sys.path:
        sys.path.insert(0, str(pipeline))
    import screen                                                      # noqa: PLC0415
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
              interface_threshold_A: float = INTERFACE_MAX_A, label: str = None,
              **topology_kwargs) -> ChainResult:
    """Pass one cohort through convergence, the filter ring, the interface QC and the topology
    screen, in that order and no other (see the module docstring for why the order is not free).

    Only Kv2.1 has a published pore-VSD interface cutoff, so for Nav1.5 and Cav1.2 the chain is the
    other three steps and the counts say the interface step was not run. That is a gap in what we
    have measured, not a decision that those channels need no interface check."""
    chan = get_channel(channel)
    label = label or f"{chan.name} cohort"
    distances = distances.copy()
    distances.attrs["dataset_label"] = label
    counts = {"channel": chan.name}

    kept, counts["convergence"] = apply_convergence(distances, convergence)
    kept, counts["filter_ring"] = apply_filter_ring(kept, chan, threshold_A=ring_threshold_A)

    if chan.interface_qc:
        before = kept["trajectory"].nunique()
        kept = apply_kv21_interface_qc(kept, threshold_A=interface_threshold_A, trajectory_level=True)
        counts["interface_qc"] = {
            "trajectories_rejected": before - kept["trajectory"].nunique(),
            "models_rejected": int(kept.attrs.get("kv21_interface_qc_rejected_rows", 0)),
            "threshold_A": float(interface_threshold_A)}
    else:
        counts["interface_qc"] = {"trajectories_rejected": 0, "models_rejected": 0,
                                  "not_run": f"no published interface cutoff for {chan.name}"}

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
