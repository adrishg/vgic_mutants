#!/usr/bin/env python3
"""Direction of the predicted Kv2.1 L403A mutation effect at the pore-VSD interface.

This is the repository-native version of the analysis behind the main Kv2.1 result:
default AlphaFold2 moves the L403A pore-VSD interface in the direction opposite to the
experimental 8SD3-to-8SDA change, and targeted masking removes that directional bias
without recovering the magnitude, the graded profile along S6, or the subunit asymmetry.

Two coordinate families are evaluated with identical machinery so that the specificity
claim is a like-for-like comparison:

* ``pore_vsd`` - the six same-subunit landmarks connecting K420, E423 and K427 to N179
  and V182 (model numbering K422/E425/K429 and N181/V184). Each structure contributes
  the mean of its four subunit-resolved distances. This is the region the mask targets.
* ``s6_ring`` - the six S6 cross-pore rings at V398, I401, A402, L403, I405 and V409
  (model 400/403/404/405/407/411). Each structure contributes the mean of the six
  intersubunit distances for that position. The mask does not target this region, so it
  is the negative control.

Estimand, matching the Methods: retained snapshots are reduced to a within-trajectory
median, AF2 model parameterizations are weighted equally within a seed, and seeds are
weighted equally. The predicted mutation effect is the L403A-minus-WT difference of
those seed-balanced means within one protocol. The experimental target is the
8SDA-minus-8SD3 change in the same coordinate. Uncertainty on the error reduction uses
whole-seed bootstrap resampling over seeds common to all four cohorts.

All four cohorts are read from the publication filter
(``_all_ok_rmsd_3A_structural_interface_qc``) so that vanilla and masked, WT and mutant,
pass through the same convergence, structural-integrity and pore-VSD interface QC.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from shared import data_access  # noqa: E402
from shared.distribution_statistics import parse_trajectory_metadata  # noqa: E402
from shared.plotting import KV21_PALETTE  # noqa: E402
from shared.seed_block_statistics import (  # noqa: E402
    reduce_trajectory_values,
    seed_model_stratum_weights,
)

BASE_SEED = 20261118
BOOTSTRAP_REPLICATES = 2000
DEFAULT_OUTPUT = ROOT / "analysis/statistics_revision/kv21_landmark_direction"
EXPERIMENTAL_DIR = ROOT / "kv21/experimental"

# Kv2.1 model numbering = displayed numbering + 2.
MODEL_OFFSET = 2

COHORTS = {
    ("WT", "vanilla"): "kv21/dataDistances/26-02-11_Kv2.1_wt_vanillaAF2_distances"
                       "_all_ok_rmsd_3A_structural_interface_qc.csv",
    ("WT", "masked"): "kv21/dataDistances/26-11-12_Kv2.1_wt_masked_v3_50pct_AF2_distances"
                      "_all_ok_rmsd_3A_structural_interface_qc.csv",
    ("L403A", "vanilla"): "kv21/dataDistances/26-02-11_Kv2.1_l403a_vanillaAF2_distances"
                          "_all_ok_rmsd_3A_structural_interface_qc.csv",
    ("L403A", "masked"): "kv21/dataDistances/26-11-12_Kv2.1_l403a_masked_v3_50pct_AF2_distances"
                         "_all_ok_rmsd_3A_structural_interface_qc.csv",
    ("F412L", "vanilla"): "kv21/dataDistances/26-02-11_Kv2.1_f412l_vanillaAF2_distances"
                          "_all_ok_rmsd_3A_structural_interface_qc.csv",
    ("F412L", "masked"): "kv21/dataDistances/26-11-12_Kv2.1_f412l_masked_v3_50pct_AF2_distances"
                         "_all_ok_rmsd_3A_structural_interface_qc.csv",
}

# F412L enters only the within-background displacement summary; the mutation-effect
# contrast needs a mutation-matched experimental reference, which F412L does not have.
MUTATION_EFFECT_CONDITIONS = ("WT", "L403A")

# (model residue A, model residue B, displayed label) - same subunit.
PORE_VSD_LANDMARKS = [
    (422, 181, "K420-N179"), (425, 181, "E423-N179"), (429, 181, "K427-N179"),
    (422, 184, "K420-V182"), (425, 184, "E423-V182"), (429, 184, "K427-V182"),
]
# (model residue, displayed label) - cross-pore ring, all six intersubunit pairs.
S6_RING_POSITIONS = [
    (400, "V398"), (403, "I401"), (404, "A402"),
    (405, "L403"), (407, "I405"), (411, "V409"),
]

CHAINS = "ABCD"


# --------------------------------------------------------------------------- experimental


def experimental_ca(pdb_id: str) -> dict[tuple[str, int], np.ndarray]:
    """Read chain/residue-keyed Cα coordinates from a deposited structure."""
    path = EXPERIMENTAL_DIR / f"{pdb_id}.pdb"
    coordinates: dict[tuple[str, int], np.ndarray] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                key = (line[21], int(line[22:26]))
                coordinates[key] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    if not coordinates:
        raise ValueError(f"No Cα atoms parsed from {path}")
    return coordinates


def experimental_pore_vsd(coordinates) -> dict[str, float]:
    """Mean same-subunit landmark distance over the four protomers."""
    values = {}
    for first, second, label in PORE_VSD_LANDMARKS:
        per_chain = [
            float(np.linalg.norm(
                coordinates[(chain, first - MODEL_OFFSET)]
                - coordinates[(chain, second - MODEL_OFFSET)]
            ))
            for chain in CHAINS
        ]
        values[label] = float(np.mean(per_chain))
    return values


def experimental_s6_ring(coordinates) -> dict[str, float]:
    """Mean of the six intersubunit distances in one cross-pore ring."""
    values = {}
    for residue, label in S6_RING_POSITIONS:
        points = [coordinates[(chain, residue - MODEL_OFFSET)] for chain in CHAINS]
        pairs = [
            float(np.linalg.norm(a - b)) for a, b in itertools.combinations(points, 2)
        ]
        values[label] = float(np.mean(pairs))
    return values


# --------------------------------------------------------------------------- predictions


def pore_vsd_columns(frame: pd.DataFrame) -> dict[str, list[str]]:
    """Same-subunit landmark columns, four per landmark."""
    residue_names = {181: "ASN181", 184: "VAL184", 422: "LYS422", 425: "GLU425", 429: "LYS429"}
    selected = {}
    for first, second, label in PORE_VSD_LANDMARKS:
        columns = []
        for chain in CHAINS:
            forward = f"CA_CA_{chain}_{residue_names[first]}_CA-{chain}_{residue_names[second]}_CA"
            reverse = f"CA_CA_{chain}_{residue_names[second]}_CA-{chain}_{residue_names[first]}_CA"
            if forward in frame.columns:
                columns.append(forward)
            elif reverse in frame.columns:
                columns.append(reverse)
        if len(columns) != len(CHAINS):
            raise KeyError(f"{label}: expected 4 same-subunit columns, found {len(columns)}")
        selected[label] = columns
    return selected


def s6_ring_columns(frame: pd.DataFrame) -> dict[str, list[str]]:
    """Cross-subunit ring columns, six per position."""
    import re

    pattern = re.compile(r"^CA_CA_([A-D])_([A-Z]{3})(\d+)_CA-([A-D])_([A-Z]{3})(\d+)_CA$")
    by_residue: dict[int, list[str]] = {}
    for column in frame.columns:
        match = pattern.match(str(column))
        if match and match.group(3) == match.group(6) and match.group(1) != match.group(4):
            by_residue.setdefault(int(match.group(3)), []).append(column)
    selected = {}
    for residue, label in S6_RING_POSITIONS:
        columns = by_residue.get(residue, [])
        if len(columns) != 6:
            raise KeyError(f"{label}: expected 6 intersubunit columns, found {len(columns)}")
        selected[label] = sorted(columns)
    return selected


def structure_level_values(frame: pd.DataFrame, columns: dict[str, list[str]]) -> pd.DataFrame:
    """One value per structure per coordinate: the mean over its constituent distances."""
    result = frame[["seed", "model_number"]].copy()
    for label, group in columns.items():
        numeric = frame[group].apply(pd.to_numeric, errors="coerce")
        result[label] = numeric.mean(axis=1)
    return result


def seed_balanced_mean(values: pd.DataFrame, label: str, seeds=None) -> float:
    """Trajectory median, then equal model mass within seed, then equal seed mass."""
    frame = values if seeds is None else values[values["seed"].isin(seeds)]
    reduced = reduce_trajectory_values(frame, label, reduction="median")
    if reduced.empty:
        return float("nan")
    weights = seed_model_stratum_weights(reduced)
    return float(np.average(reduced[label].to_numpy(float), weights=weights))


def cohort_vector(values: pd.DataFrame, labels: list[str], seeds=None) -> np.ndarray:
    return np.array([seed_balanced_mean(values, label, seeds) for label in labels])


# --------------------------------------------------------------------------- inference


def cosine_similarity(predicted: np.ndarray, experimental: np.ndarray) -> float:
    denominator = math.sqrt(float(predicted @ predicted) * float(experimental @ experimental))
    if denominator == 0:
        return float("nan")
    return float(predicted @ experimental) / denominator


def analyze(family: str, cohorts: dict, labels: list[str], experimental_delta: np.ndarray,
            rng: np.random.Generator, replicates: int):
    """Per-landmark table plus the summary row for one coordinate family."""
    predicted = {
        protocol: cohort_vector(cohorts[("L403A", protocol)][family], labels)
        - cohort_vector(cohorts[("WT", protocol)][family], labels)
        for protocol in ("vanilla", "masked")
    }

    rows = []
    for index, label in enumerate(labels):
        vanilla_error = abs(predicted["vanilla"][index] - experimental_delta[index])
        masked_error = abs(predicted["masked"][index] - experimental_delta[index])
        rows.append({
            "coordinate_family": family,
            "coordinate": label,
            "experimental_change_A": experimental_delta[index],
            "predicted_change_vanilla_A": predicted["vanilla"][index],
            "predicted_change_masked_A": predicted["masked"][index],
            "absolute_error_vanilla_A": vanilla_error,
            "absolute_error_masked_A": masked_error,
            "masking_closer_to_experiment": bool(masked_error < vanilla_error),
            "sign_agrees_vanilla": bool(
                np.sign(predicted["vanilla"][index]) == np.sign(experimental_delta[index])
            ),
            "sign_agrees_masked": bool(
                np.sign(predicted["masked"][index]) == np.sign(experimental_delta[index])
            ),
        })
    per_coordinate = pd.DataFrame(rows)

    contrast_keys = [
        key for key in cohorts if key[0] in MUTATION_EFFECT_CONDITIONS
    ]
    common = sorted(
        set.intersection(*(
            set(cohorts[key][family]["seed"].unique()) for key in contrast_keys
        ))
    )
    observed = (
        per_coordinate["absolute_error_masked_A"].mean()
        - per_coordinate["absolute_error_vanilla_A"].mean()
    )
    differences = []
    for _ in range(replicates):
        draw = rng.choice(common, size=len(common), replace=True)
        vanilla = (cohort_vector(cohorts[("L403A", "vanilla")][family], labels, draw)
                   - cohort_vector(cohorts[("WT", "vanilla")][family], labels, draw))
        masked = (cohort_vector(cohorts[("L403A", "masked")][family], labels, draw)
                  - cohort_vector(cohorts[("WT", "masked")][family], labels, draw))
        if not (np.isfinite(vanilla).all() and np.isfinite(masked).all()):
            continue
        differences.append(
            float(np.mean(np.abs(masked - experimental_delta)))
            - float(np.mean(np.abs(vanilla - experimental_delta)))
        )
    low, high = np.percentile(differences, [2.5, 97.5]) if differences else (np.nan, np.nan)

    summary = {
        "coordinate_family": family,
        "n_coordinates": len(labels),
        "n_common_seeds": len(common),
        "mean_experimental_change_A": float(experimental_delta.mean()),
        "mean_predicted_change_vanilla_A": float(predicted["vanilla"].mean()),
        "mean_predicted_change_masked_A": float(predicted["masked"].mean()),
        "cosine_similarity_vanilla": cosine_similarity(predicted["vanilla"], experimental_delta),
        "cosine_similarity_masked": cosine_similarity(predicted["masked"], experimental_delta),
        "mean_absolute_error_vanilla_A": float(per_coordinate["absolute_error_vanilla_A"].mean()),
        "mean_absolute_error_masked_A": float(per_coordinate["absolute_error_masked_A"].mean()),
        "mean_absolute_error_difference_A": float(observed),
        "CI_low_difference_A": float(low),
        "CI_high_difference_A": float(high),
        "pearson_r_vanilla_vs_experimental_profile": float(
            np.corrcoef(predicted["vanilla"], experimental_delta)[0, 1]
        ),
        "pearson_r_masked_vs_experimental_profile": float(
            np.corrcoef(predicted["masked"], experimental_delta)[0, 1]
        ),
        "n_coordinates_improved_by_masking": int(
            per_coordinate["masking_closer_to_experiment"].sum()
        ),
        "sign_agreement_vanilla": int(per_coordinate["sign_agrees_vanilla"].sum()),
        "sign_agreement_masked": int(per_coordinate["sign_agrees_masked"].sum()),
        "bootstrap_replicates": len(differences),
        "estimand": (
            "within-trajectory median; AF2 models weighted equally within seed; "
            "seeds weighted equally; whole-seed bootstrap over seeds common to all four cohorts"
        ),
        "resampling_unit": "whole recorded nominal seed label",
        "random_seed_pairing_status": (
            "not verified; matching numeric labels are nominal design keys"
        ),
    }
    return per_coordinate, summary, predicted


def ordered_interface_rmse(frame: pd.DataFrame, sd3, sda) -> pd.DataFrame:
    """Per-structure RMSE of the ordered four E423-N179 distances to each reference.

    Chain identity is arbitrary in a predicted homotetramer, so the four subunit
    distances are sorted before comparison and the experimental references are sorted
    the same way. This is the coordinate that carries the asymmetry claim: 8SDA has two
    near-wild-type interfaces and two elongated ones, which an ordered vector preserves
    and a single summary statistic does not.
    """
    columns = [f"CA_CA_{chain}_GLU425_CA-{chain}_ASN181_CA" for chain in CHAINS]
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise KeyError(f"Ordered interface RMSE requires {missing}")
    values = np.sort(frame[columns].apply(pd.to_numeric, errors="coerce").to_numpy(float), axis=1)
    reference_sd3 = np.sort([
        float(np.linalg.norm(sd3[(chain, 425 - MODEL_OFFSET)] - sd3[(chain, 181 - MODEL_OFFSET)]))
        for chain in CHAINS
    ])
    reference_sda = np.sort([
        float(np.linalg.norm(sda[(chain, 425 - MODEL_OFFSET)] - sda[(chain, 181 - MODEL_OFFSET)]))
        for chain in CHAINS
    ])
    result = frame[["seed", "model_number"]].copy()
    result["rmse_to_8SD3_A"] = np.sqrt(((values - reference_sd3) ** 2).mean(axis=1))
    result["rmse_to_8SDA_A"] = np.sqrt(((values - reference_sda) ** 2).mean(axis=1))
    result["closer_to_8SDA"] = (
        result["rmse_to_8SDA_A"] < result["rmse_to_8SD3_A"]
    ).astype(float)
    return result


def within_background_displacement(cohorts, family: str, labels: list[str]) -> pd.DataFrame:
    """Masked-minus-vanilla displacement within each sequence background.

    This is the protocol response, not a mutation effect: it says how far the same mask
    moved each background's own ensemble, which is what the roughly two-fold larger
    mutant response is based on.
    """
    rows = []
    for condition in sorted({key[0] for key in cohorts}):
        vanilla = cohorts[(condition, "vanilla")][family]
        masked = cohorts[(condition, "masked")][family]
        common = sorted(set(vanilla["seed"].unique()) & set(masked["seed"].unique()))
        displacement = (cohort_vector(masked, labels, common)
                        - cohort_vector(vanilla, labels, common))
        rows.append({
            "coordinate_family": family,
            "sequence_background": condition,
            "n_common_seeds": len(common),
            "mean_masked_minus_vanilla_A": float(displacement.mean()),
            "min_masked_minus_vanilla_A": float(displacement.min()),
            "max_masked_minus_vanilla_A": float(displacement.max()),
        })
    return pd.DataFrame(rows)


def make_figure(per_coordinate: pd.DataFrame, output: Path) -> None:
    """Per-landmark experimental target against both predicted mutation effects."""
    families = list(per_coordinate["coordinate_family"].unique())
    fig, axes = plt.subplots(1, len(families), figsize=(6.0 * len(families), 4.2), sharey=False)
    axes = np.atleast_1d(axes)
    titles = {
        "pore_vsd": "Pore–VSD landmarks (mask target)",
        "s6_ring": "S6 cross-pore rings (not targeted)",
    }
    for axis, family in zip(axes, families):
        block = per_coordinate[per_coordinate["coordinate_family"].eq(family)]
        positions = np.arange(len(block))
        axis.axhline(0.0, color="#999999", linewidth=0.8, zorder=0)
        axis.bar(positions - 0.26, block["experimental_change_A"], width=0.24,
                 color=KV21_PALETTE["experimental"], label="Experimental (8SDA − 8SD3)")
        axis.bar(positions, block["predicted_change_vanilla_A"], width=0.24,
                 color=KV21_PALETTE["L403A_VAN"], edgecolor="#555555", linewidth=0.4,
                 label="Predicted, vanilla")
        axis.bar(positions + 0.26, block["predicted_change_masked_A"], width=0.24,
                 color=KV21_PALETTE["L403A_HM"], label="Predicted, masked")
        axis.set_xticks(positions)
        axis.set_xticklabels(block["coordinate"], rotation=45, ha="right")
        axis.set_ylabel("L403A − WT change (Å)")
        axis.set_title(titles.get(family, family))
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(output.with_suffix(f".{suffix}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-replicates", type=int, default=BOOTSTRAP_REPLICATES)
    parser.add_argument("--random-seed", type=int, default=BASE_SEED)
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=True)

    sd3, sda = experimental_ca("8SD3"), experimental_ca("8SDA")
    experimental = {
        "pore_vsd": (experimental_pore_vsd(sd3), experimental_pore_vsd(sda)),
        "s6_ring": (experimental_s6_ring(sd3), experimental_s6_ring(sda)),
    }

    cohorts: dict[tuple[str, str], dict[str, pd.DataFrame]] = {}
    raw_frames: dict[tuple[str, str], pd.DataFrame] = {}
    provenance = []
    for key, relative in COHORTS.items():
        path = data_access.resolve(relative)
        frame = parse_trajectory_metadata(pd.read_csv(path, low_memory=False))
        raw_frames[key] = frame
        cohorts[key] = {
            "pore_vsd": structure_level_values(frame, pore_vsd_columns(frame)),
            "s6_ring": structure_level_values(frame, s6_ring_columns(frame)),
        }
        provenance.append({
            "condition": key[0], "protocol": key[1], "path": str(path),
            "retained_snapshots": int(len(frame)),
            "retained_trajectories": int(frame.groupby(["seed", "model_number"]).ngroups),
            "n_seed_labels": int(frame["seed"].nunique()),
        })
        print(f"{key[0]:>5} {key[1]:<8} {len(frame):>6} snapshots  {path}")

    rng = np.random.default_rng(arguments.random_seed)
    tables, summaries = [], []
    for family in ("pore_vsd", "s6_ring"):
        sd3_values, sda_values = experimental[family]
        labels = list(sd3_values)
        delta = np.array([sda_values[label] - sd3_values[label] for label in labels])
        per_coordinate, summary, _ = analyze(
            family, cohorts, labels, delta, rng, arguments.bootstrap_replicates
        )
        tables.append(per_coordinate)
        summaries.append(summary)

    displacement = pd.concat(
        [within_background_displacement(cohorts, "pore_vsd", list(experimental["pore_vsd"][0]))],
        ignore_index=True,
    )
    displacement.to_csv(
        arguments.output / "kv21_within_background_displacement.csv", index=False
    )

    ordered_rows = []
    for (condition, protocol), frame in raw_frames.items():
        if condition not in MUTATION_EFFECT_CONDITIONS:
            continue
        table = ordered_interface_rmse(frame, sd3, sda)
        ordered_rows.append({
            "sequence_background": condition,
            "protocol": protocol,
            "seed_balanced_rmse_to_8SDA_A": seed_balanced_mean(table, "rmse_to_8SDA_A"),
            "seed_balanced_rmse_to_8SD3_A": seed_balanced_mean(table, "rmse_to_8SD3_A"),
            "seed_balanced_fraction_closer_to_8SDA": seed_balanced_mean(table, "closer_to_8SDA"),
        })
    ordered = pd.DataFrame(ordered_rows)
    ordered.to_csv(arguments.output / "kv21_ordered_interface_rmse.csv", index=False)

    per_coordinate = pd.concat(tables, ignore_index=True)
    summary_frame = pd.DataFrame(summaries)
    per_coordinate.to_csv(arguments.output / "kv21_landmark_direction_per_coordinate.csv",
                          index=False)
    summary_frame.to_csv(arguments.output / "kv21_landmark_direction_summary.csv", index=False)
    pd.DataFrame(provenance).to_csv(arguments.output / "kv21_landmark_direction_cohorts.csv",
                                    index=False)
    make_figure(per_coordinate, arguments.output / "kv21_landmark_direction")

    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print("\n", per_coordinate.round(3).to_string(index=False), sep="")
        print("\n", summary_frame.round(3).to_string(index=False), sep="")
        print("\n", displacement.round(3).to_string(index=False), sep="")
        print("\n", ordered.round(4).to_string(index=False), sep="")
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
