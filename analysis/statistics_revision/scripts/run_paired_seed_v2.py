#!/usr/bin/env python3
"""Run the publication paired-seed statistical analysis.

All outputs are additive and written below ``paired_seed_v2``.  The script
uses nominal seed labels as design keys but reports that actual RNG values are
not recoverable from the available run metadata.
"""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.statistics_revision.scripts.run_seed_block_revision import (
    FINAL_DISTANCE_PATHS,
    GATE_COLUMNS,
    L403_COLUMNS,
    add_metadata,
    derive_l403a_threshold,
    master_cohort,
    read_csv,
)
from shared.paired_seed_statistics import (
    analyze_paired_seed_distance,
    joint_seed_bootstrap_contrast,
    paired_common_seed_bootstrap,
    paired_factorial_interaction,
    paired_seed_categorical_vector_contrast,
    paired_seed_ratio_of_means,
    paired_seed_retention_contrast,
    select_paired_estimand,
)


BASE_SEED = 20260824
MASK_IDS = {
    ("kv21", "vanilla"): "unmasked",
    ("kv21", "masked"): "kv21_common",
    ("nav15", "vanilla"): "unmasked",
    ("nav15", "masked"): "nav15_standard_plus_IFM",
    ("nav15", "masked_v2"): "nav15_v2",
    ("cav12", "vanilla"): "unmasked",
    ("cav12", "masked"): "condition_specific_primary_mask",
}

G402_PARTNERS = (
    "PHE1519", "VAL1520", "ALA1521", "VAL1522", "ILE1523", "MET1524",
    "ASP1525", "ASN1526", "PHE1527", "ASP1528", "TYR1529", "LEU1530",
    "THR1531", "ARG1532", "ASP1533", "TRP1534", "SER1535",
)

# The overlap screen is intentionally explicit. Adding a new shortest-distance
# column to an input CSV must not silently change the endpoint. GLU407 and
# GLU411 are excluded because they are the covalent/adjacent local environment
# of residue 406 and are not present in the current masked companion table;
# retaining them would make the vanilla--masked screen non-comparable.
R406_SCREEN_PARTNERS = (
    "PHE1519", "VAL1520", "ALA1521", "VAL1522", "ILE1523", "MET1524",
    "ASP1525", "ASN1526", "PHE1527", "ASP1528", "TYR1529", "LEU1530",
    "THR1531", "ARG1532", "ASP1533", "TRP1534", "SER1535",
    # GLU407 and GLU411 are deliberately omitted from the common panel.
)
R406_SCREEN_COLUMNS = tuple(
    f"shortest_ARG406-{partner}" for partner in R406_SCREEN_PARTNERS
)


def save(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def plot_zero_to_four(frame: pd.DataFrame, output: Path) -> None:
    """Write the publication figure beside its authoritative source table."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    x = frame["category"].to_numpy()
    width = 0.36
    ax.bar(x - width / 2, 100 * frame["probability_A"], width, label="Vanilla", color="#4C78A8")
    ax.bar(x + width / 2, 100 * frame["probability_B"], width, label="Masked", color="#E45756")
    ax.set(
        xlabel="Interfaces above experiment-anchored 12.8 Å cutoff",
        ylabel="Seed-balanced protocol sampling frequency (%)",
        xticks=x,
    )
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(figure_dir / "l403a_zero_to_four_subunits.png", dpi=300)
    fig.savefig(figure_dir / "l403a_zero_to_four_subunits.pdf")
    plt.close(fig)


def trajectory_summary(frame: pd.DataFrame, value: str, reduction: str) -> pd.DataFrame:
    work = add_metadata(frame) if not {"seed", "model_number"}.issubset(frame) else frame.copy()
    work[value] = pd.to_numeric(work[value], errors="coerce")
    work = work.dropna(subset=[value])
    grouped = work.groupby(["seed", "model_number"])[value]
    if reduction == "median":
        result = grouped.median()
    elif reduction == "mean":
        result = grouped.mean()
    elif reduction in {"earliest", "latest"}:
        ordered = work.sort_values(["seed", "model_number", "recycle_number", "pdb_file"])
        selected = ordered.groupby(["seed", "model_number"], as_index=False).head(1) \
            if reduction == "earliest" else ordered.groupby(["seed", "model_number"], as_index=False).tail(1)
        return selected[["seed", "model_number", value]].reset_index(drop=True)
    else:
        raise ValueError(reduction)
    return result.rename(value).reset_index()


def whole_seed_mean_ci(
    trajectory_frame: pd.DataFrame,
    value: str,
    *,
    n_bootstrap: int,
    random_seed: int,
) -> tuple[float, float, float, int]:
    """Estimate a mean across seeds and bootstrap complete seed blocks.

    Each contributing trajectory is reduced before entry. Seeds receive equal
    weight; trajectories with no defined conditional value do not contribute.
    """
    seed_values = (
        trajectory_frame.groupby("seed", sort=True)[value].mean().dropna()
    )
    if seed_values.empty:
        return np.nan, np.nan, np.nan, 0
    values = seed_values.to_numpy(float)
    rng = np.random.default_rng(random_seed)
    samples = values[
        rng.integers(0, len(values), size=(n_bootstrap, len(values)))
    ].mean(axis=1)
    low, high = np.quantile(samples, [.025, .975])
    return float(values.mean()), float(low), float(high), len(values)


def paired_scalar_rows(
    frame_a: pd.DataFrame,
    frame_b: pd.DataFrame,
    value: str,
    *,
    contrast: str,
    outcome: str,
    bootstrap: int,
    seed: int,
) -> list[dict[str, object]]:
    rows = []
    primary = joint_seed_bootstrap_contrast(
        frame_a, frame_b, value, n_bootstrap=bootstrap, random_seed=seed
    )
    rows.append({"contrast": contrast, "outcome": outcome, "estimand_id": "primary_joint_nominal_seed", **primary})
    common = paired_common_seed_bootstrap(
        frame_a, frame_b, value, n_bootstrap=bootstrap, random_seed=seed + 1
    )
    rows.append({"contrast": contrast, "outcome": outcome, "estimand_id": "common_seed", **common})
    model_a, model_b, audit = select_paired_estimand(
        frame_a, frame_b, value, estimand="common_model_seed"
    )
    model = paired_common_seed_bootstrap(
        model_a, model_b, value, n_bootstrap=bootstrap, random_seed=seed + 2
    )
    rows.append({
        "contrast": contrast, "outcome": outcome, "estimand_id": "common_model_seed",
        **audit, **model,
    })
    return rows


def nominalized_trajectory_table(
    frame: pd.DataFrame,
    value: str,
    *,
    nominal_seeds: list[int] | np.ndarray | pd.Series | None = None,
) -> pd.DataFrame:
    if nominal_seeds is None:
        nominal_seeds = sorted(pd.to_numeric(frame["seed"], errors="raise").astype(int).unique())
    summarized = trajectory_summary(frame, value, "mean").set_index(["seed", "model_number"])
    index = pd.MultiIndex.from_product(
        [list(nominal_seeds), [1, 2, 3, 4, 5]], names=["seed", "model_number"]
    )
    return summarized.reindex(index, fill_value=0).reset_index()


def g402_common_ca_columns(frame: pd.DataFrame, background: str) -> list[str]:
    """Return the atomically matched position-402 C-alpha partner panel.

    This deliberately rejects the legacy shortest-heavy columns so glycine can
    never be described as having a non-hydrogen side-chain distance.
    """
    residue = {"WT": "GLY402", "G402S": "SER402"}.get(background)
    if residue is None:
        raise ValueError("background must be WT or G402S")
    columns = [f"CA_{residue}_CA-{partner}_CA" for partner in G402_PARTNERS]
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise KeyError(f"Missing atomically matched position-402 C-alpha columns: {missing}")
    if any(column.startswith("shortest_") for column in columns):
        raise AssertionError("G402 primary comparison cannot use undocumented side-chain fallback")
    return columns


def make_inventory(master: pd.DataFrame, output: Path) -> None:
    inventory = master.groupby(["channel", "sequence", "protocol", "ensemble_id"], as_index=False).agg(
        nominal_seed_labels=("seed", "nunique"),
        retained_snapshots=("distance_final_cohort", "sum"),
        mapping_qc_snapshots=("mapping_QC_pass", "sum"),
        convergence_qc_snapshots=("all_ok", "sum"),
        structural_qc_snapshots=("all_ok_3", "sum"),
    )
    trajectory_counts = master[["ensemble_id", "seed", "model_number"]].drop_duplicates().groupby(
        "ensemble_id", as_index=False
    ).size().rename(columns={"size": "nominal_model_seed_trajectories"})
    inventory = inventory.merge(trajectory_counts, on="ensemble_id", validate="one_to_one")
    inventory["mask_id"] = [MASK_IDS.get((c, p), "unresolved") for c, p in zip(inventory.channel, inventory.protocol)]
    seed_ranges = master.groupby("ensemble_id").seed.agg(seed_id_min="min", seed_id_max="max").reset_index()
    inventory = inventory.merge(seed_ranges, on="ensemble_id", validate="one_to_one")
    inventory["actual_random_seed_status"] = "unavailable in run metadata; numeric labels are nominal design keys"
    inventory["source_csv"] = inventory.ensemble_id.map(FINAL_DISTANCE_PATHS).fillna("")
    inventory.to_csv(output / "DATASET_INVENTORY.tsv", sep="\t", index=False)

    trajectory = master.groupby(
        ["channel", "sequence", "protocol", "ensemble_id", "seed", "model_number"], as_index=False
    ).agg(
        trajectory_basename=("pdb_basename", "first"),
        mapping_qc=("mapping_QC_pass", "max"),
        convergence_qc=("all_ok", "max"),
        final_qc=("all_ok_3", "max"),
        analysis_final=("distance_final_cohort", "max"),
    )
    trajectory["mask_id"] = [MASK_IDS.get((c, p), "unresolved") for c, p in zip(trajectory.channel, trajectory.protocol)]
    trajectory["actual_random_seed"] = ""
    trajectory["actual_random_seed_status"] = "not recorded/recovered"
    trajectory["nominal_in_design"] = True
    trajectory = trajectory.rename(columns={
        "sequence": "sequence_background", "protocol": "msa_protocol",
        "seed": "seed_id", "model_number": "af2_model",
    })
    columns = [
        "channel", "sequence_background", "msa_protocol", "mask_id", "seed_id",
        "actual_random_seed", "actual_random_seed_status", "af2_model",
        "trajectory_basename", "nominal_in_design", "mapping_qc", "convergence_qc",
        "final_qc", "analysis_final",
    ]
    trajectory[columns].to_csv(output / "SEED_REGISTRY.tsv", sep="\t", index=False)

    trajectory_flow = master.groupby(
        ["channel", "sequence", "protocol", "ensemble_id", "seed", "model_number"], as_index=False
    ).agg(
        mapping_qc=("mapping_QC_pass", "max"), convergence_qc=("all_ok", "max"),
        structural_qc=("all_ok_3", "max"), interface_qc=("distance_source_membership", "max"),
        analysis_final=("distance_final_cohort", "max"),
    )
    flow = trajectory_flow.groupby(
        ["channel", "sequence", "protocol", "ensemble_id"], as_index=False
    ).agg(
        nominal_seeds=("seed", "nunique"),
        nominal_model_seed_trajectories=("model_number", "size"),
        mapping_qc_trajectories=("mapping_qc", "sum"),
        converged_trajectories=("convergence_qc", "sum"),
        final_structural_qc_trajectories=("structural_qc", "sum"),
        channel_interface_qc_trajectories=("interface_qc", "sum"),
        analysis_final_trajectories=("analysis_final", "sum"),
    )
    snapshots = master.groupby("ensemble_id", as_index=False).agg(
        retained_snapshots=("distance_final_cohort", "sum")
    )
    exclusion = master[~master.exclusion_reason.eq("included_in_final_distance_cohort")].groupby(
        "ensemble_id"
    ).exclusion_reason.agg(lambda values: values.value_counts().index[0] if len(values) else "none")
    flow = flow.merge(snapshots, on="ensemble_id", validate="one_to_one")
    flow["primary_exclusion_reason"] = flow.ensemble_id.map(exclusion).fillna("none")
    save(flow, output / "master_cohort_flow_summary.csv")


def retention(master: pd.DataFrame, output: Path, bootstrap: int) -> None:
    trajectory = master.groupby(
        ["channel", "sequence", "protocol", "seed", "model_number"], as_index=False
    ).agg(
        mapping_qc=("mapping_QC_pass", "max"),
        convergence_qc=("all_ok", "max"),
        structural_qc=("all_ok_3", "max"),
        analysis_final=("distance_final_cohort", "max"),
    )
    rows = []
    for (channel, sequence), backgrounds in trajectory.groupby(["channel", "sequence"]):
        if "vanilla" not in set(backgrounds.protocol):
            continue
        a = backgrounds[backgrounds.protocol.eq("vanilla")]
        for protocol, b in backgrounds[~backgrounds.protocol.eq("vanilla")].groupby("protocol"):
            for stage_index, stage in enumerate(("mapping_qc", "convergence_qc", "structural_qc", "analysis_final")):
                result = paired_seed_retention_contrast(
                    a, b, pass_col=stage,
                    nominal_seeds_a=sorted(a.seed.unique()),
                    nominal_seeds_b=sorted(b.seed.unique()),
                    n_bootstrap=bootstrap,
                    random_seed=BASE_SEED + 1000 + len(rows) + stage_index,
                )
                rows.append({
                    "channel": channel, "sequence_background": sequence,
                    "contrast": f"{protocol}-vanilla", "stage": stage,
                    "mask_id_B": MASK_IDS.get((channel, protocol), "unresolved"), **result,
                })
    save(pd.DataFrame(rows), output / "paired_retention_contrasts.csv")


def l403a(output: Path, bootstrap: int, permutations: int) -> None:
    threshold, wt_vector, mutant_vector = derive_l403a_threshold(output)
    frame = add_metadata(read_csv(ROOT / "kv21/dataDistances/analysis/L403A_E423_N179_all_structure_distances.csv"))
    frame[L403_COLUMNS] = frame[L403_COLUMNS].apply(pd.to_numeric, errors="coerce")
    frame["max_E423_N179_A"] = frame[L403_COLUMNS].max(axis=1)
    ordered = np.sort(frame[L403_COLUMNS].to_numpy(float), axis=1)
    frame["RMSE_to_8SD3_A"] = np.sqrt(np.mean((ordered - wt_vector) ** 2, axis=1))
    frame["RMSE_to_8SDA_A"] = np.sqrt(np.mean((ordered - mutant_vector) ** 2, axis=1))
    frame["closer_to_8SDA"] = (frame.RMSE_to_8SDA_A < frame.RMSE_to_8SD3_A).astype(float)
    frame["shifted_subunits"] = frame[L403_COLUMNS].ge(threshold).sum(axis=1)
    frame["any_shifted"] = frame.shifted_subunits.ge(1).astype(float)
    parts = {name: part.copy() for name, part in frame.groupby("condition")}

    rows = []
    definitions = [
        ("maximum_E423_N179_CA_distance_A", "max_E423_N179_A", "median"),
        ("ordered_vector_RMSE_to_8SD3_A", "RMSE_to_8SD3_A", "median"),
        ("ordered_vector_RMSE_to_8SDA_A", "RMSE_to_8SDA_A", "median"),
        ("fraction_closer_to_8SDA_ordered_vector", "closer_to_8SDA", "mean"),
        ("fraction_any_subunit_at_experiment_anchored_cutoff", "any_shifted", "mean"),
    ]
    for index, (outcome, column, reduction) in enumerate(definitions):
        a, b = (trajectory_summary(parts[p], column, reduction) for p in ("vanilla", "masked"))
        rows.extend(paired_scalar_rows(
            a, b, column, contrast="L403A masked - vanilla", outcome=outcome,
            bootstrap=bootstrap, seed=BASE_SEED + 2000 + index * 10,
        ))
        if column == "any_shifted":
            ratio = paired_seed_ratio_of_means(
                a, b, column, n_bootstrap=bootstrap, random_seed=BASE_SEED + 2060
            )
            rows.append({
                "contrast": "L403A masked / vanilla", "outcome": outcome,
                "estimand_id": "ratio_secondary", **ratio,
            })
    for reduction in ("earliest", "latest"):
        a, b = (trajectory_summary(parts[p], "any_shifted", reduction) for p in ("vanilla", "masked"))
        rows.extend(paired_scalar_rows(
            a, b, "any_shifted", contrast="L403A masked - vanilla",
            outcome=f"any_shifted_{reduction}_retained_recycle",
            bootstrap=bootstrap, seed=BASE_SEED + 2070 + (reduction == "latest"),
        ))
    save(pd.DataFrame(rows), output / "l403a_focal_paired_seed.csv")

    w1_rows = []
    a = trajectory_summary(parts["vanilla"], "max_E423_N179_A", "median")
    b = trajectory_summary(parts["masked"], "max_E423_N179_A", "median")
    for index, estimand in enumerate(("primary_joint_nominal_seed", "common_seed", "common_model_seed")):
        w1_rows.append(analyze_paired_seed_distance(
            a, b, "max_E423_N179_A", estimand=estimand,
            n_bootstrap=bootstrap, n_permutations=permutations,
            random_seed=BASE_SEED + 2100 + index * 10,
        ))
    save(pd.DataFrame(w1_rows), output / "l403a_w1_paired_seed.csv")

    sensitivity_rows = []
    for cutoff in np.arange(11.5, 14.5001, .25):
        work = frame.assign(_target=frame[L403_COLUMNS].ge(cutoff).any(axis=1).astype(float))
        a, b = (trajectory_summary(work[work.condition.eq(p)], "_target", "mean") for p in ("vanilla", "masked"))
        result = joint_seed_bootstrap_contrast(
            a, b, "_target", n_bootstrap=bootstrap,
            random_seed=BASE_SEED + 2200 + int(round(cutoff * 100)),
        )
        sensitivity_rows.append({
            "cutoff_A": cutoff,
            "cutoff_role": "prediction-independent experiment-anchored cutoff" if abs(cutoff - threshold) < .126 else "sensitivity grid",
            **result,
        })
    save(pd.DataFrame(sensitivity_rows), output / "l403a_threshold_sensitivity.csv")

    category_summary, category_rows = paired_seed_categorical_vector_contrast(
        parts["vanilla"], parts["masked"], category_col="shifted_subunits",
        categories=[0, 1, 2, 3, 4], n_bootstrap=bootstrap, random_seed=BASE_SEED + 2300,
    )
    for key, value in category_summary.items():
        category_rows[key] = value
    category_rows["chain_mapping_note"] = "sorted/count summary; predicted chain labels are not mapped to experimental chain identities"
    save(category_rows, output / "l403a_zero_to_four_subunits.csv")
    plot_zero_to_four(category_rows, output)

    yields = []
    for protocol in ("vanilla", "masked"):
        nominal = nominalized_trajectory_table(parts[protocol], "any_shifted")
        nominal["protocol"] = protocol
        yields.append(nominal)
    a, b = yields
    yield_result = paired_scalar_rows(
        a, b, "any_shifted", contrast="L403A masked - vanilla",
        outcome="QC_adjusted_fraction_retained_snapshots_with_any_shifted_subunit_per_nominal_trajectory",
        bootstrap=bootstrap, seed=BASE_SEED + 2400,
    )
    save(pd.DataFrame(yield_result), output / "usable_target_geometry_yields_l403a.csv")

    four = {}
    four_binary = {}
    for background in ("WT", "L403A"):
        for protocol in ("vanilla", "masked"):
            if background == "L403A":
                source = parts[protocol]
            else:
                source = add_metadata(read_csv(ROOT / FINAL_DISTANCE_PATHS[f"kv21|WT|{protocol}"]))
                source[L403_COLUMNS] = source[L403_COLUMNS].apply(pd.to_numeric, errors="coerce")
                source["max_E423_N179_A"] = source[L403_COLUMNS].max(axis=1)
            source["any_elongated_interface"] = (
                source[L403_COLUMNS].ge(threshold).any(axis=1).astype(float)
            )
            key = f"{'wt' if background == 'WT' else 'mutant'}_{protocol}"
            four[key] = trajectory_summary(
                source, "max_E423_N179_A", "median"
            )
            four_binary[key] = trajectory_summary(
                source, "any_elongated_interface", "mean"
            )
    interaction = paired_factorial_interaction(
        four, "max_E423_N179_A", statistic="mean", n_bootstrap=bootstrap,
        random_seed=BASE_SEED + 2500,
    )
    interaction.update({
        "mutant": "L403A", "outcome": "maximum E423-N179 C-alpha distance",
        "estimand": "mean complete-seed [(L403A masked-L403A vanilla)-(WT masked-WT vanilla)]",
    })
    binary_interaction = paired_factorial_interaction(
        four_binary, "any_elongated_interface", statistic="mean",
        n_bootstrap=bootstrap, random_seed=BASE_SEED + 2600,
    )
    binary_interaction.update({
        "mutant": "L403A",
        "outcome": "fraction of retained snapshots with any elongated interface",
        "estimand": (
            "mean complete-seed [(L403A masked-L403A vanilla)-"
            "(WT masked-WT vanilla)]; binary snapshot threshold"
        ),
        "threshold_A": threshold,
        "binary_rule": "1 if any of four E423-N179 C-alpha distances >= threshold",
    })
    interaction["threshold_A"] = np.nan
    interaction["binary_rule"] = "not applicable; continuous trajectory median"
    save(
        pd.DataFrame([interaction, binary_interaction]),
        output / "kv21_interactions_paired_seed.csv",
    )


def f412l(output: Path, bootstrap: int) -> None:
    path = ROOT / "kv21/dataRMSD/analysis/comparison_v5/f412l_pocket_D_paper_nexus_shortest_contacts_long_v5.csv"
    long = read_csv(path)
    value = "Shortest heavy-atom distance (Å)"
    wide = long.pivot_table(index=["Protocol", "pdb_file"], columns="Contact", values=value, aggfunc="first").reset_index()
    rename = {}
    for column in wide:
        if "L412–L316" in str(column): rename[column] = "L412_L316_A"
        elif "L412–L329" in str(column): rename[column] = "L412_L329_A"
        elif "L412–L403" in str(column): rename[column] = "L412_L403_A"
    wide = add_metadata(wide.rename(columns=rename))
    parts = {name.lower(): part.copy() for name, part in wide.groupby("Protocol")}
    rows = []
    for metric_index, column in enumerate(rename.values()):
        a, b = (trajectory_summary(parts[p], column, "median") for p in ("vanilla", "masked"))
        rows.extend(paired_scalar_rows(
            a, b, column, contrast="F412L masked - vanilla", outcome=f"continuous_{column}",
            bootstrap=bootstrap, seed=BASE_SEED + 3000 + metric_index * 100,
        ))
        for threshold in (3.5, 4.0, 4.5):
            threshold_parts = {}
            for protocol in ("vanilla", "masked"):
                work = parts[protocol].assign(_prox=pd.to_numeric(parts[protocol][column], errors="coerce").le(threshold).astype(float))
                threshold_parts[protocol] = trajectory_summary(work, "_prox", "mean")
            result = paired_scalar_rows(
                threshold_parts["vanilla"], threshold_parts["masked"], "_prox",
                contrast="F412L masked - vanilla",
                outcome=f"fraction_{column}_within_{threshold:.1f}A",
                bootstrap=bootstrap, seed=BASE_SEED + 3100 + metric_index * 20 + int(threshold * 2),
            )
            for row in result:
                row["threshold_A"] = threshold
            rows.extend(result)
        for overlap in (1.8, 2.0, 2.2):
            overlap_parts = {}
            for protocol in ("vanilla", "masked"):
                work = parts[protocol].assign(_overlap=pd.to_numeric(parts[protocol][column], errors="coerce").lt(overlap).astype(float))
                overlap_parts[protocol] = trajectory_summary(work, "_overlap", "mean")
            rows.extend(paired_scalar_rows(
                overlap_parts["vanilla"], overlap_parts["masked"], "_overlap",
                contrast="F412L masked - vanilla", outcome=f"severe_overlap_{column}_below_{overlap:.1f}A",
                bootstrap=bootstrap, seed=BASE_SEED + 3200 + metric_index * 20 + int(overlap * 10),
            ))
    save(pd.DataFrame(rows), output / "f412l_contacts_paired_seed.csv")

    yields = []
    for column in rename.values():
        condition = {}
        for protocol in ("vanilla", "masked"):
            work = parts[protocol].assign(_target=pd.to_numeric(parts[protocol][column], errors="coerce").le(4).astype(float))
            condition[protocol] = nominalized_trajectory_table(work, "_target")
        yields.extend(paired_scalar_rows(
            condition["vanilla"], condition["masked"], "_target",
            contrast="F412L masked - vanilla", outcome=f"QC_adjusted_within_4A_yield_{column}",
            bootstrap=bootstrap, seed=BASE_SEED + 3300 + len(yields),
        ))
    save(pd.DataFrame(yields), output / "usable_target_geometry_yields_f412l.csv")


def nav15(output: Path, bootstrap: int, permutations: int) -> None:
    tables = {}
    for ensemble, relative in FINAL_DISTANCE_PATHS.items():
        if not ensemble.startswith("nav15|"):
            continue
        frame = add_metadata(read_csv(ROOT / relative))
        residue = "PHE" if "|WT|" in ensemble else "GLN"
        frame["motif_receptor_A"] = frame[[
            f"CA_{residue}1170_CA-ASN1343_CA", f"CA_{residue}1170_CA-ASN1449_CA"
        ]].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        if set(GATE_COLUMNS).issubset(frame.columns):
            frame["gate_span_A"] = frame[GATE_COLUMNS].apply(pd.to_numeric, errors="coerce").max(axis=1)
        tables[ensemble] = frame
    specs = [
        ("nav15|WT|vanilla", "nav15|QQQ|vanilla", "vanilla QQQ - vanilla WT", "unmasked"),
        ("nav15|WT|vanilla", "nav15|WT|masked", "WT nav15_standard - vanilla", "nav15_standard"),
        ("nav15|WT|vanilla", "nav15|WT|masked_v2", "WT nav15_v2 - vanilla", "nav15_v2"),
        ("nav15|QQQ|vanilla", "nav15|QQQ|masked", "QQQ nav15_standard_plus_IFM - vanilla", "nav15_standard_plus_IFM"),
        ("nav15|QQQ|vanilla", "nav15|QQQ|masked_v2", "QQQ nav15_v2 - vanilla", "nav15_v2"),
    ]
    rows = []
    for index, (key_a, key_b, contrast, mask_id) in enumerate(specs):
        for outcome in ("motif_receptor_A", "gate_span_A"):
            if outcome not in tables[key_a] or outcome not in tables[key_b]:
                rows.append({
                    "contrast": contrast, "outcome": outcome, "estimand_id": "not_estimable",
                    "mask_id": mask_id, "status": "required gate coordinates absent from one registered source table",
                })
                continue
            a, b = (trajectory_summary(tables[key], outcome, "median") for key in (key_a, key_b))
            scalar = paired_scalar_rows(
                a, b, outcome, contrast=contrast, outcome=outcome,
                bootstrap=bootstrap, seed=BASE_SEED + 4000 + index * 100,
            )
            for row in scalar:
                row["mask_id"] = mask_id
            rows.extend(scalar)
            if outcome == "motif_receptor_A":
                w1 = analyze_paired_seed_distance(
                    a, b, outcome, estimand="primary_joint_nominal_seed",
                    n_bootstrap=bootstrap, n_permutations=permutations,
                    random_seed=BASE_SEED + 4050 + index * 100,
                )
                rows.append({"contrast": contrast, "outcome": "motif_receptor_W1", "estimand_id": "primary_joint_nominal_seed", "mask_id": mask_id, **w1})
    save(pd.DataFrame(rows), output / "nav15_focal_paired_seed.csv")


def g402s(output: Path, bootstrap: int, permutations: int) -> None:
    frames = {}
    source_frames = []
    for background in ("WT", "G402S"):
        for protocol in ("vanilla", "masked"):
            frame = add_metadata(read_csv(ROOT / FINAL_DISTANCE_PATHS[f"cav12|{background}|{protocol}"]))
            residue = "GLY402" if background == "WT" else "SER402"
            ca_columns = g402_common_ca_columns(frame, background)
            numeric = frame[ca_columns].apply(pd.to_numeric, errors="coerce")
            normalized = [column.replace(residue, "POSITION402") for column in ca_columns]
            numeric.columns = normalized
            for column in normalized:
                frame[column] = numeric[column]
            frames[(background, protocol)] = frame

            source = frame[["pdb_file", "seed", "model_number", "recycle_number", *normalized]].melt(
                id_vars=["pdb_file", "seed", "model_number", "recycle_number"],
                value_vars=normalized,
                var_name="source_column",
                value_name="position402_partner_CA_distance_A",
            )
            source["partner"] = source.source_column.str.extract(r"-(.+)_CA$")[0]
            source["normalized_metric"] = source["source_column"]
            source["source_column"] = source["source_column"].str.replace(
                "POSITION402", residue, regex=False
            )
            source["sequence"] = background
            source["protocol"] = protocol
            source["mask_id"] = MASK_IDS[("cav12", protocol)]
            source["metric_definition"] = (
                "position-402 C-alpha to partner C-alpha; identical atom names in WT and G402S"
            )
            source_frames.append(source)

    partner_rows = []
    partner_specs = [
        ("WT", "vanilla", "G402S", "vanilla", "vanilla G402S - vanilla WT", "primary sequence contrast"),
        ("WT", "masked", "G402S", "masked", "masked G402S - masked WT", "protocol-specific context"),
        ("G402S", "vanilla", "G402S", "masked", "G402S masked - G402S vanilla", "within-G402S masking effect"),
    ]
    for partner_index, partner in enumerate(G402_PARTNERS):
        value = f"CA_POSITION402_CA-{partner}_CA"
        for spec_index, (seq_a, protocol_a, seq_b, protocol_b, contrast, role) in enumerate(partner_specs):
            a = trajectory_summary(frames[(seq_a, protocol_a)], value, "median")
            b = trajectory_summary(frames[(seq_b, protocol_b)], value, "median")
            result = paired_scalar_rows(
                a,
                b,
                value,
                contrast=contrast,
                outcome="matched_position402_CA_distance_A",
                bootstrap=bootstrap,
                seed=BASE_SEED + 5300 + partner_index * 20 + spec_index * 3,
            )
            for row in result:
                row.update({
                    "partner": partner,
                    "comparison_role": role,
                    "metric_definition": (
                        "position-402 C-alpha to partner C-alpha; identical atom names in WT and G402S"
                    ),
                    "resampling_unit": "whole recorded nominal seed label",
                    "actual_random_seed_pairing_status": (
                        "not verified; matching numeric labels are nominal design keys"
                    ),
                })
            partner_rows.extend(result)
    save(pd.DataFrame(partner_rows), output / "g402s_matched_ca_partner_statistics.csv")
    save(
        pd.concat(source_frames, ignore_index=True),
        output / "source_data/g402s_matched_ca_1519_1535.csv.gz",
    )

    # No nearest-partner category is emitted.  A minimum over 17 candidates is
    # difficult to interpret as a contact and obscures the matched WT-mutant
    # comparison.  The complete 17-partner distance panel above is the current
    # source result; categorical outputs use only the declared coordinate rule.


def g406r(output: Path, bootstrap: int) -> None:
    def read_g406r(protocol: str) -> pd.DataFrame:
        """Load the current R406 panel, including its companion export."""
        source = ROOT / FINAL_DISTANCE_PATHS[f"cav12|G406R|{protocol}"]
        frame = add_metadata(read_csv(source))
        companion = source.with_name(source.stem + "_extra.csv")
        if any(column.startswith("shortest_ARG406-") for column in frame.columns):
            extra = None
        else:
            try:
                extra = add_metadata(read_csv(companion))
            except Exception:  # companion is optional for exports already carrying R406
                extra = None
        if extra is not None and "pdb_file" in extra.columns:
            duplicate = [column for column in extra.columns if column in frame.columns and column != "pdb_file"]
            extra = extra.drop(columns=duplicate)
            frame = frame.merge(extra, on="pdb_file", how="left", validate="one_to_one")
        return frame

    parts = {
        protocol: read_g406r(protocol)
        for protocol in ("vanilla", "masked")
    }
    rows, sensitivity, yields, denominator_rows, column_audit_rows = [], [], [], [], []
    reconciliation_rows = []
    centered = list(R406_SCREEN_COLUMNS)
    for protocol, frame in parts.items():
        observed = {column for column in frame if column.startswith("shortest_ARG406-")}
        missing = sorted(set(centered) - observed)
        unregistered = sorted(observed - set(centered))
        if missing:
            raise KeyError(
                f"{protocol} R406 overlap-screen column mismatch; "
                f"missing={missing}, unregistered={unregistered}"
            )
        numeric = frame[centered].apply(pd.to_numeric, errors="coerce")
        if not np.isfinite(numeric.to_numpy(float)).all():
            bad = numeric.columns[~np.isfinite(numeric.to_numpy(float)).all(axis=0)].tolist()
            raise ValueError(f"Nonfinite R406 overlap-screen distances in {protocol}: {bad}")
        for column, partner in zip(centered, R406_SCREEN_PARTNERS):
            values = numeric[column]
            column_audit_rows.append({
                "protocol": protocol,
                "focal_residue": "ARG406",
                "partner_residue": partner,
                "source_column": column,
                "screen_rule": "fails if distance < 2.0 A",
                "retained_snapshots": len(values),
                "missing_or_nonfinite_values": int((~np.isfinite(values)).sum()),
                "minimum_observed_A": float(values.min()),
                "maximum_observed_A": float(values.max()),
                "snapshots_below_2A": int(values.lt(2).sum()),
                "fraction_below_2A": float(values.lt(2).mean()),
                "residue_pair_verification": "verified from explicit CSV column name",
                "closest_atom_names": "not recorded in source table",
                "atom_pair_verification": (
                    "unavailable: generator and per-snapshot closest-atom identities are absent"
                ),
            })
    for overlap_cutoff in (1.8, 2.0, 2.2):
        derived = {}
        for protocol, frame in parts.items():
            work = frame.copy()
            distances = work[centered].apply(pd.to_numeric, errors="coerce")
            work["overlap_pass"] = (~distances.lt(overlap_cutoff).any(axis=1)).astype(float)
            for partner in ("ASP1528", "ASP1533"):
                proximity = pd.to_numeric(work[f"shortest_ARG406-{partner}"], errors="coerce").le(4)
                work[f"conditional_{partner}"] = np.where(work.overlap_pass.eq(1), proximity.astype(float), np.nan)
                work[f"unconditional_{partner}"] = (work.overlap_pass.eq(1) & proximity).astype(float)
            derived[protocol] = work
        for outcome in ("overlap_pass", "conditional_ASP1528", "conditional_ASP1533", "unconditional_ASP1528", "unconditional_ASP1533"):
            a, b = (trajectory_summary(derived[p], outcome, "mean") for p in ("vanilla", "masked"))
            result = paired_scalar_rows(
                a, b, outcome, contrast="G406R masked - vanilla",
                outcome=outcome, bootstrap=bootstrap,
                seed=BASE_SEED + 6000 + int(overlap_cutoff * 100) + len(rows),
            )
            for row in result:
                row["overlap_threshold_A"] = overlap_cutoff
                row["terminology"] = "passed the R406-centered overlap screen"
                row["resampling_unit"] = "whole recorded nominal seed label"
                row["actual_random_seed_pairing_status"] = (
                    "not verified; matching numeric labels are nominal design keys"
                )
            (rows if overlap_cutoff == 2.0 else sensitivity).extend(result)
        # Any overlap-pass snapshot per trajectory, kept distinct from snapshot fraction.
        for protocol in ("vanilla", "masked"):
            trajectory_any = derived[protocol].groupby(["seed", "model_number"], as_index=False).overlap_pass.max()
            derived[protocol + "_any"] = trajectory_any
        any_rows = paired_scalar_rows(
            derived["vanilla_any"], derived["masked_any"], "overlap_pass",
            contrast="G406R masked - vanilla", outcome="fraction_nominal_trajectories_with_any_overlap_pass_snapshot",
            bootstrap=bootstrap, seed=BASE_SEED + 6100 + int(overlap_cutoff * 100),
        )
        for row in any_rows:
            row["overlap_threshold_A"] = overlap_cutoff
            row["terminology"] = "passed the R406-centered overlap screen"
            row["resampling_unit"] = "whole recorded nominal seed label"
            row["actual_random_seed_pairing_status"] = (
                "not verified; matching numeric labels are nominal design keys"
            )
        (rows if overlap_cutoff == 2.0 else sensitivity).extend(any_rows)
        if overlap_cutoff == 2.0:
            for protocol, frame in derived.items():
                if protocol.endswith("_any"):
                    continue
                pass_mask = frame["overlap_pass"].eq(1)
                trajectory_any = frame.groupby(["seed", "model_number"]).overlap_pass.max()
                base = {
                    "protocol": protocol,
                    "overlap_threshold_A": overlap_cutoff,
                    "total_snapshots": len(frame),
                    "passing_snapshots": int(pass_mask.sum()),
                    "unconditional_snapshot_pass_fraction": float(pass_mask.mean()),
                    "total_nominal_trajectories": int(trajectory_any.size),
                    "trajectories_with_any_passing_snapshot": int(trajectory_any.sum()),
                    "trajectory_any_pass_fraction": float(trajectory_any.mean()),
                    "snapshot_pass_denominator_definition": "all retained QC-qualified snapshots",
                    "trajectory_pass_denominator_definition": (
                        "all retained nominal seed x AF2-model trajectories"
                    ),
                }
                denominator_rows.append({**base, "partner": "screen_only"})
                for partner in ("ASP1528", "ASP1533"):
                    contact = pd.to_numeric(
                        frame[f"shortest_ARG406-{partner}"], errors="raise"
                    ).le(4)
                    numerator = int((pass_mask & contact).sum())
                    denominator_rows.append({
                        **base,
                        "partner": partner,
                        "passing_and_within_4A_snapshots": numerator,
                        "conditional_contact_denominator_passing_snapshots": int(pass_mask.sum()),
                        "conditional_contact_fraction": (
                            numerator / int(pass_mask.sum()) if pass_mask.any() else np.nan
                        ),
                        "unconditional_contact_denominator_all_snapshots": len(frame),
                        "unconditional_pass_and_contact_fraction": numerator / len(frame),
                        "contact_rule": "distance <= 4.0 A",
                    })
            for partner in ("ASP1528", "ASP1533"):
                nominal = {}
                for protocol in ("vanilla", "masked"):
                    nominal[protocol] = nominalized_trajectory_table(derived[protocol], f"unconditional_{partner}")
                yields.extend(paired_scalar_rows(
                    nominal["vanilla"], nominal["masked"], f"unconditional_{partner}",
                    contrast="G406R masked - vanilla",
                    outcome=f"QC_adjusted_overlap_pass_and_within_4A_{partner}_yield",
                    bootstrap=bootstrap, seed=BASE_SEED + 6200 + len(yields),
                ))

                value = f"conditional_{partner}"
                protocol_values = {}
                pass_definition = "; ".join(
                    f"{column} >= {overlap_cutoff:.1f} A" for column in centered
                )
                for protocol_index, protocol in enumerate(("vanilla", "masked")):
                    frame = derived[protocol]
                    pass_mask = frame["overlap_pass"].eq(1)
                    contact_mask = pd.to_numeric(
                        frame[f"shortest_ARG406-{partner}"], errors="raise"
                    ).le(4)
                    trajectories = trajectory_summary(frame, value, "mean")
                    seed_estimate, ci_low, ci_high, n_seeds = whole_seed_mean_ci(
                        trajectories,
                        value,
                        n_bootstrap=bootstrap,
                        random_seed=(
                            BASE_SEED + 6400 + 20 * R406_SCREEN_PARTNERS.index(partner)
                            + protocol_index
                        ),
                    )
                    protocol_values[protocol] = {
                        "pooled": float((pass_mask & contact_mask).sum() / pass_mask.sum()),
                        "trajectory": float(trajectories[value].mean()),
                        "seed": seed_estimate,
                        "trajectories": trajectories,
                    }
                    reconciliation_rows.append({
                        "analysis": "coordinate_derived_R406_screen",
                        "partner": partner,
                        "protocol": protocol,
                        "pooled_passing_snapshots": int(pass_mask.sum()),
                        "pooled_contact_numerator": int((pass_mask & contact_mask).sum()),
                        "pooled_contact_denominator": int(pass_mask.sum()),
                        "pooled_conditional_fraction": protocol_values[protocol]["pooled"],
                        "conditional_trajectories": len(trajectories),
                        "total_trajectories": int(
                            frame[["seed", "model_number"]].drop_duplicates().shape[0]
                        ),
                        "per_trajectory_conditional_estimate": protocol_values[protocol]["trajectory"],
                        "contributing_seed_labels": n_seeds,
                        "seed_balanced_conditional_estimate": seed_estimate,
                        "whole_seed_CI_low": ci_low,
                        "whole_seed_CI_high": ci_high,
                        "bootstrap_replicates": bootstrap,
                        "resampling_unit": "whole recorded nominal seed label",
                        "exact_contact_field": f"shortest_ARG406-{partner}",
                        "contact_rule": "screen passes and contact distance <= 4.0 A",
                        "pass_definition": pass_definition,
                        "atom_identity_status": (
                            "source distances; closest atom names unavailable"
                        ),
                        "actual_random_seed_pairing_status": (
                            "not applicable within condition"
                        ),
                    })
                contrast = joint_seed_bootstrap_contrast(
                    protocol_values["vanilla"]["trajectories"],
                    protocol_values["masked"]["trajectories"],
                    value,
                    n_bootstrap=bootstrap,
                    random_seed=BASE_SEED + 6500 + R406_SCREEN_PARTNERS.index(partner),
                )
                reconciliation_rows.append({
                    "analysis": "coordinate_derived_R406_screen",
                    "partner": partner,
                    "protocol": "masked_minus_vanilla",
                    "pooled_passing_snapshots": np.nan,
                    "pooled_contact_numerator": np.nan,
                    "pooled_contact_denominator": np.nan,
                    "pooled_conditional_fraction": (
                        protocol_values["masked"]["pooled"]
                        - protocol_values["vanilla"]["pooled"]
                    ),
                    "conditional_trajectories": np.nan,
                    "total_trajectories": np.nan,
                    "per_trajectory_conditional_estimate": (
                        protocol_values["masked"]["trajectory"]
                        - protocol_values["vanilla"]["trajectory"]
                    ),
                    "contributing_seed_labels": int(
                        contrast["n_nominal_seed_labels_union"]
                    ),
                    "seed_balanced_conditional_estimate": contrast["estimate_B_minus_A"],
                    "whole_seed_CI_low": contrast["CI_low_B_minus_A"],
                    "whole_seed_CI_high": contrast["CI_high_B_minus_A"],
                    "bootstrap_replicates": bootstrap,
                    "resampling_unit": "whole recorded nominal seed label",
                    "exact_contact_field": f"shortest_ARG406-{partner}",
                    "contact_rule": "screen passes and contact distance <= 4.0 A",
                    "pass_definition": pass_definition,
                    "atom_identity_status": (
                        "source distances; closest atom names unavailable"
                    ),
                    "actual_random_seed_pairing_status": (
                        "matching numeric labels are nominal design keys; common RNG not verified"
                    ),
                })
    save(pd.DataFrame(rows), output / "g406r_overlap_and_contacts.csv")
    save(pd.DataFrame(sensitivity), output / "g406r_overlap_threshold_sensitivity.csv")
    save(pd.DataFrame(yields), output / "usable_target_geometry_yields_g406r.csv")
    save(pd.DataFrame(denominator_rows), output / "g406r_screen_denominators.csv")
    save(pd.DataFrame(column_audit_rows), output / "g406r_overlap_screen_columns.csv")
    save(
        pd.DataFrame(reconciliation_rows),
        output / "g406r_conditional_contact_estimand_reconciliation.csv",
    )


def run_summary(output: Path, mode: str, bootstrap: int, permutations: int) -> None:
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    report = {
        "analysis": "paired_seed_v2", "mode": mode, "date": "2026-08-24",
        "git_starting_commit": commit, "python": platform.python_version(),
        "base_random_seed": BASE_SEED, "bootstrap_replicates": bootstrap,
        "permutation_replicates": permutations,
        "primary_estimand": "marginal QC-qualified contrast with joint nominal-seed resampling",
        "pairing_provenance": "actual RNG seed values unavailable; pairing keys are nominal seed labels",
        "source_hashes": "audit/SOURCE_INPUT_HASHES.tsv",
        "source_mutation_policy": "read-only inputs; all results written under paired_seed_v2",
    }
    (output / "run_summary.json").write_text(json.dumps(report, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "analysis/statistics_revision/paired_seed_v2")
    parser.add_argument("--mode", choices=("exploratory", "publication"), default="publication")
    args = parser.parse_args()
    bootstrap, permutations = ((250, 499) if args.mode == "exploratory" else (2000, 9999))
    output = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    master = master_cohort(output)
    make_inventory(master, output)
    retention(master, output, bootstrap)
    l403a(output, bootstrap, permutations)
    f412l(output, bootstrap)
    nav15(output, bootstrap, permutations)
    g402s(output, bootstrap, permutations)
    g406r(output, bootstrap)
    run_summary(output, args.mode, bootstrap, permutations)
    print(f"paired-seed analysis completed: {output}")


if __name__ == "__main__":
    main()
