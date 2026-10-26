#!/usr/bin/env python3
"""Build publication figures, source tables, and methods metadata."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.statistics_revision.scripts.run_paired_seed_v2 import (
    G402_PARTNERS, MASK_IDS, trajectory_summary,
)
from analysis.statistics_revision.scripts.run_seed_block_revision import (
    FINAL_DISTANCE_PATHS, L403_COLUMNS, add_metadata, read_csv,
)
from scripts.ensemble_rmsf_analysis.materialize_authoritative_masks import _resolve_masks
from scripts.ensemble_rmsf_analysis.masks import compact_ranges, parse_ranges
from shared.plotting import CAV12_PALETTE, KV21_PALETTE, apply_kv21_style

OUT = ROOT / "analysis/statistics_revision/paired_seed_v2"
TABLES = ROOT / "analysis/statistics_revision/tables"
FIGURES = ROOT / "analysis/statistics_revision/figures"
SOURCE = OUT / "source_data"


def write(frame: pd.DataFrame, path: Path, *, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, sep=sep)


def figure4() -> None:
    threshold = float(pd.read_csv(OUT / "l403a_experimental_threshold_derivation.csv").derived_threshold_A.iloc[0])
    trajectory_rows, category_rows = [], []
    order = []
    for sequence in ("WT", "L403A"):
        for protocol in ("vanilla", "masked"):
            frame = add_metadata(read_csv(ROOT / FINAL_DISTANCE_PATHS[f"kv21|{sequence}|{protocol}"]))
            distances = frame[L403_COLUMNS].apply(pd.to_numeric, errors="raise")
            frame = frame.assign(
                maximum_E423_N179_CA_distance_A=distances.max(axis=1),
                elongated_interface_count=distances.ge(threshold).sum(axis=1),
            )
            label = f"{sequence} | {protocol}"
            order.append(label)
            summary = frame.groupby(["seed", "model_number"], as_index=False).agg(
                trajectory_median_maximum_distance_A=("maximum_E423_N179_CA_distance_A", "median"),
                retained_snapshots=("pdb_file", "size"),
            )
            summary.insert(0, "condition", label)
            summary["sequence"] = sequence
            summary["protocol"] = protocol
            summary["mask_id"] = MASK_IDS[("kv21", protocol)]
            summary["threshold_A"] = threshold
            trajectory_rows.append(summary)
            per_trajectory = (
                frame.groupby(["seed", "model_number", "elongated_interface_count"])
                .size().rename("snapshots").reset_index()
            )
            totals = per_trajectory.groupby(["seed", "model_number"]).snapshots.transform("sum")
            per_trajectory["trajectory_fraction"] = per_trajectory.snapshots / totals
            trajectories_index = frame[["seed", "model_number"]].drop_duplicates()
            complete_index = pd.MultiIndex.from_product(
                [
                    sorted(trajectories_index.seed.unique()),
                    sorted(trajectories_index.model_number.unique()),
                    range(5),
                ],
                names=["seed", "model_number", "elongated_interface_count"],
            )
            # Some seed-model combinations fail QC completely.  Average only
            # over retained trajectories, but include explicit zeroes for the
            # four categories that each retained trajectory did not occupy.
            retained_index = pd.MultiIndex.from_frame(trajectories_index)
            complete_index = complete_index[
                complete_index.droplevel("elongated_interface_count").isin(retained_index)
            ]
            completed = per_trajectory.set_index(
                ["seed", "model_number", "elongated_interface_count"]
            ).reindex(complete_index, fill_value=0).reset_index()
            balanced = completed.groupby("elongated_interface_count", as_index=False).trajectory_fraction.mean()
            balanced.insert(0, "condition", label)
            balanced["n_trajectories"] = frame[["seed", "model_number"]].drop_duplicates().shape[0]
            balanced["resampling_unit"] = "trajectory; each seed-model trajectory has equal weight"
            category_rows.append(balanced)
    trajectories = pd.concat(trajectory_rows, ignore_index=True)
    categories = pd.concat(category_rows, ignore_index=True)
    write(trajectories, SOURCE / "figure4_four_condition_trajectory_maxima.csv.gz")
    write(categories, SOURCE / "figure4_zero_to_four_elongated_interfaces.csv")

    apply_kv21_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)
    colors = [KV21_PALETTE["WT_VAN"], KV21_PALETTE["WT_HM"], KV21_PALETTE["L403A_VAN"], KV21_PALETTE["L403A_HM"]]
    rng = np.random.default_rng(20261006)
    for x, (label, color) in enumerate(zip(order, colors)):
        values = trajectories.loc[trajectories.condition.eq(label), "trajectory_median_maximum_distance_A"].to_numpy()
        parts = axes[0].violinplot(values, [x], widths=.75, showmedians=True, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(color); body.set_edgecolor("#4C5552"); body.set_alpha(.9)
        axes[0].scatter(x + rng.normal(0, .055, len(values)), values, s=5, alpha=.12, color="#25302B")
    axes[0].axhline(threshold, color="#7B6D86", ls="--", lw=1, label=f"Experimental cutoff ({threshold:.2f} Å)")
    axes[0].set(xticks=range(4), xticklabels=order, ylabel="Trajectory-median maximum E423–N179 Cα distance (Å)", title="A  Four-condition interface displacement")
    axes[0].tick_params(axis="x", rotation=25); axes[0].legend(frameon=False, fontsize=8)
    bottom = np.zeros(4)
    count_colors = ["#EEF2EE", "#C9DDCB", "#89B992", "#3E8C55", "#155F35"]
    for count, color in zip(range(5), count_colors):
        values = categories[categories.elongated_interface_count.eq(count)].set_index("condition").loc[order, "trajectory_fraction"].to_numpy()
        axes[1].bar(range(4), 100 * values, bottom=100 * bottom, color=color, edgecolor="white", label=str(count))
        bottom += values
    axes[1].set(xticks=range(4), xticklabels=order, ylabel="Trajectory-balanced snapshots (%)", title="B  Elongated interfaces per tetramer", ylim=(0, 100))
    axes[1].tick_params(axis="x", rotation=25); axes[1].legend(title="Count (0–4)", frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(.5, -.18))
    for ax in axes: ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Figure 4 | Kv2.1 WT and L403A under vanilla and targeted masking", fontweight="bold")
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"Figure4_Kv21_L403A_four_condition_revision.{suffix}", dpi=300 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)


def figure7() -> None:
    stats = pd.read_csv(OUT / "g402s_matched_ca_partner_statistics.csv")
    effects = stats[
        stats.contrast.eq("vanilla G402S - vanilla WT")
        & stats.estimand_id.eq("primary_joint_nominal_seed")
    ].copy()
    effect_col = "estimate_B_minus_A"
    low_col, high_col = "CI_low_B_minus_A", "CI_high_B_minus_A"
    effects["partner_label"] = effects.partner.str.replace(r"([A-Z]{3})(\d+)", r"\1 \2", regex=True)
    coordinate = pd.read_csv(OUT / "g406r_coordinate_screen_summary.csv")
    coordinate = coordinate[
        coordinate.atom_selection_rule.eq("side_chain_non_hydrogen_atoms")
    ].copy()
    screen = coordinate[coordinate.partner_residue.eq("screen_only")].copy()
    contacts = coordinate[
        coordinate.partner_residue.isin(["ASP1528", "ASP1533"])
    ].copy()

    apply_kv21_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5), constrained_layout=True)
    y = np.arange(len(effects))
    axes[0].errorbar(
        effects[effect_col], y,
        xerr=[effects[effect_col] - effects[low_col], effects[high_col] - effects[effect_col]],
        fmt="o", color=CAV12_PALETTE["G402S_HM"], ecolor="#7892B7", capsize=2,
    )
    axes[0].axvline(0, color="#666", lw=.8, ls="--")
    axes[0].set(yticks=y, yticklabels=effects.partner_label, xlabel="G402S − WT matched Cα distance (Å)", title="D  Position 402: vanilla sequence effect")
    axes[0].invert_yaxis()
    x = np.arange(2); protocols = ["vanilla", "masked"]
    values = [100 * float(screen.loc[screen.protocol.eq(p), "unconditional_pass_fraction"].iloc[0]) for p in protocols]
    bars = axes[1].bar(x, values, color=[CAV12_PALETTE["G406R_VAN"], CAV12_PALETTE["G406R_HM"]], width=.58)
    axes[1].bar_label(bars, fmt="%.1f%%", padding=3)
    axes[1].set(xticks=x, xticklabels=protocols, ylabel="All retained snapshots passing (%)", ylim=(0, 100), title="E  Passed explicit R406 side-chain overlap screen")
    for xi, protocol in enumerate(protocols):
        rows = contacts[contacts.protocol.eq(protocol)].set_index("partner_residue")
        text = "\n".join(
            f"{partner.replace('ASP', 'D')}: {100*rows.loc[partner, 'conditional_contact_fraction']:.1f}% "
            f"({int(rows.loc[partner, 'conditional_contact_numerator'])}/{int(rows.loc[partner, 'conditional_contact_denominator'])})"
            for partner in ("ASP1528", "ASP1533")
        )
        axes[1].text(xi, values[xi] + 10, text, ha="center", va="bottom", fontsize=8)
    axes[1].text(.5, -0.16, "All 19 distances use explicit side-chain non-hydrogen atoms. D1528/D1533 contacts use ≤4 Å and are conditional on passing.", transform=axes[1].transAxes, ha="center", fontsize=8)
    for ax in axes: ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Figure 7 | CaV1.2 mutation-site analyses", fontweight="bold")
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"Figure7_Cav12_matched_G402S_R406_screen_revision.{suffix}", dpi=300 if suffix == "png" else None)
    plt.close(fig)
    write(effects, SOURCE / "figure7D_G402S_matched_CA_effects.csv")
    write(pd.concat([screen, contacts], ignore_index=True), SOURCE / "figure7E_R406_screen_and_contacts.csv")


def cohort_manifest() -> None:
    master = pd.read_csv(OUT / "master_structure_cohort.csv")
    mask_lookup = {
        ("kv21", "WT", "masked"): "kv21_common", ("kv21", "L403A", "masked"): "kv21_common", ("kv21", "F412L", "masked"): "kv21_common",
        ("cav12", "WT", "masked"): "cav12_wt_common", ("cav12", "G402S", "masked"): "cav12_g402s", ("cav12", "G406R", "masked"): "cav12_g406r",
        ("nav15", "WT", "masked"): "nav15_standard", ("nav15", "WT", "masked_v2"): "nav15_v2", ("nav15", "WT", "masked_v2_noIFM"): "nav15_v2_noIFM",
        ("nav15", "QQQ", "masked"): "nav15_standard_plus_IFM", ("nav15", "QQQ", "masked_v2"): "nav15_v2",
    }
    master["condition"] = master.sequence + " | " + master.protocol
    master["mask"] = ["unmasked" if p == "vanilla" else mask_lookup.get((c, s, p), "unresolved") for c, s, p in zip(master.channel, master.sequence, master.protocol)]
    master["qc_status"] = np.where(master.distance_final_cohort, "included", "excluded")
    master["inclusion_exclusion_reason"] = master.exclusion_reason
    analyses = []
    specs = {
        "Figure4_L403A": master.channel.eq("kv21") & master.sequence.isin(["WT", "L403A"]),
        "Figure5_F412L": master.channel.eq("kv21") & master.sequence.isin(["WT", "F412L"]),
        "Figure6_Nav15": master.channel.eq("nav15"),
        "Figure7_G402S": master.channel.eq("cav12") & master.sequence.isin(["WT", "G402S"]),
        "Figure7_R406": master.channel.eq("cav12") & master.sequence.eq("G406R"),
        "cross_channel_MolProbity": pd.Series(True, index=master.index),
    }
    keep = ["analysis_id", "channel", "condition", "sequence", "protocol", "mask", "seed", "model_number", "recycle_number", "pdb_basename", "qc_status", "distance_final_cohort", "inclusion_exclusion_reason"]
    for analysis_id, selector in specs.items():
        part = master.loc[selector].copy(); part["analysis_id"] = analysis_id; analyses.append(part[keep])
    write(pd.concat(analyses, ignore_index=True), OUT / "analysis_cohort_manifest.tsv.gz", sep="\t")


def mask_tables() -> None:
    path = ROOT / "scripts/ensemble_rmsf_analysis/config/authoritative_mask_definitions.yaml"
    config = yaml.safe_load(path.read_text())
    masks = _resolve_masks(config["mask_definitions"])
    rows = []
    for channel, datasets in config["datasets"].items():
        for dataset, meta in datasets.items():
            name = meta["mask"]; positions = masks[name]
            rows.append({
                "channel": channel, "dataset": dataset, "mask_definition": name,
                "retained_directly_masked_columns": len(positions),
                "resolved_ranges": compact_ranges(positions),
                "definition_validation": "passed expected-count and controlled-difference checks",
                "production_A3M_equivalence": "not verifiable: production A3Ms are absent from repository",
            })
    write(pd.DataFrame(rows), TABLES / "mask_definition_equivalence.tsv", sep="\t")
    noifm = sorted(masks["nav15_v2_noIFM"])
    write(pd.DataFrame({"retained_masked_raw_query_column_1_based": noifm}), TABLES / "nav15_v2_noIFM_exact_retained_columns.tsv", sep="\t")
    assert set(noifm) == masks["nav15_v2"] - parse_ranges("1164-1176")


def figure6_table() -> None:
    rows = [
        ("6A", "WT and QQQ structural examples", "as labeled", "none", "structure depiction; no quantitative reference marker"),
        ("6B", "WT versus QQQ vanilla", "unmasked", "6UZ3, 8VYJ, 8VYK, 7DTC, 8T6L where coordinates exist", "sequence contrast"),
        ("6C", "WT vanilla versus masked", "nav15_standard", "none", "masking effect; original WT mask lacks gate columns"),
        ("6D", "WT nav15_standard; QQQ nav15_standard_plus_IFM", "condition-specific", "same canonical marker map", "protocol context, not a common mask"),
        ("6E", "QQQ vanilla", "unmasked", "none", "motif-gate correlation"),
        ("6F", "QQQ vanilla versus masked", "nav15_standard_plus_IFM", "7FBS gate marker allowed; 7FBS motif marker excluded", "7FBS linker/motif coordinates unresolved"),
    ]
    table = pd.DataFrame(rows, columns=["panel", "condition", "mask", "experimental_reference_markers", "note"])
    write(table, TABLES / "Figure6_mask_and_reference_markers.tsv", sep="\t")


def methods_manifest() -> None:
    rows = [
        ("Figure4", "Figure4_Kv21_L403A_four_condition_revision", 0, "descriptive equal-trajectory weighting", "seed and AF2 model parsed from filename"),
        ("Figure4 source", "source_data/figure4_four_condition_trajectory_maxima.csv.gz", 0, "one median per retained seed-model trajectory", "QC-empty trajectories are absent and counted in the cohort manifest"),
        ("Figure4 0-4 source", "source_data/figure4_zero_to_four_elongated_interfaces.csv", 0, "equal weight per retained seed-model trajectory", "five mutually exclusive categories sum to one per condition"),
        ("Figure4 binary interaction", "kv21_interactions_paired_seed.csv", 2000, "complete four-cell nominal seed label", "common RNG pairing is not independently verified"),
        ("Figure5", "Figure5_F412L_mutation_site_cleaned", 0, "retained snapshots shown descriptively", "WT-mutant and masking effects remain separately labeled"),
        ("Figure5 source", "source_data/figure5_F412L_mutation_site_distances.csv.gz", 0, "retained snapshots", "chain identity is retained in partner labels"),
        ("Figure6", "Figure6_mask_and_reference_markers.tsv", 0, "deterministic panel metadata", "7FBS motif marker excluded because coordinates are unresolved"),
        ("Figure7D", "g402s_matched_ca_partner_statistics.csv", 2000, "whole nominal seed label", "nominal labels do not prove common RNG pairing"),
        ("Figure7D source", "source_data/g402s_matched_ca_1519_1535.csv.gz", 0, "retained QC snapshots with seed/model/recycle identifiers", "all 17 matched C-alpha pairs retained"),
        ("Figure7E", "g406r_coordinate_screen_statistics.csv", 2000, "whole nominal seed label", "explicit side-chain non-hydrogen atom selection; nominal labels do not prove common RNG pairing"),
        ("Figure7E denominators", "g406r_coordinate_screen_summary.csv", 0, "descriptive retained snapshots/trajectories", "all coordinate-derived numerators and denominators printed"),
        ("Figure7E atom pairs", "source_data/g406r_coordinate_heavy_atom_pair_distances.csv.gz", 0, "one row per snapshot x atom rule x residue pair", "winning atom names stored for all 19 pairs"),
        ("Historical R406 sensitivity", "g406r_conditional_contact_estimand_reconciliation.csv", 2000, "whole nominal seed label", "historical atom identities unavailable; pooled and hierarchical estimands shown separately"),
        ("Cross-channel steric QC", "cross_channel_clashscore_summary.csv", 3000, "whole nominal seed label", "channel repository_QC cohorts"),
        ("Supplementary Figure S8", "Figure_S8_WT_First100_split_sampling_and_retention", 0, "20 nominal labels x five AF2 models; seed-level breadth", "filename is legacy; visible labels state 20 seeds / 100 model-seeds"),
        ("Supplementary Figure S8 source", "docs/tables/supplementary_figure_s8/*.csv", 0, "seed summaries and retained model-seed trajectories", "no bootstrap interval displayed in this descriptive figure"),
        ("Cohort manifest", "analysis_cohort_manifest.tsv.gz", 0, "one row per analysis x structure snapshot", "QC inclusion and exclusion retained"),
        ("Mask definitions", "mask_definition_equivalence.tsv", 0, "deterministic mask-definition comparison", "production A3Ms unavailable"),
    ]
    write(pd.DataFrame(rows, columns=["table_or_figure", "source_artifact", "bootstrap_replicates", "resampling_unit", "pairing_or_scope_caveat"]), OUT / "OUTPUT_METHODS.tsv", sep="\t")


def supplementary_manifest() -> None:
    rows = []
    for figure_dir in sorted((ROOT / "docs/figures").glob("supplementary_figure_s*")):
        supplement = figure_dir.name.replace("supplementary_figure_", "").upper()
        table_dir = ROOT / "docs/tables" / figure_dir.name
        figures = sorted(path.relative_to(ROOT).as_posix() for path in figure_dir.iterdir() if path.suffix.lower() in {".png", ".pdf"})
        sources = sorted(path.relative_to(ROOT).as_posix() for path in table_dir.iterdir()) if table_dir.is_dir() else []
        rows.append({
            "supplement": supplement,
            "figure_files": ";".join(figures),
            "source_data_files": ";".join(sources),
            "figure_present": bool(figures),
            "source_data_present": bool(sources),
            "status": "complete" if figures and sources else "missing figure or source data",
        })
    write(pd.DataFrame(rows), OUT / "SUPPLEMENTARY_OUTPUTS.tsv", sep="\t")


def main() -> None:
    for directory in (OUT, TABLES, FIGURES, SOURCE): directory.mkdir(parents=True, exist_ok=True)
    figure4(); figure7(); cohort_manifest(); mask_tables(); figure6_table(); methods_manifest(); supplementary_manifest()
    print("publication outputs completed")


if __name__ == "__main__":
    main()
