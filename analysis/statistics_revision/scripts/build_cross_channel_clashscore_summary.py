#!/usr/bin/env python3
"""Build the general cross-channel global clashscore result."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from shared.plotting import CAV12_PALETTE, KV21_PALETTE, NAV15_PALETTE, apply_kv21_style


CHANNELS = {
    "cav12": ("Cav1.2", ("WT", "G402S", "G406R"), CAV12_PALETTE),
    "kv21": ("Kv2.1", ("WT", "L403A", "F412L"), KV21_PALETTE),
    "nav15": ("Nav1.5", ("WT", "QQQ"), NAV15_PALETTE),
}


def build_table() -> pd.DataFrame:
    rows = []
    for channel, (display, sequences, _) in CHANNELS.items():
        contrasts = pd.read_csv(
            ROOT / channel / "dataExtra/molprobity/tables/contrasts.csv"
        )
        selected = contrasts[
            contrasts.cohort.eq("repository_QC")
            & contrasts.metric.eq("clashscore")
        ]
        for sequence in sequences:
            row = selected[
                selected.A.eq(f"{sequence} | vanilla")
                # 50% Kv2.1/CaV1.2 cohorts are labelled ``masked v3``.
                & selected.B.isin([f"{sequence} | masked", f"{sequence} | masked v3"])
            ].iloc[0]
            rows.append({
                "channel": display,
                "sequence": sequence,
                "vanilla_seed_balanced_mean": row.seed_mean_A,
                "masked_seed_balanced_mean": row.seed_mean_B,
                "masked_minus_vanilla": row.seed_mean_difference,
                "CI_low": row.CI_low,
                "CI_high": row.CI_high,
                "cohort": "repository_QC",
                "metric": "MolProbity clashscore",
                "software": "Phenix MolProbity",
                "software_version": "not recorded in source table",
                "invocation_options": "not recorded in source table",
                "hydrogen_treatment": "not recorded in source table",
                "analyzed_cohort": (
                    f"{display} {sequence}: structures retained by the channel-specific "
                    "repository_QC manifest"
                ),
                "source_table": f"{channel}/dataExtra/molprobity/tables/contrasts.csv",
            })
    return pd.DataFrame(rows)


def plot(table: pd.DataFrame, output: Path) -> None:
    apply_kv21_style()
    fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    for axis, (channel, (_, _, palette)) in zip(axes, CHANNELS.items()):
        display = CHANNELS[channel][0]
        part = table[table.channel.eq(display)].reset_index(drop=True)
        x = range(len(part))
        width = .36
        vanilla = axis.bar(
            [value - width / 2 for value in x],
            part.vanilla_seed_balanced_mean, width,
            color=[palette[f"{sequence}_VAN"] for sequence in part.sequence],
            edgecolor="#495057", linewidth=.5, label="Vanilla",
        )
        masked = axis.bar(
            [value + width / 2 for value in x],
            part.masked_seed_balanced_mean, width,
            color=[palette[f"{sequence}_HM"] for sequence in part.sequence],
            edgecolor="#495057", linewidth=.5, label="Masked",
        )
        axis.bar_label(vanilla, fmt="%.2f", fontsize=7, padding=2)
        axis.bar_label(masked, fmt="%.2f", fontsize=7, padding=2)
        # Leave enough headroom that the legend cannot obscure the absolute
        # value labels for the high-clashscore Kv2.1 ensembles.
        axis.set(xticks=list(x), xticklabels=part.sequence, title=display, ylim=(0, 60))
        axis.legend(frameon=False, fontsize=8, loc="upper left")
    axes[0].set_ylabel("Seed-balanced mean clashscore")
    fig.suptitle("Global steric QC across channels", fontweight="bold")
    fig.text(
        .5, .02,
        "Absolute means use each channel's repository QC cohort. Clashscore is whole-model; "
        "it is distinct from the R406-centered <2 Å screen.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=(0, .08, 1, .93))
    fig.savefig(output / "cross_channel_clashscore_absolute_means.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(output / "cross_channel_clashscore_absolute_means.pdf",
                bbox_inches="tight")
    plt.close(fig)


def write_report(table: pd.DataFrame, path: Path) -> None:
    lines = [
        "# Cross-channel global steric QC",
        "",
        "This is a general result across Cav1.2, Kv2.1, and Nav1.5. It is not part of the Cav1.2 mutation-site result.",
        "",
        "The values below are seed-balanced mean MolProbity clashscores in each channel's repository QC cohort. Clashscore counts serious all-atom overlaps per 1,000 atoms, so absolute means are shown beside masked-minus-vanilla differences.",
        "",
        "| Channel | Sequence | Vanilla mean | Masked mean | Difference | 95% CI |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in table.itertuples():
        lines.append(
            f"| {row.channel} | {row.sequence} | {row.vanilla_seed_balanced_mean:.2f} | "
            f"{row.masked_seed_balanced_mean:.2f} | {row.masked_minus_vanilla:+.2f} | "
            f"[{row.CI_low:+.2f}, {row.CI_high:+.2f}] |"
        )
    lines += [
        "",
        "Masking was associated with higher global clashscore in most ensembles, with the largest change in Kv2.1 F412L. The L403A interval included zero. Absolute levels differed across channels, so the within-ensemble contrasts are the relevant protocol comparisons.",
        "",
        "The Cav1.2 R406-centered overlap screen is a separate local endpoint. Passing that screen means none of the selected R406-centered distances was below 2 Å; it does not validate the entire pocket. MolProbity clashscore summarizes serious overlaps across the whole model and does not identify the residue responsible.",
        "",
        "Source tables: `<channel>/dataExtra/molprobity/tables/contrasts.csv`.",
        "",
        "The source tables identify the program as Phenix MolProbity but do not record the software version, command-line/options, or whether hydrogens were added, retained, or optimized. Those fields are therefore reported as `not recorded` rather than reconstructed by assumption.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    table_dir = ROOT / "analysis/statistics_revision/tables"
    figure_dir = ROOT / "analysis/statistics_revision/figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    table = build_table()
    table.to_csv(table_dir / "cross_channel_clashscore_summary.csv", index=False)
    plot(table, figure_dir)
    write_report(table, ROOT / "docs/results/CROSS_CHANNEL_STERIC_QC.md")


if __name__ == "__main__":
    main()
