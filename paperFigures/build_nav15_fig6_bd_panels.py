#!/usr/bin/env python3
"""Manuscript Fig. 6B and 6D panels: IFM/QQQ side-chain distances to the receptor.

6B: shortest side-chain non-hydrogen distance from the central motif residue
    (F1486 in WT, Q1486 in QQQ; model 1170) to N1659 and to N1765 (model 1343
    and 1449) in the WT and QQQ vanilla ensembles, with the WT experimental
    structures.
6D: the central residue-N1659 distance in WT and QQQ under the vanilla,
    standard-mask (v1) and restricted-mask (v2) protocols.

Every ensemble is its final-QC cohort (the rows of the All-OK-3 distance table)
and every distance, model or experimental, uses the same atoms: all
non-hydrogen side-chain atoms (``shared/side_chain_supplements.py``).
Outputs: ``nav15/dataDistances/analysis/figures/fig6{b,d}_*.{png,pdf}`` and a
summary table in ``nav15/dataDistances/analysis/tables/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from analysis.statistics_revision.scripts.run_seed_block_revision import (  # noqa: E402
    FINAL_DISTANCE_PATHS,
    read_csv,
)
from shared.plotting import NAV15_PALETTE, experimental_reference_style  # noqa: E402
from shared.side_chain_supplements import (  # noqa: E402
    shortest_side_chain_distance,
    side_chain_atoms,
)

ENSEMBLES = {
    "WT | vanilla": ("nav15|WT|vanilla", "PHE", NAV15_PALETTE["WT_VAN"]),
    "QQQ | vanilla": ("nav15|QQQ|vanilla", "GLN", NAV15_PALETTE["QQQ_VAN"]),
    "WT | standard mask": ("nav15|WT|masked", "PHE", NAV15_PALETTE["WT_HM"]),
    "QQQ | standard mask": ("nav15|QQQ|masked", "GLN", NAV15_PALETTE["QQQ_MASKED"]),
    "WT | restricted mask": ("nav15|WT|masked_v2", "PHE", NAV15_PALETTE["WT_MASKED_V2"]),
    "QQQ | restricted mask": ("nav15|QQQ|masked_v2", "GLN", NAV15_PALETTE["QQQ_MASKED_V2"]),
}
# Human 8VYJ/8VYK/7DTC use canonical numbering; rat 6UZ3/8T6L are shifted by +2.
EXPERIMENTAL = {"8VYJ": 0, "8VYK": 0, "7DTC": 0, "8T6L": 2, "6UZ3": 2}
N1659, N1765 = "Central residue–N1659", "Central residue–N1765"
YLABEL = "Central F/Q1486–N1659\nshortest side-chain distance (Å)"


def load() -> pd.DataFrame:
    frames = []
    for label, (key, residue, _) in ENSEMBLES.items():
        table = read_csv(ROOT / FINAL_DISTANCE_PATHS[key])
        frames.append(pd.DataFrame({
            "Ensemble": label,
            N1659: pd.to_numeric(table[f"shortest_{residue}1170-ASN1343"]),
            N1765: pd.to_numeric(table[f"shortest_{residue}1170-ASN1449"]),
        }))
    return pd.concat(frames, ignore_index=True)


def experimental() -> pd.DataFrame:
    rows = []
    for pdb_id, offset in EXPERIMENTAL.items():
        residues = side_chain_atoms(ROOT / f"nav15/experimental/{pdb_id}.pdb")
        motif = ("PHE", 1486 + offset)
        rows.append({
            "Structure": pdb_id,
            N1659: shortest_side_chain_distance(residues, motif, ("ASN", 1659 + offset)),
            N1765: shortest_side_chain_distance(residues, motif, ("ASN", 1765 + offset)),
        })
    return pd.DataFrame(rows)


def reference_handles(references: pd.DataFrame) -> list:
    handles = []
    for k, pdb_id in enumerate(references["Structure"]):
        style = experimental_reference_style(pdb_id, k)
        handles.append(Line2D([0], [0], linestyle="", marker=style["marker"], markerfacecolor="white",
                              markeredgecolor=style["color"], markeredgewidth=1.5, markersize=8,
                              label=f"{pdb_id} | WT IFM reference"))
    return handles


def panel_b(data: pd.DataFrame, references: pd.DataFrame, figures: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    handles = []
    for label in ("WT | vanilla", "QQQ | vanilla"):
        part = data[data["Ensemble"].eq(label)]
        color = ENSEMBLES[label][2]
        ax.scatter(part[N1659], part[N1765], s=9, color=color, alpha=0.45, linewidths=0, rasterized=True)
        handles.append(Line2D([0], [0], linestyle="", marker="o", color=color, markersize=7, label=label))
    for k, row in references.iterrows():
        style = experimental_reference_style(row["Structure"], k)
        ax.scatter(row[N1659], row[N1765], marker=style["marker"], s=70, facecolor="white",
                   edgecolor=style["color"], linewidth=1.5, zorder=5)
    ax.set_xlabel("Central F/Q1486–N1659\nshortest side-chain distance (Å)")
    ax.set_ylabel("Central F/Q1486–N1765\nshortest side-chain distance (Å)")
    ax.grid(color="#EEE9F2", linewidth=0.45, alpha=0.65)
    ax.legend(handles=handles + reference_handles(references), frameon=False, fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1))
    sns.despine(ax=ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figures / f"fig6b_nav15_ifm_side_chain_vanilla.{ext}", dpi=600,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_d(data: pd.DataFrame, references: pd.DataFrame, figures: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.9), sharey=True)
    for ax, protocol in zip(axes, ("vanilla", "standard mask", "restricted mask")):
        order = [f"WT | {protocol}", f"QQQ | {protocol}"]
        part = data[data["Ensemble"].isin(order)]
        sns.violinplot(data=part, x="Ensemble", y=N1659, order=order, hue="Ensemble", hue_order=order,
                       palette={label: ENSEMBLES[label][2] for label in order}, inner="quart", cut=0,
                       density_norm="width", linewidth=0.8, legend=False, ax=ax)
        for k, row in references.iterrows():
            style = experimental_reference_style(row["Structure"], k)
            ax.scatter(-0.25 + 0.125 * k, row[N1659], marker=style["marker"], s=38, facecolor="white",
                       edgecolor=style["color"], linewidth=1.2, zorder=5)
        ax.set_xticks(range(2))
        ax.set_xticklabels(["WT", "QQQ"])
        ax.set_xlabel("")
        ax.set_title(protocol.capitalize() if protocol != "vanilla" else "Vanilla", fontsize=11)
        ax.set_ylabel(YLABEL if ax is axes[0] else "")
        ax.grid(axis="y", color="#EEE9F2", linewidth=0.45, alpha=0.65)
        sns.despine(ax=ax)
    axes[-1].legend(handles=reference_handles(references), frameon=False, fontsize=8,
                    loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figures / f"fig6d_nav15_ifm_side_chain_three_protocols.{ext}", dpi=600,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    figures = ROOT / "nav15/dataDistances/analysis/figures"
    tables = ROOT / "nav15/dataDistances/analysis/tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    data, references = load(), experimental()
    panel_b(data, references, figures)
    panel_d(data, references, figures)
    summary = data.groupby("Ensemble", sort=False).agg(
        n=(N1659, "size"),
        N1659_median=(N1659, "median"), N1659_p05=(N1659, lambda s: s.quantile(0.05)),
        N1659_within_6A=(N1659, lambda s: float((s <= 6).mean())),
        N1765_median=(N1765, "median"),
    ).reset_index()
    summary.to_csv(tables / "fig6bd_nav15_ifm_side_chain_summary.csv", index=False)
    references.to_csv(tables / "fig6bd_nav15_ifm_side_chain_experimental.csv", index=False)
    print(summary.round(3).to_string(index=False))
    print(references.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
