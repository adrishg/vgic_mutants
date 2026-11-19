#!/usr/bin/env python3
"""Supplementary Fig. S6: seed-block W1 and normalized-W1 heatmaps (WT vs mutant).

Reads ``seed_block/full_panel/all_distance_seed_block_comparison_summary.csv``
(written by ``run_seed_block_distance_panel.py``) and plots, for each mutant,
the median seed-balanced W1 across shared coordinates (top) and the median W1
normalized by the pooled weighted IQR (bottom). Columns separate the primary
vanilla sequence contrast from the protocol-specific masked contrast. Distances
are reported to one decimal place, normalized W1 to two.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
FULL = ROOT / "analysis/statistics_revision/seed_block/full_panel"
ROWS = [
    ("Cav1.2 G402S", "cav_wt_g402s_van", "cav_wt_g402s_mask"),
    ("Cav1.2 G406R", "cav_wt_g406r_van", "cav_wt_g406r_mask"),
    ("Kv2.1 L403A", "kv_wt_l403a_van", "kv_wt_l403a_mask"),
    ("Kv2.1 F412L", "kv_wt_f412l_van", "kv_wt_f412l_mask"),
    ("Nav1.5 QQQ", "nav_wt_qqq_van", "nav_wt_qqq_maskv2"),
]
COLUMNS = ["Vanilla\n(primary)", "Masked\n(protocol specific)"]


def main() -> None:
    summary = pd.read_csv(FULL / "all_distance_seed_block_comparison_summary.csv").set_index("comparison_id")
    w1 = np.array([[summary.loc[v, "median_seed_balanced_W1_A"], summary.loc[m, "median_seed_balanced_W1_A"]] for _, v, m in ROWS])
    nw1 = np.array([[summary.loc[v, "median_normalized_W1"], summary.loc[m, "median_normalized_W1"]] for _, v, m in ROWS])
    counts = np.array([[summary.loc[v, "number_distances"], summary.loc[m, "number_distances"]] for _, v, m in ROWS])

    fig, axes = plt.subplots(2, 1, figsize=(5.2, 7.6))
    for ax, values, title, fmt, cmap in (
        (axes[0], w1, "Median seed-balanced W1 (Å)", "{:.1f}", "Greens"),
        (axes[1], nw1, "Median W1 / pooled weighted IQR", "{:.2f}", "Purples"),
    ):
        image = ax.imshow(values, cmap=cmap, aspect="auto", vmin=0)
        for (i, j), value in np.ndenumerate(values):
            dark = value > 0.6 * np.nanmax(values)
            ax.text(j, i, fmt.format(value), ha="center", va="center", fontsize=10,
                    color="white" if dark else "#222")
        ax.set_xticks(range(2), COLUMNS)
        ax.set_yticks(range(len(ROWS)), [r[0] for r in ROWS])
        ax.set_title(title, fontweight="bold")
        fig.colorbar(image, ax=ax, fraction=0.05, pad=0.03)
    fig.tight_layout()
    out = ROOT / "analysis/statistics_revision/figures"
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"Figure_S6_seed_block_W1_heatmaps.{ext}", dpi=600, bbox_inches="tight", facecolor="white")
    table = pd.DataFrame({
        "mutant": [r[0] for r in ROWS] * 2,
        "contrast": ["vanilla"] * len(ROWS) + ["masked"] * len(ROWS),
        "comparison_id": [r[1] for r in ROWS] + [r[2] for r in ROWS],
        "median_seed_balanced_W1_A": np.r_[w1[:, 0], w1[:, 1]],
        "median_normalized_W1": np.r_[nw1[:, 0], nw1[:, 1]],
        "number_shared_coordinates": np.r_[counts[:, 0], counts[:, 1]],
    })
    table.to_csv(ROOT / "analysis/statistics_revision/tables/Figure_S6_seed_block_W1_heatmaps.csv", index=False)
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
