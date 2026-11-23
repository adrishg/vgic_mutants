#!/usr/bin/env python3
"""Manuscript Fig. 7C/7D panels for the 50% Cav1.2 Timothy-syndrome ensembles.

7C: shortest side-chain non-hydrogen distance from R406 to T1531, D1533, D1528
    and F1527 in G406R vanilla (left half) and masked (right half) ensembles.
    Both protocols use the side-chain tables of 2026-11-23
    (``shared/side_chain_supplements.py``); the original vanilla export used
    terminal atoms only and the masked export every non-hydrogen atom.
7D: S402 Cα distance to I1523, V1520, M1524 and D1528 in G402S vanilla and
    masked ensembles (Cα, because glycine has no side chain to compare).

Inputs are the All-OK-3 filtered tables (vanilla 2026 exports and the
registered 50% masked cohorts). WT experimental markers (8WE6, 8HLP, 8FD7) come
from ``cav12_experimental_shortest_distances`` and are drawn in 7D only: the
experimental structures carry glycine at 406, which has no side-chain atom.
Outputs: ``cav12/dataDistances/analysis/figures/fig7{c,d}_*.{png,pdf}`` and the
plotted values in ``cav12/dataDistances/analysis/tables/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shared import data_access  # noqa: E402
from shared.dataset_selection import filtered_masked_distance_path  # noqa: E402
from shared.mutation_contact_analysis import cav12_experimental_shortest_distances  # noqa: E402
from shared.plotting import CAV12_PALETTE, experimental_reference_style  # noqa: E402
from shared.side_chain_supplements import apply_side_chain_supplement  # noqa: E402

VANILLA = {
    "G406R": "cav12/dataDistances/26-07-25_Cav1.2_g406r_vanillaAF2_distances_all_ok_rmsd_3A.csv",
    "G402S": "cav12/dataDistances/26-02-10_Cav12_g402s_vanillaAF2_distances_all_ok_rmsd_3A.csv",
}
PANELS = {
    "G406R": dict(
        site=406, residue="ARG406", label="R406",
        partners=["THR1531", "ASP1533", "ASP1528", "PHE1527"],
        column="shortest_ARG406-{partner}", exp_column="Shortest distance (Å)",
        ylabel="Shortest side-chain heavy-atom distance (Å)",
        stem="fig7c_g406r_r406_shortest_distances", experimental=False,
        guides=((4.0, "--", "4 Å proximity guide"), (2.0, ":", "2 Å overlap flag")),
    ),
    "G402S": dict(
        site=402, residue="SER402", label="S402",
        partners=["ILE1523", "VAL1520", "MET1524", "ASP1528"],
        column="CA_SER402_CA-{partner}_CA", exp_column="Cα distance (Å)",
        ylabel="Cα distance (Å)", stem="fig7d_g402s_s402_ca_distances", experimental=True,
        guides=(),
    ),
}
ONE = {"THR": "T", "ASP": "D", "PHE": "F", "ILE": "I", "VAL": "V", "MET": "M"}


def short(partner: str) -> str:
    return ONE[partner[:3]] + partner[3:]


def load(mutant: str, protocol: str) -> pd.DataFrame:
    if protocol == "vanilla":
        requested = VANILLA[mutant]
    else:
        requested = filtered_masked_distance_path(ROOT, "cav12", mutant)
    path = data_access.resolve(requested)
    frame = pd.read_csv(path, low_memory=False)
    extra = Path(str(path).replace("_all_ok_rmsd_3A", "_extra_all_ok_rmsd_3A"))
    if extra.is_file():
        companion = pd.read_csv(extra, low_memory=False)
        companion = companion[[c for c in companion if c == "pdb_file" or c not in frame]]
        frame = frame.merge(companion, on="pdb_file", how="left")
    return apply_side_chain_supplement(frame, Path(str(requested)).name)


def build(mutant: str, spec: dict) -> pd.DataFrame:
    rows = []
    for protocol in ("vanilla", "masked"):
        frame = load(mutant, protocol)
        for partner in spec["partners"]:
            values = pd.to_numeric(frame[spec["column"].format(partner=partner)], errors="coerce").dropna()
            rows.append(pd.DataFrame({
                "Mutant": mutant, "Protocol": protocol, "Partner": short(partner),
                "Distance (Å)": values.to_numpy(),
            }))
    table = pd.concat(rows, ignore_index=True)
    if spec["experimental"]:
        experimental = cav12_experimental_shortest_distances(ROOT, spec["site"], spec["partners"])
        experimental["Partner"] = experimental["Partner"].map(short)
    else:
        experimental = pd.DataFrame(columns=["Structure", "Partner", spec["exp_column"]])
    order = [short(p) for p in spec["partners"]]
    colors = {"vanilla": CAV12_PALETTE[f"{mutant}_VAN"], "masked": CAV12_PALETTE[f"{mutant}_HM"]}

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    sns.violinplot(data=table, x="Partner", y="Distance (Å)", hue="Protocol", order=order,
                   hue_order=["vanilla", "masked"], split=True, inner="quart", cut=0,
                   density_norm="width", palette=colors, linewidth=0.8, ax=ax)
    if ax.get_legend() is not None:
        ax.get_legend().remove()
    handles = [Patch(facecolor=colors[p], edgecolor="#444", label=f"{mutant} | {p}") for p in ("vanilla", "masked")]
    for k, (reference, part) in enumerate(experimental.groupby("Structure", sort=False)):
        style = experimental_reference_style(reference, k)
        for i, partner in enumerate(order):
            values = part.loc[part["Partner"].eq(partner), spec["exp_column"]]
            ax.scatter(np.full(len(values), i + 0.28), values, marker=style["marker"], s=60,
                       facecolor="white", edgecolor=style["color"], linewidth=1.3, zorder=5)
        handles.append(Line2D([0], [0], linestyle="", marker=style["marker"], markerfacecolor="white",
                              markeredgecolor=style["color"], markersize=8, label=f"{reference} | WT"))
    for y, style, label in spec["guides"]:
        ax.axhline(y, color="#71637A", linestyle=style, linewidth=0.9)
        handles.append(Line2D([0], [0], color="#71637A", linestyle=style, label=label))
    ax.set_xlabel(f"Partner of {spec['label']}")
    ax.set_ylabel(spec["ylabel"])
    ax.set_title(f"Cav1.2 | {mutant} | {spec['label']}-centered distances", fontweight="bold")
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    sns.despine(ax=ax)
    fig.tight_layout()

    figures = ROOT / "cav12/dataDistances/analysis/figures"
    tables = ROOT / "cav12/dataDistances/analysis/tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(figures / f"{spec['stem']}.{ext}", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    summary = table.groupby(["Protocol", "Partner"])["Distance (Å)"].agg(
        n="size", median="median",
        within_4A=lambda s: float((s <= 4).mean()), below_2A=lambda s: float((s < 2).mean()),
    ).reset_index()
    summary.to_csv(tables / f"{spec['stem']}_summary.csv", index=False)
    experimental.to_csv(tables / f"{spec['stem']}_experimental.csv", index=False)
    return summary


def main() -> None:
    for mutant, spec in PANELS.items():
        print(mutant)
        print(build(mutant, spec).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
