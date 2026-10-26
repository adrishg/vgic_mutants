"""Atom-matched backbone panels and explicitly separate glycine references."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns

from shared.plotting import CAV12_PALETTE, RMSD_REFERENCE_STYLES, format_distance_alias


def matched_ca_table(comparisons, site, mutant_residue, partners):
    """Require the same Cα endpoints for WT and mutant; never use shortest_*.

    Keep snapshot identifiers and source columns for auditing. Missing columns
    or nonfinite values fail explicitly rather than changing the plotted cohort.
    """
    records = []
    for protocol, (wt, mutant) in comparisons.items():
        for sequence, residue, frame in (
            ("WT", f"GLY{site}", wt),
            (mutant_residue, mutant_residue, mutant),
        ):
            for partner in partners:
                column = f"CA_{residue}_CA-{partner}_CA"
                if column not in frame:
                    raise KeyError(f"Missing atom-matched Cα column: {column}")
                values = pd.to_numeric(frame[column], errors="raise")
                if not np.isfinite(values).all():
                    raise ValueError(f"Nonfinite atom-matched Cα values: {column}")
                records.append(pd.DataFrame({
                    "Protocol": protocol, "Sequence": sequence,
                    "Partner": partner, "Cα distance (Å)": values.to_numpy(),
                    "pdb_file": frame["pdb_file"].to_numpy(),
                    "Source column": column,
                    "Selection": frame.attrs.get("actual_selection", "unspecified"),
                    "Source file": frame.attrs.get("source_path", "unspecified"),
                }))
    return pd.concat(records, ignore_index=True)


def plot_matched_ca(comparisons, site, mutant_residue, mutant_label, partners,
                    experimental_distances):
    """Compare matched backbone distances separately within each protocol."""
    table = matched_ca_table(comparisons, site, mutant_residue, partners)
    table["Sequence"] = table["Sequence"].replace({mutant_residue: mutant_label})
    labels = {p: format_distance_alias(p) for p in partners}
    table["Partner label"] = table.Partner.map(labels)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    for ax, protocol in zip(axes, ("vanilla", "masked")):
        suffix = "VAN" if protocol == "vanilla" else "HM"
        sns.violinplot(
            data=table[table.Protocol.eq(protocol)], x="Partner label",
            y="Cα distance (Å)", hue="Sequence", hue_order=["WT", mutant_label],
            order=list(labels.values()), split=True, cut=0, inner="quartile",
            palette={s: CAV12_PALETTE[f"{s}_{suffix}"] for s in ["WT", mutant_label]},
            linewidth=.7, ax=ax,
        )
        for index, reference in enumerate(("8WE6", "8HLP", "8FD7")):
            style = RMSD_REFERENCE_STYLES[reference]
            for _, row in experimental_distances.query("Structure == @reference").iterrows():
                if row.Partner in partners:
                    ax.scatter(partners.index(row.Partner) + (index - 1) * .12,
                               row["Cα distance (Å)"], s=32, marker=style["marker"],
                               facecolor="white", edgecolor=style["color"], zorder=5)
        ax.set(title=protocol.capitalize(), xlabel="Partner Cα", ylabel="Cα–Cα distance (Å)")
        ax.legend(title="Sequence", frameon=False)
        sns.despine(ax=ax)
    fig.suptitle(f"Cav1.2 | {mutant_label} | matched backbone geometry", fontweight="bold")
    handles = [Line2D([], [], linestyle="", marker=RMSD_REFERENCE_STYLES[r]["marker"],
                      markerfacecolor="white", markeredgecolor=RMSD_REFERENCE_STYLES[r]["color"],
                      label=f"{r} WT Cα–Cα") for r in ("8WE6", "8HLP", "8FD7")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .07), ncol=3, frameon=False)
    fig.text(.5, .025, f"Position {site} Cα to the same partner Cα in WT and mutant. "
             "Violins describe retained snapshots; they are not independent replicates.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .15, 1, .94))
    return table, fig


def plot_g402s_paired_revision(root):
    """Plot all matched Cα pairs, separating sequence and masking effects."""
    folder = root / "analysis/statistics_revision/paired_seed_v2"
    all_effects = pd.read_csv(folder / "g402s_matched_ca_partner_statistics.csv")
    all_effects = all_effects[
        all_effects.estimand_id.eq("primary_joint_nominal_seed")
    ].copy()
    partner_order = [
        "PHE1519", "VAL1520", "ALA1521", "VAL1522", "ILE1523", "MET1524",
        "ASP1525", "ASN1526", "PHE1527", "ASP1528", "TYR1529", "LEU1530",
        "THR1531", "ARG1532", "ASP1533", "TRP1534", "SER1535",
    ]
    specs = [
        ("vanilla G402S - vanilla WT", "A. Vanilla WT → G402S\nprimary sequence contrast", CAV12_PALETTE["G402S_VAN"]),
        ("G402S masked - G402S vanilla", "B. G402S vanilla → masked\nmasking effect", CAV12_PALETTE["G402S_HM"]),
    ]
    selected = all_effects[all_effects.contrast.isin([item[0] for item in specs])].copy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 7), sharey=True)
    y = np.arange(len(partner_order))
    for axis, (contrast, title, color) in zip(axes, specs):
        part = selected[selected.contrast.eq(contrast)].set_index("partner").loc[partner_order]
        values = part.estimate_B_minus_A.to_numpy(float)
        axis.errorbar(
            values, y,
            xerr=[values - part.CI_low_B_minus_A.to_numpy(float),
                  part.CI_high_B_minus_A.to_numpy(float) - values],
            fmt="o", capsize=2, color=color, ecolor="#8292A6",
        )
        axis.axvline(0, color="gray", lw=.8, ls="--")
        axis.set(yticks=y, yticklabels=[format_distance_alias(p) for p in partner_order],
                 xlabel="Matched position-402 Cα distance change (Å)", title=title)
        axis.invert_yaxis()
        sns.despine(ax=axis)
    fig.suptitle("Cav1.2 | G402S | paired-seed backbone revision", fontweight="bold")
    fig.text(.5, .03, "Each point uses the same Cα atom pair in WT and G402S. "
             "The left panel estimates the sequence effect; the right panel estimates masking within G402S.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .14, 1, .92))
    return selected, pd.read_csv(folder / "source_data/g402s_matched_ca_1519_1535.csv.gz"), fig


def plot_g406r_screen_and_contacts(root):
    """Put the R406 overlap-screen yield beside conditional contacts."""
    folder = root / "analysis/statistics_revision/paired_seed_v2"
    results = pd.read_csv(folder / "g406r_coordinate_screen_statistics.csv")
    denominators = pd.read_csv(folder / "g406r_coordinate_screen_summary.csv")
    results = results[
        results.atom_selection_rule.eq("side_chain_non_hydrogen_atoms")
    ].copy()
    denominators = denominators[
        denominators.atom_selection_rule.eq("side_chain_non_hydrogen_atoms")
    ].copy()
    primary = results[results.estimand_id.eq("primary_joint_nominal_seed")].set_index(
        "outcome"
    )
    required = [
        "unconditional_snapshot_screen_pass_fraction",
        "conditional_ASP1528_contact_fraction",
        "conditional_ASP1533_contact_fraction",
    ]
    missing = sorted(set(required) - set(primary.index))
    if missing:
        raise KeyError(f"Missing R406 screen outcomes: {missing}")

    fig, (screen_ax, contact_ax) = plt.subplots(1, 2, figsize=(10.5, 5.4), sharey=True)
    colors = [CAV12_PALETTE["G406R_VAN"], CAV12_PALETTE["G406R_HM"]]
    screen_rows = denominators[denominators.partner_residue.eq("screen_only")].set_index("protocol")
    screen_values = 100 * screen_rows.loc[["vanilla", "masked"], "unconditional_pass_fraction"].to_numpy(float)
    screen_ax.bar([0, 1], screen_values, color=colors, width=.62)
    for x, value in enumerate(screen_values):
        screen_ax.text(x, value + 2, f"{value:.1f}%", ha="center", fontweight="bold")
    screen_ax.set(xticks=[0, 1], xticklabels=["vanilla", "masked"],
                  ylabel="Fraction of retained snapshots (%)",
                  title="A. Passed explicit R406 side-chain\noverlap screen")
    screen_ax.text(.5, .60, "−31.7 percentage points\nwhole-seed 95% interval: −33.2 to −30.3",
                   transform=screen_ax.transAxes, ha="center", va="center", fontsize=9)

    x = np.arange(2)
    width = .36
    contact_rows = denominators[denominators.partner_residue.isin(["ASP1528", "ASP1533"])]
    contact_values = np.array([
        contact_rows[contact_rows.protocol.eq(protocol)].set_index("partner_residue")
        .loc[["ASP1528", "ASP1533"], "conditional_contact_fraction"].to_numpy(float)
        for protocol in ("vanilla", "masked")
    ]) * 100
    for offset, values, label, color in zip(
        (-width/2, width/2), contact_values, ("vanilla", "masked"), colors
    ):
        contact_ax.bar(x + offset, values, width, label=label, color=color)
        protocol = label
        protocol_rows = contact_rows[contact_rows.protocol.eq(protocol)].set_index("partner_residue")
        for xx, value, partner in zip(x + offset, values, ("ASP1528", "ASP1533")):
            numerator = int(protocol_rows.loc[partner, "conditional_contact_numerator"])
            denominator = int(protocol_rows.loc[partner, "conditional_contact_denominator"])
            contact_ax.text(xx, value + 1.2, f"{value:.1f}%\n({numerator}/{denominator})", ha="center", fontsize=8)
    contact_ax.set(xticks=x, xticklabels=["R406–D1528", "R406–D1533"],
                   title="B. Within-4 Å proximity\namong screen-passing snapshots")
    contact_ax.legend(title="Protocol", frameon=False)
    for axis in (screen_ax, contact_ax):
        axis.set_ylim(0, 82)
        sns.despine(ax=axis)
    fig.suptitle("Cav1.2 | G406R | overlap screen and conditional proximity",
                 fontweight="bold")
    fig.text(.5, .025,
             "Passing means none of the 19 explicit side-chain-heavy R406-centered distances was <2 Å. "
             "It does not validate the entire pocket. Conditional contact fractions use only pass snapshots.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, .91))
    return denominators, fig
