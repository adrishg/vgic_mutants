"""Build the publication-facing notebooks for the current 50% cohorts."""

from __future__ import annotations

from pathlib import Path
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]


def md(source):
    return nbf.v4.new_markdown_cell(source.strip())


def py(source):
    return nbf.v4.new_code_cell(source.strip())


def write(path, cells):
    nb = nbf.v4.new_notebook()
    nb.metadata.update({
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3"},
    })
    nb.cells = cells
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(nbf.writes(nb), encoding="utf-8")


KV_LOAD = r'''
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT))
from shared.dataset_selection import (
    distance_csv_options, kv21_active_masked_distance_path,
    load_selected_distance_csv,
)

CONDITIONS = ("WT", "L403A", "F412L")
PROTOCOL = "masked_v3_50pct"
SELECTION = "all_ok_3"

def load_condition(condition):
    original = kv21_active_masked_distance_path(ROOT, condition)
    options = distance_csv_options(ROOT, original, "kv21", condition, PROTOCOL)
    source = str(options[SELECTION])
    assert "50pct" in source and "vanilla" not in source.lower(), source
    frame = load_selected_distance_csv(
        f"Kv2.1 {condition} | 50% masked", options, SELECTION,
        fallback_to_all=False,
    ).copy()
    assert frame.attrs["actual_selection"] == SELECTION
    frame["Condition"] = condition
    frame["Protocol"] = PROTOCOL
    return frame, source

datasets, source_paths = {}, {}
for condition in CONDITIONS:
    datasets[condition], source_paths[condition] = load_condition(condition)
assert all("50pct" in p and "vanilla" not in p.lower() for p in source_paths.values())
for condition, path in source_paths.items():
    print(f"{condition}: {path} | n={len(datasets[condition]):,}")
'''

KV_ANALYSIS = r'''
DISTANCE_COLUMNS = [f"CA_CA_{chain}_GLU425_CA-{chain}_ASN181_CA" for chain in "ABCD"]
assert not any(
    [column for column in DISTANCE_COLUMNS if column not in frame.columns]
    for frame in datasets.values()
)
SHIFT_THRESHOLD_A = 12.841337

def metrics(frame):
    result = frame[["pdb_file", "Condition", "Protocol", *DISTANCE_COLUMNS]].copy()
    result[DISTANCE_COLUMNS] = result[DISTANCE_COLUMNS].apply(pd.to_numeric, errors="coerce")
    result = result.dropna(subset=DISTANCE_COLUMNS).copy()
    result["max_S6_span_A"] = result[DISTANCE_COLUMNS].max(axis=1)
    result["median_S6_span_A"] = result[DISTANCE_COLUMNS].median(axis=1)
    result["shifted_subunits"] = (result[DISTANCE_COLUMNS] >= SHIFT_THRESHOLD_A).sum(axis=1)
    result["any_shifted"] = result["shifted_subunits"].ge(1)
    result["exactly_two_shifted"] = result["shifted_subunits"].eq(2)
    result["all_four_shifted"] = result["shifted_subunits"].eq(4)
    parsed = result["pdb_file"].str.extract(r"seed_(?P<seed>\d+).*?\.r(?P<recycle>\d+)\.pdb$")
    result["seed"] = pd.to_numeric(parsed["seed"], errors="coerce")
    result["recycle"] = pd.to_numeric(parsed["recycle"], errors="coerce")
    return result

long = pd.concat([metrics(frame) for frame in datasets.values()], ignore_index=True)
summary = (
    long.groupby("Condition", sort=False)
    .agg(
        n_models=("pdb_file", "size"), n_seeds=("seed", "nunique"),
        median_max_S6_span_A=("max_S6_span_A", "median"),
        q025_max_S6_span_A=("max_S6_span_A", lambda x: x.quantile(0.025)),
        q975_max_S6_span_A=("max_S6_span_A", lambda x: x.quantile(0.975)),
        fraction_any_shifted=("any_shifted", "mean"),
        fraction_exactly_two=("exactly_two_shifted", "mean"),
        fraction_all_four=("all_four_shifted", "mean"),
    ).reset_index()
)
summary
'''

KV_CONTRAST = r'''
rng = np.random.default_rng(20261113)

def bootstrap_difference(left, right, column, n_boot=5000):
    left = pd.to_numeric(left[column], errors="coerce").dropna().to_numpy(float)
    right = pd.to_numeric(right[column], errors="coerce").dropna().to_numpy(float)
    observed = np.median(left) - np.median(right)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        draws[i] = (
            np.median(rng.choice(left, len(left), replace=True))
            - np.median(rng.choice(right, len(right), replace=True))
        )
    return observed, np.quantile(draws, [0.025, 0.975])

rows = []
for mutant in ("L403A", "F412L"):
    for metric in ("max_S6_span_A", "median_S6_span_A", "shifted_subunits"):
        observed, interval = bootstrap_difference(
            long[long.Condition == mutant], long[long.Condition == "WT"], metric
        )
        rows.append({
            "contrast": f"{mutant} - WT", "metric": metric,
            "difference_median": observed, "ci95_low": interval[0],
            "ci95_high": interval[1],
        })
contrast_table = pd.DataFrame(rows)
contrast_table
'''

KV_OUTPUT = r'''
out_fig = ROOT / "docs" / "figures" / "kv21_50pct_s6_shift"
out_tab = ROOT / "docs" / "tables" / "kv21_50pct_s6_shift"
out_fig.mkdir(parents=True, exist_ok=True)
out_tab.mkdir(parents=True, exist_ok=True)
palette = {"WT": "#4C78A8", "L403A": "#E45756", "F412L": "#59A14F"}
order = list(CONDITIONS)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
for i, condition in enumerate(order):
    values = long.loc[long.Condition == condition, "max_S6_span_A"]
    sample = values.sample(min(700, len(values)), random_state=7)
    axes[0].violinplot(values, positions=[i], showmedians=True, widths=0.78)
    axes[0].scatter(np.full(len(sample), i), sample, s=3, alpha=0.08, color=palette[condition])
axes[0].axhline(SHIFT_THRESHOLD_A, color="#444444", linestyle="--", linewidth=1.1)
axes[0].set_xticks(range(len(order)), order)
axes[0].set_ylabel("Maximum same-chain Cα span (Å)")
axes[0].set_title("S6 span")
axes[0].grid(axis="y", alpha=0.25)

for condition in order:
    values = long.loc[long.Condition == condition, "shifted_subunits"]
    counts = values.value_counts().reindex(range(5), fill_value=0) / len(values)
    axes[1].plot(range(5), counts.values, marker="o", color=palette[condition], label=condition)
axes[1].set_xticks(range(5), ["0", "1", "2", "3", "4"])
axes[1].set_xlabel("Number of shifted subunits")
axes[1].set_ylabel("Fraction of structures")
axes[1].set_title("Coordinated S6 shift")
axes[1].legend(frameon=False)
axes[1].grid(axis="y", alpha=0.25)

for i, condition in enumerate(order):
    values = long.loc[long.Condition == condition, DISTANCE_COLUMNS].to_numpy(float).ravel()
    axes[2].boxplot(values, positions=[i], widths=0.55, showfliers=False,
                    patch_artist=True,
                    boxprops={"facecolor": palette[condition], "alpha": 0.55})
axes[2].axhline(SHIFT_THRESHOLD_A, color="#444444", linestyle="--", linewidth=1.1)
axes[2].set_xticks(range(len(order)), order)
axes[2].set_ylabel("Per-subunit Cα span (Å)")
axes[2].set_title("Subunit-resolved spans")
axes[2].grid(axis="y", alpha=0.25)
fig.suptitle("Kv2.1 | 50% masked S6-shift comparison", fontsize=15, fontweight="semibold")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(out_fig / "kv21_50pct_s6_shift_comparison.png", dpi=300, bbox_inches="tight")
plt.show()

summary.to_csv(out_tab / "cohort_summary.csv", index=False)
contrast_table.to_csv(out_tab / "bootstrap_contrasts.csv", index=False)
long.to_csv(out_tab / "per_structure_s6_shift_metrics.csv", index=False)
pd.DataFrame({
    "condition": list(source_paths), "source_path": list(source_paths.values()),
    "protocol": PROTOCOL, "selection": SELECTION,
}).to_csv(out_tab / "source_manifest.csv", index=False)
print(f"Wrote {out_fig / 'kv21_50pct_s6_shift_comparison.png'}")
'''

CAV_LOAD = r'''
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT))
from shared.cav12_coupling_analysis import CAV12_GATE_COLUMNS, gate_openness
from shared.dataset_selection import (
    cav12_active_masked_distance_path, distance_csv_options,
    load_selected_distance_csv,
)

CONDITIONS = ("WT", "G402S", "G406R")
PROTOCOL = "masked_v3_50pct"
SELECTION = "all_ok_3"

def load_condition(condition):
    original = cav12_active_masked_distance_path(ROOT, condition)
    options = distance_csv_options(ROOT, original, "cav12", condition, PROTOCOL)
    source = str(options[SELECTION])
    assert "50pct" in source and "vanilla" not in source.lower(), source
    frame = load_selected_distance_csv(
        f"CaV1.2 {condition} | 50% masked", options, SELECTION,
        fallback_to_all=False,
    ).copy()
    assert frame.attrs["actual_selection"] == SELECTION
    frame["Condition"] = condition
    frame["Protocol"] = PROTOCOL
    return frame, source

datasets, source_paths = {}, {}
for condition in CONDITIONS:
    datasets[condition], source_paths[condition] = load_condition(condition)
assert all("50pct" in p and "vanilla" not in p.lower() for p in source_paths.values())
assert all(column in datasets["WT"].columns for column in CAV12_GATE_COLUMNS)
for condition, path in source_paths.items():
    print(f"{condition}: {path} | n={len(datasets[condition]):,}")
'''

CAV_ANALYSIS = r'''
cohort_rows, gate_rows = [], []
for condition in CONDITIONS:
    frame = datasets[condition]
    parsed = frame["pdb_file"].str.extract(r"seed_(?P<seed>\d+).*?\.r(?P<recycle>\d+)\.pdb$")
    gate_rows.append(pd.DataFrame({
        "Condition": condition, "gate_openness_A": gate_openness(frame),
        "seed": pd.to_numeric(parsed["seed"], errors="coerce"),
        "recycle": pd.to_numeric(parsed["recycle"], errors="coerce"),
    }))
    cohort_rows.append({
        "Condition": condition, "n_models": len(frame),
        "n_seeds": pd.to_numeric(parsed["seed"], errors="coerce").nunique(),
        "n_columns": len(frame.columns), "source_path": source_paths[condition],
        "protocol": PROTOCOL, "selection": SELECTION,
    })
cohort_summary = pd.DataFrame(cohort_rows)
gate_long = pd.concat(gate_rows, ignore_index=True).dropna(subset=["gate_openness_A"])
gate_summary = (
    gate_long.groupby("Condition", sort=False)["gate_openness_A"]
    .agg(n="size", median="median", q025=lambda x: x.quantile(0.025),
         q975=lambda x: x.quantile(0.975)).reset_index()
)
cohort_summary, gate_summary
'''

CAV_CONTACTS = r'''
CONTACTS = {
    "G402S": {
        "G402/S402–I1523": ("shortest_GLY402-ILE1523", "shortest_SER402-ILE1523"),
        "G402/S402–M1524": ("shortest_GLY402-MET1524", "shortest_SER402-MET1524"),
    },
    "G406R": {
        "G406/R406–D1528": ("shortest_GLY406-ASP1528", "shortest_ARG406-ASP1528"),
        "G406/R406–D1533": ("shortest_GLY406-ASP1533", "shortest_ARG406-ASP1533"),
    },
}
contact_rows = []
for condition, definitions in CONTACTS.items():
    for contact, (wt_column, mutant_column) in definitions.items():
        for label, column in (("WT", wt_column), (condition, mutant_column)):
            assert column in datasets[label].columns, (label, column)
            values = pd.to_numeric(datasets[label][column], errors="coerce")
            contact_rows.append({
                "Condition": label, "Mutation_context": condition,
                "Contact": contact, "distance_column": column,
                "n_models": values.notna().sum(),
                "median_distance_A": values.median(),
                "contact_fraction_le_4A": (values <= 4.0).mean(),
                "severe_overlap_fraction_lt_2A": (values < 2.0).mean(),
            })
contact_summary = pd.DataFrame(contact_rows)
contact_summary
'''

CAV_FEATURES = r'''
LOCAL_FEATURES = [
    "CA_LEU401_CA-ILE1186_CA", "CA_LEU401_CA-VAL753_CA",
    "CA_LEU401_CA-ASP1533_CA", "CA_GLU407_CA-ARG1532_CA",
    "CA_PHE408_CA-THR1531_CA", "CA_LEU401_CA-LEU749_CA",
    "CA_LEU749_CA-VAL1182_CA", "CA_LEU749_CA-ILE1516_CA",
    "CA_VAL1182_CA-ILE1516_CA",
]
LOCAL_FEATURES = [column for column in LOCAL_FEATURES if column in datasets["WT"].columns]
feature_rows = []
for mutant in CONDITIONS[1:]:
    for column in LOCAL_FEATURES:
        wt_values = pd.to_numeric(datasets["WT"][column], errors="coerce")
        mutant_values = pd.to_numeric(datasets[mutant][column], errors="coerce")
        feature_rows.append({
            "Mutant": mutant, "feature": column,
            "WT_median_A": wt_values.median(),
            "mutant_median_A": mutant_values.median(),
            "median_difference_A": mutant_values.median() - wt_values.median(),
            "WT_q025_A": wt_values.quantile(0.025), "WT_q975_A": wt_values.quantile(0.975),
            "mutant_q025_A": mutant_values.quantile(0.025),
            "mutant_q975_A": mutant_values.quantile(0.975),
        })
feature_changes = pd.DataFrame(feature_rows)
feature_changes["absolute_median_difference_A"] = feature_changes["median_difference_A"].abs()
feature_changes = feature_changes.sort_values(
    ["Mutant", "absolute_median_difference_A"], ascending=[True, False]
)
feature_changes
'''

CAV_OUTPUT = r'''
out_fig = ROOT / "docs" / "figures" / "cav12_50pct_mutant_comparison"
out_tab = ROOT / "docs" / "tables" / "cav12_50pct_mutant_comparison"
out_fig.mkdir(parents=True, exist_ok=True)
out_tab.mkdir(parents=True, exist_ok=True)
palette = {"WT": "#4C78A8", "G402S": "#E45756", "G406R": "#F28E2B"}
order = list(CONDITIONS)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for i, condition in enumerate(order):
    values = gate_long.loc[gate_long.Condition == condition, "gate_openness_A"].to_numpy(float)
    sample = np.random.default_rng(17 + i).choice(values, min(600, len(values)), replace=False)
    axes[0].violinplot(values, positions=[i], showmedians=True, widths=0.78)
    axes[0].scatter(np.full(len(sample), i), sample, s=3, alpha=0.08, color=palette[condition])
axes[0].set_xticks(range(len(order)), order)
axes[0].set_ylabel("Maximum intracellular-gate Cα span (Å)")
axes[0].set_title("Gate openness")
axes[0].grid(axis="y", alpha=0.25)

bar = contact_summary.copy()
bar["label"] = bar["Mutation_context"] + "\n" + bar["Condition"]
x = np.arange(len(bar))
axes[1].bar(x - 0.18, bar["contact_fraction_le_4A"], width=0.36,
            color=[palette[value] for value in bar["Condition"]])
axes[1].bar(x + 0.18, bar["severe_overlap_fraction_lt_2A"], width=0.36,
            color="#777777", alpha=0.7)
axes[1].set_xticks(x, bar["label"], rotation=45, ha="right")
axes[1].set_ylabel("Fraction of structures")
axes[1].set_title("Mutation-local contacts")
axes[1].legend(["≤4 Å contact", "<2 Å overlap"], frameon=False, fontsize=8)
axes[1].grid(axis="y", alpha=0.25)

heat = feature_changes.pivot(index="feature", columns="Mutant", values="median_difference_A")
heat = heat.reindex(columns=list(CONDITIONS[1:])).dropna(how="all")
im = axes[2].imshow(heat.to_numpy(float), aspect="auto", cmap="coolwarm", vmin=-1.5, vmax=1.5)
axes[2].set_yticks(range(len(heat)), [x.replace("CA_", "").replace("_CA", "") for x in heat.index], fontsize=7)
axes[2].set_xticks(range(len(heat.columns)), heat.columns)
axes[2].set_title("Median local-distance change vs WT")
axes[2].set_xlabel("Mutant")
fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04, label="Δ median distance (Å)")
fig.suptitle("CaV1.2 | 50% masked WT and mutant ensembles", fontsize=15, fontweight="semibold")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(out_fig / "cav12_50pct_gate_contacts_and_local_changes.png", dpi=300, bbox_inches="tight")
plt.show()

cohort_summary.to_csv(out_tab / "cohort_summary.csv", index=False)
gate_summary.to_csv(out_tab / "gate_summary.csv", index=False)
contact_summary.to_csv(out_tab / "mutation_local_contacts.csv", index=False)
feature_changes.to_csv(out_tab / "local_feature_changes_vs_wt.csv", index=False)
gate_long.to_csv(out_tab / "per_structure_gate_openness.csv", index=False)
pd.DataFrame({
    "condition": list(source_paths), "source_path": list(source_paths.values()),
    "protocol": PROTOCOL, "selection": SELECTION,
}).to_csv(out_tab / "source_manifest.csv", index=False)
print(f"Wrote {out_fig / 'cav12_50pct_gate_contacts_and_local_changes.png'}")
'''

KV_CELLS = [
    md("""
    # Kv2.1 50% masked S6-shift comparison

    This is the canonical comparison for the current release. It uses only the corrected masked_v3_50pct cohorts and convergence-qualified all_ok_3 tables for WT, L403A, and F412L. No vanilla cohort, legacy date, or prior release is loaded.

    The primary question is whether the L403A S6 displacement is reproduced in the current 50% ensemble, with F412L retained as a mutation-specific comparison.
    """),
    py(KV_LOAD), py(KV_ANALYSIS), py(KV_CONTRAST), py(KV_OUTPUT),
    md("""
    ## Interpretation boundary

    The thresholded counts describe the same-chain S6-span signature used for the Kv2.1 comparison. They do not by themselves establish a full coordinated pore rotation; the per-subunit distributions and the structural reference analysis should be considered together.
    """),
]

CAV_CELLS = [
    md("""
    # CaV1.2 50% masked mutant comparison

    This notebook analyzes the current corrected masked_v3_50pct release with convergence-qualified all_ok_3 tables for WT, G402S, and G406R.

    The analysis is descriptive: global intracellular-gate openness, mutation-local S6-interface contacts, and local Cα-distance changes relative to WT.
    """),
    py(CAV_LOAD), py(CAV_ANALYSIS), py(CAV_CONTACTS), py(CAV_FEATURES), py(CAV_OUTPUT),
    md("""
    ## Interpretation boundary

    A contact fraction uses the shortest heavy-atom distance and the stated 4 Å threshold; the <2 Å fraction is reported as a severe-overlap diagnostic. Differences in ensemble medians identify candidate local repacking signatures, not a causal mechanistic effect by themselves.
    """),
]

if __name__ == "__main__":
    write(ROOT / "paperFigures" / "Kv21_50pct_S6_shift_comparison.ipynb", KV_CELLS)
    write(ROOT / "paperFigures" / "CaV12_50pct_mutant_comparison.ipynb", CAV_CELLS)
    print("Built 50% comparison notebooks.")

