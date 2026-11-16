"""Generate the matched Kv2.1 WT-versus-F412L masking notebook."""

from pathlib import Path

import nbformat as nbf


HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "Supplementary_Figure_S9_Kv21_L403A_Masked_WT_vs_Mutant.ipynb"
OUT = HERE / "Supplementary_Figure_S10_Kv21_F412L_Masked_WT_vs_Mutant.ipynb"

nb = nbf.read(TEMPLATE, as_version=4)
for cell in nb.cells:
    cell.source = (
        cell.source
        .replace("supplementary_figure_s9", "supplementary_figure_s10")
        .replace("Figure_S9", "Figure_S10")
        .replace("S9", "S10")
        .replace("L403A", "F412L")
        .replace("l403a", "f412l")
    )
    if cell.cell_type == "code":
        cell.execution_count = None
        cell.outputs = []

nb.cells[0].source = r"""# Supplementary Figure S10 — does F412L alter access to the L403A-like E423–N179 shift?

This specificity control asks whether the elongated E423–N179 pore–VSD/S6 coordinate is confined to L403A or is also sampled by the other Kv2.1 mutant, F412L. The primary comparison is **WT masked versus F412L masked**, with WT and F412L vanilla ensembles retained as the baseline context.

The threshold (12.841 Å) remains the experimentally derived L403A-like definition: the midpoint between the longest 8SD3 WT subunit and the shortest elongated 8SDA L403A subunit. It is used here as a common structural coordinate, not as an F412L-specific experimental target. Because chain labels in a homotetramer can rotate, subunit distances are ranked within each structure before the asymmetry pattern is compared.

Inference uses complete input seeds as the independent units. Recycles are reduced within each seed–AF2-model trajectory, available AF2 models receive equal weight within a seed, and seeds receive equal weight between groups."""

# S10 still uses the L403A experimental structures to define the common
# E423–N179 threshold and reference markers.
nb.cells[3].source = nb.cells[3].source.replace(
    "f412l_experimental_threshold_derivation.csv",
    "l403a_experimental_threshold_derivation.csv",
)
nb.cells[7].source = nb.cells[7].source.replace("8SDA | F412L", "8SDA | L403A")
nb.cells[7].source = (
    nb.cells[7].source
    .replace("ax_d.text(2.04, 13.7,", "ax_d.text(2.04, 41.0,")
    .replace(
        "ax_d.set_ylim(0, 15.5)",
        "ax_d.set_ylim(0, max(44, float(occupancy.ci_high.max()) * 100 * 1.15))",
    )
)

nb.cells[8].source = """## Take-home message

F412L is evaluated against the same current convergence-qualified E423–N179 coordinate as WT. The exact seed-balanced shifted-state prevalence and uncertainty intervals are written to `e423_n179_group_summary.csv`; no legacy percentages are hard-coded into this notebook.

The masked-versus-WT contrast, the masking response, and the masking-by-F412L interaction are reported in `e423_n179_seed_bootstrap_contrasts.csv`. This keeps the specificity control tied to the current 50% inputs rather than to earlier dataset denominators.

As in the L403A control, the threshold-positive population is dominated by four-subunit elongation rather than the two-subunit asymmetry of experimental 8SDA. E423–N179 elongation is therefore not unique to L403A and should be interpreted as a shared Kv2.1 conformational coordinate whose baseline accessibility is especially high in F412L—not as a mutation-specific reconstruction of 8SDA."""

nb.metadata.kernelspec = {
    "display_name": "Python 3", "language": "python", "name": "python3",
}
nbf.write(nb, OUT)
print(OUT)
