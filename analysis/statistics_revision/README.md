# Publication statistics for voltage-gated ion-channel ensembles

This directory contains the executable statistical analyses used for the manuscript figures and
supplementary figures. The workflows cover CaV1.2, Kv2.1, and NaV1.5 AlphaFold2 ensembles and keep
structural geometry, quality-control retention, and prediction-run variability as separate estimands.

## Statistical design

Retained recycle snapshots are reduced within each seed–model trajectory. Available AlphaFold2 model
parameterizations receive equal weight within a seed, and seeds receive equal weight between
conditions. Continuous outcomes use trajectory medians; categorical outcomes use explicitly defined
trajectory fractions or retained-recycle summaries. Bootstrap and permutation procedures resample
whole seeds or complete paired seed labels, as specified by each analysis.

Snapshot-level percentages describe structures sampled by the prediction protocol. They are not
thermodynamic occupancies, and recycle snapshots are not treated as independent biological
replicates.

## Publication workflows

Run from the repository root:

```bash
python analysis/statistics_revision/scripts/run_paired_seed_v2.py --mode publication
python analysis/statistics_revision/scripts/run_paired_seed_full_panel.py --mode publication
MPLCONFIGDIR=/tmp/vgci_mpl python analysis/statistics_revision/scripts/run_paired_seed_breadth_rmsf.py --mode publication
python analysis/statistics_revision/scripts/run_paired_seed_reduced_depth.py
python analysis/statistics_revision/scripts/run_qc_adjusted_target_yields.py --mode publication
python analysis/statistics_revision/scripts/build_publication_outputs.py
```

The generated publication tables, source data, and figures are written below
`analysis/statistics_revision/paired_seed_v2/`, `analysis/statistics_revision/tables/`, and
`analysis/statistics_revision/figures/`. The output manifest records the resampling unit, mask
definition, and scope of each figure and source table.

## Primary analyses

- Kv2.1 L403A: maximum E423–N179 Cα distance and the distribution of elongated interfaces across
  WT and L403A under vanilla and targeted masking.
- Kv2.1 F412L: mutation-site contact distributions with chain-resolved partner labels.
- CaV1.2 G402S: matched position-402 Cα distance effects across the declared partner panel.
- CaV1.2 G406R: coordinate-derived non-hydrogen side-chain overlap screening and conditional
  contacts.
- NaV1.5 WT and IFM→QQQ: regional experimental-reference RMSD and IFM/receptor coupling analyses.
- Cross-channel steric quality control: seed-balanced MolProbity clashscore contrasts.

The complete distance panel and regional RMSF analyses are retained as effect-size and uncertainty
support for the main and supplementary figures. The reduced-depth analysis is a prespecified
sensitivity used for Supplementary Figure S8; it uses 20 nominal seed labels and five AF2 model
parameterizations per label.

## Data and reproducibility

Large distance, RMSD, RMSF, and cohort tables are resolved through `shared/data_access.py` from the
working tree, a local cache, or the versioned companion dataset. The repository stores the analysis
code, experimental reference structures, manifests, registries, notebooks, and publication figure
source-data maps.

The exact mask definitions and registered datasets are documented in `docs/MASK_REGISTRY.tsv` and
`docs/FIGURE_INDEX.md`. The current scientific limitations are listed in `docs/KNOWN_ISSUES.md`.
