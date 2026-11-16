# Figure index

## Publication main panels

- Figure 4: `analysis/statistics_revision/figures/Figure4_Kv21_L403A_four_condition.{png,pdf}`. Includes WT vanilla, WT masked, L403A vanilla, L403A masked, and the 0–4 elongated-interface distribution. Source data are in `analysis/statistics_revision/paired_seed_v2/source_data/figure4_*`.
- Figure 5 F412L panel: `analysis/statistics_revision/figures/Figure5_F412L_mutation_site_cleaned.{png,pdf}`. Partner labels retain chain identity. Source data are in `paired_seed_v2/source_data/figure5_*`.
- Figure 6 mask/reference metadata: `analysis/statistics_revision/tables/Figure6_mask_and_reference_markers.tsv`.
- Figure 7: `analysis/statistics_revision/figures/Figure7_Cav12_matched_G402S_R406_screen.{png,pdf}`. Panel D uses all 17 atom-matched position-402 Cα pairs; panel E uses the coordinate-derived side-chain-heavy R406 screen and reports conditional contact numerators and denominators.
- Cross-channel global clashscore: `analysis/statistics_revision/figures/cross_channel_clashscore_absolute_means.{png,pdf}`.

The shared sequence/protocol palette is listed in `palette.txt` and implemented
in `shared/plotting.py`.

The presentation notebooks live directly inside each channel directory. Generated figures and
tables are stored below the corresponding `dataDistances/`, `dataRMSD/`, or `dataRMSF/` analysis
directory.

## CaV1.2 experimental RMSD

Primary notebooks:

- `cav12/Cav12_WT_experimental_RMSD.ipynb`
- `cav12/Cav12_G402S_experimental_RMSD.ipynb`
- `cav12/Cav12_G406R_experimental_RMSD.ipynb`

The notebooks report retained-model counts, measurement completeness, overall and stable-core RMSD,
regional pore/S6/interface distributions, core-aligned versus locally aligned RMSD, protocol median
differences, and robust effect-size summaries. Generated outputs are under
`cav12/dataRMSD/analysis/<condition>/`.

Experimental references are 8HLP, 8WE6, and 8FD7. These are WT/state references rather than
mutation-matched G402S or G406R structures.

## Nav1.5 experimental RMSD

Primary notebooks:

- `nav15/Nav15_WT_experimental_RMSD.ipynb`
- `nav15/Nav15_QQQ_experimental_RMSD.ipynb`

The notebooks report retained-model counts, completeness, overall/core and regional RMSD,
core-aligned versus locally aligned measurements, masked-minus-vanilla median differences, and
effect-size summaries. WT and QQQ outputs are under `nav15/dataRMSD/analysis/wt/` and
`nav15/dataRMSD/analysis/qqq/`.

The main-text comparison emphasizes 8VYJ as the native full-length open reference and 7FBS as the
engineered QQQ-open reference. The broader experimental panel remains available as supplemental
state context.

## Kv2.1 experimental RMSD

Primary notebooks:

- `kv21/Kv21_WT_experimental_RMSD.ipynb`
- `kv21/Kv21_L403A_experimental_RMSD.ipynb`
- `kv21/Kv21_F412L_experimental_RMSD.ipynb`
- `kv21/Kv21_WT_mutants_RMSD_comparison.ipynb`

8SD3 is the mutation-matched WT reference and 8SDA is the L403A reference. For F412L, 8SDA is used
only as an L403A-like state comparator because no mutation-matched experimental F412L structure is
available.

## Figure outputs

The executed notebooks are the reproducible source for the publication figures. Rendered images and
source tables are generated from those notebooks and are intentionally kept out of Git when they are
available in the companion data repository.
