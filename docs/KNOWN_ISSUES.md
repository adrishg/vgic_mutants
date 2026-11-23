# Known issues

- The exact production FASTA and A3M files are not distributed with the repository. Query identity,
  starting/masked A3M hashes, and homolog-row correspondence therefore require the recorded source
  files in the project data store.
- AlphaFold random-number values are not recorded in the available run metadata. `seed_###` labels
  are nominal design keys and are not interpreted as verified RNG equality across prediction runs.
- Kv2.1 experimental files 8SD3/8SDA are rat KCNB1 (P15387), while the exact modeled FASTA is not
  distributed here. The modeled construct organism is therefore not inferred from the experimental
  reference.
- The registered WT `nav15_standard` focal table lacks the six gate-span columns; its motif-separation contrast is estimable, but the corresponding gate-span contrast is not.
- F412L side-chain shortest-distance outcomes have no atom-matched WT analogue. Mutation-specific interactions are limited to common C-alpha outcomes when available.
- G402 and G406 have no non-hydrogen side-chain atoms. Position-402 sequence comparisons therefore
  use explicit C-alpha-to-C-alpha distances; position-406 local screening uses coordinate-derived
  non-hydrogen side-chain atoms for the arginine substitution.
- MolProbity source tables do not record the invoking software version, options, or hydrogen
  treatment. Those run fields remain `not recorded` rather than being reconstructed by assumption.
- RESOLVED 2026-11-21. The Kv2.1 50% masked ensembles analyzed up to 2026-11-20 were generated from
  a 41-position subset of the registered 99-position `kv21_common` design, not from the full design.
  Only 4 of the 48 positions in the 370-417 pore-helix/S6 block were masked, so the pore domain
  received 6.7% residue coverage against 36.2% in the earlier 100% cohort. This affected all three
  sequence backgrounds (the mutants added their mutated position at 100%, giving 42); vanilla
  cohorts were unaffected. The replacement run at the full 99 positions has finished for all three
  backgrounds and is what the dataset now publishes: `26-11-20` for WT and F412L, `26-11-19` for
  L403A, under the `masked_common99_50pct` cohort name. Its applied positions, read from every
  `*.a3m.audit.json`, match the registered design position for position (256-260, 276-280, 288-328,
  370-417), giving full 48/48 pore coverage. The 41-position cohort has been withdrawn from the
  companion dataset in full: the distance and convergence tables, the coordinate shards and the
  Kv2.1 50% coordinate archive. The replacement coordinates were published on 2026-11-22 as one
  models.tar.zst per sequence background under coordinates/Kv2.1/<variant>/masked_common99_50pct/,
  packed on Hive and uploaded from there. L403A holds 5995 members rather than 6000 because two
  models met the recycle tolerance early and wrote fewer intermediates; both have their final
  structures and the cohort manifest names them. Because each mutated position falls inside the
  99-position design, the mutant cohorts apply it deterministically and total 99 rather than 100.
- The NaV1.5 QQQ standard cohort was masked at 296 positions, not the 306 its design specifies: the
  293-position standard set plus the three substituted motif positions 1169-1171, without the rest
  of the 1164-1176 block. Positions 1164-1168 and 1172-1176 are covered by nothing else, because the
  standard set jumps 1157 -> 1215. WT masked under both designs, and QQQ restricted, match their
  definitions.
- Mask sizes in `docs/MASK_REGISTRY.tsv` are read from the per-A3M audit JSON files, which record
  what was applied. Where the applied and registered sizes differ, the registry gives both and the
  manuscript's Supplementary Table S1 reports the applied size.

## 2026-11-22: stale 50pct_reanalysis_2026-11-12 figures removed from HuggingFace

All 23 image files (5 PDF + 18 PNG) under
`data/Kv2.1/extra/figures/50pct_reanalysis_2026-11-12/` and
`data/Kv2.1/rmsf/figures/50pct_reanalysis_2026-11-12/` were deleted (HF commit
`d10d3caa`). The dataset now holds no PNG or PDF of any kind: it carries models,
archives and tables only. Dataset file count 18,342 -> 18,319.

The `*_trajectory_medians.csv.gz` files in
`data/Kv2.1/extra/conformation_analysis/50pct_reanalysis_2026-11-12/` were kept for
now. The three masked ones were verified to come from the withdrawn 41-column
`masked_v3_50pct` ensembles, not from `masked_common99_50pct`: recomputing per-trajectory
medians from the published common99 raw distance tables reproduces only 5.7% (WT), 8.6%
(L403A) and 5.7% (F412L) of cells exactly, and zero of 500 trajectories match in full.
For reference, when only the convergence filter differs and the underlying ensemble is
the same, 54% of cells still match exactly. The folder's own `run_manifest.json` records
`"masked_policy": "corrected_data_pending_final_structural_qc_and_clashscore"`.

The three vanilla medians files carry `"vanilla_policy": "unchanged_valid"` in the same
manifest, and the vanilla ensembles were never regenerated. They could not be checked as
cleanly, because only convergence-filtered vanilla distance tables are published, with no
raw table to recompute from; against the filtered tables they reproduce 10-13% of cells,
below the 54% filter-only baseline. Unresolved, and not a reason to delete them.

The three masked medians files were then deleted as well (HF commit `d33a0df5`, 18,319 ->
18,316), on the grounds that the 41-column ensemble they describe is not the one the paper
uses. The three vanilla medians remain. No masked `.gz` exists anywhere in the dataset now.

Still present and still derived from the withdrawn ensembles: 15 `*_masked_*` RMSF tables
under `data/Kv2.1/rmsf/analysis/50pct_reanalysis_2026-11-12/` (paired_rmsf, mask_classes,
top_rmsf_changes, absolute_change_allocation and experimental_differences, one set per
background), the clashscore and reconciliation tables under
`data/Kv2.1/extra/conformation_analysis/50pct_reanalysis_2026-11-12/`, and
`kv21_a3m_mask_summary.csv`, which claims 99 masked positions for ensembles that in fact
carried 41 or 42. These have not been removed or relabelled; the decision is open.

## 2026-11-22: superseded Hive cohorts archived to zip and the originals deleted

35 superseded masked cohorts across the three channels were zipped (Slurm job `25885704`)
and their loose files deleted (job `25898013`). 169,265 files, 505.6 GB of source, now
123.9 GB of `.zip` sitting beside each cohort's former location; 381.7 GB freed.

Each cohort was checked three ways before deletion: entry count equal to the source,
`unzip -t` CRC over every member, and a sha256 sample of 25 files per cohort compared
against the extracted member (830 files, 0 mismatches). The delete job re-checked the entry
count against the live directory immediately before removing it, and refused any active 50%
cohort by name.

Unaffected and still unzipped: the 11 active 50% cohorts (Kv2.1 `masked_common99_50pct`
x3, Cav1.2 `masked_v3_50pct` x4, Nav1.5 `masked_v1_50pct` and `masked_v2_50pct` x4) and all
10 vanilla cohorts. To restore an archived cohort:
`cd /quobyte/yarovoygrp/ahgz/vgic_mutants/<channel>/<background> && unzip <cohort>.zip`.

