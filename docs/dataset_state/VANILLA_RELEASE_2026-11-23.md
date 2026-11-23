# Vanilla control release to HuggingFace, 2026-11-23

The ten unmasked (vanilla) control cohorts were published to
`adrishgz/vgic-mutant-structural-ensembles`, with coordinates, AlphaFold run metadata, and
the non-redundant distance tables. The dataset went from 18,316 to 18,380 files.

Everything was packed, verified and uploaded from Hive. No model data was copied to the
laptop.

## What went up

63 files, 44.1 GB, in Slurm job 25919654 (51 minutes, 63/63 uploaded, 0 failed).

| cohort | models | members | run_metadata | flavor |
|---|---|---|---|---|
| Kv2.1 WT | 1.78 GB | 6000 | 4.36 GB | multimer_v3 |
| Kv2.1 L403A | 1.78 GB | 5994 | 4.34 GB | multimer_v3 |
| Kv2.1 F412L | 1.78 GB | 6000 | 4.37 GB | multimer_v3 |
| CaV1.2 WT, G402S, G406R, G490R | 1.24 GB each | 6000 each | 2.50 GB each | ptm |
| CaV1.2 E1115K | 1.24 GB | 6000 | 2.35 GB | ptm |
| NaV1.5 WT, QQQ | 1.16 GB each | 6000 each | 2.19 GB each | ptm |

Coordinates total 43.67 GB on the Hub. Layout is
`coordinates/<channel>/<background>/vanilla/{models,run_metadata}.tar.zst`, with the
per-cohort manifests under `metadata/coordinate_pdb_release/<cohort_id>/` and an aggregate
manifest at `metadata/coordinate_pdb_release_vanilla/release_manifest.json`.

Run metadata is larger than the coordinates: 29.8 GB against 13.9 GB. The per-model scores
JSON dominates it, roughly 39 MB each, almost all of that the predicted aligned error matrix.

## Verification

Packing job 25906413, verification job 25907426: 10/10 PASS, 59,994 members streamed and
checked against the manifests by sha256 and size, 0 mismatches, 0 missing, 0 extra. Sampled
structures are sane per channel (Kv2.1 ABCD / 2400 CA / residues 1-600; CaV1.2 chain A /
1685; NaV1.5 chain A / 1572), with no missing END records. After upload, each archive's
sha256 was checked against the sha256 the Hub recorded (job 25920796).

## Which distance tables were published, and why

13 of the 35 vanilla tables on Hive. The rule, applied per cohort greedily from the richest
table down: omit a table when it is byte-identical to a kept table, or when both its columns
and its `pdb_file` rows are contained in a single kept table. Row-filtered variants
(`*_all_filtered.csv`) are kept, because which rows the filter removed cannot be recovered
from the unfiltered table. Job 25919590; the decision for all 35 is in
`vanilla_hf_release_20261122/table_shortlist.json` on Hive.

A first version of this rule was wrong and its output was discarded. It dropped any table
whose columns were covered by the *union* of the others, which threw out the richest table in
every cohort, and it let equal-column pairs drop each other circularly. It reported 2 keepers.

## Things recorded rather than repaired

- **Kv2.1 L403A has 5994 members, not 6000.** Seed 154 model 5 holds recycle intermediates
  r0-r5, no final structure, and no scores JSON (499 instead of 500). That run died partway.
  It is not the recycle early-stopping seen elsewhere, where finals are present.
- **CaV1.2 E1115K has no vanilla distance table at all.** Its coordinates are released
  without one.
- **CaV1.2 G490R has only the early 392-column table**, with no `_all` equivalent.
- **Kv2.1 WT needs two tables.** Its `_all` table has 1075 columns while L403A and F412L have
  1411. The file named `26-02-11_Kv2.1_wt_vanillaAF2test_distances.csv` is not a discarded
  test run: it carries 336 chain-B/C/D intra-subunit columns the `_all` table lacks, measured
  on the same models, with shared cells agreeing exactly. WT-against-mutant comparisons on
  those 336 columns cannot be made from `_all` alone.
- **`rmsd_convergence.csv` exists for only two vanilla cohorts**, Kv2.1 WT and CaV1.2 WT. The
  other eight were never run through that step. If the vanilla arm needs the same convergence
  treatment as the masked arm, that is computation still to run.

## Still open

The 30 stale files under `data/*/50pct_reanalysis_2026-11-12/` have not been touched: 15
masked RMSF tables, `kv21_a3m_mask_summary.csv`, 7 clashscore/reconciliation tables, and
`run_manifest.json`. Three vanilla `trajectory_medians.csv.gz` live in that same tree and
would need rehoming if it is removed wholesale.
