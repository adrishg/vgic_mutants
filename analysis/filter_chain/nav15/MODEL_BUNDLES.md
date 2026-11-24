# Where the complete Nav1.5 model sets live

Established 2026-11-23. Nothing was published, changed or re-uploaded to produce this note.

## On Hugging Face — all six cohorts are already public

Dataset: **`adrishgz/vgic-mutant-structural-ensembles`**
(<https://huggingface.co/datasets/adrishgz/vgic-mutant-structural-ensembles>)

| model set | condition | path on the Hub | layout | size |
|---|---|---|---|---|
| vanilla | WT | `coordinates/NaV1.5/WT/vanilla/models.tar.zst` | one package | 1.16 GB |
| vanilla | QQQ | `coordinates/NaV1.5/QQQ/vanilla/models.tar.zst` | one package | 1.16 GB |
| masked v1 50% | WT | `coordinates/NaV1.5/WT/masked/v1_50pct/shard-00{0..4}.tar.zst` | 5 shards | — |
| masked v1 50% | QQQ | `coordinates/NaV1.5/QQQ/masked/v1_50pct/shard-00{0..4}.tar.zst` | 5 shards | — |
| masked v2 50% | WT | `coordinates/NaV1.5/WT/masked/v2_50pct/shard-00{0..4}.tar.zst` | 5 shards | — |
| masked v2 50% | QQQ | `coordinates/NaV1.5/QQQ/masked/v2_50pct/shard-00{0..4}.tar.zst` | 5 shards | — |

Each cohort is 6,000 PDBs: 500 predictions × eleven recycles (`.r0`–`.r10`) plus each
prediction's final structure written again without a `.rN` suffix. Chain A, 1,572 residues, rat
numbering.

Alongside them:

- `coordinates/NaV1.5/{WT,QQQ}/vanilla/run_metadata.tar.zst` — 2.19 GB each, the per-model scores
  JSON (almost all of it the predicted aligned error matrix). The masked cohorts have no published
  run metadata.
- `archives/NaV1.5_50pct_coordinates.zip` — the older bundle of the 50% coordinates.
- `metadata/coordinate_pdb_release/nav15_{wt,qqq}_vanilla/` — per-cohort manifests, one line per
  model with sha256 and size. The vanilla packages were verified end to end: every member streamed
  and hashed against the manifest, and each archive's sha256 checked against the one the Hub
  recorded.
- Distance tables under `data/Nav1.5/distances/`: the two vanilla
  `26-07-27_*_distances_all_all.csv` and, for each masked cohort, the raw table, its
  convergence-filtered companion (`_all_ok_rmsd_3A`) and the IFM extras.

### The two gaps

1. **The vanilla convergence tables are not published.** `data/Nav1.5/convergence/` holds the four
   masked cohorts only. The vanilla ones exist on Hive (see below) and would need uploading before
   anyone outside could reproduce which recycle each vanilla prediction converged at. Our verdict
   CSVs report the representative per prediction, so they partly stand in for this, but they are
   not the convergence table.
2. **The masked cohorts have no run metadata on the Hub** — no pLDDT/PAE JSON. Only vanilla does.

Nothing else is missing: all six model sets are complete and public.

## On Hive — the working copies

Root: `/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5`

| model set | condition | models | distance table | convergence |
|---|---|---|---|---|
| vanilla | WT | `wt/vanilla/` (6,000 PDB, 14 GB) | `dataDistances/26-07-27_Nav15_wt_vanillaAF2_distances_all_all.csv` | `rmsd_convergence_filtering/Nav15_rmsd_convergence_wt_vanilla/all_ok_models.csv` |
| vanilla | QQQ | `qqq/vanilla/` (6,000 PDB, 14 GB) | `dataDistances/26-07-27_Nav15_qqq_vanilla_AF2_distances_all_all.csv` | `rmsd_convergence_filtering/Nav15_rmsd_convergence_qqq_vanilla/all_ok_models.csv` |
| masked v1 50% | WT | `wt/masked_v1_50pct/models/` (6,000 PDB, 13 GB) | `dataDistances/26-11-10_Nav15_wt_masked_v1_50pct_AF2_distances.csv` | `Nav15_rmsd_convergence_wt_masked_v1_50pct/all_ok_models.csv` |
| masked v1 50% | QQQ | `qqq/masked_v1_50pct/models/` (6,000 PDB, 13 GB) | `dataDistances/26-11-10_Nav15_qqq_masked_v1_50pct_AF2_distances.csv` | `Nav15_rmsd_convergence_qqq_masked_v1_50pct/all_ok_models.csv` |
| masked v2 50% | WT | `wt/masked_v2_50pct/models/` (6,000 PDB, 13 GB) | `dataDistances/26-11-10_Nav15_wt_masked_v2_50pct_AF2_distances.csv` | `Nav15_rmsd_convergence_wt_masked_v2_50pct/all_ok_models.csv` |
| masked v2 50% | QQQ | `qqq/masked_v2_50pct/models/` (6,000 PDB, 13 GB) | `dataDistances/26-11-10_Nav15_qqq_masked_v2_50pct_AF2_distances.csv` | `Nav15_rmsd_convergence_qqq_masked_v2_50pct/all_ok_models.csv` |

The vanilla model folders are the ColabFold output directories themselves, not a `models/`
subfolder, and they also hold the a3m, the score JSONs and the coverage/PAE/pLDDT plots.

Chain-chain results of this run live in
`/quobyte/yarovoygrp/ahgz/vgic_mutants/analysis/filter_chain/`: per cohort an
`*_interface.csv` and an `*_topology_verdicts.csv`, and per dataset the
`*_filter_verdicts.csv` that is the deliverable.

## What would have to change to point a collaborator at the Hub

Nothing, for the models. The three model sets are complete and public, and the masked distance
tables are there in both raw and convergence-filtered form. The two gaps above — vanilla
convergence tables, masked run metadata — are uploads that have not been made, not data that is
missing.

The verdict CSVs produced here are **not** on the Hub. Whether they go up, and under what path, is
a decision that has not been taken.
