# What changed between the November Nav1.5 set and this one

Short answer: **one model, out of two thousand.** And the two protocols that changed by nothing at
all changed by nothing at all — that is a measured result, not an absence of checking.

## What was actually sent before

Not reconstructed from memory. Read from the package that was built and shipped:

`/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5/share/Nav15_50pct_filtered_package_2026-11-11/`
(and its `.zip` with a `.sha256` beside it, 2026-11-11).

Its builder, `Nav1.5/package_nav15_50pct.slurm`, copies for each condition exactly the `pdb_file`
column of

`Nav15_rmsd_convergence_<condition>/earliest_converged_one_per_seed_model.csv`

into `models/<condition>/`, and nothing else. So the old criterion was **RMSD convergence alone,
one representative per trajectory**: the earliest recycle from which every later
successive-recycle transition stays under 3 Å. 500 models per condition. The package's own README
says so too, and warns against using `all_ok_models.csv` as the model set.

**The package holds the four masked cohorts only — `wt`/`qqq` × `masked_v1_50pct`/`masked_v2_50pct`.
There is no vanilla in it.** The vanilla protocol has never been sent, so for `nav15_vanilla` there
is nothing to compare against; its 999 passing predictions are all new to the collaborator.

Four provenance checks per condition, all sixteen passed, so every difference below is the filter
chain's doing and not a convergence step that moved underneath:

- the models shipped in `models/<condition>/` are exactly the rows of the selection CSV;
- the package's `all_ok_models.csv` is byte-identical to the one on Hive today;
- so is its `earliest_converged_one_per_seed_model.csv`;
- and that selection is exactly the rows `all_ok_models.csv` marks `earliest_converged_selected`,
  which is what the chain reads.

## The comparison

| dataset | condition | November package | new chain | in both | removed by the new chain | added by the new chain |
|---|---|---|---|---|---|---|
| masked v1 | WT | 500 | 500 | 500 | 0 | 0 |
| masked v1 | QQQ | 500 | 500 | 500 | 0 | 0 |
| masked v2 | WT | 500 | 500 | 500 | 0 | 0 |
| masked v2 | QQQ | 500 | **499** | 499 | **1** | 0 |
| vanilla | WT, QQQ | — never sent — | 500, 499 | — | — | — |

The one model that moved:

```
nav15_qqq_7fbs_masked_v2_50pct_unrelaxed_rank_047_alphafold2_ptm_model_3_seed_077.r0.pdb
removed by: topology screen (box_corners; wiring and s6_on_own_s5 both pass)
```

Per-model tables, every model named with its status under both criteria:
`nav15_masked_v1_old_vs_new.csv`, `nav15_masked_v2_old_vs_new.csv`, summary in
`nav15_old_vs_new_summary.csv`. Columns: `in_november_package`, `in_new_filter_chain`, `change`
(`unchanged_kept` / `removed_by_new_chain` / `added_by_new_chain`), plus `removed_by` and
`removed_because`.

## The recycle-0 question, answered as measured rather than as expected

The brief for this comparison expected recycle 0 to show up as models the new chain keeps that the
old criterion dropped. **It does not, and the reason is worth stating, because it changes what the
recycle-0 fix means.**

The November selection already chose recycle 0 as the representative for almost every trajectory:
498/500 (WT v1), 499/500 (QQQ v1), 500/500 (WT v2), 499/500 (QQQ v2). The new chain selects the
same recycles, because it runs the same convergence step. So the recycle-0 count is unchanged on
both sides and the `added_by_new_chain` column is zero everywhere. There is no recycle-0 delta to
separate from the topology delta; they are not two additive effects on this data.

What the recycle-0 change actually did is make the topology screen *apply* to these sets. The
screen's `[screen] first_recycle` was 1, meaning recycle 0 was never judged. With essentially every
Nav1.5 representative at recycle 0, a screen starting at recycle 1 would have had one to three
models per cohort to look at and would have reported the set screened. Both removals this run — the
masked v2 one above and the vanilla QQQ one — are recycle-0 representatives, so both are exactly
the cases that setting would have missed.

In other words: recycle 0 is not worth any models here, but without judging it the topology filter
is worth nothing either.

## How to read "nothing changed"

Three of the four previously-sent cohorts are identical, model for model. That is a real finding
about these ensembles and not a failure to test: the same chain removes 297 of 500 Kv2.1 WT
trajectories, so it is not a chain that never rejects. On Nav1.5 the selectivity-filter ring
rejected nothing (it has rejected nothing in any cohort of any channel, ever), the pore–VSD
interface rejected nothing, and the topology screen rejected one.

And treat that one as weak. It failed `box_corners` alone, which is the topology check that
transfers worst to a pseudo-tetramer: its 45–135° window was calibrated on homotetramers and
already fails 2 of 19 experimental Nav1.5 structures. One `box_corners` failure sits close to that
calibration edge.

So anything the collaborator has already run on the November masked sets stands, except for
whatever that single QQQ v2 trajectory contributed. The substantive new material is the vanilla
arm, which they have never had.

## One other thing that went out

`Nav1.5/forDiego6.tar.gz` (19.4 GB, 2026-11-13, owned by `dlopezma`) holds **fourteen complete,
unfiltered masked model trees** — every recycle `.r0`–`.r10` plus the final structure, 6,000 files
each. Listed in full (Slurm job 25933818; the listing is at
`analysis/filter_chain/forDiego6_listing_summary.txt` on Hive):

| background | protocols |
|---|---|
| WT | `masked`, `masked_ifm`, `masked_test` (1,200 files), `masked_v1_50pct`, `masked_v1new`, `masked_v2`, `masked_v2_50pct`, `masked_v2_noIFM`, `masked_v2new`, `masked_v3` |
| QQQ | `masked`, `masked_v1_50pct`, `masked_v1new`, `masked_v2`, `masked_v2_50pct`, `masked_v2new` |

So the collaborator holds the **complete unfiltered 50% ensembles** — `wt`/`qqq` ×
`masked_v1_50pct`/`masked_v2_50pct`, all 6,000 models each — as well as the filtered 500-per-
condition package. The filtered package is the selection on top of exactly those model trees.

**There is no vanilla anywhere in it**, which is consistent with the November package: the vanilla
arm has not been sent in any form.

`forDiego4_new.tar.gz` (988 KB) and `forDiego5_new.tar.gz` (685 MB), both 2026-10-16 and both at
the repository root rather than under `Nav1.5/`, have not been examined. They predate every cohort
compared here.
