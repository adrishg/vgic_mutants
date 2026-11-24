# Nav1.5 filter-chain results: which predictions passed, and why the others did not

Three CSVs, one per model set. They are **not versioned in this repository** - its `.gitignore`
keeps data out by policy - so the authoritative copies are
`/quobyte/yarovoygrp/ahgz/vgic_mutants/analysis/filter_chain/nav15_*_filter_verdicts.csv` on Hive,
and the ones beside this file are a working copy.

Three CSVs, one per model set:

| file | protocol | conditions | predictions |
|------|----------|------------|-------------|
| `nav15_vanilla_filter_verdicts.csv` | unmasked control | WT, QQQ | 500 each |
| `nav15_masked_v1_filter_verdicts.csv` | 50% MSA masking, design v1 | WT, QQQ | 500 each |
| `nav15_masked_v2_filter_verdicts.csv` | 50% MSA masking, design v2 | WT, QQQ | 500 each |

QQQ is the inactivation-removed mutant: the IFM motif replaced by three glutamines.

One row per **prediction** — one (model, seed) ColabFold trajectory, its eleven recycles taken
together. Every check but the last judges the trajectory as a whole, so one bad recycle removes the
whole trajectory, including the recycles that look fine.

A note on numbering before anything else: **these models are rat Nav1.5 (UniProt P15389), not
human.** The model's own numbering maps onto P15389 in three blocks — +0 through model residue 461,
+197 from 462 to 871, +318 from 872 to 1572 — with no mismatches. Everything in these CSVs is in
the **model's own numbering**; nothing here has been converted.

## The checks, in the order they run

| # | Check | What it tests | Cutoff |
|---|-------|---------------|--------|
| 1 | **convergence** (all-OK-3) | From which recycle on does every later successive-recycle transition stay under 3 Å? That recycle is the trajectory's *representative*; it and everything after it are *retained*. A trajectory with no such recycle is removed. | 3 Å |
| 2 | **selectivity-filter ring** | Do the four filter residues (Nav1.5 DEKA: D373, E704, K1103, A1395) stay close together in every retained recycle? Measured as the six pairwise distances. | 20 Å |
| 3 | **pore–VSD interface** | Does the lower S6 stay against its own S4–S5 linker in every retained recycle? Closest Cα–Cα approach per repeat. Calibrated on fourteen experimental structures, where this contact runs 4.8–9.2 Å and the subunit's own linker is always the nearest one. A model whose S6 has let go of that linker has come apart at the one place the pore is held. | 12 Å |
| 4 | **topology screen** | Three reference-free checks on the representative model: `wiring` (the helices connect in the right order), `box_corners` (the four repeats sit at the corners of a sensible box), `s6_on_own_s5` (each S6 packs against its own S5). | pass/fail |

Order matters in one direction only. Convergence must run first, because every later check is
defined on the recycles it retains, and because convergence itself is defined on a *contiguous*
recycle series — a check that deleted recycles out of the middle would silently redefine it. The
topology screen must run last, because it is the only per-model check: run earlier it would delete
exactly the bad recycles that checks 2 and 3 exist to catch, and those checks would then find
almost nothing. Checks 2 and 3 are absolute thresholds, so their order relative to each other
changes no result.

Recycle 0 is judged like any other. The convergence step often picks it as the representative, and
a screen that skipped it would let exactly those models through unchecked.

## The columns

| column | meaning |
|--------|---------|
| `dataset`, `condition`, `channel` | which model set, which condition (WT / QQQ), which channel |
| `prediction` | the trajectory, named by the file without its `.rN` suffix |
| `model`, `seed` | AlphaFold network (1–5) and random seed |
| `recycles_in_distance_table` | recycles measured for this prediction |
| `recycles_retained_after_convergence` | recycles from the representative onward; 0 if it never converged |
| `representative_model_file`, `representative_recycle` | the single model checks 4 judges and the analysis reports |
| `convergence_result` | `pass` / `fail` |
| `selectivity_filter_ring_max_A`, `..._cutoff_A`, `..._result` | widest of the six DEKA distances over the retained recycles, the cutoff, and the verdict |
| `pore_vsd_s6_linker_max_A`, `..._cutoff_A`, `..._result` | widest S6-to-own-linker approach over the retained recycles, the cutoff, and the verdict |
| `pore_vsd_published_max_A`, `..._cutoff_A`, `..._result` | a second, Kv2.1-only form of check 3. **Always `not_applicable` for Nav1.5**: the pairing it measures is a contact in a domain-swapped homotetramer and nowhere else. |
| `topology_wiring_result`, `topology_box_corners_result`, `topology_s6_on_own_s5_result`, `topology_result` | the three reference-free checks on the representative, and their conjunction |
| `status` | `passed` or `removed` |
| `removed_by` | the first check that removed it, empty when it passed |
| `removed_because` | one sentence saying what that means |

`not_applicable` means the check did not run — because the channel has no such cutoff, or because
the prediction had already been removed before a measurement could be made. It never means "passed".

Every `*_result` column reports what was measured, independently of whether the prediction had
already been removed earlier. `removed_by` is the one column that says which check actually took
it: the first one it failed, in the order above.

## What this run found (2026-11-23)

| model set | condition | predictions | passed | removed | by what |
|---|---|---|---|---|---|
| vanilla | WT | 500 | 500 | 0 | — |
| vanilla | QQQ | 500 | 499 | 1 | topology (`box_corners`) |
| masked v1 | WT | 500 | 500 | 0 | — |
| masked v1 | QQQ | 500 | 500 | 0 | — |
| masked v2 | WT | 500 | 500 | 0 | — |
| masked v2 | QQQ | 500 | 499 | 1 | topology (`box_corners`) |

Two removals in 3,000 predictions. Nothing was removed by convergence, by the selectivity-filter
ring or by the pore–VSD interface. Read that as the chain having no discrimination on these
ensembles rather than as a certificate: the ring has never fired anywhere, and the widest
S6-to-own-linker approach across all six cohorts stays inside the range of the experimental
structures the cutoff was calibrated on (medians near 6.5 Å, maxima near 7.6 Å, against a 12 Å
cutoff and a 9.2 Å widest experimental contact).

**Treat the two topology removals as weak.** Both failed `box_corners` alone, and `box_corners` is
the check that transfers least well to a pseudo-tetramer: its 45–135° window was set on
homotetramers, and on experimental Nav1.5 structures it reads roughly 70/75/90/125° and already
fails 2 of 19. A single `box_corners` failure is close to that calibration edge. Both removals are
QQQ; neither is WT.

**Nearly every representative is recycle 0** — 1000/1000 vanilla, 997/1000 masked v1, 999/1000
masked v2. That is why the screen had to be made to judge recycle 0: a screen starting at recycle 1
would have left effectively the entire Nav1.5 set unjudged while reporting that it had screened it.

## What to expect in these numbers

The chain has little to do on Nav1.5. The DEKA ring has never rejected a single trajectory in any
cohort of any channel — it is a tripwire that has not fired, and should be read as one, not as a
filter that cleaned anything. The pore–VSD measure sits inside the experimental envelope across the
whole set. Say what the filters did not reject as plainly as what they did; a check that removes
nothing on this data is a check with no discrimination here, not evidence of a clean ensemble.

## What changed since the set sent in November

`WHAT_CHANGED_SINCE_NOVEMBER.md` beside this file, with per-model tables
(`nav15_masked_v1_old_vs_new.csv`, `nav15_masked_v2_old_vs_new.csv`,
`nav15_old_vs_new_summary.csv`). In one line: of the four masked cohorts sent in November, three are
identical model for model and the fourth lost exactly one model. The vanilla arm was never sent.

## Where the models themselves are

See `MODEL_BUNDLES.md` beside this file.
