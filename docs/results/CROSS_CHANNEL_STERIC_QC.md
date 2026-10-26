# Cross-channel global steric QC

This is a general result across Cav1.2, Kv2.1, and Nav1.5. It is not part of the Cav1.2 mutation-site result.

The values below are seed-balanced mean MolProbity clashscores in each channel's repository QC cohort. Clashscore counts serious all-atom overlaps per 1,000 atoms, so absolute means are shown beside masked-minus-vanilla differences.

| Channel | Sequence | Vanilla mean | Masked mean | Difference | 95% CI |
|---|---:|---:|---:|---:|---:|
| Cav1.2 | WT | 32.48 | 33.60 | +1.12 | [+1.05, +1.19] |
| Cav1.2 | G402S | 32.41 | 33.72 | +1.31 | [+1.24, +1.39] |
| Cav1.2 | G406R | 32.74 | 33.70 | +0.96 | [+0.90, +1.03] |
| Kv2.1 | WT | 45.88 | 47.95 | +2.08 | [+1.11, +3.07] |
| Kv2.1 | L403A | 46.09 | 47.00 | +0.91 | [-0.07, +1.84] |
| Kv2.1 | F412L | 45.38 | 49.20 | +3.82 | [+3.03, +4.64] |
| Nav1.5 | WT | 29.99 | 32.16 | +2.17 | [+2.09, +2.24] |
| Nav1.5 | QQQ | 29.99 | 32.15 | +2.16 | [+2.08, +2.23] |

Masking was associated with higher global clashscore in most ensembles, with the largest change in Kv2.1 F412L. The L403A interval included zero. Absolute levels differed across channels, so the within-ensemble contrasts are the relevant protocol comparisons.

The Cav1.2 R406-centered overlap screen is a separate local endpoint. Passing the current screen means none of 19 coordinate-derived side-chain-heavy distances was below 2 Å; it does not validate the entire pocket. MolProbity clashscore summarizes serious overlaps across the whole model and does not identify the residue responsible.

Source tables: `<channel>/dataExtra/molprobity/tables/contrasts.csv`.

The historical exports identify the program as Phenix MolProbity but do not record the software version, command-line/options, or whether hydrogens were added, retained, or optimized. Those fields are therefore reported as `not recorded` in the source table rather than reconstructed by assumption.
