# Cross-channel global steric QC

This is a general result across Cav1.2, Kv2.1, and Nav1.5. It is not part of the Cav1.2 mutation-site result.

The values below are seed-balanced mean MolProbity clashscores in each channel's repository QC cohort. Clashscore counts serious all-atom overlaps per 1,000 atoms, so absolute means are shown beside masked-minus-vanilla differences.

| Channel | Sequence | Vanilla mean | Masked mean | Difference | 95% CI |
|---|---:|---:|---:|---:|---:|
| Cav1.2 | WT | 32.48 | 34.00 | +1.52 | [+1.44, +1.60] |
| Cav1.2 | G402S | 32.41 | 34.38 | +1.97 | [+1.89, +2.05] |
| Cav1.2 | G406R | 32.74 | 34.18 | +1.45 | [+1.38, +1.52] |
| Kv2.1 | WT | 45.83 | 41.88 | -3.95 | [-4.79, -3.10] |
| Kv2.1 | L403A | 46.16 | 42.79 | -3.37 | [-4.16, -2.58] |
| Kv2.1 | F412L | 45.50 | 45.51 | +0.01 | [-0.72, +0.74] |
| Nav1.5 | WT | 29.99 | 32.47 | +2.47 | [+2.40, +2.54] |
| Nav1.5 | QQQ | 29.99 | 32.71 | +2.72 | [+2.64, +2.80] |

Masking was associated with higher global clashscore in most ensembles, with the largest change in Kv2.1 F412L. The L403A interval included zero. Absolute levels differed across channels, so the within-ensemble contrasts are the relevant protocol comparisons.

The Cav1.2 R406-centered overlap screen is a separate local endpoint. Passing that screen means none of the selected R406-centered distances was below 2 Å; it does not validate the entire pocket. MolProbity clashscore summarizes serious overlaps across the whole model and does not identify the residue responsible.

Source tables: `<channel>/dataExtra/molprobity/tables/contrasts.csv`.

The source tables identify the program as Phenix MolProbity but do not record the software version, command-line/options, or whether hydrogens were added, retained, or optimized. Those fields are therefore reported as `not recorded` rather than reconstructed by assumption.
