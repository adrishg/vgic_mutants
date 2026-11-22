# Figure audit, 2026-11-21

Scope: experimental-reference marker consistency across the whole paper, a run
of every paper-figure notebook, and a pass over the draft's legends. Two fixes
are already committed; everything else below is open.

## 1. Experimental markers

The canonical assignment lives in `shared/plotting.py`
(`RMSD_REFERENCE_STYLES`, `NAV15_EXPERIMENTAL_STYLES`) and is reached through
`experimental_reference_style()`. Shape and color encode the same structure
redundantly, so a panel stays readable in grayscale.

| channel | reference | marker | color |
|---|---|---|---|
| Kv2.1 | 8SD3 WT | circle | `#E57373` |
| Kv2.1 | 8SDA L403A | square | `#F48FB1` |
| Kv2.1 | 9O10 / 9O11 / 9O12 / 9O13 | diamond / up / down / plus | raspberry / coral / orange / peach |
| Nav1.5 | 7FBS QQQ open | circle | `#D55E00` |
| Nav1.5 | 6UZ3 WT inactivated | square | `#0072B2` |
| Nav1.5 | 8VYJ open I | diamond | `#009E73` |
| Nav1.5 | 8VYK expanded open II | up triangle | `#CC79A7` |
| Nav1.5 | 7DTC E1784K | down triangle | `#E69F00` |
| Nav1.5 | 8T6L toxin-bound | plus | `#8C6D31` |
| Cav1.2 | 8HLP / 8HMA / 8HMB / 8WEA | circle / square / diamond / up | red / pink / raspberry / coral |
| Cav1.2 | 8WE9 / 8WE8 / 8WE7 / 8WE6 | down / plus / X / star | orange / peach / apricot / amber |
| Cav1.2 | 8FD7 / 8EOG | left / right | ochre / cream |

Markers repeat across channels, which is fine: no panel mixes two channels'
references, and the cross-channel distance figure resolves per channel.

**Fixed.** The Nav1.5 pore-size landscape panel assigned markers by row
position rather than by PDB, so all six disagreed with the aspect-ratio panel
beside it in the same figure — 6UZ3 was a circle in one and a square in the
other — and drew every reference in one dark lilac instead of its own color.
Both panels now read the canonical table, and
`tests/test_experimental_marker_consistency.py` renders them and compares the
drawn marker paths so this cannot return silently.

**Fixed.** The 8SDA marker in the Kv2.1 contact comparison was written out as
`'s'` instead of looked up. It was right, but only by hand.

**Open, needs your call.** Figure S7 marks the experimental *change*
(8SD3 → 8SDA) with a black diamond. That is a delta, not a structure, so it
does not break the one-marker-per-reference rule — but a diamond already means
9O10 in the Kv2.1 distance panels, and S9 and S10 use the canonical circle and
square for those same two PDBs. A reader going S7 → S9 → S10 sees the same two
structures three different ways. Changing it is a design decision, so I left it.

Everything else resolves canonically. Line markers on ensemble traces
(MolProbity, retention curves) are not reference markers and were left alone.

## 2. Notebook run

All 10 `paperFigures/` notebooks plus four channel notebooks were executed.

- **Pass (4):** `CaV12_50pct_mutant_comparison`, `CrossChannel_WT_RMSF`,
  `Supplementary_Figure_S7_Kv21_L403A_Experimental_Signatures`,
  `Cav12_WT_vanilla_vs_masked_v3new`.
- **Fail (9):** every one on a missing data file. No code errors anywhere.

The Cav1.2 notebooks were failing for a different reason and now pass: the
published dataset spells the channel `CaV1.2` and nests each bucket under a
`50pct/` level, which `hf_path` did not build, so every Cav1.2 fetch returned
404 although the files were there. Fixed in `shared/data_access.py`.

The nine remaining failures need two file families published:

1. `*_all_ok_rmsd_3A_structural_interface_qc.csv` for Kv2.1, vanilla and the
   99-column masked 50% set, all three backgrounds. Blocks seven notebooks, including both
   Figure 2B panels, the S6 shift comparison, S8, S9 and S10. These exist
   nowhere yet and have to be built. Hive holds the filter-chain session's
   `<cohort>_chain_qc.csv`, which is a different artifact from the published
   chain's `_structural_interface_qc.csv` that the notebooks resolve; build
   them with `scripts/filtering/build_kv21_structural_interface_qc_50pct.py`.
   A correct build gives 203/500 WT, 393/500 L403A and 277/500 F412L retained
   at 99 columns, with the G377 ring rejecting none in any cohort.
2. Nav1.5 distances at all. `data/Nav1.5/` on Hugging Face holds only `rmsd`:
   no `distances`, no `convergence`, no `rmsf`. Blocks the Figure 5 pore panels
   and the pore-shape analysis.

## 3. Draft legends

| figure | state |
|---|---|
| 2, 3, 4 | title only — no panel descriptions at all |
| 5 | panels A, B, C described; **D is not** |
| S3, S5 | no panel letters used, though both figures have parts |
| S4 | **legend missing entirely** (S3 jumps to S5) |
| S7 | describes A, B, C, D, (E1, E2), F–H — but the builder lays out 3 + 1 + 4 panels, so the bottom row is four panels, not an E split in two plus F–H. Either the letters or the layout is wrong. **G is never described on its own.** |
| S9, S10 | **legends missing entirely**, though both notebooks exist |

Also open from earlier: Figure 6 relayout, uncited panels, the data-availability
placeholder, and the missing Acknowledgments / Funding / Author Contributions /
Competing Interests sections.

## 4. Other drift found on the way

`scripts/generate_kv21_rmsd_comparison_notebook.py` no longer matches
`kv21/Kv21_WT_mutants_RMSD_comparison.ipynb`. The generator builds an 8SDA
contact panel and reads data through `VGCI_RECONCILIATION_DATA_ROOT`; the
committed notebook has neither. Regenerating would be a visible change to the
notebook, so it was not done here. Decide which one is current.

`analysis/active_masked_cohorts.tsv` has no Cav1.2 E1115K row, although E1115K
is part of the 50% set.
