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
- The Kv2.1 50% masked ensembles analyzed for the paper were generated from a 41-position subset of
  the registered 99-position `kv21_common` design, not from the full design. Only 4 of the 48
  positions in the 370-417 pore-helix/S6 block were masked, so the pore domain received 6.7% residue
  coverage against 36.2% in the earlier 100% cohort. This affects all three sequence backgrounds
  (the mutants add their mutated position at 100%, giving 42). Vanilla cohorts are unaffected. The
  applied positions are recorded in `docs/MASK_REGISTRY.tsv`; a replacement run at the full 99
  positions is in progress for L403A.
- The NaV1.5 QQQ standard cohort was masked at 296 positions, not the 306 its design specifies: the
  293-position standard set plus the three substituted motif positions 1169-1171, without the rest
  of the 1164-1176 block. Positions 1164-1168 and 1172-1176 are covered by nothing else, because the
  standard set jumps 1157 -> 1215. WT masked under both designs, and QQQ restricted, match their
  definitions.
- Mask sizes in `docs/MASK_REGISTRY.tsv` are read from the per-A3M audit JSON files, which record
  what was applied. Where the applied and registered sizes differ, the registry gives both and the
  manuscript's Supplementary Table S1 reports the applied size.

