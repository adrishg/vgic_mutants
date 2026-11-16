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
