# A3M provenance

The production A3M files are not distributed with this repository. The
analysis release therefore records the mask definitions and their numbering
conventions, while the prediction ensembles and derived tables are deposited
in the companion data repository.

Programmatic comparison with the original alignments requires an A3M tree
containing the relative paths in `mask_definitions.yaml`. With that tree, run:

```bash
python -m scripts.ensemble_rmsf_analysis.extract_masks \
  --a3m-root /path/to/production/vgic_mutants
```

The RMSF analysis uses direct-mask annotations from
`authoritative_mask_definitions.yaml` in 1-based raw AlphaFold query/model
numbering. The mask registry is the reproducibility record used by the
published analyses; it does not substitute for the original A3M files.

The command exits with an error when the A3M files are absent. When the
production tree is supplied, its output compares extracted position sets with
the registered mask definitions.
