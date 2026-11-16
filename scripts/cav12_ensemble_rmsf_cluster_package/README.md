# CaV1.2 ensemble RMSF implementation

This directory contains the reusable alignment and ensemble-profile code for
the current 50% CaV1.2 cohorts. It is analysis code only; prediction jobs,
array launchers, and intermediate shard assembly are not part of the public
workflow.

The implementation aligns each single-chain model once to the 8WE6 frame and
preserves the complete raw residue axis (1–1685). Ensemble RMSF is calculated
from aligned Cα coordinates relative to the ensemble mean. Experimental
structures are kept as separate references and are not added as ensemble
members.

The current configuration includes WT, G402S, and G406R, each with vanilla and
masked 50% cohorts.

## Reusable commands

Run from the repository root with the current 50% configuration:

```bash
python scripts/cav12_ensemble_rmsf_cluster_package/inspect_cav12_inputs.py \
  --config scripts/cav12_ensemble_rmsf_cluster_package/config/cav12_hive_50pct.yaml \
  --output-dir outputs/cav12_rmsf/input_inventory

python scripts/cav12_ensemble_rmsf_cluster_package/prepare_cav12_references.py \
  --config scripts/cav12_ensemble_rmsf_cluster_package/config/cav12_hive_50pct.yaml \
  --output-dir outputs/cav12_rmsf/references
```

The alignment and profile entry points accept an explicit manifest,
references, output directory, and deterministic task slice. They can therefore
be run serially on a local machine or distributed by an external workflow
system without changing the scientific calculation. The final publication
profiles are also reproducible from aligned arrays with
`scripts/ensemble_rmsf_analysis/recompute_all_ok3_local.py`.

## Scientific details

- Alignment core: 112–381, 524–759, 901–1145, 1160–1197, 1239–1323, and
  1370–1534.
- The mutation-centered IS6 window 382–426 is excluded from fitting but kept in
  the output profile.
- Reference mapping and residue annotations are written with SHA-256 metadata.
- Missing atoms, failed alignments, and incomplete residue coverage are
  reported explicitly rather than silently imputed.

The public analysis notebooks and downstream summaries use the generated
profiles, masks, and QC manifests; they do not require the prediction or
cluster orchestration layer.
