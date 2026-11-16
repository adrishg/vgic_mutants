# Kv2.1 ensemble RMSF implementation

This directory contains the reusable alignment and ensemble-profile code for
the current 50% Kv2.1 cohorts. It is analysis code only; prediction jobs,
array launchers, and intermediate shard assembly are not part of the public
workflow.

The implementation aligns each tetramer once to the 8SD3 frame, infers the
cyclic chain correspondence, and retains the complete raw axis (residues
1–600 in chains A–D). Ensemble RMSF is calculated from aligned Cα coordinates
relative to the ensemble mean. 8SD3 and 8SDA are evaluated as separate
experimental references.

The current configuration includes WT, L403A, and F412L, each with vanilla and
masked 50% cohorts. The mask annotations are the registered target set in
`docs/MASK_REGISTRY.tsv` and the corresponding current input manifests.

## Reusable commands

Run from the repository root with the current 50% configuration:

```bash
python scripts/kv21_ensemble_rmsf_cluster_package/inspect_kv21_inputs.py \
  --config scripts/kv21_ensemble_rmsf_cluster_package/config/kv21_hive_50pct.yaml \
  --output-dir outputs/kv21_rmsf/input_inventory

python scripts/kv21_ensemble_rmsf_cluster_package/prepare_kv21_references.py \
  --config scripts/kv21_ensemble_rmsf_cluster_package/config/kv21_hive_50pct.yaml \
  --output-dir outputs/kv21_rmsf/references
```

The alignment and profile entry points accept explicit inputs and a
deterministic task slice. They can be run serially or parallelized by any
external workflow system without changing the scientific calculation. The
publication profiles can also be recalculated locally from aligned arrays
with `scripts/ensemble_rmsf_analysis/recompute_all_ok3_local.py`.

## Scientific details

- Alignment core: 184–208, 231–247, and 262–277 in the four chains.
- Chain correspondence is restricted to cyclic, adjacency-preserving mappings;
  chains are not aligned independently.
- RMSF is positional dispersion among independent AlphaFold outputs, not an MD
  fluctuation.
- Missing Cα coordinates, failed mappings, and insufficient coverage are
  recorded in metadata and are not silently filled.
- Structural/interface QC is applied downstream through the registered
  allowlists used by the distance analyses.
