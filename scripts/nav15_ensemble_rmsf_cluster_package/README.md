# NaV1.5 ensemble RMSF implementation

This directory contains the reusable alignment and ensemble-profile code for
the current 50% NaV1.5 cohorts. It is analysis code only; prediction jobs,
array launchers, and intermediate shard assembly are not part of the public
workflow.

The implementation aligns each single-chain model once to the 6UZ3 frame and
preserves the complete raw residue axis (1–1572). Ensemble RMSF is calculated
from aligned Cα coordinates relative to the ensemble mean. Experimental
structures are evaluated separately and are not inserted into the prediction
ensemble.

The current configuration includes WT and QQQ, each with vanilla, masked v1,
and masked v2 50% cohorts. The IFM-excluded control is not part of the active
analysis.

## Reusable commands

Run from the repository root with the current 50% configuration:

```bash
python scripts/nav15_ensemble_rmsf_cluster_package/inspect_nav15_inputs.py \
  --config scripts/nav15_ensemble_rmsf_cluster_package/config/nav15_hive_50pct.yaml \
  --output-dir outputs/nav15_rmsf/input_inventory

python scripts/nav15_ensemble_rmsf_cluster_package/prepare_nav15_references.py \
  --config scripts/nav15_ensemble_rmsf_cluster_package/config/nav15_hive_50pct.yaml \
  --output-dir outputs/nav15_rmsf/references
```

The alignment and profile entry points accept explicit inputs and a
deterministic task slice. They can be run serially or parallelized by any
external workflow system without changing the scientific calculation. The
publication profiles can also be recalculated locally from aligned arrays
with `scripts/ensemble_rmsf_analysis/recompute_all_ok3_local.py`.

## Scientific details

- Alignment core: 132–282, 359–426, 505–665, 690–749, 870–1112,
  1128–1164, 1200–1405, and 1429–1466.
- The IFM motif is annotated at raw positions 1169–1171.
- All raw positions remain available for RMSF and regional summaries.
- Missing Cα coordinates, failed mappings, and insufficient coverage are
  recorded in metadata and are not silently filled.
