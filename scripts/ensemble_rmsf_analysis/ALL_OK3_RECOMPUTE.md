# Local all-OK-3 RMSF recomputation

Final-QC RMSF profiles can be recomputed locally from the aligned coordinate
arrays without reopening or realigning the original structures:

```bash
python scripts/ensemble_rmsf_analysis/recompute_all_ok3_local.py --channel all
```

The calculation streams bounded coordinate chunks, joins selection keys by
dataset and PDB basename, and validates the RMSF formula against the complete
all-model profile before writing each final profile. Kv2.1 uses the stricter
structural/interface QC allowlists; NaV1.5 and CaV1.2 use their all-OK-3
selection manifests.

The RMSF definition is

```text
RMSF_i = sqrt(E[||x_i||²] - ||E[x_i]||²)
```

where `x_i` is the aligned Cα coordinate at residue `i`, and the expectation
is taken over retained prediction snapshots. Downstream seed-level summaries
do not treat recycle snapshots as independent biological replicates.

The notebooks prefer these current profiles when available:

- `kv21_all_ok_3_symmetry_averaged_profiles.csv`
- `nav15_all_ok_3_per_residue_profiles.csv`
- `cav12_all_ok_3_per_residue_profiles.csv`

Selection manifests are stored with the channel analyses:

- `kv21/dataRMSF/qc/kv21_all_ok3_selection_manifest.csv`
- `nav15/dataRMSF/qc/nav15_all_ok3_selection_manifest.csv`
- `cav12/dataRMSF/qc/cav12_all_ok3_selection_manifest.csv`

The aligned arrays and current profiles are distributed with the companion
dataset. If they are not present in the checkout, set `VGIC_DATA_ROOT` to the
downloaded dataset directory or allow `shared.data_access` to resolve them
from the versioned data source.
