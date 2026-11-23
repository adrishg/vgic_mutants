#!/bin/bash
#SBATCH --job-name=21_conf6_compact
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/compact-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/compact-%j.err
# Laptop-sized copies of the two largest merged tables: only the rows and columns
# kv21/run_l403a_conformational_validation.py reads. The full tables stay on Hive.
set -euo pipefail
export PATH=/quobyte/yarovoygrp/ahgz/mamba_envs/bin:$PATH
python - <<'PY'
import pandas as pd
out = "/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122"
use = ["structure_id", "source_type", "condition", "protocol", "dataset", "seed", "model_number",
       "recycle", "canonical_subunit", "is_final_model", "paper_residue", "phi_deg", "psi_deg",
       "incoming_i5_O_N_distance_A", "incoming_i4_O_N_distance_A", "incoming_i5_O_H_N_angle_deg",
       "incoming_i4_O_H_N_angle_deg", "incoming_geometric_pi_like", "geometric_pi_like"]
first = True
rows = 0
for chunk in pd.read_csv(f"{out}/kv21_l403a_conformational_metrics_residue_level.csv", usecols=use,
                         chunksize=200000, low_memory=False):
    chunk = chunk[chunk.paper_residue.between(407, 411)][use]
    chunk.to_csv(f"{out}/kv21_l403a_conformational_metrics_residue_level_pi_window.csv",
                 mode="w" if first else "a", header=first, index=False)
    first = False
    rows += len(chunk)
print("pi-window residue rows", rows)
qcols = ["source_type", "condition", "protocol", "category", "status", "detail", "structure_id",
         "frame_orientation_score"]
first = True
for chunk in pd.read_csv(f"{out}/kv21_l403a_conformational_metrics_qc.csv", usecols=qcols,
                         chunksize=200000, low_memory=False):
    chunk[qcols].to_csv(f"{out}/kv21_l403a_conformational_metrics_qc_compact.csv",
                        mode="w" if first else "a", header=first, index=False)
    first = False
PY
ls -l /quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/*.csv
