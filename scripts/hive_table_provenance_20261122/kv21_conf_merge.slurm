#!/bin/bash
#SBATCH --job-name=21_conf6_merge
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/merge-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/merge-%j.err
# Merge the per-structure parts of the six current cohorts into the pipeline tables.
set -euo pipefail
unset CONFIG
R=/quobyte/yarovoygrp/ahgz/vgic_mutants
P=$R/scripts/kv21_conformation_pipeline_v2
OUT=$R/Kv2.1/conformation_c99/six_cohorts_20261122
export PATH=/quobyte/yarovoygrp/ahgz/mamba_envs/bin:$PATH
python "$P/scripts/analyze_kv21_l403a_conformation.py" \
  --repo-root "$R" --model-list "$OUT/manifest_six_cohorts.csv" \
  --experimental-wt "$R/Kv2.1/experimental/8SD3.pdb" \
  --experimental-l403a "$R/Kv2.1/experimental/8SDA.pdb" \
  --config "$P/config/kv21_l403a_conformation.yaml" \
  --output-dir "$OUT" --merge-only
python - <<'PY'
import json, pandas as pd
out = "/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122"
qc = pd.read_csv(f"{out}/kv21_l403a_conformational_metrics_qc.csv", usecols=["dataset", "category", "status"], low_memory=False)
print(qc[qc.category.isin(["structure", "merge"])].groupby(["dataset", "category", "status"]).size().to_string())
chain = pd.read_csv(f"{out}/kv21_l403a_conformational_metrics_chain_resolved.csv", usecols=["dataset", "structure_id", "source_path"], low_memory=False)
print(chain.groupby("dataset").structure_id.nunique().to_string())
print(chain.source_path.str.replace(r"/(shard_[^/]+/)?[^/]+$", "", regex=True).value_counts().to_string())
PY
