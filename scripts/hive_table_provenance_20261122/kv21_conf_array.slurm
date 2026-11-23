#!/bin/bash
#SBATCH --job-name=21_conf6
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --array=0-99%25
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/conf-%A_%a.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/conf-%A_%a.err
# Conformation metrics for whatever the prepare step could not link (WT and F412L
# masked_common99_50pct). The config path is passed explicitly: the 2026-11-20 WT
# run inherited a CONFIG variable from the submitting shell and every part failed.
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
  --output-dir "$OUT" --task-id "$SLURM_ARRAY_TASK_ID" --task-count 100 \
  --workers 1 --resume
