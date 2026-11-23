#!/bin/bash
#SBATCH --job-name=navcav_ok3_rmsf
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/analysis/all_ok3_rmsf_20261122/logs/navcav_ok3_rmsf-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/analysis/all_ok3_rmsf_20261122/logs/navcav_ok3_rmsf-%j.err
# Final-QC (All-OK-3) ensemble RMSF for Nav1.5 and Cav1.2, vanilla and current 50%
# masked cohorts together, straight from the 2026-11-11 aligned arrays on Hive.
# Replaces profiles whose vanilla rows were pasted in from an older profile file.
set -euo pipefail
R=/quobyte/yarovoygrp/ahgz/vgic_mutants
W=$R/analysis/all_ok3_rmsf_20261122
eval "$(micromamba shell hook --shell bash)"
micromamba activate bioadri
for pair in nav15:Nav1.5 cav12:Cav1.2; do
  ch=${pair%%:*}; SRC=$R/${pair##*:}/ensemble_rmsf_results_50pct
  mkdir -p "$W/$ch/dataRMSF/profiles"
  ln -sfn "$SRC/merged" "$W/$ch/dataRMSF/merged"
  ln -sfn "$SRC/references" "$W/$ch/dataRMSF/references"
  cp "$SRC/profiles/${ch}_all_models_per_residue_profiles.csv" "$W/$ch/dataRMSF/profiles/"
  cat "$SRC/merged/${ch}_merge_summary.json"
done
python "$W/code/recompute_all_ok3_local.py" --channel nav15 --repo-root "$W" \
  --dataset-alias wt_masked_v1_50pct=wt_masked --dataset-alias wt_masked_v2_50pct=wt_masked_v2 \
  --dataset-alias qqq_masked_v1_50pct=qqq_masked --dataset-alias qqq_masked_v2_50pct=qqq_masked_v2
python "$W/code/recompute_all_ok3_local.py" --channel cav12 --repo-root "$W" \
  --dataset-alias wt_masked_v3_50pct=wt_masked --dataset-alias g402s_masked_v3_50pct=g402s_masked \
  --dataset-alias g406r_masked_v3_50pct=g406r_masked --ignore-dataset e1115k_masked_v3_50pct
ls -l "$W"/*/dataRMSF/profiles
