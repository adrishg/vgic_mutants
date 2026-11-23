#!/bin/bash
#SBATCH --job-name=st23_stats
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-6
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/step23-%A_%a.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/step23-%A_%a.err
# The topology screen and steps 2 and 3 of Brandon's pipeline for one dataset. screen.py runs
# first: step 2 refuses to start without its results, because the config has a [screen]. Then
# region_dv.py cuts the step 1 DVs into the five regions and applies [screen], including our
# filter chain's omit file, and stats.py fits
# log(value) ~ condition * recycle + network with a condition:network random effect and runs the
# Dunnett contrasts against the control. All three read what step 1 wrote for that dataset only, so
# tasks share no file. Submit this with a dependency on the step 1 array.
set -euo pipefail

R=/quobyte/yarovoygrp/ahgz/vgic_mutants
P=$R/scripts/filter_chain/filter_stats/pipeline
D=$R/analysis/stats_inputs/datasets
MAMBA=/quobyte/yarovoygrp/ahgz/apps/micromamba/bin/micromamba
PREFIX=/quobyte/yarovoygrp/ahgz/apps/vgic-stats
export MAMBA_ROOT_PREFIX=/quobyte/yarovoygrp/ahgz/.mamba_envs
export TMPDIR=/quobyte/yarovoygrp/ahgz/.tmp
mkdir -p "$TMPDIR" "$R/scripts/filter_chain/logs"

DATASETS=(kv21_vanilla kv21_masked99 nav15_vanilla nav15_masked_v1 nav15_masked_v2 cav12_vanilla cav12_masked_v3)
NAME=${DATASETS[$SLURM_ARRAY_TASK_ID]}
CONFIG=$D/$NAME/$NAME.toml

set +u
eval "$($MAMBA shell hook --shell bash)"
micromamba activate "$PREFIX"
set -u

echo "screen and steps 2-3 on $NAME, started $(date) on $(hostname)"
python "$P/screen.py" --config "$CONFIG"
python "$P/region_dv.py" --config "$CONFIG"
python "$P/stats.py" --config "$CONFIG"
echo "finished $(date)"
