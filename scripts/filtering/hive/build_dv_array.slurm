#!/bin/bash
#SBATCH --job-name=st1_build
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-6
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/step1-%A_%a.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/step1-%A_%a.err
# Step 1 of Brandon's pipeline (build_dv.py) on all seven datasets, one array task each. Each task
# reads its two or three ColabFold archives whole, numbers every model against UniProt, superposes
# on the fit set and writes coords.npy and dv.csv to that dataset's own output_dir, so no two tasks
# write the same file. 64G is for the coordinate block the step holds in memory; the step is
# single-threaded Python over gemmi, and the four CPUs are only for numpy's BLAS.
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

set +u
eval "$($MAMBA shell hook --shell bash)"
micromamba activate "$PREFIX"
set -u

echo "step 1 on $NAME, started $(date) on $(hostname)"
python "$P/build_dv.py" --config "$D/$NAME/$NAME.toml"
echo "finished $(date)"
