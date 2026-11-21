#!/bin/bash
#SBATCH --job-name=fc_judge
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=0-9%5
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/judge-%A_%a.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/judge-%A_%a.err
# One array task per cohort: judge that cohort's convergence representatives with the topology
# screen. CPU only on "high", so it does not touch the gres/gpu cap the AF2 arrays use. Each task
# reads its own cohort and writes its own verdict file - one writer per file, as Quobyte needs.
set -euo pipefail

S=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain
MAMBA=/quobyte/yarovoygrp/ahgz/apps/micromamba/bin/micromamba
COLABFOLD_ENV=/quobyte/yarovoygrp/ahgz/apps/localcolabfold/colabfold-conda
export XDG_CACHE_HOME=/quobyte/yarovoygrp/ahgz/.cache
export TMPDIR=/quobyte/yarovoygrp/ahgz/.tmp
export MAMBA_ROOT_PREFIX=/quobyte/yarovoygrp/ahgz/.mamba_envs
mkdir -p "$TMPDIR" "$S/logs"

set +u
if ! command -v module >/dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
module load gcc/13.2.0
eval "$($MAMBA shell hook --shell bash)"
micromamba activate "$COLABFOLD_ENV"
set -u

echo "cohort index ${SLURM_ARRAY_TASK_ID} on $(hostname), started $(date)"
python "$S/hive_cohorts.py" judge --index "${SLURM_ARRAY_TASK_ID}"
echo "finished $(date)"
