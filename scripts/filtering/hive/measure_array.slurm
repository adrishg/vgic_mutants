#!/bin/bash
#SBATCH --job-name=fc_measure
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=0-11%5
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/measure-%A_%a.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/measure-%A_%a.err
# One array task per cohort: measure the pore-VSD interface of every model that cohort generated,
# all eleven recycles, not only the convergence representatives - a trajectory fails this check on
# any of its retained recycles, so all of them have to be read. That is about 5,500 PDBs a cohort
# against roughly 500 for the judge step, hence the longer wall time. CPU only on "high", so it
# does not touch the gres/gpu cap the AF2 arrays use, and each task writes its own file.
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
python "$S/hive_cohorts.py" measure --index "${SLURM_ARRAY_TASK_ID}"
echo "finished $(date)"
