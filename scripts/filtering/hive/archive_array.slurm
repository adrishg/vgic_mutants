#!/bin/bash
#SBATCH --job-name=fc_archive
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-17%6
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/archive-%A_%a.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/archive-%A_%a.err
# One array task per condition: pack its saved-recycle PDBs into the single .tar.gz that step 1 of
# Brandon's pipeline reads. Each task writes its own archive, which is what Quobyte wants, and
# writes through a .part so an interrupted task never leaves something that looks finished. Six at
# a time because this is read-heavy on the shared filesystem, not because of a CPU limit.
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

echo "condition index ${SLURM_ARRAY_TASK_ID} on $(hostname), started $(date)"
python "$S/archive_cohorts.py" pack --index "${SLURM_ARRAY_TASK_ID}"
echo "finished $(date)"
