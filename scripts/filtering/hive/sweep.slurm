#!/bin/bash
#SBATCH --job-name=fc_sweep
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/sweep-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/sweep-%j.err
# Every cohort through the four steps, and then through eleven filtering setups, after the measure
# and judge arrays have written their tables. Nothing here reads a model: the sweep is eleven passes
# over four CSVs per cohort, which is why trying another cutoff is cheap. The memory is for the
# distance tables, a few hundred MB each and read whole by pandas.
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

echo "started $(date) on $(hostname)"
python "$S/hive_cohorts.py" chain --all
python "$S/hive_cohorts.py" sweep --all
echo "finished $(date)"
