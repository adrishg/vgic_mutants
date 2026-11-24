#!/bin/bash
#SBATCH --job-name=fc_verdicts
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/verdicts-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/verdicts-%j.err
# The per-prediction pass/fail table, one CSV per dataset. Nothing here reads a model: it walks the
# distance table, the convergence table, the interface measures and the topology verdicts, all of
# which earlier steps wrote. The memory is for the distance tables, a few hundred MB each and read
# whole by pandas. Pass the cohorts as arguments, or --all.
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

echo "started $(date) on $(hostname): verdicts $*"
python "$S/hive_cohorts.py" verdicts "$@"
echo "finished $(date)"
