#!/bin/bash
#SBATCH --job-name=fc_compare
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/compare-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/compare-%j.err
# The November package against the chain's current survivors. Reads CSVs and hashes four pairs of
# convergence tables; the only directory walk is the package's four models/ folders, 500 files each.
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
python "$S/compare_nav15_old_vs_new.py" "$@"
echo "--- top-level entries of Nav1.5/forDiego6.tar.gz, to see what else went out"
tar -tzf /quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5/forDiego6.tar.gz 2>/dev/null | head -40 || true
echo "finished $(date)"
