#!/bin/bash
#SBATCH --job-name=fc_statsenv
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/statsenv-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/filter_chain/logs/statsenv-%j.err
# The environment Brandon's pipeline needs, which is not the ColabFold one: that has Python 3.10 and
# neither gemmi (step 1 reads structures with it) nor statsmodels (the mixed model). His
# requirements.txt asks for 3.11 or later. Solving and downloading an environment is real work, so
# it happens here and not on a login node.
set -euo pipefail

MAMBA=/quobyte/yarovoygrp/ahgz/apps/micromamba/bin/micromamba
PREFIX=/quobyte/yarovoygrp/ahgz/apps/vgic-stats
export MAMBA_ROOT_PREFIX=/quobyte/yarovoygrp/ahgz/.mamba_envs
export TMPDIR=/quobyte/yarovoygrp/ahgz/.tmp
mkdir -p "$TMPDIR"

"$MAMBA" create --yes --prefix "$PREFIX" -c conda-forge \
    python=3.12 numpy scipy pandas statsmodels gemmi biopython pytest

set +u
eval "$($MAMBA shell hook --shell bash)"
micromamba activate "$PREFIX"
set -u
python - <<'PY'
import gemmi, numpy, pandas, scipy, statsmodels, Bio, sys
print("python", sys.version.split()[0])
for m in (numpy, scipy, pandas, statsmodels, gemmi, Bio):
    print(f"  {m.__name__} {getattr(m, '__version__', '?')}")
PY
echo "environment at $PREFIX"
