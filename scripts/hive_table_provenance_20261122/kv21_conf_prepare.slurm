#!/bin/bash
#SBATCH --job-name=21_conf6_prep
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/prep-%j.out
#SBATCH --error=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122/logs/prep-%j.err
# Manifest for the six current Kv2.1 cohorts + links to the parts already computed.
set -euo pipefail
OUT=/quobyte/yarovoygrp/ahgz/vgic_mutants/Kv2.1/conformation_c99/six_cohorts_20261122
export PATH=/quobyte/yarovoygrp/ahgz/mamba_envs/bin:$PATH
python "$OUT/code/kv21_conf_prepare.py" "$OUT"
