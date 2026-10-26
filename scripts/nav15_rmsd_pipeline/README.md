# Nav1.5 recycle-RMSD convergence

This directory contains the modular cluster launchers for the four new Nav1.5
ensembles:

```text
WT  / masked_v1new
WT  / masked_v2new
QQQ / masked_v1new
QQQ / masked_v2new
```

The calculation follows the established Nav1.5 workflow recorded in the
existing `run_config.json` files:

- `af2_recycle_rmsd_convergence.py`
- four worker processes per run
- stable-core Cα blocks: 132–282, 359–426, 505–665, 690–749,
  870–1112, 1128–1164, 1200–1405, and 1429–1466
- minimum aligned coverage: 0.90
- configured successive-transition threshold: 3.0 Å
- sensitivity thresholds: 3, 4, 5, 6, and 8 Å

Each wrapper runs one ensemble only. They do not submit Slurm arrays, start
SSH agents, or launch the other three calculations.

## Cluster use

Copy this directory to the cluster scripts directory, then from the copied
directory run one wrapper at a time:

```bash
./run_wt_masked_v1new.sh
./run_wt_masked_v2new.sh
./run_qqq_masked_v1new.sh
./run_qqq_masked_v2new.sh
```

For conservative cluster execution, submit the four jobs as an `afterok`
chain so only one four-worker run is active at a time:

```bash
./submit_nav15_rmsd_chain.sh
squeue --me
```

The Slurm resource assumptions match the established cluster package
(`high`, `genome-center-grp`, four CPUs, 8 GB RAM). Adjust them in
`job_nav15_rmsd.sbatch` if the cluster policy has changed.

The output directories are:

```text
/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5/Nav15_rmsd_convergence_wt_masked_v1new
/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5/Nav15_rmsd_convergence_wt_masked_v2new
/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5/Nav15_rmsd_convergence_qqq_masked_v1new
/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5/Nav15_rmsd_convergence_qqq_masked_v2new
```

The wrappers stop if the input directory is missing, if the convergence
program is not on `PATH`, or if the output directory already exists. To rerun
an existing calculation, set `ALLOW_EXISTING_OUTPUT=1` explicitly after
checking that the previous result is incomplete.

