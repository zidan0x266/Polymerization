#!/bin/bash
# LAMMPS runs of the comparison as a SLURM job array, one node per run. From the run root:
#   python3 $REPO/examples/example_5/cases.py md > jobs_md.txt
#   mkdir -p logs && sbatch --array=1-$(wc -l < jobs_md.txt) $REPO/examples/example_5/slurm_md.sbatch
# Fill in the partition, account and tasks per node below; REPO, LMP (LAMMPS with PKG_PLUGIN), PLUGIN
# (polyreactplugin.so built against that LAMMPS) and MPIRUN (srun, or ibrun on TACC systems).
#SBATCH --job-name=poly-md
#SBATCH --partition=PARTITION
#SBATCH --account=ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
set -euo pipefail
: "${REPO:?set REPO to the repository}" "${LMP:?set LMP to the LAMMPS executable}"
export MPIRUN=${MPIRUN:-srun}
bash "$REPO/examples/example_5/run_job.sh" jobs_md.txt
