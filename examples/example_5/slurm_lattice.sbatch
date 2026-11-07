#!/bin/bash
# 3dbfm and statistical runs of the comparison as a SLURM job array, one core per run. From the run
# root (an empty directory on the scratch file system):
#   python3 $REPO/examples/example_5/cases.py lattice > jobs_lattice.txt
#   mkdir -p logs && sbatch --array=1-$(wc -l < jobs_lattice.txt) $REPO/examples/example_5/slurm_lattice.sbatch
# Fill in the partition and account below, and REPO (and POLY_BUILD if the build directory is not
# $REPO/build) with --export or here.
#SBATCH --job-name=poly-lattice
#SBATCH --partition=PARTITION
#SBATCH --account=ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
set -euo pipefail
: "${REPO:?set REPO to the repository}"
module list 2>/dev/null || true
bash "$REPO/examples/example_5/run_job.sh" jobs_lattice.txt
