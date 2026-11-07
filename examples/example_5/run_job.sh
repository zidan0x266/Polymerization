#!/bin/bash
# Run line N of a job list in the run root (the current directory, or RUN_ROOT):
#   run_job.sh jobs.txt N        N counts from 1; default: $SLURM_ARRAY_TASK_ID
# Paths can be set in the environment: POLY_BUILD (build directory of the toolkit, with the plugin),
# LMP (LAMMPS executable), PLUGIN, MPIRUN (e.g. mpirun, srun, ibrun).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
list=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
n=${2:-${SLURM_ARRAY_TASK_ID:?give the line number}}
export EXAMPLE=$here
export POLY_BUILD=${POLY_BUILD:-$repo/build}
export LMP=${LMP:-lmp}
export PLUGIN=${PLUGIN:-$POLY_BUILD/src/lammps/polyreactplugin.so}
export MPIRUN=${MPIRUN:-mpirun}
line=$(sed -n "${n}p" "$list")
[ -n "$line" ] || { echo "no line $n in $list" >&2; exit 1; }
cd "${RUN_ROOT:-$PWD}"
eval "$line"
