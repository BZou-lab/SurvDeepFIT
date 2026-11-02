#!/usr/bin/env bash
# Set resources and the array range with sbatch command-line options (see README).
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch [resources] run_simulation.slurm CONFIG OUTPUT [extra R options]" >&2
  exit 2
fi

# SLURM copies scripts into a spool directory; use the submission directory.
simulation_root="${SIMULATION_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
config="$1"
output="$2"
shift 2

# Load the site's R module/environment before submission, or set RSCRIPT explicitly.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"

"${RSCRIPT:-Rscript}" "$simulation_root/run_simulation.R" \
  --config "$config" --output "$output" \
  --rep "${SLURM_ARRAY_TASK_ID:?Submit this script as a job array}" \
  --threads "$OMP_NUM_THREADS" "$@"
