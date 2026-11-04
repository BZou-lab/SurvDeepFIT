#!/usr/bin/env bash
# Submit from numerical_experiments_v2; additional arguments go to the diagnostic.
#SBATCH --job-name=rsf-diagnostic
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
simulation_root="${SIMULATION_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
project_root="$(cd "$simulation_root/../../../.." && pwd)"
export R_LIBS_USER="$project_root/.r-library/rhel9-r-4.5.0"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
diagnostic_arguments=("$@")
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  diagnostic_arguments+=(--rep "$SLURM_ARRAY_TASK_ID")
fi
"${RSCRIPT:-/nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript}" --vanilla \
  "$simulation_root/run_rsf_diagnostic.R" \
  --config "$simulation_root/configs/rsf_diagnostic.yaml" --threads "$OMP_NUM_THREADS" "${diagnostic_arguments[@]}"
