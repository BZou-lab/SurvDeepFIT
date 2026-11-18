#!/usr/bin/env bash
#SBATCH --job-name=surv-scenario-summary
#SBATCH --account=rc_bzou_pi
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# Usage: sbatch [log paths] summarize_scenario.slurm SCRIPT_ROOT STUDY_RUNS SCENARIO OUTPUT
set -euo pipefail
[[ $# -eq 4 ]] || exit 2
export R_LIBS_USER=/proj/bzou_lab/users/Thomas/SurvDeepFIT/.r-library/rhel9-r-4.5.0
unset R_LIBS
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript --vanilla \
  "$1/summarize_scenario.R" --input "$2" --scenario "$3" --output "$4"
