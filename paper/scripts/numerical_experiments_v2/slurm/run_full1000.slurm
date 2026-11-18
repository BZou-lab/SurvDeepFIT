#!/usr/bin/env bash
#SBATCH --job-name=surv-full1000-effect08
#SBATCH --account=rc_bzou_pi
#SBATCH --partition=general
#SBATCH --array=1-12000%250
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --constraint=amd&avx512
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/full1000-%A_%a.out
#SBATCH --error=logs/full1000-%A_%a.err

# 2 scenarios x 2 training sizes x 3 correlations x 1,000 repetitions.
# Each task runs one dataset through all five methods, including full/selected
# prediction evaluation. The 24-hour limit applies separately to each task.
# Usage: sbatch [log overrides] run_full1000.slurm SIMULATION_ROOT STUDY_OUTPUT
# Use an absolute frozen source path and a new absolute output directory.
# Create the log directory before submission; SLURM opens logs before this script.
set -euo pipefail
if [[ $# -ne 2 || "$1" != /* || "$2" != /* ]]; then
  echo "Usage: sbatch run_full1000.slurm ABSOLUTE_SIMULATION_ROOT ABSOLUTE_STUDY_OUTPUT" >&2
  exit 2
fi
simulation_root="$1"
study_output="$2"
[[ "${SLURM_CPUS_PER_TASK:?}" -eq 1 ]] || exit 2
[[ "${SLURM_ARRAY_TASK_ID:?}" =~ ^[1-9][0-9]*$ ]] || exit 2
(( SLURM_ARRAY_TASK_ID <= 12000 )) || exit 2

export R_LIBS_USER=/proj/bzou_lab/users/Thomas/SurvDeepFIT/.r-library/rhel9-r-4.5.0
unset R_LIBS
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export RSCRIPT=/nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript

resource_root="$(dirname "$study_output")/task_resources"
mkdir -p "$resource_root"
resource_stem="$resource_root/task${SLURM_ARRAY_TASK_ID}"
hostname > "$resource_stem.host.txt"
lscpu > "$resource_stem.cpu.txt"

# Preserve the original deepTL build used by the successful 0.8 pilot.
"$RSCRIPT" --vanilla -e '
  stopifnot(normalizePath(dirname(find.package("deepTL"))) == normalizePath(Sys.getenv("R_LIBS_USER")),
            as.character(utils::packageVersion("deepTL")) == "1.0.0",
            !"workers" %in% names(formals(deepTL::ensemble_dnnet)))
  cat("Original deepTL:", find.package("deepTL"), "\n")
'

# The existing dispatcher maps the array ID to scenario, n, rho, and repetition.
# With 1,000 repetitions, it covers all 12,000 combinations exactly once.
/usr/bin/time -v -o "$resource_stem.time.txt" \
  bash "$simulation_root/slurm/run_test.slurm" "$simulation_root" "$study_output" 1000
