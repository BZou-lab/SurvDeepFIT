#!/usr/bin/env bash
#SBATCH --job-name=surv-test20
#SBATCH --array=1-240%200
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --constraint=amd&avx512
#SBATCH --nodes=1
#SBATCH --time=24:00:00
# One saved dataset shared by all methods for each setting/repetition.
set -euo pipefail
simulation_root="$1"
study_output="$2"
[[ "${SLURM_CPUS_PER_TASK:?}" -eq 1 ]] || exit 2
export R_LIBS_USER=/proj/bzou_lab/users/Thomas/SurvDeepFIT/.r-library/rhel9-r-4.5.0
unset R_LIBS
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export RSCRIPT=/nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript
resource_root="$(dirname "$study_output")/task_resources"
mkdir -p "$resource_root"
resource_stem="$resource_root/task${SLURM_ARRAY_TASK_ID:?}"
hostname > "$resource_stem.host.txt"
lscpu > "$resource_stem.cpu.txt"
# Guard against accidentally selecting the experimental deepTL installation.
"$RSCRIPT" --vanilla -e '
  stopifnot(normalizePath(dirname(find.package("deepTL"))) == normalizePath(Sys.getenv("R_LIBS_USER")),
            as.character(utils::packageVersion("deepTL")) == "1.0.0",
            !"workers" %in% names(formals(deepTL::ensemble_dnnet)))
  cat("Original deepTL:", find.package("deepTL"), "\n")
'
/usr/bin/time -v -o "$resource_stem.time.txt" \
  bash "$simulation_root/slurm/run_test.slurm" "$simulation_root" "$study_output" 20
