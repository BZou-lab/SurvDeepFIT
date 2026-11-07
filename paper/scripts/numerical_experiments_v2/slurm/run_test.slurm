#!/usr/bin/env bash
# One dataset and all methods per task, with repetitions supplied at submission.
set -euo pipefail
simulation_root="$1"
study_output="$2"
repetitions="${3:-5}"
[[ "$repetitions" =~ ^[1-9][0-9]*$ ]] || exit 2
export R_LIBS_USER=/proj/bzou_lab/users/Thomas/SurvDeepFIT/.r-library/rhel9-r-4.5.0
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export RSCRIPT=/nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript
task_index=$((SLURM_ARRAY_TASK_ID - 1))
if (( task_index < 0 || task_index >= 12 * repetitions )); then exit 2; fi
repetition=$((task_index % repetitions + 1))
setting_index=$((task_index / repetitions))
training_sizes=(1000 3000)
correlations=(0 0.3 0.5)
scenario=$((setting_index / 6 + 1))
training_n="${training_sizes[$((setting_index % 2))]}"
rho="${correlations[$((setting_index / 2 % 3))]}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
cd "$simulation_root"
exec "$RSCRIPT" --vanilla run_simulation.R \
  --config "configs/scenario${scenario}.yaml" --output "$study_output" \
  --repetitions "$repetitions" --rep "$repetition" --n "$training_n" --rho "$rho" \
  --threads "$OMP_NUM_THREADS"
