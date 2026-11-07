#!/usr/bin/env bash
set -euo pipefail
project_root="$1"
benchmark_root="$2"
run_root="$3"
workers="$4"
[[ -f "$run_root/VALIDATED" ]] || { echo "Setup validation has not passed" >&2; exit 1; }
[[ "$workers" =~ ^(1|2|4)$ ]] || exit 2
[[ "$SLURM_CPUS_PER_TASK" -eq "$workers" ]] || exit 2
index=$((SLURM_ARRAY_TASK_ID - 1))
(( index >= 0 && index < 12 )) || exit 2
repetition=$((index % 3 + 1))
setting=$((index / 3))
scenario=$((setting / 2 + 1))
sizes=(1000 3000)
training_n="${sizes[$((setting % 2))]}"
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export PARALLEL_R_LIBRARY="$project_root/.r-library-parallel/rhel9-r-4.5.0"
export R_LIBS_USER="$PARALLEL_R_LIBRARY:$project_root/.r-library/rhel9-r-4.5.0"
unset R_LIBS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SURVDEEPFIT_WORKERS="$workers"
export BENCHMARK_TASK_DIR="$run_root/tasks/w${workers}/task${SLURM_ARRAY_TASK_ID}"
output="$run_root/runs/w${workers}"
printf -v rep_directory 'rep%04d' "$repetition"
export BENCHMARK_RESULT_FILE="$output/scenario${scenario}/n${training_n}/rho0.3/$rep_directory/SurvDeepFIT/result.rds"
mkdir -p "$BENCHMARK_TASK_DIR"
cd "$benchmark_root"
/usr/bin/time -v -o "$BENCHMARK_TASK_DIR/resources.txt" \
  /nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript --vanilla run_parallel_benchmark.R \
  --config "configs/scenario${scenario}.yaml" --output "$output" \
  --repetitions 3 --rep "$repetition" --n "$training_n" --rho 0.3 \
  --methods SurvDeepFIT --threads 1
