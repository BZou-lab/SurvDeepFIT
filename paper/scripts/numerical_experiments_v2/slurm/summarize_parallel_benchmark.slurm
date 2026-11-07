#!/usr/bin/env bash
set -euo pipefail
project_root="$1"
benchmark_root="$2"
run_root="$3"
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export R_LIBS_USER="$project_root/.r-library/rhel9-r-4.5.0"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset R_LIBS
/nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript --vanilla \
  "$benchmark_root/summarize_parallel_benchmark.R" "$run_root"
