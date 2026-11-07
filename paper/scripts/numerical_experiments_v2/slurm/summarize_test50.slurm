#!/usr/bin/env bash
set -euo pipefail
simulation_root="$1"
study_output="$2"
export R_LIBS_USER=/proj/bzou_lab/users/Thomas/SurvDeepFIT/.r-library/rhel9-r-4.5.0
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export RSCRIPT=/nas/longleaf/rhel9/apps/r/4.5.0/bin/Rscript
cd "$simulation_root"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${RSCRIPT:?}" --vanilla summarize_results.R --input "$study_output" --output "$study_output/summary"
"${RSCRIPT:?}" --vanilla make_feature_selection_tables.R --input "$study_output/summary" --output "$study_output/summary"
