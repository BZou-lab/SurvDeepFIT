#!/usr/bin/env bash
set -euo pipefail
project_root="$1"
benchmark_root="$2"
run_root="$3"
original_library="$project_root/.r-library/rhel9-r-4.5.0"
parallel_library="$project_root/.r-library-parallel/rhel9-r-4.5.0"
export LD_LIBRARY_PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/lib64:${LD_LIBRARY_PATH:-}"
export PATH="/nas/longleaf/rhel9/apps/gcc/12.2.0/bin:/nas/longleaf/rhel9/apps/r/4.5.0/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# R CMD INSTALL also launches R; prevent personal startup files from replacing
# the project library with a personal library for another R version.
export R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null
unset R_LIBS
export R_LIBS_USER="$original_library"
mkdir -p "$parallel_library"
exec 9>"$parallel_library/.install.lock"
flock -n 9
# Compile a disposable copy; preserve the archived source snapshot.
build_root="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/deeptl-parallel-build.XXXXXX")"
trap 'rm -rf "$build_root"' EXIT
cp -R "$run_root/source/deepTL-parallel" "$build_root/deepTL"
chmod -R u+w "$build_root/deepTL"
R CMD INSTALL --preclean --library="$parallel_library" "$build_root/deepTL"
Rscript --vanilla "$benchmark_root/validate_parallel_benchmark.R" original "$run_root/original_smoke.rds"
export R_LIBS_USER="$parallel_library:$original_library"
Rscript --vanilla "$benchmark_root/validate_parallel_benchmark.R" parallel "$run_root/original_smoke.rds"
touch "$run_root/VALIDATED"
