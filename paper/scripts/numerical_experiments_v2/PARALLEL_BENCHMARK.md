This benchmark keeps the scenario YAML files and scientific workflow unchanged.
It runs only SurvDeepFIT: scenarios 1 and 2, n = 1000 and 3000, rho = 0.3,
repetitions 1 through 3, and 1/2/4 ensemble workers (36 executions).
Numerical libraries use one thread per process in all configurations.

The original package_sources/deepTL and .r-library are preserved. Experimental
source is package_sources/deepTL-parallel; installation is isolated in
.r-library-parallel/rhel9-r-4.5.0. The package name remains deepTL, with version
1.0.0.9001. R library ordering selects the experimental installation only for
these benchmark jobs. Existing dependencies are reused without upgrades.

The experimental package defaults to the original sequential algorithm. The
benchmark explicitly enables per-member L'Ecuyer RNG streams with identical
bootstrap samples and training streams for 1/2/4 workers. This seeded mode is
reproducible across worker counts, but does not reproduce the original seed
sequence. Network fits run in forked processes on Linux; aggregation retains
the original ordering and selection rules. Permutation evaluation is unchanged.

Setup compiles the package on a compute node and runs small implementation
checks against the original installed package and across worker counts. The
smoke checks use small synthetic inputs; their settings never reach the full
benchmark. Benchmark arrays depend on successful setup and validation.

Each submission preserves a source snapshot, source checksums, task manifest,
submission record, logs, and outputs in its own results/parallel-benchmark-*
directory. Arrays request 1/2/4 CPUs respectively, 16 GB per job, and 24 hours.
Each array allows up to four simultaneous tasks. These are runtime comparisons,
not a controlled fixed-total-CPU throughput experiment; compare allocated
core-hours alongside elapsed time. CPU models and placement may affect timings;
the task metadata records the host and GNU time records CPU time and process
maximum RSS. SLURM accounting should be consulted for job resource use.

tasks/wN/taskI/ensembles.csv records each ensemble's training and total ensemble
time. Task metadata records total workflow time; result.rds records selection
and prediction times. Original scientific checkpoints are stored under runs/wN.
The summary job runs after the arrays finish, including failures. It writes
timings.csv, agreement.csv, and runtime_summary.csv. Agreement checks compare
datasets, feature selection, RNG metadata, and prediction metrics exactly,
excluding execution times. Saved scientific outputs should be examined before
any production use. Summary fails if a task is missing, failed, or disagrees.

To rerun summary manually, use the project R environment and run:

    Rscript summarize_parallel_benchmark.R /absolute/path/to/submission/directory

No original SEER source, settings, package installation, or results are changed.
