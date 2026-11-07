# Wrapper intentionally lives beside run_simulation.R, whose --file root lookup
# therefore resolves to the same directory. Scientific workflow is reused intact.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
benchmark_root <- dirname(normalizePath(script))
workers <- as.integer(Sys.getenv("SURVDEEPFIT_WORKERS"))
task_directory <- Sys.getenv("BENCHMARK_TASK_DIR")
stopifnot(workers %in% c(1L, 2L, 4L), nzchar(task_directory))
dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)
expected_library <- normalizePath(Sys.getenv("PARALLEL_R_LIBRARY"))
stopifnot(identical(normalizePath(dirname(find.package("deepTL"))), expected_library),
          as.character(utils::packageVersion("deepTL")) == "1.0.0.9001",
          all(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")) == "1"))
options(deepTL.ensemble.workers = workers, deepTL.ensemble.reproducible = TRUE,
        deepTL.ensemble.timing_file = file.path(task_directory, "ensembles.csv"))
metadata <- list(workers = workers, numerical_threads = 1L,
  package_path = find.package("deepTL"), package_version = utils::packageVersion("deepTL"),
  arguments = commandArgs(), session = utils::sessionInfo(),
  started_utc = format(Sys.time(), tz = "UTC"), host = Sys.info()[["nodename"]],
  slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
  slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID"))
saveRDS(metadata, file.path(task_directory, "metadata.rds"))
started <- proc.time()[["elapsed"]]
source(file.path(benchmark_root, "run_simulation.R"))
metadata$total_seconds <- proc.time()[["elapsed"]] - started
metadata$finished_utc <- format(Sys.time(), tz = "UTC")
saveRDS(metadata, file.path(task_directory, "metadata.rds"))
result <- readRDS(Sys.getenv("BENCHMARK_RESULT_FILE"))
if (result$selection_status == "failed" || is.null(result$predictions) ||
    any(result$predictions$status == "failed"))
  stop("Benchmark scientific stage failed; inspect result.rds and task log")
cat("Benchmark completed with", workers, "workers in", metadata$total_seconds, "seconds\n")
