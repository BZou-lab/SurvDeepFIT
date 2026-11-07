args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
run_root <- normalizePath(args[1])
rows <- list()
comparison <- list()
for (task in 1:12) {
  index <- task - 1L
  repetition <- index %% 3L + 1L
  setting <- index %/% 3L
  scenario <- setting %/% 2L + 1L
  n <- c(1000L, 3000L)[setting %% 2L + 1L]
  scientific <- vector("list", 3)
  dataset <- vector("list", 3)
  for (j in seq_along(c(1L, 2L, 4L))) {
    workers <- c(1L, 2L, 4L)[j]
    task_dir <- file.path(run_root, "tasks", paste0("w", workers), paste0("task", task))
    result_dir <- file.path(run_root, "runs", paste0("w", workers),
      paste0("scenario", scenario), paste0("n", n), "rho0.3", sprintf("rep%04d", repetition))
    result_file <- file.path(result_dir, "SurvDeepFIT/result.rds")
    metadata_file <- file.path(task_dir, "metadata.rds")
    row <- data.frame(task = task, scenario = scenario, n = n, rho = 0.3,
      repetition = repetition, workers = workers, status = "missing",
      total_seconds = NA_real_, selection_seconds = NA_real_,
      prediction_seconds = NA_real_, training_seconds = NA_real_,
      allocated_core_hours = NA_real_)
    if (file.exists(result_file) && file.exists(metadata_file)) {
      result <- readRDS(result_file)
      metadata <- readRDS(metadata_file)
      row$status <- if (result$selection_status == "failed" || is.null(result$predictions) ||
        any(result$predictions$status == "failed")) "failed" else "ok"
      row$total_seconds <- if (is.null(metadata$total_seconds)) NA_real_ else metadata$total_seconds
      row$selection_seconds <- result$selection_elapsed_seconds
      if (!is.null(result$predictions)) row$prediction_seconds <- sum(result$predictions$elapsed_seconds)
      timing_file <- file.path(task_dir, "ensembles.csv")
      if (file.exists(timing_file)) row$training_seconds <- sum(read.csv(timing_file)$training_seconds)
      row$allocated_core_hours <- workers * row$total_seconds / 3600
      scientific[[j]] <- result[c("selection_status", "selection_error", "selection_warnings",
        "features", "feature_metrics", "candidate_features", "important_features", "rng")]
      prediction <- result$predictions
      if (!is.null(prediction)) prediction$elapsed_seconds <- NULL
      scientific[[j]]$predictions <- prediction
      dataset[[j]] <- readRDS(file.path(result_dir, "data.rds"))$value
    }
    rows[[length(rows) + 1L]] <- row
  }
  for (j in 2:3) {
    complete <- !is.null(scientific[[1]]) && !is.null(scientific[[j]])
    comparison[[length(comparison) + 1L]] <- data.frame(task = task, scenario = scenario,
      n = n, repetition = repetition, workers = c(2L, 4L)[j - 1L], complete = complete,
      identical_data = if (complete) identical(dataset[[1]], dataset[[j]]) else NA,
      identical_outputs = if (complete) identical(scientific[[1]], scientific[[j]]) else NA)
  }
}
timings <- do.call(rbind, rows)
agreement <- do.call(rbind, comparison)
baseline <- timings$total_seconds[timings$workers == 1L]
timings$speedup <- baseline[match(timings$task, timings$task[timings$workers == 1L])] / timings$total_seconds
write.csv(timings, file.path(run_root, "timings.csv"), row.names = FALSE)
write.csv(agreement, file.path(run_root, "agreement.csv"), row.names = FALSE)
valid <- timings[timings$status == "ok" & is.finite(timings$total_seconds), ]
if (nrow(valid)) {
  summary <- aggregate(cbind(total_seconds, selection_seconds, prediction_seconds,
    training_seconds, allocated_core_hours, speedup) ~ scenario + n + workers,
    valid, mean, na.action = na.pass)
  write.csv(summary, file.path(run_root, "runtime_summary.csv"), row.names = FALSE)
  print(summary)
}
print(agreement)
if (any(timings$status != "ok") || anyNA(agreement$identical_outputs) ||
    !all(agreement$identical_data) || !all(agreement$identical_outputs))
  stop("Some tasks are missing, failed, or disagree. Inspect timings.csv and agreement.csv")
cat("All 36 executions completed; all 24 paired comparisons agree exactly.\n")
