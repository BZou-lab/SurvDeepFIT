summary_root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                     winslash = "/"))
source(file.path(summary_root, "R/workflow.R"))

summary_arguments <- function(arguments) {
  if (length(arguments) != 4 || !identical(arguments[c(1, 3)], c("--input", "--output"))) {
    stop("Usage: Rscript summarize_results.R --input RESULTS_DIRECTORY --output SUMMARY_DIRECTORY")
  }
  list(input = arguments[2], output = arguments[4])
}

record_context <- function(record) {
  data.frame(scenario = record$context$scenario, n = record$context$n_train,
             n_test = record$context$n_test, rho = record$context$rho,
             repetition = record$context$repetition, seed = record$context$base_seed,
             Method = record$method, signature = record$signature,
             selection_status = record$selection_status, selection_error = record$selection_error,
             selection_warnings = paste(record$selection_warnings, collapse = " | "))
}

summary_stat <- function(values) {
  valid <- is.finite(values)
  c(mean = if (any(valid)) mean(values[valid]) else NA_real_,
    sd = if (sum(valid) > 1) stats::sd(values[valid]) else NA_real_,
    n_valid = sum(valid), n_missing = sum(!valid))
}

group_summary <- function(data, keys, metrics) {
  groups <- split(seq_len(nrow(data)), interaction(data[, keys, drop = FALSE], drop = TRUE))
  do.call(rbind, lapply(groups, function(indices) {
    subset <- data[indices, , drop = FALSE]
    row <- subset[1, keys, drop = FALSE]
    row$n_attempted <- nrow(subset)
    for (metric in metrics) {
      statistics <- summary_stat(subset[[metric]])
      for (statistic in names(statistics)) row[[paste0(statistic, "_", metric)]] <- statistics[[statistic]]
    }
    row
  }))
}

summarize_records <- function(records) {
  if (!length(records)) stop("No result.rds files found.")
  contexts <- do.call(rbind, lapply(records, record_context))
  identity <- contexts[, c("scenario", "n", "rho", "repetition", "Method")]
  if (anyDuplicated(identity)) stop("Duplicate method/repetition results; restrict --input to one study.")
  signatures <- split(contexts$signature, contexts$scenario)
  if (any(vapply(signatures, function(x) length(unique(x)) != 1, logical(1)))) {
    stop("Different configurations/code/environments were mixed within a scenario.")
  }

  feature_metrics <- list()
  features <- list()
  predictions <- list()
  for (i in seq_along(records)) {
    record <- records[[i]]
    context <- record_context(record)
    metrics <- record$feature_metrics
    if (is.null(metrics)) {
      metrics <- data.frame(average_power = NA_real_, type1_error = NA_real_, rank_sum = NA_real_,
                            n_selected = NA_integer_, n_invalid_pvalues = NA_integer_,
                            inference_applicable = record$method != "XGBoost")
    }
    feature_metrics[[i]] <- cbind(context, metrics)

    feature <- record$features
    if (is.null(feature)) {
      feature <- data.frame(feature = record$candidate_features,
                            important = record$candidate_features %in% record$important_features,
                            score = NA_real_, pvalue = NA_real_, valid_inference = FALSE,
                            selected = NA, rejection = NA_real_, rank = NA_real_)
    }
    features[[i]] <- cbind(context[rep(1, nrow(feature)), , drop = FALSE], feature)

    prediction <- record$predictions
    if (is.null(prediction)) {
      prediction <- data.frame(feature_set = record$expected_feature_sets, cindex = NA_real_,
                               ibs = NA_real_, status = "selection_failed", error = record$selection_error,
                               warnings = "", elapsed_seconds = NA_real_, ibs_status = "not_run",
                               ibs_error = "", ibs_warnings = "", ibs_start = NA_real_, ibs_end = NA_real_)
    }
    predictions[[i]] <- cbind(context[rep(1, nrow(prediction)), , drop = FALSE], prediction)
  }
  feature_metrics <- do.call(rbind, feature_metrics)
  features <- do.call(rbind, features)
  predictions <- do.call(rbind, predictions)
  features$selection_frequency <- as.numeric(features$selected)
  keys <- c("scenario", "n", "rho", "Method")
  coverage <- do.call(rbind, lapply(unique(contexts$scenario), function(scenario) {
    scenario_records <- records[contexts$scenario == scenario]
    design <- scenario_records[[1]]$design
    grid <- expand.grid(n = design$training_sizes, rho = design$correlations, Method = design$methods,
                        stringsAsFactors = FALSE)
    grid$scenario <- scenario
    grid$n_expected <- max(vapply(scenario_records, function(record) record$design$repetitions, numeric(1)))
    grid$n_recorded <- vapply(seq_len(nrow(grid)), function(i) {
      sum(contexts$scenario == scenario & contexts$n == grid$n[i] &
            contexts$rho == grid$rho[i] & contexts$Method == grid$Method[i])
    }, integer(1))
    grid$n_selection_failed <- vapply(seq_len(nrow(grid)), function(i) {
      sum(contexts$scenario == scenario & contexts$n == grid$n[i] &
            contexts$rho == grid$rho[i] & contexts$Method == grid$Method[i] & contexts$selection_status == "failed")
    }, integer(1))
    grid$n_missing_repetitions <- grid$n_expected - grid$n_recorded
    grid[, c(keys, "n_expected", "n_recorded", "n_missing_repetitions", "n_selection_failed")]
  }))
  list(coverage = coverage,
       feature_recovery = group_summary(feature_metrics, keys, c("average_power", "type1_error", "rank_sum")),
       feature_rates = group_summary(features, c(keys, "feature", "important"),
                                     c("rejection", "selection_frequency", "rank")),
       prediction_summary = group_summary(predictions, c(keys, "feature_set"), c("cindex", "ibs")),
       feature_results = features, prediction_results = predictions, status = feature_metrics)
}

write_summary <- function(summary, output) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  atomic_save(summary, file.path(output, "summary.rds"))
  for (name in names(summary)) {
    utils::write.csv(summary[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
  }
}

if (sys.nframe() == 0L) {
  arguments <- summary_arguments(commandArgs(trailingOnly = TRUE))
  files <- list.files(arguments$input, pattern = "^result\\.rds$", recursive = TRUE, full.names = TRUE)
  records <- lapply(files, readRDS)
  write_summary(summarize_records(records), arguments$output)
  message("Summarized ", length(records), " method/repetition results.")
}
