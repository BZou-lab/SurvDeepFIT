script_argument <- grep("^--file=", commandArgs(), value = TRUE)
source_file <- if (sys.nframe() > 0L) sys.frame(1)$ofile else NULL
diagnostic_root <- dirname(normalizePath(if (!is.null(source_file)) source_file else
  sub("^--file=", "", script_argument[1])))
for (file in c("workflow.R", "generate_data.R", "fit_models.R", "evaluate_features.R")) {
  source(file.path(diagnostic_root, "R", file))
}

diagnostic_options <- function(arguments) {
  if (!length(arguments) || length(arguments) %% 2L ||
      any(!startsWith(arguments[seq(1, length(arguments), 2)], "--"))) {
    stop("Use --config FILE [--output DIR] [--rep ID] [--design signal|global_null|both] [--threads N].")
  }
  keys <- sub("^--", "", arguments[seq(1, length(arguments), 2)])
  if (anyDuplicated(keys) || any(!keys %in% c("config", "output", "rep", "design", "threads"))) {
    stop("Duplicate or unknown diagnostic option.")
  }
  stats::setNames(as.list(arguments[seq(2, length(arguments), 2)]), keys)
}

initialize_diagnostic_metadata <- function(path, metadata) {
  # A shared run manifest is immutable. Serialize creation/checking across nodes.
  # Rewriting an existing manifest from every array task caused startup read races.
  lock <- paste0(path, ".lock")
  attempts <- 0L
  while (!dir.create(lock, showWarnings = FALSE)) {
    attempts <- attempts + 1L
    if (attempts >= 200L) stop("Provenance is locked: ", lock,
      ". Check running initializers before removing a stale lock.")
    Sys.sleep(0.1)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(path)) {
    existing <- readRDS(path)
    if (!identical(existing$signature, metadata$signature)) {
      stop("Output configuration changed; choose a new directory.")
    }
    return(invisible(existing))
  }
  atomic_save(metadata, path)
  invisible(metadata)
}

global_null_data <- function(data, config) {
  # Recover the exponential event draw used by the unchanged Gompertz generator.
  # Reuse covariates, censoring, split, and event uniforms, setting only log-risk to zero.
  gamma <- config$baseline$gamma
  lambda <- config$baseline$lambda
  exponential_draw <- (lambda / gamma) * exp(data$log_risk) * expm1(gamma * data$true_time)
  true_time <- log1p((gamma / lambda) * exponential_draw) / gamma
  nt <- nrow(data$train)
  observed_time <- pmin(true_time, data$censor_time)
  status <- as.integer(true_time <= data$censor_time)
  data$train$time <- observed_time[seq_len(nt)]
  data$train$status <- status[seq_len(nt)]
  data$test$time <- observed_time[nt + seq_len(nrow(data$test))]
  data$test$status <- status[nt + seq_len(nrow(data$test))]
  data$true_time <- true_time
  data$log_risk[] <- 0
  data$important_features <- character()
  data$censoring_train <- mean(data$train$status == 0)
  data$censoring_test <- mean(data$test$status == 0)
  data
}

diagnostic_metrics <- function(data, inference, table_name, config, arm, runtime) {
  features <- data$feature_names
  table <- inference[[table_name]]
  stopifnot(!anyDuplicated(rownames(table)), setequal(rownames(table), features))
  p <- stats::setNames(table[features, "pvalue"], features)
  vimp <- inference$vmp[features]
  result <- evaluate_features("RSF", features, data$important_features, p, vimp, config)
  truth <- result$features$important
  # A finite p-value alone does not establish that resampling inference was estimable.
  available <- rowSums(is.finite(inference$vmpS[features, , drop = FALSE]))
  se <- if (table_name == "var.jk.sel.Z") inference$se.jk.Z[features] else inference$se.Z[features]
  estimable <- available == config$models$rsf$subsamples & is.finite(vimp) &
    is.finite(se) & se > 0
  result$features$valid_inference <- result$features$valid_inference & estimable
  result$features$selected[!result$features$valid_inference] <- NA
  selected <- result$features$selected
  ranks <- rank(-vimp, ties.method = "average", na.last = "keep")
  n1 <- sum(truth)
  n0 <- sum(!truth)
  ranking_auc <- if (n1 && n0 && all(is.finite(ranks))) {
    1 - (sum(ranks[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  } else NA_real_
  complete <- all(result$features$valid_inference)
  metrics <- data.frame(arm = arm,
    power = if (n1 && all(result$features$valid_inference[truth])) mean(selected[truth]) else NA_real_,
    type1_error = if (n0 && all(result$features$valid_inference[!truth])) mean(selected[!truth]) else NA_real_,
    n_selected = if (complete) sum(selected) else NA_integer_,
    n_selected_null = if (all(result$features$valid_inference[!truth])) sum(selected[!truth]) else NA_integer_,
    ranking_auc = ranking_auc,
    vimp_rank_sum = if (n1 && all(is.finite(ranks))) sum(ranks[truth]) else NA_real_,
    pvalue_rank_sum = if (n1) result$metrics$rank_sum else NA_real_,
    runtime_seconds = runtime, n_invalid = sum(!result$features$valid_inference),
    min_available_subsamples = min(available), n_negative_vimp = sum(vimp < 0, na.rm = TRUE),
    n_zero_pvalues = sum(p == 0, na.rm = TRUE),
    n_pvalue_ci_disagreements = sum((p < config$selection$alpha) != (table[features, "lower"] > 0), na.rm = TRUE))
  result$features$arm <- arm
  result$features$full_data_vimp_rank <- ranks
  result$features$available_subsamples <- available
  result$features$se_sub <- inference$se.Z[features]
  result$features$se_jackknife <- inference$se.jk.Z[features]
  result$features$ci_lower <- table[features, "lower"]
  result$features$ci_upper <- table[features, "upper"]
  result$features$rejection <- as.numeric(selected)
  list(metrics = metrics, features = result$features)
}

run_diagnostic_fit <- function(data, config, task_state, importance = NULL) {
  # NULL explicitly restores the original grow/subsample path for the legacy arm,
  # even when the study configuration requests permutation importance.
  config$models$rsf$importance <- importance
  started <- proc.time()[["elapsed"]]
  fitted <- with_rng(stage_rng(task_state, 40L), fit_model("RSF", data$train, data$feature_names, config))
  # Identical calls/RNG streams to legacy inference, retaining raw estimates for auditing.
  sampled <- with_rng(stage_rng(task_state, 41L),
    randomForestSRC::subsample(fitted, B = config$models$rsf$subsamples))
  inference <- randomForestSRC::extract.subsample(sampled, alpha = config$selection$alpha, raw = TRUE)
  elapsed <- proc.time()[["elapsed"]] - started
  list(inference = inference, elapsed_seconds = elapsed,
    settings = list(importance = sampled$rf$forest$importance, block.size = sampled$rf$block.size,
      subratio = sampled$subratio, ntree = fitted$ntree, mtry = fitted$mtry,
      nodesize = fitted$nodesize, splitrule = fitted$splitrule, use.uno = !is.null(fitted$forest$uno.weights)))
}

summarize_diagnostic <- function(output, signature) {
  paths <- list.files(output, "^diagnostic\\.rds$", full.names = TRUE, recursive = TRUE)
  records <- lapply(paths, readRDS)
  if (any(vapply(records, function(x) !identical(x$signature, signature), logical(1)))) {
    stop("Mixed diagnostic signatures; choose a separate output directory.")
  }
  metrics <- do.call(rbind, lapply(records, `[[`, "metrics"))
  features <- do.call(rbind, lapply(records, `[[`, "features"))
  groups <- split(seq_len(nrow(metrics)), interaction(metrics$design, metrics$arm, drop = TRUE))
  measured <- c("power", "type1_error", "n_selected", "n_selected_null", "ranking_auc", "vimp_rank_sum",
                "pvalue_rank_sum", "runtime_seconds", "n_invalid", "min_available_subsamples")
  summary <- do.call(rbind, lapply(groups, function(ii) {
    z <- metrics[ii, ]
    row <- z[1, c("design", "arm")]
    row$n_repetitions <- nrow(z)
    for (metric in measured) {
      values <- z[[metric]][is.finite(z[[metric]])]
      row[[metric]] <- if (length(values)) mean(values) else NA_real_
      row[[paste0(metric, "_mcse")]] <- if (length(values) > 1) sd(values) / sqrt(length(values)) else NA_real_
      row[[paste0(metric, "_n_valid")]] <- length(values)
    }
    row
  }))
  utils::write.csv(metrics, file.path(output, "repetition_metrics.csv"), row.names = FALSE)
  utils::write.csv(features, file.path(output, "feature_results.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
  print(summary[, c("design", "arm", "n_repetitions", "power", "type1_error", "n_selected",
                    "n_selected_null", "ranking_auc", "runtime_seconds")], row.names = FALSE)
  invisible(summary)
}

main_diagnostic <- function(arguments) {
  opt <- diagnostic_options(arguments)
  if (is.null(opt$config)) stop("Supply --config configs/rsf_diagnostic.yaml.")
  require_package("yaml")
  diagnostic <- yaml::read_yaml(opt$config)
  config <- read_config(file.path(dirname(normalizePath(opt$config)), diagnostic$study_config))
  stopifnot(config$scenario == 1, diagnostic$alternative_table == "var.jk.sel.Z",
            diagnostic$alpha == 0.05)
  config$selection$alpha <- diagnostic$alpha
  threads <- if (is.null(opt$threads)) diagnostic$threads else as.numeric(opt$threads)
  if (!is.finite(threads) || threads < 1 || threads != floor(threads)) stop("Invalid threads.")
  options(rf.cores = threads)
  Sys.setenv(OMP_NUM_THREADS = threads, OPENBLAS_NUM_THREADS = threads, MKL_NUM_THREADS = threads)
  require_package("randomForestSRC")
  require_package("survival")
  repetitions <- if (is.null(opt$rep)) seq_len(diagnostic$repetitions) else as.numeric(opt$rep)
  if (any(!is.finite(repetitions) | repetitions < 1 | repetitions > diagnostic$repetitions |
          repetitions != floor(repetitions))) stop("Invalid repetition ID.")
  designs <- unlist(diagnostic$designs)
  if (!is.null(opt$design) && opt$design != "both") designs <- opt$design
  if (any(!designs %in% c("signal", "global_null"))) stop("Unknown design.")
  settings <- expand.grid(n_train = config$training_sizes, rho = config$correlations)
  setting_index <- which(settings$n_train == diagnostic$n_train & settings$rho == diagnostic$rho)
  stopifnot(length(setting_index) == 1L)
  output <- if (is.null(opt$output)) file.path(diagnostic_root, "results/rsf-diagnostic") else opt$output
  provenance <- run_provenance(diagnostic_root, config)
  extra <- tools::md5sum(c(file.path(diagnostic_root, "run_rsf_diagnostic.R"), normalizePath(opt$config)))
  signature <- object_hash(list(provenance$signature, extra, diagnostic, threads))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  metadata <- list(signature = signature, diagnostic = diagnostic, provenance = provenance,
                   diagnostic_md5 = extra, software = software_versions(), threads = threads)
  metadata_path <- file.path(output, "provenance.rds")
  initialize_diagnostic_metadata(metadata_path, metadata)
  for (repetition in repetitions) {
    state <- task_rng(config, setting_index, repetition)
    signal <- generate_data(config, diagnostic$n_train, diagnostic$rho, stage_rng(state, 0L))
    for (design in designs) {
      directory <- file.path(output, design, sprintf("rep%04d", repetition))
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      record_path <- file.path(directory, "diagnostic.rds")
      if (!is.null(read_checkpoint(record_path, signature))) next
      lock <- file.path(directory, ".lock")
      if (!dir.create(lock, showWarnings = FALSE)) stop("Diagnostic task already locked: ", directory)
      tryCatch({
        data <- if (design == "signal") signal else global_null_data(signal, config)
        atomic_save(list(signature = signature, data = data, task_rng = state), file.path(directory, "data.rds"))
        fits <- list()
        for (importance in c("legacy", "permute")) {
          path <- file.path(directory, paste0(importance, ".rds"))
          cached <- read_checkpoint(path, signature)
          if (is.null(cached)) {
            message("RSF diagnostic: ", design, "; repetition ", repetition, "; ", importance)
            stage <- capture_stage(run_diagnostic_fit(data, config, state,
              if (importance == "legacy") NULL else importance))
            atomic_save(list(signature = signature, stage = stage), path)
            if (!stage$ok) stop(stage$error)
            cached <- list(stage = stage)
          }
          if (!cached$stage$ok) stop("Cached fit failed: ", cached$stage$error)
          fits[[importance]] <- cached$stage$value
          message("  elapsed ", round(cached$stage$value$elapsed_seconds, 1), " seconds")
        }
        arms <- list(original_anti_np = c("legacy", "var.sel"),
                     permutation_np = c("permute", "var.sel"),
                     permutation_jackknife = c("permute", "var.jk.sel.Z"),
                     permutation_normal = c("permute", "var.sel.Z"))
        results <- lapply(names(arms), function(arm) {
          spec <- arms[[arm]]
          fit <- fits[[spec[1]]]
          diagnostic_metrics(data, fit$inference, spec[2], config, arm, fit$elapsed_seconds)
        })
        metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
        features <- do.call(rbind, lapply(results, `[[`, "features"))
        for (name in c("metrics", "features")) {
          z <- get(name)
          z$design <- design
          z$repetition <- repetition
          assign(name, z)
        }
        atomic_save(list(signature = signature, metrics = metrics, features = features), record_path)
        summarize_diagnostic(output, signature)
      }, finally = unlink(lock))
    }
  }
  # Refresh summaries even when all requested repetitions were already checkpointed.
  summarize_diagnostic(output, signature)
}

if (sys.nframe() == 0L) main_diagnostic(commandArgs(trailingOnly = TRUE))
