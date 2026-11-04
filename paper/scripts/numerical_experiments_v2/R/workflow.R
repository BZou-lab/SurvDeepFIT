method_names <- c("SurvDeepFIT", "Cox", "AFT", "RSF", "XGBoost")

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package is unavailable: ", package,
         ". Configure the R environment before running; no packages are installed automatically.")
  }
}

read_config <- function(path) {
  require_package("yaml")
  config <- yaml::read_yaml(path, eval.expr = FALSE)
  config$training_sizes <- as.numeric(unlist(config$training_sizes))
  config$correlations <- as.numeric(unlist(config$correlations))
  config$methods <- unlist(config$methods)
  validate_config(config)
  config
}

validate_config <- function(config) {
  positive <- c(config$training_sizes, config$repetitions, config$block_size,
                config$baseline$gamma, config$baseline$lambda,
                config$censoring$shape, config$censoring$scale)
  if (!config$scenario %in% 1:2 || any(!is.finite(positive)) || any(positive <= 0) ||
      any(config$correlations < 0 | config$correlations >= 1) ||
      !length(config$methods) || any(!config$methods %in% method_names)) {
    stop("Invalid scenario, dimensions, correlations, distribution parameters, or methods.")
  }
  integers <- c(config$training_sizes, config$repetitions, config$block_size, config$base_seed)
  if (any(integers != floor(integers)) || config$block_size < 2 ||
      config$base_seed < 1 || config$base_seed > .Machine$integer.max ||
      any(config$training_sizes * config$test_fraction < 2) ||
      anyDuplicated(config$training_sizes) || anyDuplicated(config$correlations) ||
      anyDuplicated(config$methods)) {
    stop("Use unique settings, integer dimensions/seeds, block_size >= 2, and at least two test subjects.")
  }
  if (!isTRUE(config$selection$alpha > 0 && config$selection$alpha < 1) ||
      config$selection$xgboost_top < 1 || config$selection$xgboost_top > 5 * config$block_size ||
      !config$selection$rsf_table %in% c("var.sel", "var.sel.Z", "var.jk.sel.Z") ||
      config$selection$survdeepfit_pvalue != "importance_pval_x") {
    stop("Invalid selection settings. SurvDeepFIT inference is restricted to the documented legacy p-value.")
  }
  counts <- c(config$models$survdeepfit$folds, config$models$survdeepfit$permutations,
              config$models$survdeepfit$ensemble_size, config$models$survdeepfit$epochs,
              config$models$survdeepfit$batch_size, config$models$rsf$ntree,
              config$models$rsf$subsamples, config$models$xgboost$nrounds,
              config$models$xgboost$max_depth, config$evaluation$ibs_points,
              config$selection$xgboost_top)
  if (any(!is.finite(counts)) || any(counts < 1 | counts != floor(counts)) ||
      config$models$survdeepfit$folds < 2 || config$models$survdeepfit$permutations < 2 ||
      min(config$training_sizes) / config$models$survdeepfit$folds < 2 ||
      config$models$rsf$subsamples < 2 || config$models$aft$distribution != "weibull") {
    stop("Invalid model counts, cross-fitting sizes, or AFT distribution.")
  }
  quantiles <- unlist(config$evaluation$ibs_event_quantiles)
  if (length(quantiles) != 2 || any(quantiles <= 0 | quantiles >= 1) ||
      quantiles[1] >= quantiles[2] || config$evaluation$ibs_points < 2 ||
      config$evaluation$ibs_grid_source != "training" ||
      config$evaluation$censoring_estimate != "test_km" ||
      any(!unlist(config$evaluation$feature_sets) %in% c("full", "selected"))) {
    stop("Invalid IBS grid or prediction settings; use training quantiles and test-set censoring KM.")
  }
  invisible(config)
}

parse_options <- function(arguments) {
  options <- list()
  i <- 1L
  while (i <= length(arguments)) {
    key <- sub("^--", "", arguments[i])
    if (!startsWith(arguments[i], "--")) stop("Expected --option, received: ", arguments[i])
    if (key %in% c("help", "validate-only")) {
      options[[key]] <- TRUE
    } else {
      i <- i + 1L
      if (i > length(arguments) || startsWith(arguments[i], "--")) stop("Missing value for --", key)
      options[[key]] <- arguments[i]
    }
    i <- i + 1L
  }
  allowed <- c("help", "validate-only", "config", "output", "rep", "repetitions",
               "n", "rho", "seed", "methods", "threads")
  if (any(!names(options) %in% allowed)) stop("Unknown option: ", paste(setdiff(names(options), allowed), collapse = ", "))
  options
}

with_rng <- function(state, expression) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  assign(".Random.seed", state, envir = .GlobalEnv)
  force(expression)
}

task_rng <- function(config, setting_index, repetition) {
  # Separate streams for tasks; method/stage substreams stay inside each task stream.
  with_rng_seed <- function() {
    set.seed(config$base_seed, kind = "L'Ecuyer-CMRG")
    state <- .Random.seed
    n_settings <- length(config$training_sizes) * length(config$correlations)
    offset <- 2L * ((repetition - 1L) * n_settings + setting_index - 1L) + config$scenario - 1L
    if (offset > 0) for (i in seq_len(offset)) state <- parallel::nextRNGStream(state)
    state
  }
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- .Random.seed
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  with_rng_seed()
}

stage_rng <- function(task_state, index) {
  state <- task_state
  if (index > 0) for (i in seq_len(index)) state <- parallel::nextRNGSubStream(state)
  state
}

capture_stage <- function(expression) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  error <- ""
  value <- tryCatch(withCallingHandlers(force(expression), warning = function(warning) {
    warnings <<- c(warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(condition) {
    error <<- conditionMessage(condition)
    NULL
  })
  convergence_problem <- any(grepl(
    "converg|ran out of iterations|infinite|singular|numerically", warnings, ignore.case = TRUE))
  if (!nzchar(error) && convergence_problem) error <- "Fit rejected because of a convergence/numerical warning."
  ok <- !nzchar(error)
  list(ok = ok, status = if (!ok) "failed" else if (length(warnings)) "warning" else "ok",
       value = if (ok) value else NULL, error = error, warnings = unique(warnings),
       elapsed_seconds = proc.time()[["elapsed"]] - started)
}

atomic_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("checkpoint-", tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(temporary))
  saveRDS(object, temporary)
  # Windows cannot replace an existing file by rename. Keep a recovery copy
  # through that replacement; POSIX/Longleaf uses a single atomic rename.
  previous <- paste0(path, ".previous")
  if (.Platform$OS.type == "windows" && file.exists(path)) {
    if (file.exists(previous)) unlink(previous)
    if (!file.rename(path, previous)) stop("Could not preserve previous checkpoint: ", path)
  }
  if (!file.rename(temporary, path)) {
    if (file.exists(previous) && !file.exists(path)) file.rename(previous, path)
    stop("Could not replace checkpoint: ", path)
  }
  if (file.exists(previous)) unlink(previous)
  invisible(path)
}

object_hash <- function(object) {
  temporary <- tempfile(fileext = ".rds")
  on.exit(unlink(temporary))
  saveRDS(object, temporary, version = 2)
  unname(tools::md5sum(temporary))
}

run_provenance <- function(root, config) {
  repository <- normalizePath(file.path(root, "../../.."), winslash = "/")
  code <- c(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE),
            file.path(root, c("run_simulation.R", "summarize_results.R")),
            file.path(repository, "R/source_permfit_survival_20230610.r"))
  hashes <- tools::md5sum(code)
  names(hashes) <- substring(names(hashes), nchar(repository) + 2L)
  effective_config <- config
  effective_config$repetitions <- NULL
  git <- function(arguments) {
    tryCatch(system2("git", c("-C", shQuote(repository), arguments), stdout = TRUE, stderr = FALSE),
             error = function(e) "unavailable")
  }
  versions <- software_versions()$packages
  list(signature = object_hash(list(config = effective_config, code = hashes,
                                    R_version = R.version.string, packages = versions)),
       config = config, code_md5 = hashes, git_commit = git(c("rev-parse", "HEAD")),
       git_status = git(c("status", "--short")), created_utc = format(Sys.time(), tz = "UTC"),
       R_version = R.version.string)
}

software_versions <- function() {
  packages <- c("survival", "MASS", "yaml", "deepTL", "Hmisc", "timeROC", "randomForestSRC", "xgboost")
  versions <- vapply(packages, function(package) {
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
  }, character(1))
  list(packages = versions, session = utils::sessionInfo())
}

read_checkpoint <- function(path, signature) {
  if (!file.exists(path) && file.exists(paste0(path, ".previous"))) path <- paste0(path, ".previous")
  if (!file.exists(path)) return(NULL)
  checkpoint <- readRDS(path)
  if (!identical(checkpoint$signature, signature)) {
    stop("Configuration/code changed for existing output: ", path,
         ". Use a new output directory; existing results were preserved.")
  }
  checkpoint
}

run_method <- function(method, data, config, task_state, directory, context,
                       signature, core_file, threads) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  selection_path <- file.path(directory, "selection.rds")
  selection <- read_checkpoint(selection_path, signature)
  method_offset <- 10L * match(method, method_names)
  if (is.null(selection) || !selection$stage$ok) {
    stage <- capture_stage(fit_feature_selection(
      method, data, config, stage_rng(task_state, method_offset),
      stage_rng(task_state, method_offset + 1L), core_file, threads))
    selection <- list(signature = signature, stage = stage)
    atomic_save(selection, selection_path)
  }

  predictions <- list()
  if (selection$stage$ok) {
    selected <- selection$stage$value$features$feature[selection$stage$value$features$selected %in% TRUE]
    for (feature_set in unlist(config$evaluation$feature_sets)) {
      path <- file.path(directory, paste0("prediction_", feature_set, ".rds"))
      checkpoint <- read_checkpoint(path, signature)
      if (is.null(checkpoint) || !checkpoint$stage$ok ||
          (!is.null(checkpoint$stage$value) && checkpoint$stage$value$ibs_status == "failed")) {
        features <- if (feature_set == "full") data$feature_names else selected
        stage <- capture_stage({
          if (!length(features)) stop("No selected features; no prediction fallback is specified in the manuscript.")
          model <- if (feature_set == "full") selection$stage$value$model else NULL
          if (is.null(model)) {
            model <- with_rng(stage_rng(task_state, method_offset + if (feature_set == "full") 2L else 3L),
                              fit_model(method, data$train, features, config, threads))
          }
          times <- evaluation_times(data$train, config)
          evaluate_prediction(method, model, data$train, data$test, features, times)
        })
        checkpoint <- list(signature = signature, stage = stage)
        atomic_save(checkpoint, path)
      }
      stage <- checkpoint$stage
      row <- if (stage$ok) stage$value else data.frame(
        cindex = NA_real_, ibs = NA_real_, ibs_status = "not_run", ibs_error = "",
        ibs_warnings = "", ibs_start = NA_real_, ibs_end = NA_real_,
        fit_initialization = NA_character_, fit_iterations = NA_integer_, fit_max_abs_score = NA_real_)
      row$feature_set <- feature_set
      row$status <- stage$status
      row$error <- stage$error
      row$warnings <- paste(stage$warnings, collapse = " | ")
      row$elapsed_seconds <- stage$elapsed_seconds
      predictions[[feature_set]] <- row
    }
  }

  stage <- selection$stage
  record <- list(schema_version = 1L, signature = signature, context = context,
                 design = config[c("scenario", "training_sizes", "correlations", "repetitions", "methods")],
                 method = method, selection_status = stage$status, selection_error = stage$error,
                 selection_warnings = stage$warnings, selection_elapsed_seconds = stage$elapsed_seconds,
                 features = if (stage$ok) stage$value$features else NULL,
                 feature_metrics = if (stage$ok) stage$value$feature_metrics else NULL,
                 candidate_features = data$feature_names, important_features = data$important_features,
                 expected_feature_sets = unlist(config$evaluation$feature_sets),
                 predictions = if (length(predictions)) do.call(rbind, predictions) else NULL,
                 rng = list(model = stage_rng(task_state, method_offset),
                            inference = stage_rng(task_state, method_offset + 1L)),
                 software = software_versions())
  atomic_save(record, file.path(directory, "result.rds"))
  invisible(record)
}
