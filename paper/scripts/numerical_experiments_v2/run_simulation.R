script_argument <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", script_argument[1]), winslash = "/"))
for (file in c("workflow.R", "generate_data.R", "fit_models.R", "evaluate_features.R", "evaluate_prediction.R")) {
  source(file.path(root, "R", file))
}

main <- function(arguments) {
  options <- parse_options(arguments)
  if (isTRUE(options$help)) {
    cat("Rscript run_simulation.R --config FILE [--rep ID | --repetitions N]\n",
        "  [--n TRAINING_N] [--rho RHO] [--seed BASE_SEED] [--methods Cox,AFT]\n",
        "  [--threads N] [--output DIRECTORY] [--validate-only]\n")
    return(invisible(NULL))
  }
  if (is.null(options$config)) stop("Supply --config configs/scenario1.yaml or scenario2.yaml.")
  config <- read_config(options$config)
  if (!is.null(options$seed)) config$base_seed <- as.numeric(options$seed)
  if (!is.null(options$repetitions)) config$repetitions <- as.numeric(options$repetitions)
  validate_config(config)

  repetitions <- if (is.null(options$rep)) seq_len(config$repetitions) else as.numeric(options$rep)
  if (any(!is.finite(repetitions)) || any(repetitions < 1 | repetitions != floor(repetitions)) ||
      any(repetitions > config$repetitions)) stop("Repetition IDs must be integers within configured repetitions.")
  methods <- if (is.null(options$methods)) config$methods else strsplit(options$methods, ",", fixed = TRUE)[[1]]
  if (any(!methods %in% config$methods) || anyDuplicated(methods)) stop("--methods must select unique configured methods.")
  threads <- if (is.null(options$threads)) 1L else as.numeric(options$threads)
  if (!is.finite(threads) || threads < 1 || threads != floor(threads)) stop("--threads must be a positive integer.")
  options(rf.cores = threads)
  Sys.setenv(OMP_NUM_THREADS = threads, OPENBLAS_NUM_THREADS = threads, MKL_NUM_THREADS = threads)

  settings <- expand.grid(n_train = config$training_sizes, rho = config$correlations)
  settings$setting_index <- seq_len(nrow(settings))
  if (!is.null(options$n)) settings <- settings[settings$n_train == as.numeric(options$n), , drop = FALSE]
  if (!is.null(options$rho)) settings <- settings[settings$rho == as.numeric(options$rho), , drop = FALSE]
  if (!nrow(settings)) stop("--n/--rho do not match the configuration grid.")
  if (isTRUE(options$`validate-only`)) {
    cat("Configuration valid;", nrow(settings), "settings;", length(repetitions), "repetitions; methods:", methods, "\n")
    return(invisible(NULL))
  }

  require_package("survival")
  provenance <- run_provenance(root, config)
  provenance$signature <- object_hash(list(provenance$signature, threads = threads))
  output <- if (is.null(options$output)) file.path(root, "results") else options$output
  core_file <- file.path(root, "../../../R/source_permfit_survival_20230610.r")
  for (repetition in repetitions) for (i in seq_len(nrow(settings))) {
    setting <- settings[i, ]
    directory <- file.path(output, paste0("scenario", config$scenario), paste0("n", setting$n_train),
                           paste0("rho", setting$rho), sprintf("rep%04d", repetition))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    lock <- file.path(directory, ".lock")
    if (!dir.create(lock, showWarnings = FALSE)) stop("Task already locked: ", directory)
    tryCatch({
      data_path <- file.path(directory, "data.rds")
      checkpoint <- read_checkpoint(data_path, provenance$signature)
      task_state <- task_rng(config, setting$setting_index, repetition)
      if (is.null(checkpoint)) {
        data <- generate_data(config, setting$n_train, setting$rho, stage_rng(task_state, 0L))
        checkpoint <- list(signature = provenance$signature, value = data)
        atomic_save(checkpoint, data_path)
      }
      data <- checkpoint$value
      context <- list(scenario = config$scenario, n_train = nrow(data$train), n_test = nrow(data$test),
                      rho = setting$rho, repetition = repetition, base_seed = config$base_seed,
                      task_rng = task_state, censoring_train = data$censoring_train,
                      censoring_test = data$censoring_test, threads = threads)
      atomic_save(list(context = context, provenance = provenance), file.path(directory, "run.rds"))
      for (method in methods) {
        message("Scenario ", config$scenario, "; n=", setting$n_train, "; rho=", setting$rho,
                "; repetition=", repetition, "; method=", method)
        run_method(method, data, config, task_state, file.path(directory, method), context,
                   provenance$signature, core_file, threads)
      }
    }, finally = unlink(lock))
  }
}

main(commandArgs(trailingOnly = TRUE))
