# Load dependencies and inspect interfaces only; never generate data or fit models.
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                             winslash = "/"))
source(file.path(root, "R/workflow.R"))
source(file.path(root, "R/fit_models.R"))

main <- function(arguments) {
  options <- parse_options(arguments)
  if (isTRUE(options$help)) {
    cat("Rscript check_environment.R --config FILE [--methods Cox,AFT]\n")
    return(invisible(NULL))
  }
  if (is.null(options$config)) stop("Supply --config configs/scenario1.yaml or scenario2.yaml.")
  config <- read_config(options$config)
  selected <- if (is.null(options$methods)) config$methods else strsplit(options$methods, ",", fixed = TRUE)[[1]]
  if (!length(selected) || any(!selected %in% config$methods) || anyDuplicated(selected)) {
    stop("--methods must select unique configured methods.")
  }
  packages <- c("yaml", "MASS", "survival")
  if ("SurvDeepFIT" %in% selected) packages <- c(packages, "deepTL", "magrittr", "Hmisc", "timeROC")
  if ("RSF" %in% selected) packages <- c(packages, "randomForestSRC")
  if ("XGBoost" %in% selected) packages <- c(packages, "xgboost")

  checks <- lapply(unique(packages), function(package) {
    warnings <- character()
    error <- ""
    loaded <- tryCatch(withCallingHandlers({
      suppressPackageStartupMessages(loadNamespace(package))
      TRUE
    }, warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = function(condition) {
      error <<- conditionMessage(condition)
      FALSE
    })
    data.frame(package = package, available = loaded,
               version = if (loaded) as.character(utils::packageVersion(package)) else NA_character_,
               error = error, warnings = paste(unique(warnings), collapse = " | "))
  })
  checks <- do.call(rbind, checks)
  cat(R.version.string, "\nLibrary paths:\n", paste(.libPaths(), collapse = "\n"), "\n", sep = "")
  print(checks[, c("package", "available", "version")], row.names = FALSE)
  for (i in seq_len(nrow(checks))) {
    if (nzchar(checks$error[i])) cat("\n", checks$package[i], " load error: ", checks$error[i], "\n", sep = "")
    if (nzchar(checks$warnings[i])) cat("\n", checks$package[i], " warnings: ", checks$warnings[i], "\n", sep = "")
  }
  if (!all(checks$available)) {
    stop("Dependency check failed. Resolve the listed load errors in this R environment before fitting models.")
  }

  if ("SurvDeepFIT" %in% selected) {
    core <- load_legacy_inference(file.path(root, "../../../R/source_permfit_survival_20230610.r"))
    stopifnot(is.function(core$permfit_survival), is.function(deepTL::ensemble_dnnet))
    if (!methods::isClass("PermFIT")) stop("The loaded deepTL build does not register the PermFIT class.")
  }
  if ("RSF" %in% selected && !"ntree" %in% names(formals(randomForestSRC::rfsrc))) {
    stop("Unexpected randomForestSRC interface: rfsrc has no ntree argument.")
  }
  if ("XGBoost" %in% selected) {
    for (name in c("xgb.DMatrix", "xgb.train", "xgb.importance")) getExportedValue("xgboost", name)
  }
  cat("PASS: configuration, package loading, and required interfaces for", paste(selected, collapse = ", "), "\n")
  cat("No data generated or model fits performed. Package-specific training still requires the pilot.\n")
}

main(commandArgs(trailingOnly = TRUE))
