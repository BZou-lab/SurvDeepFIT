# Deterministic checks; no forest fitting or study outputs.
file_argument <- grep("^--file=", commandArgs(), value = TRUE)[1]
root <- dirname(normalizePath(sub("^--file=", "", file_argument)))
runner <- file.path(root, "run_rsf_diagnostic.R")
scope <- new.env()
scope$commandArgs <- function(...) paste0("--file=", runner)
sys.source(runner, envir = scope)
config <- read_config(file.path(root, "configs/scenario1.yaml"))
stopifnot(config$models$rsf$importance == "permute",
          config$selection$rsf_table == "var.sel")
state <- task_rng(config, 1, 1)
data <- generate_data(config, 1000, 0, stage_rng(state, 0))
null <- scope$global_null_data(data, config)
stopifnot(identical(data$train[data$feature_names], null$train[null$feature_names]),
          identical(data$censor_time, null$censor_time),
          length(null$important_features) == 0, all(null$log_risk == 0))
hazard <- function(time) config$baseline$lambda/config$baseline$gamma * expm1(config$baseline$gamma*time)
stopifnot(isTRUE(all.equal(hazard(data$true_time)*exp(data$log_risk), hazard(null$true_time))))
features <- data$feature_names
truth <- features %in% data$important_features
vimp <- setNames(as.numeric(truth), features)
tab <- data.frame(lower=vimp-.1645, mean=vimp, upper=vimp+.1645,
                  pvalue=ifelse(truth, .001, .5), row.names=features)
inference <- list(vmp=vimp, vmpS=matrix(rep(vimp,100),nrow=50,dimnames=list(features,NULL)),
                  se.Z=setNames(rep(.1,50),features), se.jk.Z=setNames(rep(.1,50),features),
                  var.sel=tab, var.sel.Z=tab, var.jk.sel.Z=tab)
a <- scope$diagnostic_metrics(data, inference, "var.sel", config, "test", 1)
stopifnot(a$metrics$power==1, a$metrics$type1_error==0,
          a$metrics$n_selected==5, a$metrics$n_selected_null==0,
          a$metrics$ranking_auc==1, a$metrics$vimp_rank_sum==15,
          a$metrics$n_invalid==0)
b <- scope$diagnostic_metrics(null, inference, "var.sel", config, "test", 1)
stopifnot(is.na(b$metrics$power), is.na(b$metrics$ranking_auc),
          b$metrics$type1_error==.1, b$metrics$n_selected_null==5)
inference$var.sel$pvalue[1] <- NA_real_
c <- scope$diagnostic_metrics(data, inference, "var.sel", config, "test", 1)
stopifnot(is.na(c$metrics$power), is.na(c$metrics$n_selected),
          c$metrics$type1_error==0, c$metrics$n_invalid==1)
local({
  directory <- tempfile("rsf-metadata-check-", tmpdir = tempdir())
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  path <- file.path(directory, "provenance.rds")
  results <- parallel::mclapply(seq_len(8), function(i) {
    scope$initialize_diagnostic_metadata(path, list(signature = "same", writer = i))
  }, mc.cores = if (.Platform$OS.type == "unix") 4L else 1L)
  stopifnot(all(vapply(results, function(x) identical(x, results[[1]]), logical(1))))
  hash <- tools::md5sum(path)
  scope$initialize_diagnostic_metadata(path, list(signature = "same", writer = 99))
  stopifnot(identical(hash, tools::md5sum(path)))
  bad <- try(scope$initialize_diagnostic_metadata(path, list(signature = "different")), silent = TRUE)
  stopifnot(inherits(bad, "try-error"), identical(hash, tools::md5sum(path)),
            !dir.exists(paste0(path, ".lock")))
})
cat("PASS: permutation/nonparametric study settings, global-null event uniforms/covariates/censoring, selection counts,\n",
    "ranking/ties, empty signal set, missing-p-value denominators, and concurrent\n",
    "immutable provenance initialization. No forests fitted.\n")
