# Small deterministic checks only: no forests, boosting, DNNs, or full simulation study.
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                             winslash = "/"))
for (file in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(file)
source(file.path(root, "summarize_results.R"))
require_package("survival")

for (file in list.files(root, "\\.R$", recursive = TRUE, full.names = TRUE)) parse(file)
config <- read_config(file.path(root, "configs/scenario1.yaml"))
pilot_config <- config
second <- read_config(file.path(root, "configs/scenario2.yaml"))
stopifnot(identical(config[names(config) != "scenario"], second[names(second) != "scenario"]))

config$training_sizes <- 120
config$block_size <- 2
config$evaluation$ibs_event_quantiles <- c(0.1, 0.6)
config$evaluation$ibs_points <- 20
validate_config(config)
seed <- task_rng(config, 1, 1)
stopifnot(!identical(seed, task_rng(config, 1, 2)),
          !identical(seed, task_rng(config, 2, 1)),
          !identical(stage_rng(seed, 1), stage_rng(seed, 2)))

for (scenario in 1:2) for (rho in config$correlations) {
  config$scenario <- scenario
  data <- generate_data(config, 120, rho, seed)
  stopifnot(identical(data, generate_data(config, 120, rho, seed)),
            nrow(data$train) == 120, nrow(data$test) == 24,
            ncol(data$train) == 12, length(data$important_features) == scenario + 4,
            all(data$train$status %in% 0:1), all(data$train$time > 0),
            identical(names(data$train), names(data$test)),
            identical(c(data$train$time, data$test$time), pmin(data$true_time, data$censor_time)))
  reconstructed <- if (scenario == 1) {
    as.numeric(as.matrix(rbind(data$train, data$test)[, data$important_features]) %*% c(0.2, 0.4, 0.6, 0.8, 1))
  } else {
    x <- rbind(data$train, data$test)
    x$Z1 + 2 * x$X1 * x$Z2 + 0.5 * x$X3^2 + x$X5 * x$X7
  }
  # Current pilot design halves the original coefficients and omits n-dependent scaling.
  stopifnot(isTRUE(all.equal(data$log_risk, reconstructed / 2)))
}

config$scenario <- 1
data <- generate_data(config, 120, 0, seed)
fixture_p <- ifelse(data$feature_names %in% data$important_features, 0.01, 0.5)
features <- evaluate_features("Cox", data$feature_names, data$important_features,
                              fixture_p, rep(0, 10), config)
stopifnot(features$metrics$average_power == 1, features$metrics$type1_error == 0,
          features$metrics$rank_sum == 15)
fixture_p[which(!features$features$important)[1]] <- NA
invalid <- evaluate_features("Cox", data$feature_names, data$important_features,
                             fixture_p, rep(0, 10), config)
stopifnot(is.na(invalid$metrics$type1_error), is.na(invalid$metrics$rank_sum))
xgb <- evaluate_features("XGBoost", data$feature_names, data$important_features,
                         rep(NA_real_, 10), seq_len(10), config)
stopifnot(is.na(xgb$metrics$average_power), is.na(xgb$metrics$type1_error),
          xgb$metrics$n_selected == 10)

# Isolate adapter checks from Weibull AFT misspecification under Gompertz times.
# Fixed quantile permutations give a small, deterministic Weibull fixture.
model_data <- data
fixture <- rbind(data$train, data$test)
index <- seq_len(nrow(fixture))
event <- stats::qweibull(((37 * index) %% nrow(fixture) + 0.5) / nrow(fixture),
                         shape = 2, scale = exp(0.2 * fixture$X1 - 0.1 * fixture$X3))
censor <- stats::qweibull(((53 * index) %% nrow(fixture) + 0.5) / nrow(fixture),
                          shape = 1.5, scale = 2.5)
fixture$time <- pmin(event, censor)
fixture$status <- as.integer(event <= censor)
model_data$train <- fixture[seq_len(nrow(data$train)), , drop = FALSE]
model_data$test <- fixture[nrow(data$train) + seq_len(nrow(data$test)), , drop = FALSE]

# Two cheap native fits verify extraction, risk direction, and probability dimensions.
for (method in c("Cox", "AFT")) {
  stage <- capture_stage(fit_feature_selection(method, model_data, config, seed, stage_rng(seed, 1), "unused"))
  if (!stage$ok) stop(method, " smoke fit failed: ", stage$error, " ", paste(stage$warnings, collapse = "; "))
  stopifnot(identical(stage$value$features$feature, data$feature_names),
            !any(c("(Intercept)", "Log(scale)") %in% stage$value$features$feature))
  prediction <- evaluate_prediction(method, stage$value$model, model_data$train, model_data$test,
                                    data$feature_names, evaluation_times(model_data$train, config))
  stopifnot(is.finite(prediction$cindex), is.finite(prediction$ibs),
            prediction$cindex >= 0, prediction$cindex <= 1,
            prediction$ibs >= 0, prediction$ibs <= 1)
}

# Exercise the pilot-sized AFT fit using the reported setting and RNG stream.
# This is one inexpensive native fit; no DNN/forest/boosting training is involved.
pilot_data <- generate_data(pilot_config, 1000, 0, task_rng(pilot_config, 1, 1))
aft <- capture_stage(fit_model("AFT", pilot_data$train, pilot_data$feature_names, pilot_config))
if (!aft$ok) stop("Pilot-sized AFT regression check failed: ", aft$error,
                  " ", paste(aft$warnings, collapse = "; "))
solver <- attr(aft$value, "aft_solver")
stopifnot(solver$initialization == "intercept_only", solver$max_abs_score < 0.001,
          all(solver$initial[pilot_data$feature_names] == 0),
          aft$value$loglik[2] >= aft$value$loglik[1])
formula <- stats::as.formula("Surv(time, status) ~ .", env = asNamespace("survival"))
aft_longer <- fit_weibull_aft(formula, pilot_data$train, survival::survreg.control(maxiter = 100))
stopifnot(isTRUE(all.equal(stats::coef(aft$value), stats::coef(aft_longer), tolerance = 1e-7)))
rejected <- capture_stage(fit_weibull_aft(formula, pilot_data$train,
                                        survival::survreg.control(maxiter = 1)))
stopifnot(!rejected$ok)
cat("PASS: pilot-sized Weibull AFT converges from intercept-only initialization; failed fits remain rejected.\n")

# Namespace-only deepTL loading must support both legacy inference and v2 predictions.
if (requireNamespace("deepTL", quietly = TRUE) && requireNamespace("Hmisc", quietly = TRUE) &&
    requireNamespace("timeROC", quietly = TRUE) && requireNamespace("magrittr", quietly = TRUE)) {
  core <- load_legacy_inference(file.path(root, "../../../R/source_permfit_survival_20230610.r"))
  check_survdeepfit_dispatch(core)
  cat("PASS: real dnnet/dnnetEnsemble prediction fixtures without training.\n")
  # Check formula scope, numerical results, and absence of global/RNG mutations.
  search_before <- search()
  original_timeROC_environment <- environment(timeROC::timeROC)
  with_rng(seed, {
    rng_before <- .Random.seed
    check_survdeepfit_inference_metrics(core)
    stopifnot(identical(rng_before, .Random.seed))
  })
  stopifnot(identical(search_before, search()),
            identical(original_timeROC_environment, environment(timeROC::timeROC)),
            identical(body(core$timeROC), body(timeROC::timeROC)),
            identical(formals(core$timeROC), formals(timeROC::timeROC)))
  cat("PASS: real legacy timeROC/IPCW/formula path; package scope, search path, and RNG preserved.\n")
} else {
  cat("SKIP: deepTL prediction fixtures need loadable deepTL, Hmisc, timeROC, and magrittr.\n")
}

# Exact fixtures: inclusive Breslow denominator, risk direction, and uncensored IBS.
train <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
probability <- cox_survival_probabilities(train, rep(0, 3), 0, c(0.5, 1, 2, 3))
stopifnot(isTRUE(all.equal(as.numeric(probability), exp(-c(0, 1/3, 1/3 + 1/2, 1/3 + 1/2)))))
stopifnot(harrell_cindex(c(1, 2, 3), c(1, 1, 1), c(3, 2, 1)) == 1,
          harrell_cindex(c(1, 2, 3), c(1, 1, 1), c(1, 2, 3)) == 0)
test <- data.frame(time = c(1, 2, 3), status = c(1, 1, 1))
perfect <- outer(test$time, c(0.5, 1.5, 2.5), ">") * 1
stopifnot(integrated_brier_score(test, perfect, c(0.5, 1.5, 2.5)) == 0,
          !capture_stage({ warning("did not converge"); 1 })$ok,
          capture_stage({ warning("review this warning"); 1 })$status == "warning")
tied_test <- data.frame(time = c(1, 1, 2), status = c(1, 0, 1))
stopifnot(isTRUE(all.equal(
  integrated_brier_score(tied_test, matrix(0.5, 3, 2), c(1, 1.5)), 5 / 24)))

# Mock expensive methods to exercise retry/resume and serialization without training.
mock <- new.env(parent = globalenv())
sys.source(file.path(root, "R/workflow.R"), mock)
mock$calls <- 0L
mock$fit_feature_selection <- function(...) {
  mock$calls <- mock$calls + 1L
  if (mock$calls == 1L) stop("deliberate smoke-test failure")
  list(model = list(), features = features$features, feature_metrics = features$metrics)
}
mock$fit_model <- function(...) list()
mock$evaluate_prediction <- function(...) data.frame(
  cindex = 0.7, ibs = 0.15, ibs_status = "ok", ibs_error = "", ibs_warnings = "",
  ibs_start = 1, ibs_end = 2, fit_initialization = NA_character_,
  fit_iterations = NA_integer_, fit_max_abs_score = NA_real_)
validation_root <- file.path(root, "results", "validation")
dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)
directory <- tempfile("smoke-", tmpdir = validation_root)
context <- list(scenario = 1, n_train = 120, n_test = 24, rho = 0, repetition = 1, base_seed = config$base_seed)
for (i in 1:3) {
  mock$run_method("Cox", data, config, seed, directory, context, "smoke", "unused", 1)
  if (i == 1) {
    failed_record <- readRDS(file.path(directory, "result.rds"))
    failed_summary <- summarize_records(list(failed_record))
    stopifnot(failed_summary$feature_recovery$n_attempted == 1,
              failed_summary$feature_recovery$n_missing_type1_error == 1,
              all(failed_summary$prediction_summary$n_missing_cindex == 1))
  }
}
stopifnot(mock$calls == 2L)
record <- readRDS(file.path(directory, "result.rds"))
stopifnot(record$selection_status == "ok", nrow(record$predictions) == 2)
failed_record$context$repetition <- 2
summary <- summarize_records(list(record, failed_record))
write_summary(summary, file.path(directory, "summary"))
stopifnot(summary$feature_recovery$mean_average_power == 1,
          summary$feature_recovery$n_attempted == 2,
          summary$feature_recovery$n_missing_average_power == 1,
          nrow(readRDS(file.path(directory, "summary/summary.rds"))$prediction_summary) == 2,
          file.exists(file.path(directory, "summary/feature_rates.csv")),
          inherits(try(read_checkpoint(file.path(directory, "selection.rds"), "changed"), silent = TRUE), "try-error"))
cat("PASS: syntax, configs, data mechanisms, RNG streams, native Cox/AFT fits, metrics, failures, retry/resume, and aggregation.\n")
cat("Temporary smoke outputs:", directory, "\n")
