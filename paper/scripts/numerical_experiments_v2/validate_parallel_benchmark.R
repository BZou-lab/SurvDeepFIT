# Small implementation checks only; production YAML files are never modified.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(normalizePath(script))
args <- commandArgs(trailingOnly = TRUE)
output <- args[2]
suppressPackageStartupMessages(library(deepTL))
set.seed(1729)
x <- matrix(rnorm(90 * 5), 90, 5, dimnames = list(NULL, paste0("X", 1:5)))
input <- importDnnetSurv(x, rexp(90), rbinom(90, 1, 0.7))
controls <- list(n.hidden = c(6, 4), activate = "relu", l1.reg = 0.0001,
  early.stop.det = 3, n.batch = 20, n.epoch = 3,
  learning.rate.adaptive = "adam", plot = FALSE)
fit <- function() {
  set.seed(814)
  model <- ensemble_dnnet(input, n.ensemble = 4, esCtrl = controls, verbose = FALSE)
  list(model = model, seed = .Random.seed,
       prediction = methods::getGeneric("predict", where = asNamespace("deepTL"))(model, x))
}
if (args[1] == "original") {
  saveRDS(fit(), output)
  quit(status = 0)
}
stopifnot(as.character(packageVersion("deepTL")) == "1.0.0.9001")
# Original default algorithm, bootstrap draws, fit and RNG advancement survive.
stopifnot(identical(fit(), readRDS(output)))
options(deepTL.ensemble.reproducible = TRUE)
models <- lapply(c(1L, 2L, 4L), function(workers) {
  options(deepTL.ensemble.workers = workers)
  fit()
})
stopifnot(identical(models[[1]], models[[2]]), identical(models[[1]], models[[3]]))
for (file in c("workflow.R", "generate_data.R", "fit_models.R", "evaluate_features.R", "evaluate_prediction.R"))
  source(file.path(root, "R", file))
core_file <- file.path(root, "../../../R/source_permfit_survival_20230610.r")
core <- load_legacy_inference(core_file)
check_survdeepfit_dispatch(core)
check_survdeepfit_inference_metrics(core)
config <- read_config(file.path(root, "configs/scenario1.yaml"))
config$models$survdeepfit$folds <- 2
config$models$survdeepfit$permutations <- 2
config$models$survdeepfit$ensemble_size <- 4
config$models$survdeepfit$epochs <- 3
config$models$survdeepfit$early_stop <- 3
config$models$survdeepfit$hidden_nodes <- c(6, 4)
data <- generate_data(config, 90, 0.3, stage_rng(task_rng(config, 3, 1), 0))
state <- task_rng(config, 3, 1)
selections <- lapply(c(1L, 2L, 4L), function(workers) {
  options(deepTL.ensemble.workers = workers)
  fit_feature_selection("SurvDeepFIT", data, config,
    stage_rng(state, 10), stage_rng(state, 11), core_file, threads = 1L)
})
stopifnot(identical(selections[[1]], selections[[2]]),
          identical(selections[[1]], selections[[3]]))
options(deepTL.ensemble.workers = 0L)
stopifnot(inherits(try(fit(), silent = TRUE), "try-error"))
cat("PASS: unchanged default matches original installed deepTL; 1/2/4 workers match\n",
    "for models, RNG state, predictions, and feature-selection outputs.\n")
