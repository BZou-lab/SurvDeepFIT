as_survival_input <- function(data, features) {
  deepTL::importDnnetSurv(x = data[, features, drop = FALSE],
                         y = data$time, e = data$status)
}

network_controls <- function(config) {
  settings <- config$models$survdeepfit
  list(n.hidden = unlist(settings$hidden_nodes), activate = settings$activation,
       l1.reg = settings$l1_penalty, early.stop.det = settings$early_stop,
       n.batch = settings$batch_size, n.epoch = settings$epochs,
       learning.rate.adaptive = "adam", plot = FALSE)
}

load_legacy_inference <- function(core_file) {
  for (package in c("deepTL", "magrittr", "survival", "Hmisc", "timeROC")) {
    require_package(package)
  }
  core <- new.env(parent = globalenv())
  core$`%>%` <- getExportedValue("magrittr", "%>%")
  for (name in c("appendArg", "importDnnetSurv", "splitDnnet")) {
    core[[name]] <- getExportedValue("deepTL", name)
  }
  core$Surv <- survival::Surv
  core$rcorr.cens <- Hmisc::rcorr.cens
  core$timeROC <- timeROC::timeROC
  sys.source(core_file, envir = core)
  core
}

fit_model <- function(method, train, features, config, threads = 1L) {
  formula <- stats::as.formula("Surv(time, status) ~ .", env = asNamespace("survival"))
  data <- train[, c("time", "status", features), drop = FALSE]

  model <- switch(method,
    Cox = survival::coxph(formula, data = data, x = TRUE, y = TRUE),
    AFT = survival::survreg(formula, data = data, dist = config$models$aft$distribution),
    RSF = {
      require_package("randomForestSRC")
      randomForestSRC::rfsrc(formula, data = data, ntree = config$models$rsf$ntree)
    },
    XGBoost = {
      require_package("xgboost")
      labels <- ifelse(train$status == 1, train$time, -train$time)
      matrix <- xgboost::xgb.DMatrix(as.matrix(train[, features, drop = FALSE]), label = labels)
      settings <- config$models$xgboost
      xgboost::xgb.train(
        params = list(objective = "survival:cox", max_depth = settings$max_depth,
                      eta = settings$eta, nthread = threads,
                      seed = sample.int(.Machine$integer.max, 1)),
        data = matrix, nrounds = settings$nrounds, verbose = 0)
    },
    SurvDeepFIT = {
      require_package("deepTL")
      deepTL::ensemble_dnnet(as_survival_input(train, features),
                            n.ensemble = config$models$survdeepfit$ensemble_size,
                            esCtrl = network_controls(config))
    },
    stop("Unknown method: ", method)
  )

  if (method %in% c("Cox", "AFT") && any(!is.finite(stats::coef(model)))) {
    stop(method, " fit has nonfinite or aliased coefficients.")
  }
  if (method == "AFT" && (!is.finite(model$scale) || model$scale <= 0)) {
    stop("AFT fit has an invalid scale.")
  }
  model
}

fit_feature_selection <- function(method, data, config, model_rng, inference_rng,
                                  core_file, threads = 1L) {
  features <- data$feature_names
  train <- data$train
  fitted <- NULL

  if (method == "SurvDeepFIT") {
    core <- load_legacy_inference(core_file)
    settings <- config$models$survdeepfit
    inference <- with_rng(inference_rng, core$permfit_survival(
      train = as_survival_input(train, features), method = "ensemble_dnnet",
      k_fold = settings$folds, n_perm = settings$permutations,
      n.ensemble = settings$ensemble_size, esCtrl = network_controls(config)))
    importance <- inference@importance
    pvalues <- stats::setNames(importance[[config$selection$survdeepfit_pvalue]],
                               importance$var_name)
    scores <- stats::setNames(importance$importance, importance$var_name)
  } else {
    fitted <- with_rng(model_rng, fit_model(method, train, features, config, threads))
    if (method %in% c("Cox", "AFT")) {
      table <- if (method == "Cox") summary(fitted)$coefficients else summary(fitted)$table
      p_column <- if (method == "Cox") "Pr(>|z|)" else "p"
      # Match predictor names; never include the AFT intercept or log(scale).
      pvalues <- stats::setNames(table[features, p_column], features)
      scores <- stats::setNames(stats::coef(fitted)[features], features)
      inference <- table
    } else if (method == "RSF") {
      inference <- with_rng(inference_rng, randomForestSRC::extract.subsample(
        randomForestSRC::subsample(fitted, B = config$models$rsf$subsamples)))
      table <- inference[[config$selection$rsf_table]]
      if (is.null(table) || !"pvalue" %in% colnames(table)) {
        stop("RSF inference table lacks pvalue: ", config$selection$rsf_table)
      }
      if (!setequal(rownames(table), features)) {
        stop("RSF importance feature names do not match the candidate predictors.")
      }
      pvalues <- stats::setNames(table[features, "pvalue"], features)
      score_column <- if ("mean" %in% names(table)) "mean" else "median"
      scores <- stats::setNames(table[features, score_column], features)
    } else {
      inference <- xgboost::xgb.importance(feature_names = features, model = fitted)
      scores <- stats::setNames(rep(0, length(features)), features)
      scores[inference$Feature] <- inference$Gain
      pvalues <- stats::setNames(rep(NA_real_, length(features)), features)
    }
  }

  if (!setequal(names(pvalues), features) || !setequal(names(scores), features)) {
    stop(method, " inference omitted or duplicated candidate predictors.")
  }
  feature_results <- evaluate_features(method, features, data$important_features,
                                       pvalues[features], scores[features], config)
  if (method != "XGBoost" && any(!feature_results$features$valid_inference)) {
    warning(method, ": some feature p-values are not estimable; affected metrics remain NA.")
  }
  list(model = fitted, inference = inference, features = feature_results$features,
       feature_metrics = feature_results$metrics)
}
