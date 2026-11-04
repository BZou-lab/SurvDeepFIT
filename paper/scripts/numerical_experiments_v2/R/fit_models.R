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

deep_prediction_generic <- function() {
  require_package("deepTL")
  generic <- methods::getGeneric("predict", where = asNamespace("deepTL"))
  if (is.null(generic)) stop("The loaded deepTL build has no S4 prediction generic.")
  generic
}

predict_deep_model <- function(model, newdata) {
  generic <- deep_prediction_generic()
  generic(model, as.matrix(newdata))
}

scoped_timeROC <- function() {
  # timeROC 0.4.1 creates an unqualified Surv() formula inside its own frame.
  # Formula lookup follows that function's lexical scope, not the legacy caller.
  # Keep package helpers visible while supplying Surv in a private child scope.
  fun <- timeROC::timeROC
  scope <- new.env(parent = environment(fun))
  scope$Surv <- survival::Surv
  environment(fun) <- scope
  fun
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
  core$timeROC <- scoped_timeROC()
  sys.source(core_file, envir = core)
  # Namespace loading alone leaves unqualified predict() resolving to the S3 generic.
  core$predict <- deep_prediction_generic()
  core
}

check_survdeepfit_dispatch <- function(core) {
  # A hand-built one-hidden-layer network tests real package methods without training.
  network <- methods::new("dnnet",
    norm = list(x.center = c(0, 0), x.scale = c(1, 1), y.center = 0, y.scale = 1),
    weight = list(diag(2), matrix(c(1, -0.5), ncol = 1)),
    bias = list(matrix(0, 1, 2), matrix(0, 1, 1)),
    loss = 0, loss.traj = 0, label = character(), model.type = "survival",
    model.spec = list(activate = "relu", pathway = FALSE))
  second <- network
  second@bias[[2]] <- matrix(1, 1, 1)
  ensemble <- methods::new("dnnetEnsemble", model.list = list(network, second),
                           model.type = "survival", model.spec = network@model.spec,
                           loss = c(0, 0), keep = c(TRUE, TRUE))
  x <- matrix(c(1, 2, 2, 1, 3, 0), ncol = 2, byrow = TRUE,
              dimnames = list(NULL, c("X1", "X2")))
  input <- deepTL::importDnnetSurv(x = x, y = c(3, 2, 1), e = c(1, 1, 1))
  expected <- c(0, 1.5, 3)
  for (method in c("dnnet", "ensemble_dnnet")) {
    model <- if (method == "dnnet") network else ensemble
    target <- expected + if (method == "dnnet") 0 else 0.5
    legacy <- core$predict_mod_permfit(model, input, method, "survival")
    direct <- predict_deep_model(model, x)
    if (!isTRUE(all.equal(as.numeric(legacy[[1]]), target)) ||
        !isTRUE(all.equal(as.numeric(direct), target))) {
      stop("deepTL prediction dispatch or values failed for ", method, ".")
    }
  }
  risk <- predict_risk("SurvDeepFIT", ensemble, data.frame(x), colnames(x))
  if (!isTRUE(all.equal(risk, expected + 0.5))) stop("V2 ensemble risk prediction failed.")
  invisible(TRUE)
}

check_survdeepfit_inference_metrics <- function(core) {
  # Exercise timeROC -> pec::ipcw -> prodlim -> formula evaluation, without training.
  time <- seq_len(24)
  status <- rep(c(1, 0, 1), 8)
  risk <- -time
  auc <- core$timedAUC(status, time, risk)
  reverse_auc <- core$timedAUC(status, time, -risk)
  input <- deepTL::importDnnetSurv(x = cbind(X1 = time, X2 = -time), y = time, e = status)
  difference <- core$log_lik_diff("survival", risk, -risk, input,
                                 y_hatcoxl = risk, y_hat0coxl = -risk)
  stopifnot(isTRUE(all.equal(auc, 1)), isTRUE(all.equal(reverse_auc, 0)),
            length(difference) == 3, all(is.finite(difference)),
            isTRUE(all.equal(unname(difference[c(1, 3)]), c(1, 1))))
  invisible(TRUE)
}

checked_aft_fit <- function(formula, data, ...) {
  stage <- capture_stage(survival::survreg(formula, data = data, dist = "weibull", ...))
  for (message in stage$warnings) warning(message, call. = FALSE)
  if (!stage$ok) stop("AFT fit rejected: ", stage$error, call. = FALSE)
  model <- stage$value
  if (any(!is.finite(stats::coef(model))) || !is.finite(model$scale) || model$scale <= 0 ||
      any(!is.finite(model$loglik)) || any(!is.finite(model$var)) || any(diag(model$var) <= 0)) {
    stop("AFT fit has invalid coefficients, scale, likelihood, or covariance.")
  }
  model
}

fit_weibull_aft <- function(formula, data, control = survival::survreg.control()) {
  null_formula <- stats::update.formula(formula, . ~ 1)
  baseline <- checked_aft_fit(null_formula, data, control = control, score = TRUE)
  matrix <- stats::model.matrix(stats::delete.response(stats::terms(formula, data = data)), data)
  initial <- stats::setNames(rep(0, ncol(matrix) + 1L), c(colnames(matrix), "Log(scale)"))
  initial["(Intercept)"] <- unname(stats::coef(baseline)["(Intercept)"])
  initial["Log(scale)"] <- log(baseline$scale)
  model <- checked_aft_fit(formula, data, init = unname(initial), control = control, score = TRUE)
  attr(model, "aft_solver") <- list(initialization = "intercept_only", initial = initial,
                                   null_iterations = baseline$iter, iterations = model$iter,
                                   max_abs_score = max(abs(model$score)), control = control)
  model
}

fit_model <- function(method, train, features, config, threads = 1L) {
  formula <- stats::as.formula("Surv(time, status) ~ .", env = asNamespace("survival"))
  data <- train[, c("time", "status", features), drop = FALSE]

  model <- switch(method,
    Cox = survival::coxph(formula, data = data, x = TRUE, y = TRUE),
    AFT = fit_weibull_aft(formula, data),
    RSF = {
      require_package("randomForestSRC")
      settings <- config$models$rsf
      if (is.null(settings$importance)) {
        # Preserve the original grow/subsample path when no RSF option is supplied.
        randomForestSRC::rfsrc(formula, data = data, ntree = settings$ntree)
      } else {
        if (!settings$importance %in% c("anti", "permute", "random")) {
          stop("RSF importance must be anti, permute, or random.")
        }
        randomForestSRC::rfsrc(formula, data = data, ntree = settings$ntree,
                              importance = settings$importance, block.size = 1L)
      }
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
  if (method == "RSF" && model$ntree != config$models$rsf$ntree) {
    stop("RSF returned a different number of trees than requested.")
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
