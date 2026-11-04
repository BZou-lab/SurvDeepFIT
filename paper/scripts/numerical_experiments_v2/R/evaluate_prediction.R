predict_risk <- function(method, model, data, features) {
  newdata <- data[, features, drop = FALSE]
  prediction <- switch(method,
    Cox = stats::predict(model, newdata = newdata, type = "lp", reference = "zero"),
    AFT = -stats::predict(model, newdata = newdata, type = "lp"),
    RSF = stats::predict(model, newdata = newdata)$predicted,
    XGBoost = stats::predict(model, xgboost::xgb.DMatrix(as.matrix(newdata)), outputmargin = TRUE),
    SurvDeepFIT = predict_deep_model(model, as.matrix(newdata))
  )
  prediction <- as.numeric(prediction)
  if (length(prediction) != nrow(data) || any(!is.finite(prediction))) {
    stop(method, " returned missing, nonfinite, or incorrectly sized risk predictions.")
  }
  prediction
}

harrell_cindex <- function(time, status, risk) {
  if (length(risk) != length(time) || any(!is.finite(risk))) {
    stop("C-index requires one finite risk prediction per test subject.")
  }
  concordance <- survival::concordance(
    survival::Surv(time, status) ~ risk, reverse = TRUE, timewt = "n")$concordance
  if (!is.finite(concordance)) stop("No estimable concordance in the test set.")
  as.numeric(concordance)
}

cox_survival_probabilities <- function(train, train_log_risk, test_log_risk, times) {
  event_times <- sort(unique(train$time[train$status == 1]))
  if (!length(event_times)) stop("Training set has no observed events.")
  center <- max(train_log_risk)
  relative_risk <- exp(train_log_risk - center)
  increments <- vapply(event_times, function(time) {
    sum(train$time == time & train$status == 1) / sum(relative_risk[train$time >= time])
  }, numeric(1))
  # Inclusive risk sets implement Breslow; include failures at the event time.
  baseline <- c(0, cumsum(increments))[findInterval(times, event_times) + 1L]
  exp(-outer(exp(test_log_risk - center), baseline))
}

predict_survival <- function(method, model, train, test, features, times) {
  if (method == "AFT") {
    lp <- -predict_risk(method, model, test, features)
    survival <- exp(-exp(outer(-lp, log(times), "+") / model$scale))
  } else if (method == "RSF") {
    prediction <- stats::predict(model, newdata = test[, features, drop = FALSE])
    indices <- findInterval(times, prediction$time.interest)
    survival <- cbind(1, prediction$survival)[, indices + 1L, drop = FALSE]
  } else {
    survival <- cox_survival_probabilities(
      train, predict_risk(method, model, train, features),
      predict_risk(method, model, test, features), times)
  }
  if (!identical(dim(survival), c(nrow(test), length(times))) ||
      any(!is.finite(survival)) || any(survival < 0 | survival > 1)) {
    stop(method, " returned invalid survival probabilities.")
  }
  survival
}

integrated_brier_score <- function(test, probabilities, times) {
  censor_fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, data = test)
  # Event weights use G(T-); survivor weights use G(t).
  censor_survival <- c(1, censor_fit$surv)
  event_g <- censor_survival[findInterval(test$time, censor_fit$time, left.open = TRUE) + 1L]
  time_g <- censor_survival[findInterval(times, censor_fit$time) + 1L]
  if (any(time_g <= 0)) stop("IBS interval extends beyond positive test-set censoring support.")
  brier <- vapply(seq_along(times), function(j) {
    event <- test$time <= times[j] & test$status == 1
    alive <- test$time > times[j]
    if (any(event_g[event] <= 0) || (any(alive) && time_g[j] <= 0)) {
      stop("IBS censoring weights are not supported over the requested interval.")
    }
    sum(probabilities[event, j]^2 / event_g[event]) / nrow(test) +
      if (any(alive)) sum((1 - probabilities[alive, j])^2 / time_g[j]) / nrow(test) else 0
  }, numeric(1))
  sum(diff(times) * (head(brier, -1) + tail(brier, -1)) / 2) / diff(range(times))
}

evaluation_times <- function(train, config) {
  interval <- stats::quantile(train$time[train$status == 1],
                              unlist(config$evaluation$ibs_event_quantiles), names = FALSE)
  if (any(!is.finite(interval)) || interval[1] <= 0 || interval[1] >= interval[2]) {
    stop("Training event times do not define a positive IBS interval.")
  }
  seq(interval[1], interval[2], length.out = config$evaluation$ibs_points)
}

evaluate_prediction <- function(method, model, train, test, features, times) {
  risk <- predict_risk(method, model, test, features)
  cindex <- harrell_cindex(test$time, test$status, risk)
  # Keep a valid C-index even when survival probabilities or IPCW cannot be estimated.
  ibs <- capture_stage({
    probabilities <- predict_survival(method, model, train, test, features, times)
    integrated_brier_score(test, probabilities, times)
  })
  solver <- attr(model, "aft_solver")
  data.frame(cindex = cindex, ibs = if (ibs$ok) ibs$value else NA_real_,
             ibs_status = ibs$status, ibs_error = ibs$error,
             ibs_warnings = paste(ibs$warnings, collapse = " | "),
             ibs_start = min(times), ibs_end = max(times),
             fit_initialization = if (is.null(solver)) NA_character_ else solver$initialization,
             fit_iterations = if (is.null(solver)) NA_integer_ else solver$iterations,
             fit_max_abs_score = if (is.null(solver)) NA_real_ else solver$max_abs_score)
}
