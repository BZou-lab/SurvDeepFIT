evaluate_features <- function(method, features, important, pvalues, scores, config) {
  truth <- features %in% important
  if (method == "XGBoost") {
    valid <- rep(FALSE, length(features))
    selected <- rep(FALSE, length(features))
    selected[head(order(-scores, features), config$selection$xgboost_top)] <- TRUE
    ranks <- rank(-scores, ties.method = "average", na.last = "keep")
    rejection <- rep(NA_real_, length(features))
  } else {
    valid <- is.finite(pvalues) & pvalues >= 0 & pvalues <= 1
    selected <- rep(NA, length(features))
    selected[valid] <- pvalues[valid] < config$selection$alpha
    ranks <- rank(replace(pvalues, !valid, NA_real_),
                  ties.method = "average", na.last = "keep")
    rejection <- as.numeric(selected)
  }

  feature_table <- data.frame(feature = features, important = truth,
                              score = as.numeric(scores), pvalue = as.numeric(pvalues),
                              valid_inference = valid, selected = selected,
                              rejection = rejection, rank = as.numeric(ranks))
  # A partial p-value vector must not silently reduce a type-I-error denominator.
  metrics <- data.frame(
    average_power = if (method != "XGBoost" && all(valid[truth])) mean(rejection[truth]) else NA_real_,
    type1_error = if (method != "XGBoost" && all(valid[!truth])) mean(rejection[!truth]) else NA_real_,
    rank_sum = if (all(is.finite(ranks))) sum(ranks[truth]) else NA_real_,
    n_selected = sum(selected %in% TRUE), n_invalid_pvalues = sum(!valid),
    inference_applicable = method != "XGBoost")
  if (method == "XGBoost") metrics$n_invalid_pvalues <- NA_integer_
  list(features = feature_table, metrics = metrics)
}
