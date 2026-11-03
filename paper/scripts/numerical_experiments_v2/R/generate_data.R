generate_data <- function(config, n_train, rho, rng_state) {
  require_package("MASS")
  with_rng(rng_state, {
    n_test <- as.integer(round(n_train * config$test_fraction))
    n <- n_train + n_test
    block_size <- config$block_size
    covariance <- (1 - rho) * diag(block_size) + rho

    continuous <- do.call(cbind, lapply(seq_len(4), function(i) {
      MASS::mvrnorm(n, rep(0, block_size), covariance)
    }))
    binary <- matrix(stats::rbinom(n * block_size, 1, 0.4), n, block_size)
    colnames(continuous) <- paste0("X", seq_len(4 * block_size))
    colnames(binary) <- paste0("Z", seq_len(block_size))
    features <- data.frame(continuous, binary, check.names = FALSE)

    important <- c("X1", paste0("X", c(1, 2, 3) * block_size + 1), "Z1")
    if (config$scenario == 1) {
      log_risk <- as.numeric(as.matrix(features[, important]) %*% c(0.1, 0.2, 0.3, 0.4, 0.5))
    } else {
      important <- c(important, "Z2")
      log_risk <- 0.5 * features$Z1 + features$X1 * features$Z2 +
        0.25 * features[[paste0("X", block_size + 1)]]^2 + 0.5 *
        features[[paste0("X", 2 * block_size + 1)]] *
        features[[paste0("X", 3 * block_size + 1)]]
    }

    # Remove the scaling factor
    # log_risk <- log_risk / sqrt(n_train / 1000)

    # Exact inversion of the manuscript's Gompertz cumulative hazard.
    # This avoids simsurv's internal set.seed(), without changing the distribution.
    gamma <- config$baseline$gamma
    lambda <- config$baseline$lambda
    log_hazard_argument <- log(-log(stats::runif(n))) + log(gamma / lambda) - log_risk
    true_time <- (pmax(log_hazard_argument, 0) +
      log1p(exp(-abs(log_hazard_argument)))) / gamma
    censor_time <- stats::rweibull(n, config$censoring$shape, config$censoring$scale)
    observed <- data.frame(time = pmin(true_time, censor_time),
                           status = as.integer(true_time <= censor_time), features)

    list(train = observed[seq_len(n_train), , drop = FALSE],
         test = observed[n_train + seq_len(n_test), , drop = FALSE],
         feature_names = names(features), important_features = important,
         true_time = true_time, censor_time = censor_time, log_risk = log_risk,
         censoring_train = mean(observed$status[seq_len(n_train)] == 0),
         censoring_test = mean(observed$status[n_train + seq_len(n_test)] == 0))
  })
}
