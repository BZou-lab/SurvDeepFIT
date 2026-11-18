# Manuscript Tables 1/2 layout, using the existing failure-aware summaries.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4, identical(args[c(1, 3)], c("--input", "--output")))
summary <- readRDS(file.path(args[2], "summary.rds"))
output <- args[4]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
methods <- c("SurvDeepFIT", "AFT", "Cox", "RSF", "XGBoost")
planned_repetitions <- max(summary$coverage$n_expected)
format_value <- function(x, digits) {
  if (length(x) != 1 || !is.finite(x)) "NA" else formatC(x, format = "f", digits = digits)
}
for (scenario in sort(unique(summary$coverage$scenario))) {
  important <- c("X1", "X11", "X21", "X31", "Z1", if (scenario == 2) "Z2")
  rows <- list()
  markdown <- c(sprintf("# Scenario %d: feature selection (%d planned repetitions)", scenario, planned_repetitions),
    "", "Cells average available repetitions. NA denotes unavailable or inapplicable inference.",
    "XGBoost has no feature-level p-values; only its rank sum is reported.",
    "See coverage.csv and feature_rates.csv/feature_recovery.csv for recorded and valid counts.",
    "The current design uses 0.8 times the original effect coefficients without sample-size scaling, and RSF permutation/nonparametric inference.", "")
  for (n in c(1000, 3000)) {
    markdown <- c(markdown, sprintf("## Training n = %d", n), "",
      paste(c("rho", "Feature / metric", methods), collapse = " | "),
      paste(rep("---", length(methods) + 2), collapse = " | "))
    for (rho in c(0, 0.3, 0.5)) for (metric in c(important, "Avg Power", "Type I Error", "Avg Rank Sum")) {
      row <- data.frame(scenario = scenario, n = n, rho = rho, metric = metric)
      for (method in methods) {
        if (metric %in% important) {
          z <- subset(summary$feature_rates, scenario == row$scenario & n == row$n & rho == row$rho & Method == method & feature == metric)
          value <- z$mean_rejection
        } else {
          z <- subset(summary$feature_recovery, scenario == row$scenario & n == row$n & rho == row$rho & Method == method)
          key <- c("Avg Power" = "mean_average_power", "Type I Error" = "mean_type1_error", "Avg Rank Sum" = "mean_rank_sum")[[metric]]
          value <- z[[key]]
        }
        row[[method]] <- if (length(value) == 1) value else NA_real_
      }
      rows[[length(rows) + 1L]] <- row
      digits <- if (metric == "Avg Rank Sum") 2 else 3
      markdown <- c(markdown, paste(c(rho, metric, vapply(row[methods], format_value, character(1), digits = digits)), collapse = " | "))
    }
    markdown <- c(markdown, "")
  }
  write.csv(do.call(rbind, rows), file.path(output, sprintf("feature_selection_scenario%d.csv", scenario)), row.names = FALSE)
  writeLines(markdown, file.path(output, sprintf("feature_selection_scenario%d.md", scenario)))
}
message("Wrote manuscript-style feature-selection tables for scenarios: ",
        paste(sort(unique(summary$coverage$scenario)), collapse = ", "), ".")
