# Summarize one completed scenario without reading the other scenario's results.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(normalizePath(script, winslash = "/"))
source(file.path(root, "summarize_results.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L || !identical(args[c(1, 3, 5)], c("--input", "--scenario", "--output"))) {
  stop("Usage: Rscript summarize_scenario.R --input STUDY_RUNS --scenario 1 --output SEPARATE_OUTPUT")
}
scenario <- as.integer(args[4])
stopifnot(length(scenario) == 1L, !is.na(scenario), scenario %in% 1:2)
input <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[6], winslash = "/", mustWork = FALSE)
if (basename(output) == "summary" || output == input ||
    startsWith(output, paste0(input, "/scenario"))) {
  stop("Use a separate output folder, e.g. STUDY_RUNS/summary_scenario1; reserve summary/ for the full job.")
}
files <- list.files(file.path(input, paste0("scenario", scenario)),
                    pattern = "^result\\.rds$", recursive = TRUE, full.names = TRUE)
message("Reading ", length(files), " saved method results for scenario ", scenario, ".")
records <- lapply(files, readRDS)
summary <- summarize_records(records)
stopifnot(all(summary$coverage$scenario == scenario))
write_summary(summary, output)
rm(records)
invisible(gc())

status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--vanilla", shQuote(file.path(root, "make_feature_selection_tables.R")),
                    "--input", shQuote(output), "--output", shQuote(output)))
if (status != 0L) stop("Feature-selection table generation failed.")

coverage <- summary$coverage
prediction <- summary$prediction_summary
selection <- summary$feature_recovery
selected <- prediction[prediction$feature_set == "selected", ]
overview <- merge(selection, selected, by = c("scenario", "n", "rho", "Method"),
                  suffixes = c("_selection", "_prediction"))
overview <- overview[order(overview$n, overview$rho, overview$Method), ]
columns <- c("scenario", "n", "rho", "Method", "mean_average_power", "mean_type1_error",
             "mean_rank_sum", "mean_cindex", "mean_ibs", "n_valid_average_power",
             "n_valid_type1_error", "n_valid_cindex", "n_valid_ibs")
write.csv(overview[, columns], file.path(output, "overview.csv"), row.names = FALSE)
report <- c(sprintf("# Scenario %d summary", scenario), "",
  paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)),
  paste("Input:", input), "",
  sprintf("Saved method results: %d / %d planned.", sum(coverage$n_recorded), sum(coverage$n_expected)),
  sprintf("Missing method/repetition records: %d. Selection failures: %d.",
          sum(coverage$n_missing_repetitions), sum(coverage$n_selection_failed)),
  sprintf("Failed prediction stages: %d. Failed IBS evaluations: %d.",
          sum(summary$prediction_results$status %in% c("failed", "selection_failed")),
          sum(summary$prediction_results$ibs_status %in% "failed")), "",
  "Means use available valid results; CSV files report valid and missing counts.",
  "XGBoost has no feature-level p-values, so its power and Type I error are unavailable.",
  "C-index and IBS below use selected-feature models; prediction_summary.csv includes full models too.", "",
  "n | rho | Method | Power | Type I error | Rank sum | C-index | IBS",
  "--- | --- | --- | --- | --- | --- | --- | ---")
number <- function(x) if (is.finite(x)) sprintf("%.4f", x) else "NA"
for (i in seq_len(nrow(overview))) {
  row <- overview[i, ]
  values <- vapply(row[c("mean_average_power", "mean_type1_error", "mean_rank_sum",
                         "mean_cindex", "mean_ibs")], number, character(1))
  report <- c(report, paste(c(row$n, row$rho, row$Method, values), collapse = " | "))
}
writeLines(report, file.path(output, "overview.md"))
if (any(coverage$n_missing_repetitions > 0)) warning("Scenario has missing repetitions; inspect coverage.csv.")
message("Wrote scenario ", scenario, " summary to ", output)
