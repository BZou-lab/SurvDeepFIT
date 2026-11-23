# Create five-method archive tables in stacked and side-by-side layouts.
# Usage: Rscript --vanilla make_feature_selection_tex_with_XGBoost.R --study STUDY_DIRECTORY
# Writes LaTeX only, without modifying the four-method tables or compiling TeX.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1] != "--study") {
  stop("Usage: Rscript make_feature_selection_tex_with_XGBoost.R --study STUDY_DIRECTORY")
}
study <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
input <- file.path(study, "runs", "summary")
output <- file.path(study, "tables_with_XGBoost")
dir.create(output, showWarnings = FALSE)
methods <- c("SurvDeepFIT", "Cox", "AFT", "RSF", "XGBoost")
labels <- c(X1 = "$X_1$", X11 = "$X_{p_0+1}$", X21 = "$X_{2p_0+1}$",
            X31 = "$X_{3p_0+1}$", Z1 = "$Z_1$", Z2 = "$Z_2$",
            "Avg Power" = "Avg Power", "Type I Error" = "Type I Error",
            "Avg Rank Sum" = "Avg Rank Sum")
header <- paste(c("$\\rho$", "\\textbf{Feature}",
                  paste0("\\textbf{", methods, "}")), collapse = " & ")
coverage <- read.csv(file.path(input, "coverage.csv"))

for (scenario in 1:2) {
  data <- read.csv(file.path(input, sprintf("feature_selection_scenario%d.csv", scenario)))
  cover <- coverage[coverage$scenario == scenario, ]
  stopifnot(all(data$scenario == scenario), nrow(cover) == 30L,
            all(cover$n_expected == 1000), all(cover$n_recorded == 1000),
            all(cover$n_missing_repetitions == 0))
  metrics <- c("X1", "X11", "X21", "X31", "Z1", if (scenario == 2) "Z2",
               "Avg Power", "Type I Error", "Avg Rank Sum")
  stopifnot(nrow(data) == 6L * length(metrics))
  caption <- paste("Feature selection performance under",
                    if (scenario == 1) "linear" else "nonlinear", "survival mechanisms")
  heading <- function(n) sprintf(
    "\\multicolumn{7}{c}{\\textbf{Scenario %d%s:} $n=%s$}",
    scenario, if (n == 1000) "a" else "b", if (n == 1000) "1{,}000" else "3{,}000")
  body <- function(sizes) {
    lines <- character()
    for (rho in c(0, 0.3, 0.5)) {
      for (metric in metrics) {
        parts <- lapply(sizes, function(n) {
          row <- data[data$n == n & data$rho == rho & data$metric == metric, ]
          stopifnot(nrow(row) == 1L)
          digits <- if (metric == "Avg Rank Sum") 2L else 3L
          values <- vapply(methods, function(method) {
            value <- row[[method]]
            if (!is.finite(value)) {
              if (method == "XGBoost") "\\textemdash{}" else "NA"
            } else sprintf(paste0("%.", digits, "f"), value)
          }, character(1))
          rho_cell <- if (metric == "X1") as.character(rho) else ""
          if (length(sizes) == 2L && n == 3000) rho_cell <- paste0("\\hspace*{5pt}", rho_cell)
          c(rho_cell, labels[[metric]], values)
        })
        if (metric == "Avg Power") lines <- c(lines, "\\hdashline[0.5pt/1.5pt]")
        lines <- c(lines, paste0(paste(unlist(parts), collapse = " & "), " \\\\"))
      }
      lines <- c(lines, "\\hline")
    }
    lines
  }

  for (layout in c("stacked", "side_by_side")) {
    side_by_side <- layout == "side_by_side"
    suffix <- if (side_by_side) "_side_by_side" else ""
    filename <- sprintf("feature_selection_scenario%d%s.tex", scenario, suffix)
    label <- sprintf("tab:simulation-feature-selection-scenario%d-with-xgboost%s",
                     scenario, if (side_by_side) "-side-by-side" else "")
    spec <- if (side_by_side) "llccccc@{\\hspace{10pt}}|llccccc" else "llccccc"
    tex <- c(
      paste0("% LaTeX fragment: include with \\input{tables_with_XGBoost/", filename, "}."),
      "% Required document preamble (load once; before \\begin{document}):",
      "% \\usepackage{booktabs}", "% \\usepackage{arydshln}",
      "% XGBoost has no feature-level p-values; em dashes denote inapplicable inference.",
      "\\begin{table*}[htbp]", "\\centering", paste0("\\caption{", caption, "}"),
      paste0("\\label{", label, "}"), "\\vspace{4pt}",
      "\\begingroup", "\\scriptsize", "\\sffamily",
      "\\setlength{\\tabcolsep}{4pt}", "\\renewcommand{\\arraystretch}{1.03}",
      paste0("\\begin{tabular}{", spec, "}"))
    if (side_by_side) {
      tex <- c(tex, "\\midrule", paste0(heading(1000), " & ", heading(3000), " \\\\"),
        "\\hline", paste0(header, " & \\hspace*{5pt}", header, " \\\\"),
        "\\hline", body(c(1000, 3000)))
    } else {
      for (n in c(1000, 3000)) {
        if (tail(tex, 1L) == "\\hline") tex[length(tex)] <- "\\midrule"
        else tex <- c(tex, "\\midrule")
        tex <- c(tex, paste0(heading(n), " \\\\"), "\\hline",
                 paste0(header, " \\\\"), "\\hline", body(n))
      }
    }
    tex <- c(tex, "\\end{tabular}", "\\endgroup", "\\end{table*}")
    writeLines(tex, file.path(output, filename))
    message("Wrote ", file.path(output, filename))
  }
}
