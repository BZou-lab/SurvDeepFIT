# Alternative layout only: leaves the existing stacked tables unchanged.
# Usage: Rscript --vanilla make_feature_selection_tex_side_by_side.R --study STUDY_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1] != "--study") {
  stop("Usage: Rscript make_feature_selection_tex_side_by_side.R --study STUDY_DIRECTORY")
}
study <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
input <- file.path(study, "runs", "summary")
output <- file.path(study, "tables")
dir.create(output, showWarnings = FALSE)
methods <- c("SurvDeepFIT", "Cox", "AFT", "RSF")
labels <- c(X1 = "$X_1$", X11 = "$X_{p_0+1}$", X21 = "$X_{2p_0+1}$",
            X31 = "$X_{3p_0+1}$", Z1 = "$Z_1$", Z2 = "$Z_2$",
            "Avg Power" = "Avg Power", "Type I Error" = "Type I Error",
            "Avg Rank Sum" = "Avg Rank Sum")
header <- paste(c("$\\rho$", "\\textbf{Feature}",
                  paste0("\\textbf{", methods, "}")), collapse = " & ")
right_header <- paste0("\\hspace*{5pt}", header)
coverage <- read.csv(file.path(input, "coverage.csv"))

for (scenario in 1:2) {
  data <- read.csv(file.path(input, sprintf("feature_selection_scenario%d.csv", scenario)))
  cover <- coverage[coverage$scenario == scenario, ]
  stopifnot(all(data$scenario == scenario), all(cover$n_expected == 1000),
            all(cover$n_recorded == 1000), all(cover$n_missing_repetitions == 0))
  metrics <- c("X1", "X11", "X21", "X31", "Z1", if (scenario == 2) "Z2",
               "Avg Power", "Type I Error", "Avg Rank Sum")
  stopifnot(nrow(data) == 6L * length(metrics))
  caption <- paste("Feature selection performance under",
                    if (scenario == 1) "linear" else "nonlinear", "survival mechanisms")
  filename <- sprintf("feature_selection_scenario%d_side_by_side.tex", scenario)
  tex <- c(
    paste0("% LaTeX fragment: include with \\input{tables/", filename, "}."),
    "% Required document preamble (load once; before \\begin{document}):",
    "% \\usepackage{booktabs}", "% \\usepackage{arydshln}",
    "% Alternative layout; unique labels permit comparison with the stacked tables.",
    "\\begin{table*}[htbp]", "\\centering",
    paste0("\\caption{", caption, "}"),
    sprintf("\\label{tab:simulation-feature-selection-scenario%d-side-by-side}", scenario),
    "\\vspace{4pt}", "\\begingroup", "\\scriptsize", "\\sffamily",
    "\\setlength{\\tabcolsep}{4pt}", "\\renewcommand{\\arraystretch}{1.03}",
    "\\begin{tabular}{llcccc@{\\hspace{10pt}}|llcccc}", "\\midrule",
    sprintf(paste0("\\multicolumn{6}{c}{\\textbf{Scenario %da:} $n=1{,}000$} & ",
                   "\\multicolumn{6}{c}{\\textbf{Scenario %db:} $n=3{,}000$} \\\\"),
            scenario, scenario),
    "\\hline", paste0(header, " & ", right_header, " \\\\"), "\\hline")
  for (rho in c(0, 0.3, 0.5)) {
    for (metric in metrics) {
      halves <- lapply(c(1000, 3000), function(n) {
        row <- data[data$n == n & data$rho == rho & data$metric == metric, ]
        stopifnot(nrow(row) == 1L)
        digits <- if (metric == "Avg Rank Sum") 2L else 3L
        values <- vapply(methods, function(method) {
          x <- row[[method]]
          if (!is.finite(x)) "NA" else sprintf(paste0("%.", digits, "f"), x)
        }, character(1))
        rho_cell <- if (metric == "X1") as.character(rho) else ""
        if (n == 3000) rho_cell <- paste0("\\hspace*{5pt}", rho_cell)
        c(rho_cell, labels[[metric]], values)
      })
      if (metric == "Avg Power") tex <- c(tex, "\\hdashline[0.5pt/1.5pt]")
      tex <- c(tex, paste0(paste(unlist(halves), collapse = " & "), " \\\\"))
    }
    tex <- c(tex, "\\hline")
  }
  tex <- c(tex, "\\end{tabular}", "\\endgroup", "\\end{table*}")
  writeLines(tex, file.path(output, filename))
  message("Wrote ", file.path(output, filename))
}
