# Manuscript pages 6-7 layout, using the final full-precision summary CSVs.
# Usage: Rscript --vanilla make_feature_selection_tex.R --study STUDY_DIRECTORY
# Writes LaTeX fragments only; never invokes a TeX compiler.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1] != "--study") {
  stop("Usage: Rscript make_feature_selection_tex.R --study STUDY_DIRECTORY")
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
coverage <- read.csv(file.path(input, "coverage.csv"))

for (scenario in 1:2) {
  data <- read.csv(file.path(input, sprintf("feature_selection_scenario%d.csv", scenario)))
  cover <- subset(coverage, scenario == data$scenario[1])
  stopifnot(all(data$scenario == scenario), all(cover$n_expected == 1000),
            all(cover$n_recorded == 1000), all(cover$n_missing_repetitions == 0))
  metrics <- c("X1", "X11", "X21", "X31", "Z1", if (scenario == 2) "Z2",
               "Avg Power", "Type I Error", "Avg Rank Sum")
  stopifnot(nrow(data) == 6L * length(metrics))
  caption <- if (scenario == 1) {
    "Feature selection performance under linear survival mechanisms"
  } else {
    "Feature selection performance under nonlinear survival mechanisms"
  }
  tex <- c(
    "% LaTeX fragment: include with \\input{tables/feature_selection_scenarioN.tex}.",
    "% Required document preamble (load once; before \\begin{document}):",
    "% \\usepackage{booktabs}",
    "% \\usepackage{arydshln}",
    "% Use table* for a full-width table in a two-column manuscript.",
    "\\begin{table*}[htbp]", "\\centering",
    paste0("\\caption{", caption, "}"),
    sprintf("\\label{tab:simulation-feature-selection-scenario%d}", scenario),
    "\\vspace{4pt}",
    "\\begingroup", "\\scriptsize", "\\sffamily",
    "\\setlength{\\tabcolsep}{4pt}", "\\renewcommand{\\arraystretch}{1.03}",
    "\\begin{tabular}{llcccc}")
  for (n in c(1000, 3000)) {
    letter <- if (n == 1000) "a" else "b"
    if (tail(tex, 1L) == "\\hline") {
      tex[length(tex)] <- "\\midrule"
    } else if (tail(tex, 1L) != "\\midrule") {
      tex <- c(tex, "\\midrule")
    }
    tex <- c(tex,
      sprintf("\\multicolumn{6}{c}{\\textbf{Scenario %d%s:} $n=%s$} \\\\",
              scenario, letter, if (n == 1000) "1{,}000" else "3{,}000"),
      "\\hline",
      paste0("$\\rho$ & \\textbf{Feature} & ",
        paste(paste0("\\textbf{", methods, "}"), collapse = " & "), " \\\\"),
      "\\hline")
    for (rho in c(0, 0.3, 0.5)) {
      for (metric in metrics) {
        row <- data[data$n == n & data$rho == rho & data$metric == metric, ]
        stopifnot(nrow(row) == 1L)
        digits <- if (metric == "Avg Rank Sum") 2L else 3L
        values <- vapply(methods, function(method) {
          x <- row[[method]]
          if (!is.finite(x)) "NA" else sprintf(paste0("%.", digits, "f"), x)
        }, character(1))
        if (metric == "Avg Power") tex <- c(tex, "\\hdashline[0.5pt/1.5pt]")
        tex <- c(tex,
          paste0(paste(c(if (metric == "X1") as.character(rho) else "",
                         labels[[metric]], values), collapse = " & "), " \\\\"))
      }
      tex <- c(tex, "\\hline")
    }
  }
  tex <- c(tex, "\\end{tabular}", "\\endgroup", "\\end{table*}")
  tex <- gsub("scenarioN", paste0("scenario", scenario), tex, fixed = TRUE)
  file <- file.path(output, sprintf("feature_selection_scenario%d.tex", scenario))
  writeLines(tex, file)
  message("Wrote ", file)
}
