# Full-study selected-feature mean +/- SD figures, matching test20-concurrency200.
# Usage: Rscript --vanilla plot_prediction_figures_with_XGBoost.R --study /absolute/study
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1] != "--study") {
  stop("Usage: Rscript plot_prediction_figures_with_XGBoost.R --study STUDY_DIRECTORY")
}
study_dir <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
figure_dir <- file.path(study_dir, "figures_with_XGBoost")
method_levels <- c("SurvDeepFIT", "Cox", "AFT", "RSF", "XGBoost")
method_colors <- c(SurvDeepFIT = "#D55E00", Cox = "#56B4E9", AFT = "#009E73",
                   RSF = "#E69F00", XGBoost = "#CC79A7")
raw <- read.csv(file.path(study_dir, "runs/summary/prediction_results.csv"))
selected <- subset(raw, feature_set == "selected" & Method %in% method_levels)
expected <- expand.grid(scenario = 1:2, n = c(1000, 3000), rho = c(0, 0.3, 0.5),
                        Method = method_levels, repetition = 1:1000,
                        stringsAsFactors = FALSE)
key <- function(d) do.call(paste, c(d[c("scenario", "n", "rho", "Method", "repetition")], sep = "/"))
stopifnot(nrow(selected) == nrow(expected), !anyDuplicated(key(selected)),
          setequal(key(selected), key(expected)))
reference <- subset(read.csv(file.path(study_dir, "runs/summary/prediction_summary.csv")),
                    feature_set == "selected" & Method %in% method_levels)
group_key <- function(d) do.call(paste, c(d[c("scenario", "n", "rho", "Method")], sep = "/"))

make_figure <- function(metric, label, stem) {
  groups <- split(selected, interaction(selected$scenario, selected$n, selected$rho,
                                       selected$Method, drop = TRUE))
  df <- do.call(rbind, lapply(groups, function(d) {
    stopifnot(nrow(d) == 1000L)
    valid <- d$status %in% c("ok", "warning") & is.finite(d[[metric]])
    if (metric == "ibs") valid <- valid & d$ibs_status %in% c("ok", "warning")
    values <- d[[metric]][valid]
    stopifnot(length(values) > 1L, all(values >= 0 & values <= 1))
    data.frame(d[1, c("scenario", "n", "rho", "Method")],
               Mean = mean(values), SD = sd(values), repetitions = length(values),
               planned_repetitions = nrow(d), missing_repetitions = sum(!valid))
  }))
  rownames(df) <- NULL
  matched <- match(group_key(df), group_key(reference))
  stopifnot(!anyNA(matched), nrow(df) == 12L * length(method_levels),
            isTRUE(all.equal(df$Mean, reference[[paste0("mean_", metric)]][matched], tolerance = 1e-12)),
            isTRUE(all.equal(df$SD, reference[[paste0("sd_", metric)]][matched], tolerance = 1e-12)),
            all(df$repetitions == reference[[paste0("n_valid_", metric)]][matched]))
  df$Method <- factor(df$Method, levels = method_levels)
  df$Panel <- factor(paste0("Scenario ", df$scenario, " | n = ", format(df$n, big.mark = ",", trim = TRUE)),
                     levels = c("Scenario 1 | n = 1,000", "Scenario 1 | n = 3,000",
                                "Scenario 2 | n = 1,000", "Scenario 2 | n = 3,000"))
  df$Correlation <- factor(paste0("\u03c1 = ", df$rho),
                           levels = paste0("\u03c1 = ", c(0, 0.3, 0.5)))
  y_limits <- if (metric == "cindex") {
    c(min(0.5, floor(min(df$Mean - df$SD) * 20) / 20),
      min(1, ceiling(max(df$Mean + df$SD) * 20) / 20))
  } else c(0, ceiling(max(df$Mean + df$SD) * 50) / 50)
  p <- ggplot(df, aes(Method, Mean, fill = Method)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15) +
    scale_fill_manual(values = method_colors, breaks = method_levels, name = "Method") +
    facet_grid(Panel ~ Correlation) +
    coord_cartesian(ylim = y_limits) +
    labs(x = NULL, y = label, title = NULL, subtitle = NULL, caption = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 11),
          legend.key.size = grid::unit(4, "mm"), plot.title = element_blank(),
          plot.subtitle = element_blank(), plot.caption = element_blank(),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(size = 11))
  dir.create(figure_dir, showWarnings = FALSE)
  write.csv(df, file.path(figure_dir, paste0(stem, "_data.csv")), row.names = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), p,
         width = 12, height = 11, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), p,
         width = 12, height = 11, units = "in", dpi = 300, bg = "white")
  message(label, ": saved PDF, PNG, and plotted data; valid repetitions per bar: ",
          min(df$repetitions), "-", max(df$repetitions), ".")
}

make_figure("cindex", "C-index", "simulation_CIndex_selected")
make_figure("ibs", "IBS", "simulation_IBS_selected")
