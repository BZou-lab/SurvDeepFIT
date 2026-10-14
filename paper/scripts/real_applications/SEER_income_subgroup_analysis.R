# Full-cohort SEER-HCC income subgroup analysis

library(survival)
library(ggplot2)
library(dplyr)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script_arg) == 1L) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)))
} else {
  getwd() # For interactive use, start in the repository root.
}
cohort <- readRDS(file.path(root, "data", "SEER_HCC_2018_2021_raw.rds"))
required <- c("adjusted_household_income", "OS_MONTHS", "OS_STATUS",
              "Race", "Grade", "T_Stage", "Age", "time_to_trt_days")
if (!is.data.frame(cohort) || !all(required %in% names(cohort))) {
  stop("The raw cohort must contain: ", paste(required, collapse = ", "))
}
income <- as.character(cohort$adjusted_household_income)
if (anyNA(income) || any(!income %in% c("<$80,000", "$80,000+"))) {
  stop("Every subject must have income classified as <$80,000 or $80,000+.")
}
cohort$income_group <- factor(ifelse(income == "$80,000+", "High", "Low"),
                              levels = c("Low", "High"))
group_counts <- table(cohort$income_group)
if (any(group_counts < 2L)) stop("Each income group must contain at least two subjects.")
numeric_vars <- c("OS_MONTHS", "OS_STATUS", "Age", "time_to_trt_days")
if (!all(vapply(cohort[numeric_vars], is.numeric, logical(1))) ||
    any(!is.finite(as.matrix(cohort[numeric_vars]))) ||
    any(cohort$OS_MONTHS < 0) || any(cohort$Age < 0) ||
    any(cohort$time_to_trt_days < 0) || !all(cohort$OS_STATUS %in% c(0, 1))) {
  stop("Outcomes, ages, and treatment times must be complete and valid numeric values.")
}
race_levels <- c("White", "Black", "Asian or Pacific Islander")
if (anyNA(cohort$Race) || any(!cohort$Race %in% race_levels) ||
    anyNA(cohort$Grade) || any(!cohort$Grade %in% c("1", "2", "3"))) {
  stop("Race and grade must be complete and match the existing cohort categories.")
}
cohort$Race <- factor(cohort$Race, levels = race_levels)
cohort$Grade <- factor(cohort$Grade, levels = c("1", "2", "3"),
                       labels = c("I", "II", "III"))
t_stage_levels <- c("T1", "T2", "T3", "T4")
if (anyNA(cohort$T_Stage) || any(!cohort$T_Stage %in% t_stage_levels)) {
  stop("T stage must be complete and classified as T1, T2, T3, or T4.")
}
cohort$T_Stage <- factor(cohort$T_Stage, levels = t_stage_levels)

# Color and label settings for figures 
income_colors <- c(Low = "#0072B2", High = "#D55E00")
category_colors <- c("#0072B2", "#D55E00", "#009E73", "#E69F00")
income_labels <- c(Low = "Low (<$80,000)", High = "High ($80,000+)")
income_axis_labels <- setNames(
  paste0(income_labels, "\n(N = ", format(as.integer(group_counts), big.mark = ","), ")"),
  names(group_counts)
)
income_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.spacing.x = grid::unit(0.1, "cm"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_blank()
  )

# KM curves from original raw times
km_fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ income_group, data = cohort)
km_summary <- summary(km_fit, censored = TRUE)
km_data <- data.frame(
  time = km_summary$time,
  survival = km_summary$surv,
  income_group = sub("^income_group=", "", as.character(km_summary$strata))
)
km_data <- bind_rows(
  data.frame(time = 0, survival = 1, income_group = c("Low", "High")),
  km_data
)
km_data$income_group <- factor(km_data$income_group, levels = c("Low", "High"))
km_data <- arrange(km_data, income_group, time)
logrank <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ income_group, data = cohort)
logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE)
p_label <- if (logrank_p < 0.0001) "Log-rank p < 0.0001" else
  paste0("Log-rank p = ", format.pval(logrank_p, digits = 3))
km_plot <- ggplot(km_data, aes(time, survival, color = income_group)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = income_colors, labels = income_labels, name = "Income") +
  scale_x_continuous(breaks = seq(0, ceiling(max(cohort$OS_MONTHS) / 10) * 10, 10)) +
  coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = 0.05 * max(cohort$OS_MONTHS), y = 0.15,
           label = p_label, hjust = 0, size = 3.8) +
  labs(x = "Time (months)", y = "Survival probability") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) + income_theme

# Percentages use all subjects within each income group as the denominator.
make_income_stacked_bar <- function(variable, legend_title) {
  plot_data <- cohort %>%
    count(income_group, category = .data[[variable]], .drop = FALSE) %>%
    group_by(income_group) %>%
    mutate(proportion = n / sum(n), label = sprintf("%.1f%%", 100 * proportion)) %>%
    ungroup()
  ggplot(plot_data, aes(income_group, proportion, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 3.8) +
    scale_fill_manual(values = setNames(category_colors[seq_len(nlevels(cohort[[variable]]))],
                                       levels(cohort[[variable]])),
                      name = legend_title, drop = FALSE) +
    scale_x_discrete(labels = income_axis_labels) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Income", y = "Percentage within income group") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + income_theme
}

# Adjust the box widths, whiskers, and outliers
make_income_boxplot <- function(variable, y_label) {
  ggplot(cohort, aes(income_group, .data[[variable]], fill = income_group)) +
    geom_boxplot(width = 0.6, coef = 1.5, linewidth = 0.4,
                 outlier.size = 1, outlier.alpha = 0.6) +
    scale_fill_manual(values = income_colors, guide = "none") +
    scale_x_discrete(labels = income_axis_labels) +
    labs(x = "Income", y = y_label) + income_theme
}

plots <- list(
  KM = km_plot,
  race = make_income_stacked_bar("Race", "Race"),
  tumor_grade = make_income_stacked_bar("Grade", "Tumor grade"),
  t_stage = make_income_stacked_bar("T_Stage", "T stage"),
  age = make_income_boxplot("Age", "Age at diagnosis (years)"),
  time_to_treatment = make_income_boxplot(
    "time_to_trt_days", "Time from diagnosis to treatment (days)")
)

# Summarize all subjects within each income group
treatment_time_summary <- cohort %>%
  group_by(income_group) %>%
  summarise(
    mean_days = mean(time_to_trt_days),
    median_days = median(time_to_trt_days),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("Mean: %.1f days\nMedian: %.1f days", mean_days, median_days))
treatment_label_y <- max(cohort$time_to_trt_days) +
  0.11 * max(1, max(cohort$time_to_trt_days))
plots$time_to_treatment <- plots$time_to_treatment +
  geom_text(
    data = treatment_time_summary,
    aes(x = income_group, label = label), y = treatment_label_y,
    inherit.aes = FALSE, size = 3.8, fontface = "bold"
  ) +
  expand_limits(y = treatment_label_y) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

figure_dir <- file.path(root, "figures", "seer_income_analysis")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
for (name in names(plots)) {
  stem <- file.path(figure_dir, paste0("SEER_income_", name))
  # The race panel needs a little more width for the three-category legend.
  width <- if (name == "race") 7 else 6
  ggsave(paste0(stem, ".pdf"), plot = plots[[name]], width = width, height = 4.5,
         units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), plot = plots[[name]], width = width, height = 4.5,
         units = "in", dpi = 600, bg = "white")
}
message("Saved six income subgroup figures (PDF and PNG) in ", figure_dir,
        "; all ", nrow(cohort), " subjects: Low n = ", group_counts[["Low"]],
        ", High n = ", group_counts[["High"]], ".")


# Adjusted Cox model analysis ---------------------------------------------------

# Create factor variables
cohort$N_Stage <- as.factor(cohort$N_Stage)
cohort$M_Stage <- as.factor(cohort$M_Stage)
cohort$AFP_status <- as.factor(cohort$AFP_status)

# Fit a multivariable Cox proportional hazards model
Cox_fit <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ 
                   income_group + Age + Race + Sex + 
                   T_Stage + N_Stage + M_Stage + tumor_size + Grade + AFP_status, 
                 data = cohort)

# Model summary
cox_summary <- summary(Cox_fit)

# View the results
cox_results <- data.frame(
  Variable = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[, "exp(coef)"],
  Lower_95CI = cox_summary$conf.int[, "lower .95"],
  Upper_95CI = cox_summary$conf.int[, "upper .95"],
  P_value = cox_summary$coefficients[, "Pr(>|z|)"]
)

print(cox_results)

# Test proportional hazards assumption
ph_test <- cox.zph(Cox_fit)

# View results
print(ph_test)

# Fit a stratified Cox model, where strata are created for variables violating the PH assumption 
Cox_stratified <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~
    income_group + Age + Race + Sex +
    N_Stage + M_Stage + tumor_size +
    strata(T_Stage) +
    strata(Grade) +
    strata(AFP_status),
  data = cohort
)

summary(Cox_stratified)

# Recheck proportional hazards assumption
cox.zph(Cox_stratified)

# Save model estimates and PH diagnostics ---------------------------------------

table_dir <- file.path(root, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# Use the same estimate columns as the unstratified cox_results table.
# Stratification variables have separate baseline hazards, so they have no HR
# coefficients and will not appear in the stratified model's estimate table.
stratified_summary <- summary(Cox_stratified)
cox_stratified_results <- data.frame(
  Variable = rownames(stratified_summary$coefficients),
  HR = stratified_summary$coefficients[, "exp(coef)"],
  Lower_95CI = stratified_summary$conf.int[, "lower .95"],
  Upper_95CI = stratified_summary$conf.int[, "upper .95"],
  P_value = stratified_summary$coefficients[, "Pr(>|z|)"],
  row.names = NULL
)

# Retain both the term-specific and GLOBAL Schoenfeld-residual tests
ph_test_stratified <- cox.zph(Cox_stratified)
make_ph_results <- function(test) {
  data.frame(
    Variable = rownames(test$table),
    Chi_square = test$table[, "chisq"],
    DF = test$table[, "df"],
    P_value = test$table[, "p"],
    row.names = NULL
  )
}
cox_ph_results <- make_ph_results(ph_test)
cox_stratified_ph_results <- make_ph_results(ph_test_stratified)

# Save results as tables in LaTeX format 
format_cox_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
escape_cox_latex <- function(x) {
  # Escape special characters in coefficient names and table text for LaTeX.
  for (character in c("_", "%", "&", "$", "#")) {
    x <- gsub(character, paste0("\\", character), x, fixed = TRUE)
  }
  gsub("<", "\\textless{}", x, fixed = TRUE)
}
save_cox_table <- function(display_table, stem, caption, note) {
  # Write a standalone table environment, matching the project's .tex tables.
  lines <- c(
    "\\begin{table}[!htbp]", "\\centering", "\\small",
    paste0("\\caption{", escape_cox_latex(caption), "}"),
    paste0("\\begin{tabular}{", paste(rep("l", ncol(display_table)), collapse = ""), "}"),
    "\\toprule",
    paste0(paste(escape_cox_latex(names(display_table)), collapse = " & "), " \\\\"),
    "\\midrule"
  )
  for (i in seq_len(nrow(display_table))) {
    cells <- as.character(unlist(display_table[i, ], use.names = FALSE))
    lines <- c(lines, paste0(paste(escape_cox_latex(cells), collapse = " & "), " \\\\"))
  }
  writeLines(c(
    lines, "\\bottomrule", "\\end{tabular}",
    "\\begin{minipage}{0.95\\linewidth}", "\\vspace{0.5em}",
    paste0("\\footnotesize ", escape_cox_latex(note)),
    "\\end{minipage}", "\\end{table}"
  ), file.path(table_dir, paste0(stem, ".tex")))
}
# Translate model terms only for display 
cox_variable_labels <- c(
  income_group = "Income Group", Age = "Age", Race = "Race", Sex = "Sex",
  T_Stage = "T Stage", N_Stage = "N Stage", M_Stage = "M Stage",
  tumor_size = "Tumor Size", Grade = "Tumor Grade", AFP_status = "AFP Status"
)
format_cox_variable <- function(terms) {
  vapply(terms, function(term) {
    matches <- names(cox_variable_labels)[startsWith(term, names(cox_variable_labels))]
    if (!length(matches)) return(term) # Preserve GLOBAL and any unknown terms.
    variable <- matches[which.max(nchar(matches))]
    level <- substring(term, nchar(variable) + 1L)
    label <- unname(cox_variable_labels[variable])
    if (nzchar(level)) paste0(label, " (", level, ")") else label
  }, character(1), USE.NAMES = FALSE)
}
format_cox_estimates <- function(results) {
  # Three decimals retain the small per-mm tumor-size effect and its interval.
  data.frame(
    Variable = format_cox_variable(results$Variable),
    HR = sprintf("%.3f", results$HR),
    `95% CI` = sprintf("(%.3f, %.3f)", results$Lower_95CI, results$Upper_95CI),
    `P-value` = format_cox_p(results$P_value),
    check.names = FALSE
  )
}
format_ph_results <- function(results) {
  data.frame(
    Variable = format_cox_variable(results$Variable),
    `Chi-square` = sprintf("%.2f", results$Chi_square),
    DF = results$DF,
    `P-value` = format_cox_p(results$P_value),
    check.names = FALSE
  )
}

# Report the actual number of observations and events used by each model,
# including any complete-case exclusions performed by coxph.
model_note <- function(fit) {
  references <- vapply(fit$xlevels, function(levels) levels[1], character(1))
  references <- references[!startsWith(names(references), "strata(")]
  paste0(
    "N = ", fit$n, "; events = ", fit$nevent,
    ". HR: adjusted hazard ratio; CI: confidence interval. P-values are Wald tests. ",
    "Categorical reference levels: ",
    paste(paste(format_cox_variable(names(references)), references, sep = " = "), collapse = "; "),
    ". Age and tumor size HRs are per year and per mm, respectively."
  )
}
save_cox_table(
  format_cox_estimates(cox_results),
  "seer_hcc_income_cox_unstratified",
  "SEER-HCC income analysis: adjusted unstratified Cox model",
  model_note(Cox_fit)
)
save_cox_table(
  format_cox_estimates(cox_stratified_results),
  "seer_hcc_income_cox_stratified",
  "SEER-HCC income analysis: adjusted stratified Cox model",
  paste0(model_note(Cox_stratified),
         " Separate baseline hazards are used for T stage, tumor grade, and AFP status; ",
         "these stratification variables have no estimated HRs.")
)
ph_note <- paste0(
  "Schoenfeld-residual tests of the proportional hazards (PH) assumption ",
  "using the KM time transformation. Each row tests a model term jointly; ",
  "GLOBAL tests all estimated coefficients jointly. Small p-values suggest ",
  "departure from PH; large p-values do not establish that PH holds."
)
save_cox_table(
  format_ph_results(cox_ph_results),
  "seer_hcc_income_cox_unstratified_ph",
  "SEER-HCC income analysis: PH checks for the unstratified Cox model", ph_note
)
save_cox_table(
  format_ph_results(cox_stratified_ph_results),
  "seer_hcc_income_cox_stratified_ph",
  "SEER-HCC income analysis: PH checks for the stratified Cox model",
  paste0(ph_note, " Stratification variables are not tested as coefficients.")
)
message("Saved both Cox model results and their separate PH-check tables ",
        "as LaTeX in ", table_dir, ".")


