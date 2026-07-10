#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
results <- file.path(root, "results")

data <- read.csv(
  file.path(results, "gse31210_locked_score_analysis_dataset.csv"),
  check.names = FALSE,
  na.strings = c("", "NA")
)
data$gender <- factor(data$gender, levels = c("female", "male"))
data$stage_group <- factor(data$stage_group, levels = c("IA", "IB", "II"))

model <- coxph(
  Surv(RFS_time_months, RFS_event) ~ axis_score + age + stage_group + strata(gender),
  data = data
)
summary_model <- summary(model)
coefficients <- as.data.frame(summary_model$coefficients)
intervals <- as.data.frame(summary_model$conf.int)
terms <- rownames(coefficients)

results_table <- data.frame(
  model = "cox_axis_age_stage_stratified_by_gender",
  term = terms,
  n = model$n,
  events = model$nevent,
  HR = intervals[terms, "exp(coef)"],
  CI95_low = intervals[terms, "lower .95"],
  CI95_high = intervals[terms, "upper .95"],
  p = coefficients[terms, "Pr(>|z|)"],
  row.names = NULL,
  check.names = FALSE
)
write.csv(results_table, file.path(results, "gse31210_gender_stratified_sensitivity.csv"), row.names = FALSE)

ph <- cox.zph(model)
ph_table <- data.frame(
  term = rownames(ph$table),
  chisq = ph$table[, "chisq"],
  df = ph$table[, "df"],
  p = ph$table[, "p"],
  row.names = NULL
)
write.csv(ph_table, file.path(results, "gse31210_gender_stratified_sensitivity_ph.csv"), row.names = FALSE)

axis_row <- results_table[results_table$term == "axis_score", ]
report <- c(
  "# GSE31210 Gender-Stratified Cox Sensitivity Analysis",
  "",
  "- Purpose: address the original adjusted model's gender-specific proportional-hazards diagnostic by stratifying the baseline hazard on gender.",
  "- Scope: the exact eight-gene score, RFS endpoint, age adjustment, stage adjustment, eligibility criteria, and analysis dataset were unchanged.",
  paste0("- Analysis set: ", model$n, " participants; ", model$nevent, " relapses."),
  paste0("- Axis HR: ", sprintf("%.4f", axis_row$HR), "; 95% CI ", sprintf("%.4f", axis_row$CI95_low), "-", sprintf("%.4f", axis_row$CI95_high), "; p=", format(axis_row$p, digits = 4, scientific = TRUE), "."),
  "- Interpretation boundary: this is a pre-specified diagnostic-response sensitivity analysis, not a new score-development or subgroup-discovery analysis."
)
writeLines(report, file.path(results, "gse31210_gender_stratified_sensitivity.md"))
cat("Wrote GSE31210 gender-stratified sensitivity package\n")
