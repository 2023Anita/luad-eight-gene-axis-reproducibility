#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survival))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)][1]
script_path <- sub("^--file=", "", file_arg)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
processed <- file.path(root, "data", "processed")
reports <- file.path(root, "reports")
dir.create(processed, recursive = TRUE, showWarnings = FALSE)
dir.create(reports, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(processed, "formal_modeling_dataset.csv")
df <- read.csv(data_path, stringsAsFactors = FALSE)
df$OS_time_days <- as.numeric(df$OS_time_days)
df$OS_event <- as.integer(df$OS_event)
df$age <- suppressWarnings(as.numeric(df$age))
df$gender <- factor(df$gender)
df$stage_group <- factor(df$stage_group, levels = c("I", "II", "III", "IV"))
df$axis_group_median <- factor(df$axis_group_median, levels = c("low", "high"))
df$axis_group_tertile <- factor(df$axis_group_tertile, levels = c("low", "mid", "high"))
df$axis_group_quartile <- factor(df$axis_group_quartile, levels = c("Q1", "Q2", "Q3", "Q4"))
df$four_quadrant <- factor(df$four_quadrant)

surv_obj <- Surv(df$OS_time_days, df$OS_event)

fmt_p <- function(p) {
  if (is.na(p)) return("")
  formatC(p, format = "g", digits = 6)
}

extract_cox <- function(model, model_name) {
  s <- summary(model)
  coefs <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  terms <- rownames(coefs)
  data.frame(
    model = model_name,
    term = terms,
    n = model$n,
    events = model$nevent,
    HR = ci[terms, "exp(coef)"],
    CI95_low = ci[terms, "lower .95"],
    CI95_high = ci[terms, "upper .95"],
    p = coefs[terms, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}

models <- list(
  cox_axis_continuous = coxph(surv_obj ~ axis_score, data = df),
  cox_axis_median = coxph(surv_obj ~ axis_group_median, data = df),
  cox_risk_immune = coxph(surv_obj ~ risk_score + immune_score, data = df),
  cox_axis_no_stat1 = coxph(surv_obj ~ axis_score_no_STAT1, data = df),
  cox_axis_alt_immune = coxph(surv_obj ~ axis_score_alt_immune, data = df)
)

multivar_df <- df[!is.na(df$age) & df$gender != "" & !is.na(df$stage_group), ]
multivar <- coxph(Surv(OS_time_days, OS_event) ~ axis_score + age + gender + stage_group, data = multivar_df)
models$cox_axis_adjusted_age_gender_stage <- multivar

cox_rows <- do.call(rbind, Map(extract_cox, models, names(models)))
write.csv(cox_rows, file.path(processed, "survival_cox_results.csv"), row.names = FALSE)

lr_axis <- survdiff(surv_obj ~ axis_group_median, data = df)
lr_p <- pchisq(lr_axis$chisq, df = length(lr_axis$n) - 1, lower.tail = FALSE)
lr_tertile <- survdiff(surv_obj ~ axis_group_tertile, data = df)
lr_tertile_p <- pchisq(lr_tertile$chisq, df = length(lr_tertile$n) - 1, lower.tail = FALSE)
lr_quartile <- survdiff(surv_obj ~ axis_group_quartile, data = df)
lr_quartile_p <- pchisq(lr_quartile$chisq, df = length(lr_quartile$n) - 1, lower.tail = FALSE)

logrank_rows <- data.frame(
  comparison = c("axis_group_median", "axis_group_tertile", "axis_group_quartile"),
  chisq = c(lr_axis$chisq, lr_tertile$chisq, lr_quartile$chisq),
  df = c(length(lr_axis$n) - 1, length(lr_tertile$n) - 1, length(lr_quartile$n) - 1),
  p = c(lr_p, lr_tertile_p, lr_quartile_p)
)
write.csv(logrank_rows, file.path(processed, "survival_logrank_results.csv"), row.names = FALSE)

fit <- survfit(surv_obj ~ axis_group_median, data = df)
fit_summary <- summary(fit)
km_rows <- data.frame(
  time = fit_summary$time,
  n.risk = fit_summary$n.risk,
  n.event = fit_summary$n.event,
  survival = fit_summary$surv,
  std.err = fit_summary$std.err,
  lower = fit_summary$lower,
  upper = fit_summary$upper,
  strata = fit_summary$strata
)
write.csv(km_rows, file.path(processed, "km_curve_axis_group_median.csv"), row.names = FALSE)

stage_summary <- aggregate(axis_score ~ stage_group, data = df, FUN = function(x) c(n = length(x), mean = mean(x), median = median(x)))
stage_out <- data.frame(
  stage_group = stage_summary$stage_group,
  n = stage_summary$axis_score[, "n"],
  mean = stage_summary$axis_score[, "mean"],
  median = stage_summary$axis_score[, "median"]
)
write.csv(stage_out, file.path(processed, "axis_score_stage_summary_r.csv"), row.names = FALSE)

axis_cont <- cox_rows[cox_rows$model == "cox_axis_continuous" & cox_rows$term == "axis_score", ]
axis_median <- cox_rows[cox_rows$model == "cox_axis_median", ]
adjusted_axis <- cox_rows[cox_rows$model == "cox_axis_adjusted_age_gender_stage" & cox_rows$term == "axis_score", ]

report <- c(
  "# LUAD Axis Survival Analysis Results",
  "",
  paste0("- Date: ", Sys.Date()),
  "- Runtime: R survival package, local TCGA Xena LUAD starter data only.",
  "- Boundary: no additional data download and no package installation in this step.",
  "",
  "## Log-rank Tests",
  "",
  "| comparison | chisq | df | p |",
  "|---|---:|---:|---:|",
  apply(logrank_rows, 1, function(r) paste0("| ", r[["comparison"]], " | ", round(as.numeric(r[["chisq"]]), 4), " | ", r[["df"]], " | ", fmt_p(as.numeric(r[["p"]])), " |")),
  "",
  "## Cox Models: Key Terms",
  "",
  "| model | term | n | events | HR | 95% CI | p |",
  "|---|---|---:|---:|---:|---|---:|",
  paste0("| continuous axis | axis_score | ", axis_cont$n, " | ", axis_cont$events, " | ", round(axis_cont$HR, 4), " | ", round(axis_cont$CI95_low, 4), "-", round(axis_cont$CI95_high, 4), " | ", fmt_p(axis_cont$p), " |"),
  paste0("| median axis | ", axis_median$term[1], " | ", axis_median$n[1], " | ", axis_median$events[1], " | ", round(axis_median$HR[1], 4), " | ", round(axis_median$CI95_low[1], 4), "-", round(axis_median$CI95_high[1], 4), " | ", fmt_p(axis_median$p[1]), " |"),
  paste0("| adjusted axis | axis_score | ", adjusted_axis$n, " | ", adjusted_axis$events, " | ", round(adjusted_axis$HR, 4), " | ", round(adjusted_axis$CI95_low, 4), "-", round(adjusted_axis$CI95_high, 4), " | ", fmt_p(adjusted_axis$p), " |"),
  "",
  "## Stage Trend",
  "",
  "| stage | n | axis mean | axis median |",
  "|---|---:|---:|---:|",
  apply(stage_out, 1, function(r) paste0("| ", r[["stage_group"]], " | ", r[["n"]], " | ", round(as.numeric(r[["mean"]]), 4), " | ", round(as.numeric(r[["median"]]), 4), " |")),
  "",
  "## Interpretation Guardrail",
  "",
  "- These are internal TCGA-LUAD results and should not yet be written as final manuscript claims.",
  "- The next evidence step should be GEO external validation if the user approves downloading a small validation dataset.",
  "- GDC full STAR-counts remains unnecessary for the current route."
)
writeLines(report, file.path(reports, "survival_analysis_results.md"))

cat("Wrote survival analysis results\n")
