#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survival))

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1] |> sub("^--file=", "", x = _)), ".."), mustWork = TRUE)
input <- file.path(root, "processed")
output <- file.path(root, "results")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

genes <- c("CA9", "MKI67", "STAT1", "HIF1A", "CD4", "CD8A", "CD8B", "TIGIT")
expression <- read.csv(file.path(input, "gse72094_locked_score_expression.csv"), check.names = FALSE, na.strings = c("", "NA"))
clinical <- read.csv(file.path(input, "gse72094_clinical_raw.csv"), check.names = FALSE, na.strings = c("", "NA"))
df <- merge(clinical, expression, by = "sample_id", all = FALSE)

df$OS_time_days <- suppressWarnings(as.numeric(df$survival_time_in_days))
df$OS_event <- ifelse(df$vital_status == "Dead", 1L, ifelse(df$vital_status == "Alive", 0L, NA_integer_))
df$age <- suppressWarnings(as.numeric(df$age_at_diagnosis))
df$gender <- factor(df$gender, levels = c("F", "M"))
df$stage_group <- factor(
  ifelse(grepl("^1", df$stage), "I", ifelse(grepl("^2", df$stage), "II", ifelse(grepl("^3", df$stage), "III", ifelse(grepl("^4", df$stage), "IV", NA)))),
  levels = c("I", "II", "III", "IV")
)

missingness <- data.frame(
  variable = c("OS_time_days", "OS_event", "age", "gender", "stage_group", genes),
  missing_n = c(sum(is.na(df$OS_time_days) | df$OS_time_days <= 0), sum(is.na(df$OS_event)), sum(is.na(df$age)), sum(is.na(df$gender)), sum(is.na(df$stage_group)), sapply(genes, function(gene) sum(is.na(df[[gene]])))),
  total_n = nrow(df)
)
missingness$missing_pct <- missingness$missing_n / missingness$total_n * 100
write.csv(missingness, file.path(output, "gse72094_missingness.csv"), row.names = FALSE)

eligible <- df[!is.na(df$OS_time_days) & df$OS_time_days > 0 & !is.na(df$OS_event) & complete.cases(df[, genes]), c("sample_id", "OS_time_days", "OS_event", "age", "gender", "stage_group", genes)]
for (gene in genes) eligible[[gene]] <- as.numeric(scale(eligible[[gene]]))
eligible$risk_score <- rowMeans(eligible[, c("CA9", "MKI67", "STAT1", "HIF1A")])
eligible$immune_score <- rowMeans(eligible[, c("CD4", "CD8A", "CD8B", "TIGIT")])
eligible$axis_score <- eligible$risk_score - eligible$immune_score
eligible$axis_group_median <- factor(ifelse(eligible$axis_score > median(eligible$axis_score), "high", "low"), levels = c("low", "high"))
write.csv(eligible, file.path(output, "gse72094_locked_score_analysis_dataset.csv"), row.names = FALSE)

extract_cox <- function(model, model_name) {
  summary_model <- summary(model)
  coefficients <- as.data.frame(summary_model$coefficients)
  intervals <- as.data.frame(summary_model$conf.int)
  terms <- rownames(coefficients)
  data.frame(
    model = model_name,
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
}

surv <- Surv(eligible$OS_time_days, eligible$OS_event)
models <- list(
  cox_axis_continuous = coxph(surv ~ axis_score, data = eligible),
  cox_axis_median = coxph(surv ~ axis_group_median, data = eligible),
  cox_risk_immune = coxph(surv ~ risk_score + immune_score, data = eligible)
)

adjusted <- eligible[complete.cases(eligible[, c("age", "gender", "stage_group")]), ]
models$cox_axis_adjusted_age_gender_stage <- coxph(Surv(OS_time_days, OS_event) ~ axis_score + age + gender + stage_group, data = adjusted)
clinical_model <- coxph(Surv(OS_time_days, OS_event) ~ age + gender + stage_group, data = adjusted)
clinical_plus_axis_model <- coxph(Surv(OS_time_days, OS_event) ~ age + gender + stage_group + axis_score, data = adjusted)
models$cox_clinical_baseline <- clinical_model
models$cox_clinical_plus_axis <- clinical_plus_axis_model

cox_results <- do.call(rbind, Map(extract_cox, models, names(models)))
write.csv(cox_results, file.path(output, "gse72094_cox_results.csv"), row.names = FALSE)

logrank <- survdiff(surv ~ axis_group_median, data = eligible)
logrank_results <- data.frame(
  comparison = "axis_group_median",
  chisq = logrank$chisq,
  df = length(logrank$n) - 1,
  p = pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
)
write.csv(logrank_results, file.path(output, "gse72094_logrank_results.csv"), row.names = FALSE)

km <- summary(survfit(surv ~ axis_group_median, data = eligible))
km_results <- data.frame(time = km$time, n.risk = km$n.risk, n.event = km$n.event, survival = km$surv, std.err = km$std.err, lower = km$lower, upper = km$upper, strata = km$strata)
write.csv(km_results, file.path(output, "gse72094_km_curve_data.csv"), row.names = FALSE)

set.seed(20260710)
bootstrap_cindex <- function(data, formula, b = 200L) {
  apparent <- summary(coxph(formula, data = data))$concordance[1]
  optimism <- numeric(b)
  for (i in seq_len(b)) {
    sample_index <- sample.int(nrow(data), replace = TRUE)
    boot <- data[sample_index, ]
    fit <- try(coxph(formula, data = boot), silent = TRUE)
    if (inherits(fit, "try-error")) {
      optimism[i] <- NA_real_
      next
    }
    apparent_boot <- summary(fit)$concordance[1]
    test_boot <- try(concordance(fit, newdata = data)$concordance, silent = TRUE)
    optimism[i] <- if (inherits(test_boot, "try-error")) NA_real_ else apparent_boot - test_boot
  }
  valid <- optimism[!is.na(optimism)]
  data.frame(
    apparent_cindex = apparent,
    optimism_mean = mean(valid),
    optimism_adjusted_cindex = apparent - mean(valid),
    bootstrap_runs = length(valid),
    stringsAsFactors = FALSE
  )
}

cindex_clinical <- bootstrap_cindex(adjusted, Surv(OS_time_days, OS_event) ~ age + gender + stage_group)
cindex_plus_axis <- bootstrap_cindex(adjusted, Surv(OS_time_days, OS_event) ~ age + gender + stage_group + axis_score)
cindex_results <- rbind(cbind(model = "clinical_baseline", cindex_clinical), cbind(model = "clinical_plus_axis", cindex_plus_axis))
write.csv(cindex_results, file.path(output, "gse72094_bootstrap_cindex.csv"), row.names = FALSE)

lr <- anova(clinical_model, clinical_plus_axis_model, test = "Chisq")
lr_results <- data.frame(
  comparison = "clinical_baseline_vs_clinical_plus_axis",
  df = lr$Df[2],
  chisq = lr$Chisq[2],
  p = lr$`Pr(>|Chi|)`[2]
)
write.csv(lr_results, file.path(output, "gse72094_likelihood_ratio_test.csv"), row.names = FALSE)

ph <- cox.zph(models$cox_axis_adjusted_age_gender_stage)
ph_results <- data.frame(term = rownames(ph$table), chisq = ph$table[, "chisq"], df = ph$table[, "df"], p = ph$table[, "p"], row.names = NULL)
write.csv(ph_results, file.path(output, "gse72094_ph_assumption.csv"), row.names = FALSE)

report <- c(
  "# GSE72094 Locked-Score External Validation: Evidence Results",
  "",
  "- Dataset: GSE72094 / GPL15048; processed series matrix and public GEO platform metadata.",
  "- Scope: pre-specified locked-score external OS analysis; no gene selection, score refitting, or manuscript modification.",
  paste0("- Eligible analysis set: ", nrow(eligible), " participants; ", sum(eligible$OS_event), " deaths."),
  paste0("- Adjusted clinical-model set: ", nrow(adjusted), " participants; ", sum(adjusted$OS_event), " deaths."),
  "",
  "## Generated Artifacts",
  "",
  "- `gse72094_missingness.csv`: eligibility and covariate missingness.",
  "- `gse72094_locked_score_analysis_dataset.csv`: locked-score analysis dataset.",
  "- `gse72094_cox_results.csv`, `gse72094_logrank_results.csv`, and `gse72094_km_curve_data.csv`: survival evidence tables.",
  "- `gse72094_bootstrap_cindex.csv` and `gse72094_likelihood_ratio_test.csv`: clinical incremental-value analyses.",
  "- `gse72094_ph_assumption.csv`: proportional-hazards diagnostic.",
  "",
  "## Interpretation Boundary",
  "",
  "- This file documents executed evidence only. Effect interpretation and manuscript claims are deferred to S4.",
  "- Time-dependent AUC was not calculated because the required package is absent and no package installation was authorized."
)
writeLines(report, file.path(output, "gse72094_evidence_results.md"))
cat("Wrote GSE72094 locked-score evidence package\n")
