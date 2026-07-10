#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input <- file.path(root, "processed")
output <- file.path(root, "results")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

genes <- c("CA9", "MKI67", "STAT1", "HIF1A", "CD4", "CD8A", "CD8B", "TIGIT")
expression <- read.csv(file.path(input, "gse31210_locked_score_expression.csv"), check.names = FALSE, na.strings = c("", "NA"))
clinical <- read.csv(file.path(input, "gse31210_clinical_raw.csv"), check.names = FALSE, na.strings = c("", "NA"))
df <- merge(clinical, expression, by = "sample_id", all = FALSE)

df$RFS_time_months <- suppressWarnings(as.numeric(df$months_before_relapse_censor))
df$RFS_event <- ifelse(df$relapse == "relapsed", 1L, ifelse(df$relapse == "not relapsed", 0L, NA_integer_))
df$eligible_tumor <- df$tissue == "primary lung tumor"
df$eligible_prognosis <- df$exclude_for_prognosis_analysis_due_to_incomplete_resection_or_adjuvant_therapy == "none"
df$age <- suppressWarnings(as.numeric(df$age_years))
df$gender <- factor(df$gender, levels = c("female", "male"))
df$stage_group <- factor(df$pathological_stage, levels = c("IA", "IB", "II"))

missingness <- data.frame(
  variable = c("eligible_tumor", "eligible_prognosis", "RFS_time_months", "RFS_event", "age", "gender", "stage_group", genes),
  missing_n = c(
    sum(!df$eligible_tumor, na.rm = TRUE),
    sum(!df$eligible_prognosis, na.rm = TRUE),
    sum(is.na(df$RFS_time_months) | df$RFS_time_months <= 0),
    sum(is.na(df$RFS_event)),
    sum(is.na(df$age)),
    sum(is.na(df$gender)),
    sum(is.na(df$stage_group)),
    sapply(genes, function(gene) sum(is.na(df[[gene]])))
  ),
  total_n = nrow(df)
)
missingness$missing_pct <- missingness$missing_n / missingness$total_n * 100
write.csv(missingness, file.path(output, "gse31210_missingness.csv"), row.names = FALSE)

eligible <- df[
  df$eligible_tumor & df$eligible_prognosis &
    !is.na(df$RFS_time_months) & df$RFS_time_months > 0 &
    !is.na(df$RFS_event) & complete.cases(df[, genes]),
  c("sample_id", "RFS_time_months", "RFS_event", "age", "gender", "stage_group", genes)
]
for (gene in genes) eligible[[gene]] <- as.numeric(scale(eligible[[gene]]))
eligible$risk_score <- rowMeans(eligible[, c("CA9", "MKI67", "STAT1", "HIF1A")])
eligible$immune_score <- rowMeans(eligible[, c("CD4", "CD8A", "CD8B", "TIGIT")])
eligible$axis_score <- eligible$risk_score - eligible$immune_score
eligible$axis_group_median <- factor(ifelse(eligible$axis_score > median(eligible$axis_score), "high", "low"), levels = c("low", "high"))
write.csv(eligible, file.path(output, "gse31210_locked_score_analysis_dataset.csv"), row.names = FALSE)

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

surv <- Surv(eligible$RFS_time_months, eligible$RFS_event)
models <- list(
  cox_axis_continuous = coxph(surv ~ axis_score, data = eligible),
  cox_axis_median = coxph(surv ~ axis_group_median, data = eligible),
  cox_risk_immune = coxph(surv ~ risk_score + immune_score, data = eligible)
)

adjusted <- eligible[complete.cases(eligible[, c("age", "gender", "stage_group")]), ]
models$cox_axis_adjusted_age_gender_stage <- coxph(Surv(RFS_time_months, RFS_event) ~ axis_score + age + gender + stage_group, data = adjusted)
clinical_model <- coxph(Surv(RFS_time_months, RFS_event) ~ age + gender + stage_group, data = adjusted)
clinical_plus_axis_model <- coxph(Surv(RFS_time_months, RFS_event) ~ age + gender + stage_group + axis_score, data = adjusted)
models$cox_clinical_baseline <- clinical_model
models$cox_clinical_plus_axis <- clinical_plus_axis_model

cox_results <- do.call(rbind, Map(extract_cox, models, names(models)))
write.csv(cox_results, file.path(output, "gse31210_cox_results.csv"), row.names = FALSE)

logrank <- survdiff(surv ~ axis_group_median, data = eligible)
logrank_results <- data.frame(
  comparison = "axis_group_median",
  chisq = logrank$chisq,
  df = length(logrank$n) - 1,
  p = pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
)
write.csv(logrank_results, file.path(output, "gse31210_logrank_results.csv"), row.names = FALSE)

km <- summary(survfit(surv ~ axis_group_median, data = eligible))
km_results <- data.frame(time = km$time, n.risk = km$n.risk, n.event = km$n.event, survival = km$surv, std.err = km$std.err, lower = km$lower, upper = km$upper, strata = km$strata)
write.csv(km_results, file.path(output, "gse31210_km_curve_data.csv"), row.names = FALSE)

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

cindex_clinical <- bootstrap_cindex(adjusted, Surv(RFS_time_months, RFS_event) ~ age + gender + stage_group)
cindex_plus_axis <- bootstrap_cindex(adjusted, Surv(RFS_time_months, RFS_event) ~ age + gender + stage_group + axis_score)
cindex_results <- rbind(cbind(model = "clinical_baseline", cindex_clinical), cbind(model = "clinical_plus_axis", cindex_plus_axis))
write.csv(cindex_results, file.path(output, "gse31210_bootstrap_cindex.csv"), row.names = FALSE)

lr <- anova(clinical_model, clinical_plus_axis_model, test = "Chisq")
lr_results <- data.frame(
  comparison = "clinical_baseline_vs_clinical_plus_axis",
  df = lr$Df[2],
  chisq = lr$Chisq[2],
  p = lr$`Pr(>|Chi|)`[2]
)
write.csv(lr_results, file.path(output, "gse31210_likelihood_ratio_test.csv"), row.names = FALSE)

ph <- cox.zph(models$cox_axis_adjusted_age_gender_stage)
ph_results <- data.frame(term = rownames(ph$table), chisq = ph$table[, "chisq"], df = ph$table[, "df"], p = ph$table[, "p"], row.names = NULL)
write.csv(ph_results, file.path(output, "gse31210_ph_assumption.csv"), row.names = FALSE)

report <- c(
  "# GSE31210 Locked-Score External Validation: Evidence Results",
  "",
  "- Dataset: GSE31210 / GPL570; public GEO series matrix and platform metadata.",
  "- Endpoint: relapse-free survival in months; event defined as the supplied `relapsed` status.",
  "- Eligibility: primary lung tumors only; cases marked for prognosis exclusion because of incomplete resection or adjuvant therapy were excluded.",
  "- Scope: pre-specified locked-score secondary external RFS analysis; no gene selection, score refitting, or manuscript modification.",
  paste0("- Eligible analysis set: ", nrow(eligible), " participants; ", sum(eligible$RFS_event), " relapses."),
  paste0("- Adjusted clinical-model set: ", nrow(adjusted), " participants; ", sum(adjusted$RFS_event), " relapses."),
  "",
  "## Generated Artifacts",
  "",
  "- `gse31210_missingness.csv`: eligibility and covariate missingness.",
  "- `gse31210_locked_score_analysis_dataset.csv`: locked-score analysis dataset.",
  "- `gse31210_cox_results.csv`, `gse31210_logrank_results.csv`, and `gse31210_km_curve_data.csv`: RFS evidence tables.",
  "- `gse31210_bootstrap_cindex.csv` and `gse31210_likelihood_ratio_test.csv`: clinical incremental-value analyses.",
  "- `gse31210_ph_assumption.csv`: proportional-hazards diagnostic.",
  "",
  "## Interpretation Boundary",
  "",
  "- This file documents executed evidence only. Effect interpretation and manuscript claims are deferred to S4.",
  "- Time-dependent AUC was not calculated because the required package is absent and no package installation was authorized."
)
writeLines(report, file.path(output, "gse31210_evidence_results.md"))
cat("Wrote GSE31210 locked-score evidence package\n")
