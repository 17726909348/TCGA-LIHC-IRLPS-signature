# ============================================================================
# Supplementary Code 5: Cut-off Validation and Sensitivity Analysis
# ============================================================================
# Project: TCGA-LIHC IRLPS Signature for HCC Prognosis
# Purpose: Validate optimal cut-off, perform sensitivity analysis,
#          conduct bootstrap internal validation
# ============================================================================
# Author: Ruili Zhou
# Date: 2024
# R version: 4.2.1
# Required packages: survival
# ============================================================================

# ============================================================================
# 1. Environment setup
# ============================================================================

rm(list = ls())
library(survival)
set.seed(12345)

# ============================================================================
# 2. Flexible path configuration
# ============================================================================
# USER: Modify the following line to point to your data folder
# The script will look for "riskScore.txt" generated from Code 3

input_file <- "riskScore.txt"

if (!file.exists(input_file)) {
    stop(paste("Error:", input_file, "not found.",
               "\nPlease run Code 3 first or set the correct working directory."))
}

cat("Working directory:", getwd(), "\n")
cat("Input file found. Starting analysis...\n")

# ============================================================================
# 3. Read data and define cut-off
# ============================================================================

rt <- read.table(input_file, header = TRUE, sep = "\t", 
                 check.names = FALSE, row.names = 1)

current_cutoff <- 0.880  # Youden index from 3-year ROC

cat("\n========== Data Summary ==========\n")
cat("Number of patients:", nrow(rt), "\n")
cat("Number of events (death):", sum(rt$fustat == 1), "\n")
cat("Risk score range:", round(range(rt$riskScore), 3), "\n")

# ============================================================================
# 4. Youden index cut-off validation
# ============================================================================

rt$group_youden <- ifelse(rt$riskScore > current_cutoff, 1, 0)
cox_youden <- coxph(Surv(futime, fustat) ~ group_youden, data = rt, ties = "efron")
sum_youden <- summary(cox_youden)

cat("\n========== Cut-off 0.880 (Youden index) ==========\n")
cat("High-risk group (n):", sum(rt$group_youden == 1), "\n")
cat("Low-risk group (n):", sum(rt$group_youden == 0), "\n")
cat("HR:", round(sum_youden$coefficients[1, "exp(coef)"], 2))
cat(", 95% CI:", round(sum_youden$conf.int[1, "lower .95"], 2), "-", 
    round(sum_youden$conf.int[1, "upper .95"], 2))
cat(", P =", round(sum_youden$coefficients[1, "Pr(>|z|)"], 4), "\n")

# ============================================================================
# 5. Sensitivity analysis across alternative cut-offs
# ============================================================================

percentiles <- c(0.30, 0.40, 0.50, 0.60, 0.70)
alt_cutoffs <- quantile(rt$riskScore, probs = percentiles)

cat("\n========== Sensitivity Analysis ==========\n")
cat("Method\t\t\tCutoff\tHR\t95% CI\t\tP\n")
cat("--------------------------------------------------------\n")

cutoff_names <- c("30% percentile", "40% percentile", "50% percentile", 
                  "60% percentile", "70% percentile", "Median")

for (i in 1:length(percentiles)) {
    cutoff_val <- alt_cutoffs[i]
    group_name <- paste0("group_p", percentiles[i]*100)
    rt[[group_name]] <- ifelse(rt$riskScore > cutoff_val, 1, 0)
    cox_obj <- coxph(Surv(futime, fustat) ~ rt[[group_name]], data = rt, ties = "efron")
    sum_obj <- summary(cox_obj)
    
    cat(cutoff_names[i], "\t", round(cutoff_val, 3), "\t",
        round(sum_obj$coefficients[1, "exp(coef)"], 2), "\t",
        round(sum_obj$conf.int[1, "lower .95"], 2), "-", 
        round(sum_obj$conf.int[1, "upper .95"], 2), "\t",
        round(sum_obj$coefficients[1, "Pr(>|z|)"], 4), "\n")
}

# Median cut-off
median_cutoff <- median(rt$riskScore)
rt$group_median <- ifelse(rt$riskScore > median_cutoff, 1, 0)
cox_median <- coxph(Surv(futime, fustat) ~ group_median, data = rt, ties = "efron")
sum_median <- summary(cox_median)
cat("Median\t\t\t", round(median_cutoff, 3), "\t",
    round(sum_median$coefficients[1, "exp(coef)"], 2), "\t",
    round(sum_median$conf.int[1, "lower .95"], 2), "-", 
    round(sum_median$conf.int[1, "upper .95"], 2), "\t",
    round(sum_median$coefficients[1, "Pr(>|z|)"], 4), "\n")

# ============================================================================
# 6. Bootstrap optimism correction (100 iterations)
# ============================================================================

n_bootstrap <- 100
n_samples <- nrow(rt)
boot_hr <- numeric(n_bootstrap)

cat("\n========== Bootstrap Validation (100 iterations) ==========\n")

for (i in 1:n_bootstrap) {
    if (i %% 20 == 0) cat("Completed", i, "of", n_bootstrap, "\n")
    
    boot_idx <- sample(1:n_samples, size = n_samples, replace = TRUE)
    boot_data <- rt[boot_idx, ]
    boot_data$boot_group <- ifelse(boot_data$riskScore > current_cutoff, 1, 0)
    
    if (sum(boot_data$boot_group) == 0 || sum(boot_data$boot_group) == n_samples) {
        boot_hr[i] <- NA
        next
    }
    
    boot_cox <- tryCatch({
        coxph(Surv(futime, fustat) ~ boot_group, data = boot_data, ties = "efron")
    }, error = function(e) NULL)
    
    if (!is.null(boot_cox)) {
        boot_hr[i] <- summary(boot_cox)$coefficients[1, "exp(coef)"]
    } else {
        boot_hr[i] <- NA
    }
}

boot_hr_clean <- boot_hr[!is.na(boot_hr)]
original_hr <- sum_youden$coefficients[1, "exp(coef)"]
optimism <- mean(boot_hr_clean) - original_hr
corrected_hr <- original_hr - optimism
hr_ci <- quantile(boot_hr_clean, probs = c(0.025, 0.975))

cat("\n========== Bootstrap Results ==========\n")
cat("Original HR:", round(original_hr, 2), "\n")
cat("Bootstrap mean HR:", round(mean(boot_hr_clean), 2), "\n")
cat("Optimism:", round(optimism, 2), "\n")
cat("Optimism-corrected HR:", round(corrected_hr, 2), "\n")
cat("95% CI:", round(hr_ci[1], 2), "-", round(hr_ci[2], 2), "\n")

# ============================================================================
# 7. Save results
# ============================================================================

bootstrap_results <- data.frame(
    Metric = c("Original_HR", "Bootstrap_mean_HR", "Optimism", 
               "Optimism_corrected_HR", "CI_lower", "CI_upper"),
    Value = c(original_hr, mean(boot_hr_clean), optimism, 
              corrected_hr, hr_ci[1], hr_ci[2])
)
write.table(bootstrap_results, file = "bootstrap_results.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResults saved to: bootstrap_results.txt\n")

# ============================================================================
# 8. Session information
# ============================================================================
sessionInfo()