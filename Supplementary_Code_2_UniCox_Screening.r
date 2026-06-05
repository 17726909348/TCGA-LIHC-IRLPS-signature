# ============================================================================
# Supplementary Code 2: Univariate Cox Regression Screening
# ============================================================================
# Project: TCGA-LIHC IRLPS Signature for HCC Prognosis
# Purpose: Perform univariate Cox regression for each lncRNA pair,
#          screen for prognostic pairs with P < 0.01,
#          output screened data for LASSO regression
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
# The script will look for "pairTime.txt" generated from Code 1

# Option A: Set working directory (uncomment and modify)
# setwd("PATH_TO_YOUR_DATA_FOLDER")

# Option B: Let the script use current working directory
# Make sure "pairTime.txt" is in the current directory

input_file <- "pairTime.txt"

# Check if input file exists
if (!file.exists(input_file)) {
    stop(paste("Error:", input_file, "not found.",
               "\nPlease run Code 1 first or set the correct working directory."))
}

cat("Working directory:", getwd(), "\n")
cat("Input file found. Starting analysis...\n")

# ============================================================================
# 3. Read input data
# ============================================================================

rt <- read.table(input_file, header = TRUE, sep = "\t", 
                 check.names = FALSE, row.names = 1)

# Convert survival time from days to years
rt$futime <- rt$futime / 365

# ============================================================================
# 4. Parameters
# ============================================================================

pFilter <- 0.01    # Significance threshold for screening

# ============================================================================
# 5. Univariate Cox regression analysis
# ============================================================================

outTab <- data.frame()
sigGenes <- c("futime", "fustat")

cat("Starting univariate Cox regression analysis...\n")
cat("Total genes to analyze:", ncol(rt) - 2, "\n")

for (gene in colnames(rt[, 3:ncol(rt)])) {
    cox <- coxph(Surv(futime, fustat) ~ rt[, gene], data = rt, ties = "efron")
    coxSummary <- summary(cox)
    coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
    
    if (coxP < pFilter) {
        sigGenes <- c(sigGenes, gene)
        outTab <- rbind(outTab,
                        cbind(gene = gene,
                              HR = coxSummary$conf.int[, "exp(coef)"],
                              HR.95L = coxSummary$conf.int[, "lower .95"],
                              HR.95H = coxSummary$conf.int[, "upper .95"],
                              pvalue = coxP))
    }
}

cat("Significant genes found:", nrow(outTab), "\n")

# ============================================================================
# 6. Output results
# ============================================================================

write.table(outTab, file = "uniCox.txt", sep = "\t", row.names = FALSE, quote = FALSE)

surSigExp <- rt[, sigGenes]
surSigExp <- cbind(id = row.names(surSigExp), surSigExp)
write.table(surSigExp, file = "uniSigExp.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("Output saved to: uniCox.txt, uniSigExp.txt\n")

# ============================================================================
# 7. Session information
# ============================================================================
sessionInfo()