# ============================================================================
# Supplementary Code 3: LASSO and Multivariate Cox Regression
# ============================================================================
# Project: TCGA-LIHC IRLPS Signature for HCC Prognosis
# Purpose: Perform LASSO regression, build multivariate Cox model,
#          calculate risk scores, generate forest plots
# ============================================================================
# Author: Ruili Zhou
# Date: 2024
# R version: 4.2.1
# Required packages: survival, survminer, glmnet
# ============================================================================

# ============================================================================
# 1. Environment setup
# ============================================================================

rm(list = ls())
library(survival)
library(survminer)
library(glmnet)
set.seed(12345)

# ============================================================================
# 2. Flexible path configuration
# ============================================================================
# USER: Modify the following line to point to your data folder
# The script will look for "uniSigExp.txt" generated from Code 2

input_file <- "uniSigExp.txt"

if (!file.exists(input_file)) {
    stop(paste("Error:", input_file, "not found.",
               "\nPlease run Code 2 first or set the correct working directory."))
}

cat("Working directory:", getwd(), "\n")
cat("Input file found. Starting analysis...\n")

# ============================================================================
# 3. Read input data
# ============================================================================

rt <- read.table(input_file, header = TRUE, sep = "\t", 
                 check.names = FALSE, row.names = 1)

# ============================================================================
# 4. LASSO regression
# ============================================================================
# Note: glmnet default standardize = TRUE applies Z-score transformation

# Prepare data matrix
x <- as.matrix(rt[, c(3:ncol(rt))])
y <- data.matrix(Surv(rt$futime, rt$fustat))

# Fit LASSO path
fit <- glmnet(x, y, family = "cox", maxit = 1000)

# Plot LASSO path
pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

# Cross-validated LASSO
cvfit <- cv.glmnet(x, y, family = "cox", maxit = 1000)

# Plot cross-validation results
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")
dev.off()

# Extract coefficients at lambda.min
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
lassoGene <- c("futime", "fustat", lassoGene)

# Output LASSO-selected data
lassoSigExp <- rt[, lassoGene]
lassoSigExp <- cbind(id = row.names(lassoSigExp), lassoSigExp)
write.table(lassoSigExp, file = "lasso.SigExp.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)

cat("LASSO selected", length(lassoGene) - 2, "genes\n")

# ============================================================================
# 5. Multivariate Cox regression
# ============================================================================

rt <- read.table("lasso.SigExp.txt", header = TRUE, sep = "\t", 
                 check.names = FALSE, row.names = 1)

multiCox <- coxph(Surv(futime, fustat) ~ ., data = rt, ties = "efron")
multiCox <- step(multiCox, direction = "both")
multiCoxSum <- summary(multiCox)

# ============================================================================
# 6. Proportional hazards assumption test
# ============================================================================

ph_test <- cox.zph(multiCox)
print(ph_test)
cat("Global PH test P =", ph_test$global[3], "\n")

# ============================================================================
# 7. Output model parameters
# ============================================================================

outTab <- data.frame()
outTab <- cbind(
    coef = multiCoxSum$coefficients[, "coef"],
    HR = multiCoxSum$conf.int[, "exp(coef)"],
    HR.95L = multiCoxSum$conf.int[, "lower .95"],
    HR.95H = multiCoxSum$conf.int[, "upper .95"],
    pvalue = multiCoxSum$coefficients[, "Pr(>|z|)"]
)
outTab <- cbind(id = row.names(outTab), outTab)
outTab <- gsub("`", "", outTab)
write.table(outTab, file = "multi.Cox.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================================
# 8. Calculate risk scores
# ============================================================================

riskScore <- predict(multiCox, type = "risk", newdata = rt)
coxGene <- rownames(multiCoxSum$coefficients)
coxGene <- gsub("`", "", coxGene)
outCol <- c("futime", "fustat", coxGene)
riskOut <- cbind(rt[, outCol], riskScore)
riskOut <- cbind(id = rownames(riskOut), riskOut)
write.table(riskOut, file = "riskScore.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Risk scores saved to: riskScore.txt\n")

# ============================================================================
# 9. Generate forest plots (function definition)
# ============================================================================

bioForest <- function(coxFile = NULL, forestFile = NULL, forestCol = NULL) {
    
    rt <- read.table(coxFile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    gene <- rownames(rt)
    hr <- sprintf("%.3f", rt$"HR")
    hrLow <- sprintf("%.3f", rt$"HR.95L")
    hrHigh <- sprintf("%.3f", rt$"HR.95H")
    Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
    pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
    
    pdf(file = forestFile, width = 9, height = 7)
    n <- nrow(rt)
    nRow <- n + 1
    ylim <- c(1, nRow)
    layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
    
    xlim <- c(0, 3)
    par(mar = c(4, 2.5, 2, 1))
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
    text.cex <- 0.8
    text(0, n:1, gene, adj = 0, cex = text.cex)
    text(2.08 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
    text(2.08 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)
    text(3.12, n:1, Hazard.ratio, adj = 1, cex = text.cex)
    text(3.12, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
    
    par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
    xlim <- c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ylab = "", 
         xaxs = "i", xlab = "Hazard ratio")
    arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, 
           code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
    abline(v = 1, col = "black", lty = 2, lwd = 2)
    boxcolor <- ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.6)
    axis(1)
    dev.off()
}

# Generate forest plots
bioForest(coxFile = "multi.Cox.txt", forestFile = "model.multiForest.pdf", 
          forestCol = c("red", "green"))

uniRT <- read.table("uniCox.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
uniRT <- uniRT[coxGene, ]
uniRT <- cbind(id = row.names(uniRT), uniRT)
write.table(uniRT, file = "unicox.forest.txt", sep = "\t", row.names = FALSE, quote = FALSE)
bioForest(coxFile = "unicox.forest.txt", forestFile = "model.uniForest.pdf", 
          forestCol = c("red", "green"))

cat("Forest plots saved to: model.multiForest.pdf, model.uniForest.pdf\n")

# ============================================================================
# 10. Session information
# ============================================================================
sessionInfo()