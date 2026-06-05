# ============================================================================
# Supplementary Code 4: ROC Curves and Optimal Cut-off
# ============================================================================
# Project: TCGA-LIHC IRLPS Signature for HCC Prognosis
# Purpose: Plot time-dependent ROC curves, determine optimal cut-off,
#          classify patients into risk groups
# ============================================================================
# Author: Ruili Zhou
# Date: 2024
# R version: 4.2.1
# Required packages: survivalROC
# ============================================================================

# ============================================================================
# 1. Environment setup
# ============================================================================

rm(list = ls())
library(survivalROC)
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
# 3. Read risk score file
# ============================================================================

rt <- read.table(input_file, header = TRUE, sep = "\t", 
                 check.names = FALSE, row.names = 1)

cat("Number of patients:", nrow(rt), "\n")
cat("Number of events:", sum(rt$fustat), "\n")

# ============================================================================
# 4. Determine optimal cut-off (Youden index)
# ============================================================================

predictTime <- 1
roc <- survivalROC(Stime = rt$futime, 
                   status = rt$fustat, 
                   marker = rt$riskScore, 
                   predict.time = predictTime, 
                   method = "KM")

youden_index <- roc$TP - roc$FP
optimal_cutoff <- roc$cut.values[which.max(youden_index)]

cat("\nOptimal cut-off value (Youden index):", optimal_cutoff, "\n")

# ============================================================================
# 5. Plot ROC curve with optimal cut-off point
# ============================================================================

pdf(file = "ROC.cutoff.pdf", width = 5.5, height = 5.5)
plot(roc$FP, roc$TP, type = "l", xlim = c(0, 1), ylim = c(0, 1), 
     col = "black", xlab = "False positive rate", ylab = "True positive rate",
     lwd = 2, cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
polygon(x = c(0, roc$FP, 1, 0), y = c(0, roc$TP, 1, 0), 
        col = "#24B35D", border = NA)
segments(0, 0, 1, 1, lty = 2)
points(roc$FP[which.max(youden_index)], roc$TP[which.max(youden_index)], 
       pch = 20, col = "red", cex = 1.5)
text(0.85, 0.1, paste0("AUC = ", sprintf("%.3f", roc$AUC)), cex = 1.2)
text(roc$FP[which.max(youden_index)] + 0.15, 
     roc$TP[which.max(youden_index)] - 0.05, 
     paste0("Cutoff = ", sprintf("%0.3f", optimal_cutoff)))
dev.off()

# ============================================================================
# 6. Classify patients into risk groups
# ============================================================================

risk <- as.vector(ifelse(rt$riskScore > optimal_cutoff, "high", "low"))
outTab <- cbind(rt, risk)
write.table(cbind(id = rownames(outTab), outTab), 
            file = "risk.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nRisk group distribution:\n")
cat("High-risk group:", sum(risk == "high"), "\n")
cat("Low-risk group:", sum(risk == "low"), "\n")

# ============================================================================
# 7. Plot multi-year ROC curves
# ============================================================================

rocCol <- c("red", "green", "blue")
aucText <- c()

pdf(file = "ROC.multiTime.pdf", width = 6, height = 6)
par(oma = c(0.5, 1, 0, 1), font.lab = 1.5, font.axis = 1.5)

for (t in 1:3) {
    roc <- survivalROC(Stime = rt$futime, status = rt$fustat, 
                       marker = rt$riskScore, predict.time = t, method = "KM")
    if (t == 1) {
        plot(roc$FP, roc$TP, type = "l", xlim = c(0, 1), ylim = c(0, 1), 
             col = rocCol[t], xlab = "False positive rate", ylab = "True positive rate",
             lwd = 2)
        abline(0, 1)
    } else {
        lines(roc$FP, roc$TP, col = rocCol[t], lwd = 2)
    }
    aucText <- c(aucText, paste0(t, "-year (AUC = ", sprintf("%.3f", roc$AUC), ")"))
}

legend("bottomright", aucText, lwd = 2, bty = "n", col = rocCol)
dev.off()

cat("\nOutput saved to: ROC.cutoff.pdf, ROC.multiTime.pdf, risk.txt\n")

# ============================================================================
# 8. Session information
# ============================================================================
sessionInfo()