# ============================================================================
# Supplementary Code 1: Data Preprocessing
# ============================================================================
# Project: TCGA-LIHC IRLPS Signature for HCC Prognosis
# Manuscript: Observational TCGA-LIHC analysis: construction of an 
#            immune-related long non-coding RNA pair signature
# Purpose: Merge lncRNA pair expression data with clinical survival data,
#          remove normal samples, match samples, output merged data
# ============================================================================
# Author: Ruili Zhou
# Date: 2024
# R version: 4.2.1
# Required packages: limma
# ============================================================================

# ============================================================================
# 1. Environment setup with flexible path handling
# ============================================================================

# Clear environment
rm(list = ls())

# Load required packages
library(limma)

# Set random seed for reproducibility
set.seed(12345)

# ============================================================================
# 2. Flexible path configuration
# ============================================================================
# USER: Choose ONE of the following methods:

# Method 1 (Recommended): Set the project root directory
# Example: setwd("C:/Users/YourName/YourProjectFolder")
# Then place input files in this folder

# Method 2: Use RStudio project (if using .Rproj file)
# This automatically sets working directory to the project root

# Method 3: Specify file paths directly
# Modify the lines below to point to your data files

# ----------------------------------------------------------------------------
# Please modify the following variables to match your local setup:
# ----------------------------------------------------------------------------

# Option A: Set working directory to where your input files are located
# setwd("PATH_TO_YOUR_DATA_FOLDER")

# Option B: Specify full paths to input files (uncomment and modify)
# pairFile <- "PATH_TO_YOUR_INPUT_FILE/lncrnaPair.txt"
# cliFile <- "PATH_TO_YOUR_INPUT_FILE/time.txt"

# ----------------------------------------------------------------------------
# Default settings (assumes input files are in current working directory)
# ----------------------------------------------------------------------------

pairFile <- "lncrnaPair.txt"     # Immune lncRNA pair expression file
cliFile <- "time.txt"            # Clinical survival data file

# Check if input files exist
if (!file.exists(pairFile)) {
    stop(paste("Error: Input file not found:", pairFile, 
               "\nPlease set the correct working directory or file path."))
}

if (!file.exists(cliFile)) {
    stop(paste("Error: Input file not found:", cliFile,
               "\nPlease set the correct working directory or file path."))
}

cat("Working directory:", getwd(), "\n")
cat("Input files found. Starting analysis...\n")

# ============================================================================
# 3. Read lncRNA pair expression data
# ============================================================================

rt <- read.table(pairFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]

# Convert to numeric matrix
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), 
               nrow = nrow(exp), 
               dimnames = dimnames)
data <- avereps(data)

# ============================================================================
# 4. Remove normal samples (keep only tumor samples)
# ============================================================================
# TCGA sample barcode: 
#   - 01 = tumor (primary solid tumor)
#   - 11 = normal (solid tissue normal)

group <- sapply(strsplit(colnames(data), "\\-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)  # Reclassify
data <- data[, group == 0]       # Keep only tumor samples

# Simplify column names to patient IDs
colnames(data) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", 
                       "\\1\\-\\2\\-\\3", colnames(data))
data <- t(data)
data <- avereps(data)

# ============================================================================
# 5. Read clinical survival data
# ============================================================================

cli <- read.table(cliFile, sep = "\t", check.names = FALSE, 
                  header = TRUE, row.names = 1)

# ============================================================================
# 6. Match samples between expression and clinical data
# ============================================================================

sameSample <- intersect(row.names(data), row.names(cli))
data <- data[sameSample, ]
cli <- cli[sameSample, ]

# ============================================================================
# 7. Merge and output
# ============================================================================

out <- cbind(cli, data)
out <- cbind(id = row.names(out), out)
write.table(out, file = "pairTime.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("Output saved to: pairTime.txt\n")

# ============================================================================
# 8. Session information for reproducibility
# ============================================================================
sessionInfo()