markdown

# Immune-related lncRNA Pair Signature for HCC Prognosis

[![R version](https://img.shields.io/badge/R-4.2.1-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![TCGA](https://img.shields.io/badge/Data-TCGA--LIHC-orange.svg)](https://portal.gdc.cancer.gov/projects/TCGA-LIHC)

## 📖 Overview

This repository contains the complete analysis code for the manuscript:

**"Observational TCGA-LIHC analysis: construction of an immune-related long non-coding RNA pair signature for predicting prognosis and therapeutic response in hepatocellular carcinoma"**

*Running title: lncRNA signature for HCC prognosis*

**Authors**: Ruili Zhou, Fanglan Yao, Yuxi Shuai, Kuo Qi, Yan Teng

## 📊 Study Design

This study develops and validates an immune-related lncRNA pair signature (IRLPS) for predicting prognosis and therapeutic response in hepatocellular carcinoma (HCC) using the TCGA-LIHC dataset.

The key innovation is the **within-sample ranking strategy** for lncRNA expression, which effectively mitigates batch effects across different sequencing platforms and detection times.

## 🗂️ Repository Structure

TCGA-LIHC-IRLPS-signature/
│
├── Code/
│ ├── 01_data_preprocessing.R # Data preprocessing and sample matching
│ ├── 02_univariate_Cox.R # Univariate Cox regression screening
│ ├── 03_LASSO_multivariate_Cox.R # LASSO and multivariate Cox model
│ ├── 04_ROC_analysis.R # ROC curves and cut-off determination
│ └── 05_validation_sensitivity.R # Bootstrap validation and sensitivity analysis
│
├── README.md # This file
├── .gitignore # Git ignore file
└── LICENSE # MIT License
text


## 🔧 Requirements

### Software
- **R version 4.2.1** or higher

### Required R Packages

```r
# Install all required packages
packages <- c("limma", "survival", "survminer", "glmnet", "survivalROC")
install.packages(packages)

Package	Version	Purpose
limma	≥ 3.52.0	Data preprocessing and normalization
survival	≥ 3.5.0	Cox regression analysis
survminer	≥ 0.4.9	Kaplan-Meier curves (optional)
glmnet	≥ 4.1-4	LASSO regression
survivalROC	≥ 1.0.3	Time-dependent ROC curves
📥 Data Availability

The datasets used in this study are publicly available from:
Data source	Link
TCGA-LIHC	https://portal.gdc.cancer.gov/projects/TCGA-LIHC
ImmPort (immune genes)	https://www.immport.org/

Input files required:

    lncrnaPair.txt - Immune lncRNA pair expression matrix

    time.txt - Clinical survival data (time and status)

🚀 Usage Instructions
Step 1: Prepare Data

    Download TCGA-LIHC data from GDC Portal

    Obtain immune-related lncRNA pairs (see manuscript methods)

    Prepare lncrnaPair.txt and time.txt files

Step 2: Set Working Directory

Option A (RStudio):
text

Session → Set Working Directory → Choose Directory → Select your data folder

Option B (R command):
r

setwd("C:/Users/YourName/YourProjectFolder")

Step 3: Run Analysis

Run the scripts in numerical order:
r

# Run sequentially
source("Code/01_data_preprocessing.R")
source("Code/02_univariate_Cox.R")
source("Code/03_LASSO_multivariate_Cox.R")
source("Code/04_ROC_analysis.R")
source("Code/05_validation_sensitivity.R")

Or run each script individually in R console or RStudio.
📊 Analysis Workflow
text

Raw RNA-seq data (TCGA-LIHC)
        ↓
    Data preprocessing (01)
        ↓
    Immune lncRNA identification
        ↓
    lncRNA pair construction (within-sample ranking)
        ↓
    Univariate Cox screening (P < 0.01) (02)
        ↓
    LASSO regression (1000 iterations) (03)
        ↓
    Multivariate Cox model (12 IRLPs)
        ↓
    Risk score calculation
        ↓
    ROC analysis & optimal cut-off (04)
        ↓
    Bootstrap validation & sensitivity analysis (05)

📈 Key Results
Performance Metric	Value (95% CI)
1-year AUC	0.869 (0.825-0.913)
3-year AUC	0.830 (0.781-0.879)
5-year AUC	0.822 (0.771-0.873)
Optimism-corrected HR	5.62 (4.23-9.44)
📝 Reproducibility
Parameter	Setting
Random seed	12345 (fixed for all analyses)
LASSO standardization	Z-score (mean=0, SD=1) via glmnet(standardize=TRUE)
Bootstrap iterations	100 for internal validation
Ties handling	Efron method for Cox regression
Significance threshold	P < 0.01 for univariate screening
📄 Citation

If you use this code in your research, please cite our manuscript:

    Zhou R, Yao F, Shuai Y, Qi K, Teng Y. Observational TCGA-LIHC analysis: construction of an immune-related long non-coding RNA pair signature for predicting prognosis and therapeutic response in hepatocellular carcinoma. [Journal], [Year].

📧 Contact

For questions about the code or analysis, please contact:

    Kuo Qi: ldyy_qik@lzu.edu.cn

    Yan Teng: ldyy_qik@lzu.edu.cn

📜 License

This code is released under the MIT License. See LICENSE file for details.
⚠️ Disclaimer

This code is provided for research purposes only. Users should independently verify the results and ensure compliance with TCGA data usage policies. The authors assume no responsibility for any consequences arising from the use of this code.
