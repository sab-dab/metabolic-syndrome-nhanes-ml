Machine Learning Prediction of Metabolic Syndrome Using Multi-Cycle NHANES Data

This repository contains the R analysis code supporting the manuscript:

Machine Learning Prediction of Metabolic Syndrome Using Multi-Cycle NHANES Data: Quantifying the Impact of a Diagnostic Component on Model Performance

Overview

This project evaluates machine-learning prediction of metabolic syndrome (MetS) using pooled NHANES data from five survey cycles spanning 2007–2008 through 2015–2016.

The analysis compares:

Logistic regression

XGBoost

Two predictor sets are evaluated:

Primary leakage-resistant model
Excludes waist circumference because waist circumference is itself one of the diagnostic components used to define metabolic syndrome.

Sensitivity model
Includes waist circumference to quantify how incorporation of a diagnostic component affects apparent model performance.

The workflow also includes:

Multi-cycle NHANES harmonization

Adult-only eligibility filtering

Metabolic syndrome outcome derivation

Descriptive summaries

Stratified five-fold cross-validation

Leakage-safe XGBoost early stopping

Temporal validation using 2015–2016 as a held-out later cycle

Calibration analysis

SHAP-based XGBoost interpretation

Repository Structure

metabolic-syndrome-nhanes-ml/
├── README.md
├── 01_NHANES_MetS_Full_Analysis.R
├── LICENSE
└── results_final/

The results_final/ directory is created automatically when the script is run.

Data Source

The analysis uses publicly available data from the National Health and Nutrition Examination Survey (NHANES).

NHANES data are available from the U.S. Centers for Disease Control and Prevention (CDC):

https://www.cdc.gov/nchs/nhanes/

The following survey cycles are used:

2007–2008

2009–2010

2011–2012

2013–2014

2015–2016

The raw NHANES .XPT files are not included in this repository.

Required NHANES Components

For each cycle, the analysis expects files containing the following components:

Demographics

Day 1 dietary totals

Body measures

Blood pressure

Fasting glucose

Triglycerides

HDL cholesterol

The analysis script searches for standard NHANES file prefixes such as:

DEMO

DR1TOT

BMX

BPX

GLU

TRIGLY

HDL

Expected Folder Structure

Download the required NHANES .XPT files and organize them as follows:

metabolic-syndrome-nhanes-ml/
├── 01_NHANES_MetS_Full_Analysis.R
├── NHANES_dataset_07-08/
├── NHANES_dataset_09-10/
├── NHANES_dataset_11-12/
├── NHANES_dataset_13-14/
├── NHANES_dataset_15-16/
└── results_final/

Place the corresponding NHANES .XPT files inside each cycle-specific folder.

Software Requirements

The analysis was written in R and requires the following packages:

foreign
dplyr
tidyr
stringr
caret
xgboost
pROC
ggplot2
Matrix
readr

Install the required packages once before running the analysis:

install.packages(c(
  "foreign",
  "dplyr",
  "tidyr",
  "stringr",
  "caret",
  "xgboost",
  "pROC",
  "ggplot2",
  "Matrix",
  "readr"
))

Running the Analysis

Clone or download this repository.

Download the required NHANES .XPT files from the CDC NHANES website.

Place the files into the cycle-specific folders shown above.

Set the repository root as the working directory in R.

Run:

source("01_NHANES_MetS_Full_Analysis.R")

The script uses:

PARENT_DIR <- getwd()

so the working directory should be the repository root.

Metabolic Syndrome Definition

Metabolic syndrome is defined as the presence of at least three of five criteria:

Elevated waist circumference

Elevated triglycerides

Elevated fasting glucose

Reduced HDL cholesterol

Elevated blood pressure

The script derives these criteria directly from the harmonized NHANES variables and then constructs a binary MetS outcome.

Predictor Sets

Primary model

The primary predictor set excludes waist circumference:

Sex

Age

Total energy intake

Protein intake

Fiber intake

Fiber per 1,000 kcal

Protein per 1,000 kcal

Sensitivity model

The sensitivity model includes all primary predictors plus:

Waist circumference

This comparison is intended to quantify the impact of including a diagnostic component of the outcome definition among the model predictors.

Validation Strategy

Five-fold cross-validation

The analysis uses stratified five-fold outer cross-validation.

For XGBoost, the number of boosting rounds is selected using an inner training/validation split only. The outer test fold is not used for early stopping.

Temporal validation

Models are trained using:

2007–2008

2009–2010

2011–2012

2013–2014

and evaluated on:

2015–2016

This provides a temporally separated validation assessment.

Main Outputs

The script generates files including:

NHANES_combined_adults_2007_2016.rds
NHANES_combined_adults_2007_2016.csv
Descriptive_Overall.csv
Descriptive_By_Cycle.csv
CV_Fold_Level_Results.csv
CV_Model_Comparison_Summary.csv
Temporal_Validation_Results.csv
Temporal_Calibration_Data_XGBoost_No_Waist.csv
Temporal_Calibration_Data_XGBoost_With_Waist.csv
Temporal_Calibration_XGBoost_No_Waist.png
Temporal_Calibration_XGBoost_With_Waist.png
Final_XGBoost_No_Waist.rds
Final_Logistic_No_Waist.rds
XGBoost_SHAP_Feature_Importance_Corrected.csv
XGBoost_SHAP_Feature_Importance_Corrected.png
XGBoost_SHAP_Beeswarm_Corrected.png
XGBoost_SHAP_Summary_Publication.png
OOF_Predictions_XGBoost_No_Waist.csv
OOF_ROC_XGBoost_No_Waist.png
sessionInfo.txt

Reproducibility

Random seeds are set within the script for cross-validation, XGBoost tuning, temporal validation, and SHAP plotting steps.

The script also saves sessionInfo() so that package and R version information can be recorded alongside the analysis outputs.

Interpretation

The primary purpose of the study is not to introduce a new machine-learning algorithm. Instead, it evaluates how predictor selection can influence apparent model performance when a predictor is also part of the clinical outcome definition.

The waist-inclusive analysis should therefore be interpreted as a diagnostic-component-inclusive classification analysis, whereas the waist-excluded analysis is intended to provide a more independent estimate of predictive performance.

Citation

If you use this code, please cite the associated manuscript:

Dabeer S, Palitanawala H. Machine Learning Prediction of Metabolic Syndrome Using Multi-Cycle NHANES Data: Quantifying the Impact of a Diagnostic Component on Model Performance.

Journal citation and DOI will be added after publication.

Authors

Sabira Dabeer
Arizona State University

Hatim Palitanawala

Contact

For questions regarding the manuscript or analysis, please contact the corresponding author through the contact information provided 
