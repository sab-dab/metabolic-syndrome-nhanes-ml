# Machine Learning Prediction of Metabolic Syndrome Using Multi-Cycle NHANES Data

This repository contains the R analysis code supporting the manuscript:

**Machine Learning Prediction of Metabolic Syndrome Using Multi-Cycle NHANES Data: Quantifying the Impact of a Diagnostic Component on Model Performance**

## Overview

This project evaluates machine-learning prediction of metabolic syndrome (MetS) using pooled NHANES data from five survey cycles spanning 2007–2008 through 2015–2016.

The analysis compares:

- Logistic regression
- XGBoost

Two predictor sets are evaluated:

1. **Primary leakage-resistant model**  
   Excludes waist circumference because waist circumference is itself one of the diagnostic components used to define metabolic syndrome.

2. **Sensitivity model**  
   Includes waist circumference to quantify how incorporation of a diagnostic component affects apparent model performance.

The workflow also includes:

- Multi-cycle NHANES harmonization
- Adult-only eligibility filtering
- Metabolic syndrome outcome derivation
- Descriptive summaries
- Stratified five-fold cross-validation
- Leakage-safe XGBoost early stopping
- Temporal validation using 2015–2016 as a held-out later cycle
- Calibration analysis
- SHAP-based XGBoost interpretation

## Repository Structure

```text
metabolic-syndrome-nhanes-ml/
├── README.md
├── 01_NHANES_MetS_Full_Analysis.R
└── results_final/
```

The `results_final/` directory is created automatically when the analysis script is run.

## Data Source

This analysis uses publicly available data from the **National Health and Nutrition Examination Survey (NHANES)** conducted by the U.S. Centers for Disease Control and Prevention.

NHANES data are available at:

https://www.cdc.gov/nchs/nhanes/

The following survey cycles are included:

- 2007–2008
- 2009–2010
- 2011–2012
- 2013–2014
- 2015–2016

The raw NHANES `.XPT` files are **not included in this repository**.

## Required NHANES Components

For each survey cycle, the analysis uses files containing the following components:

- Demographics
- Day 1 dietary totals
- Body measures
- Blood pressure
- Fasting glucose
- Triglycerides
- HDL cholesterol

The script searches for standard NHANES file prefixes including:

- `DEMO`
- `DR1TOT`
- `BMX`
- `BPX`
- `GLU`
- `TRIGLY`
- `HDL`

## Expected Folder Structure

Download the required NHANES `.XPT` files and organize them as follows:

```text
metabolic-syndrome-nhanes-ml/
├── 01_NHANES_MetS_Full_Analysis.R
├── NHANES_dataset_07-08/
├── NHANES_dataset_09-10/
├── NHANES_dataset_11-12/
├── NHANES_dataset_13-14/
├── NHANES_dataset_15-16/
└── results_final/
```

Place the corresponding NHANES `.XPT` files inside each cycle-specific folder.

## Software Requirements

The analysis was written in R and requires the following packages:

```text
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
```

Install the required packages before running the analysis:

```r
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
```

## Running the Analysis

1. Clone or download this repository.
2. Download the required NHANES `.XPT` files from the CDC NHANES website.
3. Place the files into the cycle-specific folders shown above.
4. Set the repository root as the working directory in R.
5. Run:

```r
source("01_NHANES_MetS_Full_Analysis.R")
```

The script uses:

```r
PARENT_DIR <- getwd()
```

Therefore, the working directory should be the repository root.

## Metabolic Syndrome Definition

Metabolic syndrome is defined as the presence of at least **three of five criteria**:

- Elevated waist circumference
- Elevated triglycerides
- Elevated fasting glucose
- Reduced HDL cholesterol
- Elevated blood pressure

The analysis script derives these criteria from the harmonized NHANES variables and constructs a binary `MetS` outcome.

## Predictor Sets

### Primary Model

The primary predictor set excludes waist circumference and includes:

- Sex
- Age
- Total energy intake
- Protein intake
- Fiber intake
- Fiber per 1,000 kcal
- Protein per 1,000 kcal

Waist circumference is excluded because it is itself part of the metabolic syndrome diagnostic definition.

### Sensitivity Model

The sensitivity model includes all primary predictors plus:

- Waist circumference

This comparison evaluates how inclusion of a diagnostic component affects apparent model performance.

## Validation Strategy

### Five-Fold Cross-Validation

The analysis uses stratified five-fold outer cross-validation.

For XGBoost, the number of boosting rounds is selected using an **inner training/validation split only**. The outer test fold is not used for early stopping.

### Temporal Validation

Models are trained using NHANES cycles:

- 2007–2008
- 2009–2010
- 2011–2012
- 2013–2014

and evaluated on:

- 2015–2016

This provides a temporally separated validation assessment.

## Main Outputs

The script generates files including:

```text
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
```

## Reproducibility

Random seeds are set within the analysis script for:

- Cross-validation
- XGBoost tuning
- Temporal validation
- SHAP plotting

The script also saves `sessionInfo()` so that the R version and package environment used for the analysis can be documented.

## Interpretation

The purpose of this study is not to introduce a new machine-learning algorithm.

Instead, the analysis evaluates how predictor selection can influence apparent model performance when a predictor is also part of the clinical outcome definition.

The waist-inclusive analysis should therefore be interpreted as a **diagnostic-component-inclusive classification analysis**, while the waist-excluded analysis provides a more independent estimate of predictive performance.

## Citation

If you use this code, please cite the associated manuscript:

**Dabeer S, Palitanawala H. Machine Learning Prediction of Metabolic Syndrome Using Multi-Cycle NHANES Data: Quantifying the Impact of a Diagnostic Component on Model Performance.**

Journal citation and DOI will be added after publication.

## Authors

**Sabira Dabeer**  
Arizona State University

**Hatim Palitanawala**

## Contact

For questions regarding the manuscript or analysis, please contact the corresponding author using the contact information provided in the manuscript.
