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
