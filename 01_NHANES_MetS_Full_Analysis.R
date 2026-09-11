# ============================================================
# Machine Learning Prediction of Metabolic Syndrome
# Using Multi-Cycle NHANES Data
#
# Reproducible analysis pipeline for:
#   1) Multi-cycle NHANES harmonization
#   2) Adult-only eligibility
#   3) Metabolic syndrome (MetS) outcome derivation
#   4) Descriptive summaries
#   5) Leakage-safe stratified 5-fold cross-validation
#   6) Logistic regression vs XGBoost
#   7) Sensitivity analysis with/without waist circumference
#   8) Temporal validation (2007–2014 -> 2015–2016)
#   9) Calibration analysis
#  10) SHAP-based XGBoost interpretation
#  11) Out-of-fold ROC analysis
#
# Notes:
# - This script does NOT install packages automatically.
# - Install the required packages once before running:
#   foreign, dplyr, tidyr, stringr, caret, xgboost, pROC,
#   ggplot2, Matrix, readr
# - Run the script from the repository root directory.
# - Raw NHANES XPT files are not included in the repository.
# ============================================================

# -----------------------------
# 0. PACKAGES
# -----------------------------
library(foreign)
library(dplyr)
library(tidyr)
library(stringr)
library(caret)
library(xgboost)
library(pROC)
library(ggplot2)
library(Matrix)
library(readr)

set.seed(123)

# -----------------------------
# 1. CONFIGURATION
# -----------------------------
# The script assumes the five NHANES cycle folders listed below are
# located inside the repository root. See README.md for data download
# and folder-placement instructions.
PARENT_DIR <- getwd()

FOLDERS <- c(
  "NHANES_dataset_07-08",
  "NHANES_dataset_09-10",
  "NHANES_dataset_11-12",
  "NHANES_dataset_13-14",
  "NHANES_dataset_15-16"
)

RESULTS_DIR <- file.path(PARENT_DIR, "results_final")
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# -----------------------------
# 2. HELPER FUNCTIONS
# -----------------------------

read_xpt_safe <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    stop("Missing file: ", path)
  }
  as_tibble(read.xport(path))
}

clean_names <- function(df) {
  nm <- names(df)
  nm <- str_trim(nm)
  nm <- str_replace(nm, "^P_", "")
  names(df) <- nm
  df
}

find_file_any <- function(folder, prefixes) {
  files <- list.files(
    folder,
    pattern = "\\.xpt$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) return(NA_character_)

  base <- toupper(basename(files))

  for (p in prefixes) {
    P <- toupper(p)
    hit <- files[
      startsWith(base, paste0(P, "_")) |
      base == paste0(P, ".XPT") |
      startsWith(base, P)
    ]

    if (length(hit) > 0) return(hit[1])
  }

  NA_character_
}

pick_col <- function(df, candidates, label = "unknown") {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) > 0) return(hit[1])

  stop(
    "Could not find required column for ", label, ".\n",
    "Tried: ", paste(candidates, collapse = ", "), "\n",
    "Available columns (first 40): ",
    paste(head(names(df), 40), collapse = ", ")
  )
}

safe_row_mean <- function(df, cols) {
  out <- rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

compute_metrics <- function(y_true, y_prob, threshold = 0.5) {
  y_pred <- ifelse(y_prob >= threshold, 1, 0)

  TP <- sum(y_true == 1 & y_pred == 1)
  TN <- sum(y_true == 0 & y_pred == 0)
  FP <- sum(y_true == 0 & y_pred == 1)
  FN <- sum(y_true == 1 & y_pred == 0)

  accuracy <- (TP + TN) / (TP + TN + FP + FN)

  sensitivity <- ifelse(
    (TP + FN) > 0,
    TP / (TP + FN),
    NA_real_
  )

  specificity <- ifelse(
    (TN + FP) > 0,
    TN / (TN + FP),
    NA_real_
  )

  precision <- ifelse(
    (TP + FP) > 0,
    TP / (TP + FP),
    NA_real_
  )

  f1 <- ifelse(
    !is.na(precision) && !is.na(sensitivity) &&
      (precision + sensitivity) > 0,
    2 * precision * sensitivity / (precision + sensitivity),
    NA_real_
  )

  roc_obj <- pROC::roc(
    response = y_true,
    predictor = y_prob,
    quiet = TRUE
  )

  auc_value <- as.numeric(pROC::auc(roc_obj))

  data.frame(
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    ROC_AUC = auc_value
  )
}

make_calibration_table <- function(y_true, y_prob, bins = 10) {
  tibble(
    truth = y_true,
    pred = y_prob
  ) %>%
    mutate(bin = ntile(pred, bins)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(pred),
      observed_rate = mean(truth),
      n = n(),
      .groups = "drop"
    )
}

plot_calibration <- function(calib_df, title_text) {
  ggplot(calib_df, aes(x = mean_pred, y = observed_rate)) +
    geom_point() +
    geom_line() +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = title_text,
      x = "Mean predicted risk",
      y = "Observed MetS rate"
    )
}

# Select XGBoost nrounds using an INNER validation split only.
# The outer test fold is never used for early stopping.
select_xgb_nrounds <- function(train_df, features, seed = 123) {
  set.seed(seed)

  inner_idx <- caret::createDataPartition(
    factor(train_df$MetS),
    p = 0.80,
    list = FALSE
  )

  inner_train <- train_df[inner_idx, , drop = FALSE]
  inner_valid <- train_df[-inner_idx, , drop = FALSE]

  dtrain <- xgb.DMatrix(
    data = as.matrix(inner_train[, features, drop = FALSE]),
    label = inner_train$MetS
  )

  dvalid <- xgb.DMatrix(
    data = as.matrix(inner_valid[, features, drop = FALSE]),
    label = inner_valid$MetS
  )

  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = 0.05,
    max_depth = 4,
    subsample = 0.8,
    colsample_bytree = 0.8
  )

  fit <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 1000,
    evals = list(train = dtrain, valid = dvalid),
    early_stopping_rounds = 30,
    verbose = 0
  )

  best_nrounds <- fit$best_iteration

  if (is.null(best_nrounds) || is.na(best_nrounds)) {
    best_nrounds <- 300
  }

  best_nrounds
}

fit_xgb_fixed_rounds <- function(train_df, features, nrounds) {
  dtrain <- xgb.DMatrix(
    data = as.matrix(train_df[, features, drop = FALSE]),
    label = train_df$MetS
  )

  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = 0.05,
    max_depth = 4,
    subsample = 0.8,
    colsample_bytree = 0.8
  )

  xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )
}

run_cv_logistic <- function(model_df, folds, features, label) {
  bind_rows(lapply(seq_along(folds), function(i) {
    idx_test <- folds[[i]]

    train_df <- model_df[-idx_test, , drop = FALSE]
    test_df  <- model_df[idx_test, , drop = FALSE]

    form <- as.formula(
      paste("MetS ~", paste(features, collapse = " + "))
    )

    fit <- glm(
      form,
      data = train_df,
      family = binomial()
    )

    prob <- predict(
      fit,
      newdata = test_df,
      type = "response"
    )

    met <- compute_metrics(
      y_true = test_df$MetS,
      y_prob = prob,
      threshold = 0.5
    )

    met$Fold <- i
    met$Model <- label
    met
  }))
}

run_cv_xgb <- function(model_df, folds, features, label, base_seed = 123) {
  bind_rows(lapply(seq_along(folds), function(i) {
    idx_test <- folds[[i]]

    train_df <- model_df[-idx_test, , drop = FALSE]
    test_df  <- model_df[idx_test, , drop = FALSE]

    best_nrounds <- select_xgb_nrounds(
      train_df = train_df,
      features = features,
      seed = base_seed + i
    )

    fit <- fit_xgb_fixed_rounds(
      train_df = train_df,
      features = features,
      nrounds = best_nrounds
    )

    dtest <- xgb.DMatrix(
      data = as.matrix(test_df[, features, drop = FALSE])
    )

    prob <- predict(fit, dtest)

    met <- compute_metrics(
      y_true = test_df$MetS,
      y_prob = prob,
      threshold = 0.5
    )

    met$Fold <- i
    met$Model <- label
    met$Best_nrounds <- best_nrounds
    met
  }))
}

summarise_cv <- function(cv_results) {
  cv_results %>%
    group_by(Model) %>%
    summarise(
      Accuracy = mean(Accuracy, na.rm = TRUE),
      Sensitivity = mean(Sensitivity, na.rm = TRUE),
      Specificity = mean(Specificity, na.rm = TRUE),
      Precision = mean(Precision, na.rm = TRUE),
      F1 = mean(F1, na.rm = TRUE),
      ROC_AUC_mean = mean(ROC_AUC, na.rm = TRUE),
      ROC_AUC_SD = sd(ROC_AUC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(ROC_AUC_mean))
}

# -----------------------------
# 3. PROCESS ONE NHANES CYCLE
# -----------------------------
process_cycle_folder <- function(folder_name) {
  folder <- file.path(PARENT_DIR, folder_name)

  if (!dir.exists(folder)) {
    stop("Folder not found: ", folder)
  }

  f_demo <- find_file_any(folder, c("DEMO", "P_DEMO"))
  f_diet <- find_file_any(folder, c("DR1TOT", "P_DR1TOT"))
  f_bmx  <- find_file_any(folder, c("BMX", "P_BMX"))
  f_bpx  <- find_file_any(folder, c("BPX", "P_BPX"))
  f_glu  <- find_file_any(folder, c("GLU", "P_GLU"))
  f_trig <- find_file_any(folder, c("TRIGLY", "P_TRIGLY"))
  f_hdl  <- find_file_any(folder, c("HDL", "P_HDL"))

  demo <- clean_names(read_xpt_safe(f_demo))
  diet <- clean_names(read_xpt_safe(f_diet))
  bmx  <- clean_names(read_xpt_safe(f_bmx))
  bpx  <- clean_names(read_xpt_safe(f_bpx))
  glu  <- clean_names(read_xpt_safe(f_glu))
  trig <- clean_names(read_xpt_safe(f_trig))
  hdl  <- clean_names(read_xpt_safe(f_hdl))

  seqn_col_demo <- pick_col(demo, c("SEQN", "seqn"), "SEQN (demo)")
  seqn_col_diet <- pick_col(diet, c("SEQN", "seqn"), "SEQN (diet)")
  seqn_col_bmx  <- pick_col(bmx,  c("SEQN", "seqn"), "SEQN (bmx)")
  seqn_col_bpx  <- pick_col(bpx,  c("SEQN", "seqn"), "SEQN (bpx)")
  seqn_col_glu  <- pick_col(glu,  c("SEQN", "seqn"), "SEQN (glu)")
  seqn_col_trig <- pick_col(trig, c("SEQN", "seqn"), "SEQN (trig)")
  seqn_col_hdl  <- pick_col(hdl,  c("SEQN", "seqn"), "SEQN (hdl)")

  age_col <- pick_col(
    demo,
    c("RIDAGEYR", "RIDAGE_YR", "AGE", "Age"),
    "Age"
  )

  sex_col <- pick_col(
    demo,
    c("RIAGENDR", "SEX", "Sex"),
    "Sex"
  )

  kcal_col <- pick_col(
    diet,
    c("DR1TKCAL", "KCAL", "DR1_KCAL"),
    "Calories"
  )

  prot_col <- pick_col(
    diet,
    c("DR1TPROT", "PROT", "DR1_PROT"),
    "Protein"
  )

  fibe_col <- pick_col(
    diet,
    c("DR1TFIBE", "FIBER", "DR1_FIBE"),
    "Fiber"
  )

  waist_col <- pick_col(
    bmx,
    c("BMXWAIST", "WAIST", "WC", "BMX_WAIST"),
    "Waist"
  )

  glu_col <- pick_col(
    glu,
    c("LBXGLU", "GLU", "LBX_GLU"),
    "Glucose"
  )

  trig_col <- pick_col(
    trig,
    c("LBXTR", "TG", "TRIG", "LBX_TR"),
    "Triglycerides"
  )

  hdl_col <- pick_col(
    hdl,
    c("LBDHDD", "HDL", "LBD_HDD"),
    "HDL"
  )

  sbp_candidates <- c(
    "BPXSY1","BPXSY2","BPXSY3","BPXSY4",
    "SBP1","SBP2","SBP3","SBP4"
  )

  dbp_candidates <- c(
    "BPXDI1","BPXDI2","BPXDI3","BPXDI4",
    "DBP1","DBP2","DBP3","DBP4"
  )

  sbp_cols <- sbp_candidates[sbp_candidates %in% names(bpx)]
  dbp_cols <- dbp_candidates[dbp_candidates %in% names(bpx)]

  if (length(sbp_cols) == 0 || length(dbp_cols) == 0) {
    stop("Could not find BP columns in folder: ", folder_name)
  }

  demo_std <- demo %>%
    transmute(
      SEQN = .data[[seqn_col_demo]],
      RIDAGEYR = as.numeric(.data[[age_col]]),
      RIAGENDR = as.numeric(.data[[sex_col]])
    )

  diet_std <- diet %>%
    transmute(
      SEQN = .data[[seqn_col_diet]],
      DR1TKCAL = as.numeric(.data[[kcal_col]]),
      DR1TPROT = as.numeric(.data[[prot_col]]),
      DR1TFIBE = as.numeric(.data[[fibe_col]])
    )

  bmx_std <- bmx %>%
    transmute(
      SEQN = .data[[seqn_col_bmx]],
      BMXWAIST = as.numeric(.data[[waist_col]])
    )

  glu_std <- glu %>%
    transmute(
      SEQN = .data[[seqn_col_glu]],
      LBXGLU = as.numeric(.data[[glu_col]])
    )

  trig_std <- trig %>%
    transmute(
      SEQN = .data[[seqn_col_trig]],
      LBXTR = as.numeric(.data[[trig_col]])
    )

  hdl_std <- hdl %>%
    transmute(
      SEQN = .data[[seqn_col_hdl]],
      LBDHDD = as.numeric(.data[[hdl_col]])
    )

  bpx_std <- bpx %>%
    mutate(
      SBP = safe_row_mean(., sbp_cols),
      DBP = safe_row_mean(., dbp_cols)
    ) %>%
    transmute(
      SEQN = .data[[seqn_col_bpx]],
      SBP,
      DBP
    )

  df <- demo_std %>%
    inner_join(diet_std, by = "SEQN") %>%
    inner_join(bmx_std,  by = "SEQN") %>%
    inner_join(glu_std,  by = "SEQN") %>%
    inner_join(trig_std, by = "SEQN") %>%
    inner_join(hdl_std,  by = "SEQN") %>%
    inner_join(bpx_std,  by = "SEQN")

  # Adult-only analysis
  df <- df %>%
    filter(RIDAGEYR >= 18) %>%
    filter(DR1TKCAL > 0) %>%
    mutate(
      fiber_per_1000kcal =
        DR1TFIBE / (DR1TKCAL / 1000),
      protein_per_1000kcal =
        DR1TPROT / (DR1TKCAL / 1000)
    ) %>%
    mutate(
      across(
        c(
          fiber_per_1000kcal,
          protein_per_1000kcal
        ),
        ~ ifelse(is.finite(.x), .x, NA_real_)
      )
    ) %>%
    drop_na(
      RIDAGEYR, RIAGENDR,
      DR1TKCAL, DR1TPROT, DR1TFIBE,
      BMXWAIST, LBXGLU, LBXTR, LBDHDD,
      SBP, DBP,
      fiber_per_1000kcal,
      protein_per_1000kcal
    )

  male <- df$RIAGENDR == 1
  female <- df$RIAGENDR == 2

  df <- df %>%
    mutate(
      waist_risk = as.integer(
        (male & BMXWAIST >= 102) |
        (female & BMXWAIST >= 88)
      ),
      tg_risk = as.integer(LBXTR >= 150),
      glu_risk = as.integer(LBXGLU >= 100),
      hdl_risk = as.integer(
        (male & LBDHDD < 40) |
        (female & LBDHDD < 50)
      ),
      bp_risk = as.integer(
        SBP >= 130 | DBP >= 85
      ),
      risk_count =
        waist_risk +
        tg_risk +
        glu_risk +
        hdl_risk +
        bp_risk,
      MetS = as.integer(risk_count >= 3),
      cycle = folder_name
    )

  df
}

# -----------------------------
# 4. BUILD COMBINED DATASET
# -----------------------------
all_cycles <- lapply(FOLDERS, function(fn) {
  cat("\nProcessing:", fn, "\n")
  out <- process_cycle_folder(fn)
  cat("Rows retained:", nrow(out), "\n")
  out
})

df_all <- bind_rows(all_cycles)

cat("\n==============================\n")
cat("Combined adult dataset created\n")
cat("Total rows:", nrow(df_all), "\n")
cat("MetS prevalence:", round(mean(df_all$MetS), 4), "\n")
cat("==============================\n")

saveRDS(
  df_all,
  file.path(RESULTS_DIR, "NHANES_combined_adults_2007_2016.rds")
)

write_csv(
  df_all,
  file.path(RESULTS_DIR, "NHANES_combined_adults_2007_2016.csv")
)

# -----------------------------
# 5. DESCRIPTIVE SUMMARIES
# -----------------------------
overall_summary <- df_all %>%
  summarise(
    n = n(),
    age_mean = mean(RIDAGEYR),
    age_sd = sd(RIDAGEYR),
    age_median = median(RIDAGEYR),
    age_min = min(RIDAGEYR),
    age_max = max(RIDAGEYR),
    male_n = sum(RIAGENDR == 1),
    female_n = sum(RIAGENDR == 2),
    MetS_n = sum(MetS == 1),
    MetS_prevalence = mean(MetS)
  )

cycle_summary <- df_all %>%
  group_by(cycle) %>%
  summarise(
    n = n(),
    MetS_n = sum(MetS == 1),
    MetS_prevalence = mean(MetS),
    mean_age = mean(RIDAGEYR),
    .groups = "drop"
  )

write_csv(
  overall_summary,
  file.path(RESULTS_DIR, "Descriptive_Overall.csv")
)

write_csv(
  cycle_summary,
  file.path(RESULTS_DIR, "Descriptive_By_Cycle.csv")
)

print(overall_summary)
print(cycle_summary)

# -----------------------------
# 6. MODELING DATASET
# -----------------------------
model_df <- df_all %>%
  select(
    MetS,
    cycle,
    RIAGENDR,
    RIDAGEYR,
    BMXWAIST,
    DR1TKCAL,
    DR1TPROT,
    DR1TFIBE,
    fiber_per_1000kcal,
    protein_per_1000kcal
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(is.infinite(.x), NA_real_, .x)
    )
  ) %>%
  drop_na()

# Primary leakage-resistant feature set:
# excludes waist because waist is itself part of MetS definition
features_no_waist <- c(
  "RIAGENDR",
  "RIDAGEYR",
  "DR1TKCAL",
  "DR1TPROT",
  "DR1TFIBE",
  "fiber_per_1000kcal",
  "protein_per_1000kcal"
)

# Sensitivity / augmented feature set:
# includes waist circumference
features_with_waist <- c(
  features_no_waist,
  "BMXWAIST"
)

# -----------------------------
# 7. STRATIFIED 5-FOLD OUTER CV
# -----------------------------
set.seed(123)

folds <- caret::createFolds(
  factor(model_df$MetS),
  k = 5,
  list = TRUE,
  returnTrain = FALSE
)

# Logistic, no waist
cv_log_no_waist <- run_cv_logistic(
  model_df = model_df,
  folds = folds,
  features = features_no_waist,
  label = "Logistic: no waist"
)

# XGBoost, no waist
cv_xgb_no_waist <- run_cv_xgb(
  model_df = model_df,
  folds = folds,
  features = features_no_waist,
  label = "XGBoost: no waist",
  base_seed = 1000
)

# Logistic, with waist
cv_log_with_waist <- run_cv_logistic(
  model_df = model_df,
  folds = folds,
  features = features_with_waist,
  label = "Logistic: with waist"
)

# XGBoost, with waist
cv_xgb_with_waist <- run_cv_xgb(
  model_df = model_df,
  folds = folds,
  features = features_with_waist,
  label = "XGBoost: with waist",
  base_seed = 2000
)

all_cv <- bind_rows(
  cv_log_no_waist,
  cv_xgb_no_waist,
  cv_log_with_waist,
  cv_xgb_with_waist
)

cv_summary <- summarise_cv(all_cv)

write_csv(
  all_cv,
  file.path(RESULTS_DIR, "CV_Fold_Level_Results.csv")
)

write_csv(
  cv_summary,
  file.path(RESULTS_DIR, "CV_Model_Comparison_Summary.csv")
)

print(cv_summary)

# -----------------------------
# 8. TEMPORAL VALIDATION
#    Train: 2007–2014
#    Test:  2015–2016
# -----------------------------
temporal_train <- model_df %>%
  filter(cycle != "NHANES_dataset_15-16")

temporal_test <- model_df %>%
  filter(cycle == "NHANES_dataset_15-16")

# Logistic: no waist
form_temporal_log_no_waist <- as.formula(
  paste(
    "MetS ~",
    paste(features_no_waist, collapse = " + ")
  )
)

fit_temporal_log_no_waist <- glm(
  form_temporal_log_no_waist,
  data = temporal_train,
  family = binomial()
)

prob_temporal_log_no_waist <- predict(
  fit_temporal_log_no_waist,
  newdata = temporal_test,
  type = "response"
)

met_temporal_log_no_waist <- compute_metrics(
  temporal_test$MetS,
  prob_temporal_log_no_waist
) %>%
  mutate(Model = "Logistic: no waist")

# XGBoost: no waist
best_rounds_temporal_no_waist <- select_xgb_nrounds(
  temporal_train,
  features_no_waist,
  seed = 3001
)

fit_temporal_xgb_no_waist <- fit_xgb_fixed_rounds(
  temporal_train,
  features_no_waist,
  best_rounds_temporal_no_waist
)

dtest_temporal_no_waist <- xgb.DMatrix(
  as.matrix(
    temporal_test[, features_no_waist, drop = FALSE]
  )
)

prob_temporal_xgb_no_waist <- predict(
  fit_temporal_xgb_no_waist,
  dtest_temporal_no_waist
)

met_temporal_xgb_no_waist <- compute_metrics(
  temporal_test$MetS,
  prob_temporal_xgb_no_waist
) %>%
  mutate(Model = "XGBoost: no waist")

# Logistic: with waist
form_temporal_log_with_waist <- as.formula(
  paste(
    "MetS ~",
    paste(features_with_waist, collapse = " + ")
  )
)

fit_temporal_log_with_waist <- glm(
  form_temporal_log_with_waist,
  data = temporal_train,
  family = binomial()
)

prob_temporal_log_with_waist <- predict(
  fit_temporal_log_with_waist,
  newdata = temporal_test,
  type = "response"
)

met_temporal_log_with_waist <- compute_metrics(
  temporal_test$MetS,
  prob_temporal_log_with_waist
) %>%
  mutate(Model = "Logistic: with waist")

# XGBoost: with waist
best_rounds_temporal_with_waist <- select_xgb_nrounds(
  temporal_train,
  features_with_waist,
  seed = 3002
)

fit_temporal_xgb_with_waist <- fit_xgb_fixed_rounds(
  temporal_train,
  features_with_waist,
  best_rounds_temporal_with_waist
)

dtest_temporal_with_waist <- xgb.DMatrix(
  as.matrix(
    temporal_test[, features_with_waist, drop = FALSE]
  )
)

prob_temporal_xgb_with_waist <- predict(
  fit_temporal_xgb_with_waist,
  dtest_temporal_with_waist
)

met_temporal_xgb_with_waist <- compute_metrics(
  temporal_test$MetS,
  prob_temporal_xgb_with_waist
) %>%
  mutate(Model = "XGBoost: with waist")

temporal_results <- bind_rows(
  met_temporal_log_no_waist,
  met_temporal_xgb_no_waist,
  met_temporal_log_with_waist,
  met_temporal_xgb_with_waist
)

write_csv(
  temporal_results,
  file.path(RESULTS_DIR, "Temporal_Validation_Results.csv")
)

print(temporal_results)

# -----------------------------
# 9. TEMPORAL CALIBRATION
# -----------------------------
calib_xgb_no_waist <- make_calibration_table(
  temporal_test$MetS,
  prob_temporal_xgb_no_waist
)

p_calib_xgb_no_waist <- plot_calibration(
  calib_xgb_no_waist,
  "Temporal Calibration: XGBoost Without Waist"
)

ggsave(
  filename = file.path(
    RESULTS_DIR,
    "Temporal_Calibration_XGBoost_No_Waist.png"
  ),
  plot = p_calib_xgb_no_waist,
  width = 6,
  height = 6,
  dpi = 300
)

calib_xgb_with_waist <- make_calibration_table(
  temporal_test$MetS,
  prob_temporal_xgb_with_waist
)

p_calib_xgb_with_waist <- plot_calibration(
  calib_xgb_with_waist,
  "Temporal Calibration: XGBoost With Waist"
)

ggsave(
  filename = file.path(
    RESULTS_DIR,
    "Temporal_Calibration_XGBoost_With_Waist.png"
  ),
  plot = p_calib_xgb_with_waist,
  width = 6,
  height = 6,
  dpi = 300
)

write_csv(
  calib_xgb_no_waist,
  file.path(
    RESULTS_DIR,
    "Temporal_Calibration_Data_XGBoost_No_Waist.csv"
  )
)

write_csv(
  calib_xgb_with_waist,
  file.path(
    RESULTS_DIR,
    "Temporal_Calibration_Data_XGBoost_With_Waist.csv"
  )
)

# -----------------------------
# 10. FINAL XGBOOST MODEL
#     Primary: no-waist model
# -----------------------------
best_rounds_final <- select_xgb_nrounds(
  model_df,
  features_no_waist,
  seed = 4001
)

fit_final_xgb <- fit_xgb_fixed_rounds(
  model_df,
  features_no_waist,
  best_rounds_final
)

saveRDS(
  fit_final_xgb,
  file.path(RESULTS_DIR, "Final_XGBoost_No_Waist.rds")
)

# Final logistic model
final_formula_log <- as.formula(
  paste(
    "MetS ~",
    paste(features_no_waist, collapse = " + ")
  )
)

fit_final_logistic <- glm(
  final_formula_log,
  data = model_df,
  family = binomial()
)

saveRDS(
  fit_final_logistic,
  file.path(RESULTS_DIR, "Final_Logistic_No_Waist.rds")
)


# -----------------------------
# 11. SHAP ANALYSIS FOR FINAL XGBOOST MODEL
#     Primary model: without waist circumference
# -----------------------------
X_final <- as.matrix(
  model_df[, features_no_waist, drop = FALSE]
)

d_final <- xgb.DMatrix(data = X_final)

# ------------------------------------------------------------
# 2. Obtain SHAP contributions
# ------------------------------------------------------------

shap_matrix <- predict(
  fit_final_xgb,
  d_final,
  predcontrib = TRUE
)

shap_df <- as.data.frame(shap_matrix)

# Check column names
print(colnames(shap_df))

# ------------------------------------------------------------
# 3. Remove intercept / baseline column correctly
# ------------------------------------------------------------

intercept_cols <- colnames(shap_df)[
  grepl(
    "BIAS|Intercept",
    colnames(shap_df),
    ignore.case = TRUE
  )
]

if (length(intercept_cols) > 0) {
  shap_features <- shap_df %>%
    select(-all_of(intercept_cols))
} else {
  shap_features <- shap_df
}

cat(
  "Removed intercept column(s):",
  paste(intercept_cols, collapse = ", "),
  "\n"
)

# Sanity check
print(colnames(shap_features))

# ------------------------------------------------------------
# 4. Calculate mean absolute SHAP importance
# ------------------------------------------------------------

shap_importance <- tibble(
  Feature = colnames(shap_features),
  MeanAbsSHAP = colMeans(
    abs(shap_features),
    na.rm = TRUE
  )
) %>%
  arrange(desc(MeanAbsSHAP))

print(shap_importance)

# Save numerical results
write_csv(
  shap_importance,
  file.path(
    RESULTS_DIR,
    "XGBoost_SHAP_Feature_Importance_Corrected.csv"
  )
)

# ------------------------------------------------------------
# 5. Give variables publication-friendly labels
# ------------------------------------------------------------

shap_importance_plot <- shap_importance %>%
  mutate(
    Feature_Label = case_when(
      Feature == "RIDAGEYR" ~ "Age",
      Feature == "RIAGENDR" ~ "Sex",
      Feature == "DR1TKCAL" ~ "Total energy intake",
      Feature == "DR1TPROT" ~ "Protein intake",
      Feature == "DR1TFIBE" ~ "Fiber intake",
      Feature == "fiber_per_1000kcal" ~ "Fiber per 1,000 kcal",
      Feature == "protein_per_1000kcal" ~ "Protein per 1,000 kcal",
      TRUE ~ Feature
    )
  )

# ------------------------------------------------------------
# 6. Plot mean absolute SHAP importance
# ------------------------------------------------------------

p_shap_importance <- ggplot(
  shap_importance_plot,
  aes(
    x = reorder(Feature_Label, MeanAbsSHAP),
    y = MeanAbsSHAP
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Feature Importance in the Waist-Excluded XGBoost Model",
    x = NULL,
    y = "Mean absolute SHAP value"
  ) +
  theme_minimal(base_size = 14)

print(p_shap_importance)

ggsave(
  filename = file.path(
    RESULTS_DIR,
    "XGBoost_SHAP_Feature_Importance_Corrected.png"
  ),
  plot = p_shap_importance,
  width = 8,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 7. OPTIONAL: SHAP BEESWARM-STYLE PLOT
# Shows both magnitude AND direction of effects
# ============================================================

# Convert SHAP values to long format
shap_long <- shap_features %>%
  mutate(RowID = row_number()) %>%
  pivot_longer(
    cols = -RowID,
    names_to = "Feature",
    values_to = "SHAP"
  )

# Convert original feature values to long format
feature_long <- as.data.frame(X_final) %>%
  mutate(RowID = row_number()) %>%
  pivot_longer(
    cols = -RowID,
    names_to = "Feature",
    values_to = "FeatureValue"
  )

# Combine SHAP contributions with actual feature values
shap_plot_data <- shap_long %>%
  left_join(
    feature_long,
    by = c("RowID", "Feature")
  ) %>%
  mutate(
    Feature_Label = case_when(
      Feature == "RIDAGEYR" ~ "Age",
      Feature == "RIAGENDR" ~ "Sex",
      Feature == "DR1TKCAL" ~ "Total energy intake",
      Feature == "DR1TPROT" ~ "Protein intake",
      Feature == "DR1TFIBE" ~ "Fiber intake",
      Feature == "fiber_per_1000kcal" ~ "Fiber per 1,000 kcal",
      Feature == "protein_per_1000kcal" ~ "Protein per 1,000 kcal",
      TRUE ~ Feature
    )
  )

# Order features by mean absolute SHAP
feature_order <- shap_importance_plot %>%
  arrange(MeanAbsSHAP) %>%
  pull(Feature_Label)

shap_plot_data$Feature_Label <- factor(
  shap_plot_data$Feature_Label,
  levels = feature_order
)

# Random sample for plotting if dataset is very large
set.seed(123)

shap_plot_sample <- shap_plot_data %>%
  group_by(Feature_Label) %>%
  slice_sample(prop = 1) %>%   # randomly shuffle rows within each feature
  slice_head(n = 1500) %>%     # keep up to 1500 rows per feature
  ungroup()

p_shap_beeswarm <- ggplot(
  shap_plot_sample,
  aes(
    x = SHAP,
    y = Feature_Label,
    color = FeatureValue
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "SHAP Contributions in the Waist-Excluded XGBoost Model",
    x = "SHAP value",
    y = NULL,
    color = "Feature value"
  ) +
  theme_minimal(base_size = 14)

print(p_shap_beeswarm)

ggsave(
  filename = file.path(
    RESULTS_DIR,
    "XGBoost_SHAP_Beeswarm_Corrected.png"
  ),
  plot = p_shap_beeswarm,
  width = 9,
  height = 6,
  dpi = 300
)

# ============================================================
# PUBLICATION-QUALITY SHAP SUMMARY PLOT
# ============================================================

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Scale feature values WITHIN each feature
#    0 = lowest observed value
#    1 = highest observed value
# ------------------------------------------------------------

shap_plot_final <- shap_plot_data %>%
  group_by(Feature_Label) %>%
  mutate(
    FeatureValue_scaled =
      (FeatureValue - min(FeatureValue, na.rm = TRUE)) /
      (max(FeatureValue, na.rm = TRUE) -
         min(FeatureValue, na.rm = TRUE))
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 2. Randomly sample maximum 2000 observations per feature
# ------------------------------------------------------------

set.seed(123)

shap_plot_final <- shap_plot_final %>%
  group_by(Feature_Label) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = 2000) %>%
  ungroup()

# ------------------------------------------------------------
# 3. Preserve feature order according to SHAP importance
# ------------------------------------------------------------

feature_order <- shap_importance_plot %>%
  arrange(MeanAbsSHAP) %>%
  pull(Feature_Label)

shap_plot_final$Feature_Label <- factor(
  shap_plot_final$Feature_Label,
  levels = feature_order
)

# ------------------------------------------------------------
# 4. Create SHAP summary plot
# ------------------------------------------------------------

p_shap_summary <- ggplot(
  shap_plot_final,
  aes(
    x = SHAP,
    y = Feature_Label,
    color = FeatureValue_scaled
  )
) +
  geom_jitter(
    height = 0.18,
    width = 0,
    alpha = 0.55,
    size = 1.1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    breaks = c(0, 1),
    labels = c("Low", "High")
  ) +
  labs(
    title = "SHAP Summary Plot",
    subtitle = "Waist-excluded XGBoost model",
    x = "SHAP value",
    y = NULL,
    color = "Feature value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 13
    ),
    axis.text.y = element_text(
      size = 12
    ),
    legend.position = "right"
  )

print(p_shap_summary)

# ------------------------------------------------------------
# 5. Save high-resolution figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    RESULTS_DIR,
    "XGBoost_SHAP_Summary_Publication.png"
  ),
  plot = p_shap_summary,
  width = 9,
  height = 6,
  dpi = 600
)

# -----------------------------
# 12. OPTIONAL: OUTER-CV ROC CURVE
#     Generate out-of-fold predictions for primary XGBoost
# -----------------------------
oof_predictions <- bind_rows(lapply(seq_along(folds), function(i) {
  idx_test <- folds[[i]]

  train_df <- model_df[-idx_test, , drop = FALSE]
  test_df  <- model_df[idx_test, , drop = FALSE]

  best_nrounds <- select_xgb_nrounds(
    train_df,
    features_no_waist,
    seed = 5000 + i
  )

  fit <- fit_xgb_fixed_rounds(
    train_df,
    features_no_waist,
    best_nrounds
  )

  dtest <- xgb.DMatrix(
    as.matrix(
      test_df[, features_no_waist, drop = FALSE]
    )
  )

  prob <- predict(fit, dtest)

  tibble(
    Fold = i,
    Truth = test_df$MetS,
    Probability = prob
  )
}))

roc_oof <- pROC::roc(
  oof_predictions$Truth,
  oof_predictions$Probability,
  quiet = TRUE
)

png(
  filename = file.path(
    RESULTS_DIR,
    "OOF_ROC_XGBoost_No_Waist.png"
  ),
  width = 1800,
  height = 1600,
  res = 300
)

plot(
  roc_oof,
  main = paste0(
    "Out-of-Fold ROC: XGBoost Without Waist\nAUC = ",
    round(as.numeric(auc(roc_oof)), 3)
  )
)

abline(a = 0, b = 1, lty = 2)
dev.off()

write_csv(
  oof_predictions,
  file.path(
    RESULTS_DIR,
    "OOF_Predictions_XGBoost_No_Waist.csv"
  )
)


# -----------------------------
# 13. SAVE SESSION INFO
# -----------------------------
sink(file.path(RESULTS_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nAnalysis complete.\n")
cat("Results saved in:\n", RESULTS_DIR, "\n")
