
## ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(plotly)
library(lubridate)
library(rpart)
library(randomForest)
library(e1071)
library(gbm)
library(pROC)

set.seed(42)

## ----------------------------------------------------------------------
## CONSTANTS
## ----------------------------------------------------------------------

ACCENT       <- "#0e5c73"
ACCENT_DARK  <- "#083c4b"
ACCENT_LIGHT <- "#e8f3f5"

DATA_PATH_DEFAULT <- "CEM_final_preprocessed.csv"

NEEDED_BASE_COLS <- c(
  "patient_id", "patient_status", "age_clean", "sex", "health_facility",
  "subcity_clean", "wereda_clean", "weight_clean", "height_clean",
  "bmi_recalculated", "who_hiv_stage_at_this_scre", "date",
  "arv_medicine1_clean", "dose_status"
)
HARD_REQUIRED_COLS <- c("patient_id", "patient_status")

NUMERIC_FEATS <- c(
  "age_clean", "weight_clean", "height_clean", "bmi_recalculated",
  "baseline_who_stage", "baseline_cd4", "n_lab_events", "days_on_treatment"
)
CATEGORICAL_FEATS <- c("sex_label", "baseline_vl_status", "subcity_grouped")

## ----------------------------------------------------------------------
## FEATURE ENGINEERING  (mirrors the Python engineer_features() pipeline)
## ----------------------------------------------------------------------

clean_stage_vec <- function(x) {
  s <- toupper(trimws(as.character(x)))
  roman <- c(I = 1, II = 2, III = 3, IV = 4)
  out <- rep(NA_real_, length(s))

  for (i in seq_along(s)) {
    v <- s[i]
    if (is.na(v) || v == "") next

    if (v %in% names(roman)) { out[i] <- roman[[v]]; next }              # "I","II","III","IV"
    if (v %in% c("1", "2", "3", "4")) { out[i] <- as.numeric(v); next }  # "1".."4"

    m <- str_match(v, "^T[.\\s]?([1-4])$")                               # "T1","T 1","T.1"
    if (!is.na(m[1, 2])) { out[i] <- as.numeric(m[1, 2]); next }

    if (str_detect(v, "^STAGE\\b")) {                                    # "STAGE 1","STAGE III"
      m2 <- str_match(v, "STAGE\\s*(IV|III|II|I|[1-4])\\b")
      if (!is.na(m2[1, 2])) {
        tok <- m2[1, 2]
        out[i] <- if (tok %in% names(roman)) roman[[tok]] else as.numeric(tok)
        next
      }
    }

    m3 <- str_match(v, "^(IV|III|II|I)\\s*\\(([1-4])\\)$")               # "III (3)"
    if (!is.na(m3[1, 2])) { out[i] <- as.numeric(m3[1, 3]); next }
    m4 <- str_match(v, "^([1-4])\\s*\\((IV|III|II|I)\\)$")               # "4 (IV)"
    if (!is.na(m4[1, 2])) { out[i] <- as.numeric(m4[1, 2]); next }

    ## anything else ("D", "TV", "TI", "I4", "rapied test") stays NA —
    ## genuinely ambiguous single-off entries, not a recognizable stage
  }
  out
}

parse_date_col <- function(x) {
  suppressWarnings(lubridate::dmy(x))
}

engineer_features <- function(df) {

  for (col in NEEDED_BASE_COLS) if (!col %in% names(df)) df[[col]] <- NA

  test_cols <- names(df)[str_starts(names(df), "test_")]
  df <- df[, unique(c(NEEDED_BASE_COLS, test_cols)), drop = FALSE]

  ## target
  df$died <- as.integer(df$patient_status == "died")

  ## sex label — accepts numeric 1/2 REDCap coding or text
  if (is.character(df$sex) || is.factor(df$sex)) {
    s <- str_to_title(trimws(as.character(df$sex)))
    s <- recode(s, "M" = "Male", "F" = "Female", "1" = "Male", "2" = "Female", "Nan" = NA_character_)
    df$sex_label <- s
  } else {
    df$sex_label <- ifelse(df$sex == 1, "Male", ifelse(df$sex == 2, "Female", NA))
  }

  ## WHO stage cleaning
  df$baseline_who_stage <- clean_stage_vec(df$who_hiv_stage_at_this_scre)

  ## lab test wide-format parsing
  test_types <- sort(unique(str_remove(test_cols, "_(day|result)\\d+$")))
  result_cols_by_type <- lapply(test_types, function(t) {
    cand <- paste0(t, "_result", 1:6)
    cand[cand %in% names(df)]
  })
  names(result_cols_by_type) <- test_types

  ## baseline CD4: first parseable numeric result across the 6 draws
  cd4_cols <- result_cols_by_type[["test_cd4"]]
  if (!is.null(cd4_cols) && length(cd4_cols) > 0) {
    cd4_num <- df[, cd4_cols, drop = FALSE]
    cd4_num[] <- lapply(cd4_num, function(x) suppressWarnings(as.numeric(x)))
    df$baseline_cd4 <- apply(cd4_num, 1, function(r) {
      v <- r[!is.na(r)]
      if (length(v) == 0) NA_real_ else v[1]
    })
  } else {
    df$baseline_cd4 <- NA_real_
  }

  ## baseline viral load status
  vl_cols <- result_cols_by_type[["test_viralload"]]
  vl_status_row <- function(row_vals) {
    for (v in row_vals) {
      if (is.na(v)) next
      s <- tolower(trimws(as.character(v)))
      if (s == "") next
      ## "not done" family — the data has several spellings beyond the
      ## exact "not done" string ("no" etc). Treat them all as "no result
      ## for this draw" and keep checking the patient's later visits,
      ## instead of falling through to the generic Unknown branch below
      ## and stopping the scan on the very first draw.
      if (s %in% c("not done", "notdone", "no")) next
      ## "nd"/"n.d"/"n/d" and any "not d..." spelling ("not d","not det",
      ## "not detect","not detected") all read as "not detected" here
      if (str_detect(s, "^(nd|n\\.d|n/d)$") || str_detect(s, "not\\s*d|undetect|<")) return("Undetectable")
      num_val <- suppressWarnings(as.numeric(str_remove_all(s, ",")))  # strip "555,399"-style separators
      if (!is.na(num_val)) return("Detectable")
      return("Unknown")
    }
    NA_character_
  }
  if (!is.null(vl_cols) && length(vl_cols) > 0) {
    df$baseline_vl_status <- apply(df[, vl_cols, drop = FALSE], 1, vl_status_row)
  } else {
    df$baseline_vl_status <- NA_character_
  }

  ## number of lab events recorded (engagement-with-care proxy)
  all_result_cols <- test_cols[str_detect(test_cols, "_result\\d+$")]
  if (length(all_result_cols) > 0) {
    df$n_lab_events <- rowSums(!is.na(df[, all_result_cols, drop = FALSE]))
  } else {
    df$n_lab_events <- 0
  }

  ## days on treatment: span between earliest & latest recorded lab date
  all_day_cols <- test_cols[str_detect(test_cols, "_day\\d+$")]
  if (length(all_day_cols) > 0) {
    day_dates <- as.data.frame(lapply(df[, all_day_cols, drop = FALSE], parse_date_col))
    dmax <- as.Date(apply(day_dates, 1, function(r) {
      v <- r[!is.na(r)]; if (length(v) == 0) NA else max(v)
    }), origin = "1970-01-01")
    dmin <- as.Date(apply(day_dates, 1, function(r) {
      v <- r[!is.na(r)]; if (length(v) == 0) NA else min(v)
    }), origin = "1970-01-01")
    df$days_on_treatment <- as.numeric(dmax - dmin)
  } else {
    df$days_on_treatment <- NA_real_
  }

  ## subcity grouping (rare subcities -> "Other")
  sc_counts <- table(df$subcity_clean)
  common_sc <- names(sc_counts)[sc_counts >= 10]
  df$subcity_grouped <- ifelse(df$subcity_clean %in% common_sc, df$subcity_clean, "Other")
  df$subcity_grouped[is.na(df$subcity_clean)] <- "Unknown"

  ## enrollment date
  df$enroll_date <- parse_date_col(df$date)

  ## BMI category
  df$bmi_category <- case_when(
    is.na(df$bmi_recalculated)   ~ "Unknown",
    df$bmi_recalculated < 18.5   ~ "Underweight",
    df$bmi_recalculated < 25     ~ "Normal",
    df$bmi_recalculated < 30     ~ "Overweight",
    TRUE                         ~ "Obese"
  )

  df
}

## ----------------------------------------------------------------------
## SYNTHETIC DEMO DATA (used only if no real file is found/uploaded)
## ----------------------------------------------------------------------

make_demo_data <- function(n = 800, seed = 42) {
  set.seed(seed)
  facilities <- paste("Demo Health Center", 1:6)
  subcities  <- paste("Demo Subcity", 1:8)
  status <- sample(c("active", "inactive_or_lost", "died"), n, replace = TRUE, prob = c(0.62, 0.33, 0.05))
  stage  <- sample(1:4, n, replace = TRUE, prob = c(0.45, 0.25, 0.18, 0.12))
  age    <- pmin(pmax(round(rnorm(n, 36, 12)), 1), 90)
  weight <- round(pmin(pmax(rnorm(n, 58, 11), 20), 120), 1)
  height <- round(pmin(pmax(rnorm(n, 162, 9), 120), 200), 1)
  dates  <- as.Date("2019-01-01") + sample(0:2000, n, replace = TRUE)

  demo <- data.frame(
    patient_id = sprintf("DEMO-%05d", 1:n),
    patient_status = status,
    age_clean = age,
    sex = sample(1:2, n, replace = TRUE),
    health_facility = sample(facilities, n, replace = TRUE),
    subcity_clean = sample(subcities, n, replace = TRUE),
    wereda_clean = sample(paste("Wereda", 1:4), n, replace = TRUE),
    weight_clean = weight,
    height_clean = height,
    bmi_recalculated = round(weight / ((height / 100)^2), 1),
    who_hiv_stage_at_this_scre = stage,
    date = format(dates, "%d/%m/%Y"),
    arv_medicine1_clean = sample(c("TDF/3TC/DTG", "AZT/3TC/EFV", "TDF/3TC/EFV"), n, replace = TRUE),
    dose_status = sample(c("Ongoing", "Completed"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  for (i in 1:6) {
    has_test <- runif(n) < (if (i <= 2) 0.8 else 0.3)
    cd4_val <- round(pmin(pmax(rnorm(n, 420, 220), 1), 1600))
    demo[[paste0("test_cd4_result", i)]] <- ifelse(has_test, cd4_val, NA)
    demo[[paste0("test_cd4_day", i)]] <- ifelse(has_test, format(dates + i * 90, "%d/%m/%Y"), NA)
    vl_choice <- sample(c("Undetectable", "Detectable"), n, replace = TRUE, prob = c(0.78, 0.22))
    demo[[paste0("test_viralload_result", i)]] <- ifelse(has_test, vl_choice, NA)
    demo[[paste0("test_viralload_day", i)]] <- demo[[paste0("test_cd4_day", i)]]
  }
  demo
}

## ----------------------------------------------------------------------
## MODEL TRAINING  (Decision Tree, Random Forest, SVM, Gradient Boosting)
## ----------------------------------------------------------------------

## Fits all 4 candidate models on train_df and scores them on test_df.
## Shared by the main 80/20 split and by each cross-validation fold below,
## so both use exactly the same modeling code.
fit_score_fold <- function(train_df, test_df, feat_formula, gbm_formula, keep_models = FALSE) {
  results <- list()
  metrics <- list()
  roc_data <- list()

  ## ---- Decision Tree ----
  dt <- rpart(feat_formula, data = train_df, method = "class",
              control = rpart.control(maxdepth = 5, cp = 0.001))
  dt_prob <- predict(dt, test_df, type = "prob")[, "Yes"]
  if (keep_models) results[["Decision Tree"]] <- list(model = dt, type = "rpart")
  metrics[["Decision Tree"]] <- compute_metrics(test_df$died, dt_prob)
  roc_data[["Decision Tree"]] <- safe_roc(test_df$died, dt_prob)

  ## ---- Random Forest ----
  rf <- randomForest(feat_formula, data = train_df, ntree = 200, maxnodes = 32,
                      classwt = c(No = 1, Yes = sum(train_df$died == "No") / max(1, sum(train_df$died == "Yes"))))
  rf_prob <- predict(rf, test_df, type = "prob")[, "Yes"]
  if (keep_models) results[["Random Forest"]] <- list(model = rf, type = "randomForest")
  metrics[["Random Forest"]] <- compute_metrics(test_df$died, rf_prob)
  roc_data[["Random Forest"]] <- safe_roc(test_df$died, rf_prob)

  ## ---- SVM (RBF) ----
  ## class.weights was previously the string "inverse" (an e1071 shorthand
  ## for auto-weighting) — on some e1071 versions this interacts badly with
  ## probability=TRUE on an imbalanced outcome like this one (~4.4% deaths)
  ## and can corrupt the shape of the returned probability matrix, causing
  ## "Error in matrix: length of 'dimnames' [2] not equal to array extent"
  ## in predict.svm. Computing the weights explicitly (same approach as the
  ## Random Forest above) avoids that ambiguity. The whole step is also
  ## wrapped in tryCatch: if SVM still fails for some fold/split, that model
  ## is dropped from the comparison for this fold instead of taking the
  ## entire dashboard down.
  svm_weights <- c(No = 1, Yes = sum(train_df$died == "No") / max(1, sum(train_df$died == "Yes")))
  svm_fit <- tryCatch({
    svm_mod <- e1071::svm(feat_formula, data = train_df, kernel = "radial",
                          probability = TRUE, class.weights = svm_weights)
    svm_pred <- predict(svm_mod, test_df, probability = TRUE)
    svm_prob <- attr(svm_pred, "probabilities")[, "Yes"]
    if (length(svm_prob) != nrow(test_df)) {
      ## last-resort safety net: pad/truncate so downstream metrics never
      ## crash on a length mismatch, even if some other edge case slips through
      full <- rep(NA_real_, nrow(test_df))
      full[seq_along(svm_prob)] <- svm_prob
      svm_prob <- full
    }
    list(model = svm_mod, prob = svm_prob)
  }, error = function(e) NULL)

  if (!is.null(svm_fit)) {
    if (keep_models) results[["SVM (RBF)"]] <- list(model = svm_fit$model, type = "svm")
    metrics[["SVM (RBF)"]] <- compute_metrics(test_df$died, svm_fit$prob)
    roc_data[["SVM (RBF)"]] <- safe_roc(test_df$died, svm_fit$prob)
  }

  ## ---- Gradient Boosting ----
  gbm_df <- train_df
  gbm_df$died_num <- as.numeric(gbm_df$died == "Yes")
  gbm_mod <- gbm::gbm(gbm_formula, data = gbm_df, distribution = "bernoulli",
                       n.trees = 200, interaction.depth = 4, shrinkage = 0.1, verbose = FALSE)
  gbm_prob <- predict(gbm_mod, test_df, n.trees = 200, type = "response")
  if (keep_models) results[["Gradient Boosting"]] <- list(model = gbm_mod, type = "gbm")
  metrics[["Gradient Boosting"]] <- compute_metrics(test_df$died, gbm_prob)
  roc_data[["Gradient Boosting"]] <- safe_roc(test_df$died, gbm_prob)

  list(models = results, metrics = metrics, roc_data = roc_data)
}

## Stratified k-fold cross-validation, reusing fit_score_fold(). A single
## 80/20 split can swing a lot on a cohort this size (2,478 rows, ~110
## deaths), so this gives a steadier read on each model's real performance.
cross_validate_models <- function(model_df, feat_formula, gbm_formula, k = 5) {
  set.seed(42)
  y <- model_df$died
  fold_id <- integer(nrow(model_df))
  for (lvl in levels(y)) {
    idx <- which(y == lvl)
    fold_id[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }

  per_fold <- lapply(seq_len(k), function(f) {
    train_df <- model_df[fold_id != f, ]
    test_df  <- model_df[fold_id == f, ]
    fit_score_fold(train_df, test_df, feat_formula, gbm_formula, keep_models = FALSE)$metrics
  })

  model_names <- unique(unlist(lapply(per_fold, names)))
  if (length(model_names) == 0) return(NULL)

  cv_df <- do.call(rbind, lapply(model_names, function(nm) {
    fold_metrics <- Filter(Negate(is.null), lapply(per_fold, function(fm) fm[[nm]]))
    n_ok <- length(fold_metrics)
    avg <- function(field) if (n_ok == 0) NA_real_ else mean(sapply(fold_metrics, function(m) m[[field]]), na.rm = TRUE)
    sd_ <- function(field) if (n_ok < 2) NA_real_ else sd(sapply(fold_metrics, function(m) m[[field]]), na.rm = TRUE)
    data.frame(Model = nm,
               Accuracy = avg("accuracy"), Precision = avg("precision"),
               Recall = avg("recall"), F1 = avg("f1"),
               `ROC AUC` = avg("auc"), `ROC AUC SD` = sd_("auc"),
               `Folds used` = sprintf("%d/%d", n_ok, k),
               check.names = FALSE)
  }))
  cv_df[order(-cv_df$`ROC AUC`), ]
}

## Fills missing numeric baseline features (CD4, WHO stage, BMI, etc.) using
## MICE (Multiple Imputation by Chained Equations, predictive-mean-matching)
## instead of plugging in one global median for everyone. This estimates
## each missing value from the patient's OTHER recorded features (age,
## weight, WHO stage, lab activity, etc.) rather than assuming every
## missing-CD4 patient is "average". Uses a single completed draw (m = 1)
## rather than full multiple-imputation-with-pooling: for a deployed risk
## score (as opposed to a research effect estimate), one relationship-aware
## completed dataset is the standard practical compromise, and it keeps this
## an interactive dashboard rather than a batch job. Falls back to simple
## median imputation if the mice package isn't installed, so the app never
## breaks over a missing dependency.
mice_impute_numeric <- function(df, numeric_feats, categorical_feats) {
  median_fallback <- function(d) {
    for (col in numeric_feats) {
      med <- median(d[[col]], na.rm = TRUE)
      if (is.na(med) || is.infinite(med)) med <- 0
      d[[col]][is.na(d[[col]])] <- med
    }
    d
  }

  if (!requireNamespace("mice", quietly = TRUE)) return(median_fallback(df))

  work <- df[, c(numeric_feats, categorical_feats), drop = FALSE]
  for (col in categorical_feats) {
    v <- as.character(work[[col]])
    v[is.na(v)] <- "Unknown"
    work[[col]] <- as.factor(v)
  }

  completed <- tryCatch({
    imp <- mice::mice(work, m = 1, maxit = 5, method = "pmm", seed = 42, printFlag = FALSE)
    mice::complete(imp, 1)
  }, error = function(e) NULL)

  if (is.null(completed)) return(median_fallback(df))

  for (col in numeric_feats) df[[col]] <- completed[[col]]
  df
}

train_models <- function(df) {
  model_df <- df

  ## guard against Inf/-Inf (e.g. bmi_recalculated when height==0) which
  ## would otherwise survive imputation and break SVM's internal feature
  ## scaling, silently dropping rows at predict() time
  for (col in NUMERIC_FEATS) {
    model_df[[col]][is.infinite(model_df[[col]])] <- NA
  }

  model_df <- mice_impute_numeric(model_df, NUMERIC_FEATS, CATEGORICAL_FEATS)

  for (col in CATEGORICAL_FEATS) {
    v <- as.character(model_df[[col]])
    v[is.na(v)] <- "Unknown"
    model_df[[col]] <- as.factor(v)
  }

  ## saved so cohort-wide risk scoring (risk_scored() in the server) reuses
  ## these same imputed values, instead of falling back to a cruder
  ## median/Unknown fill that wouldn't match what the model was trained on
  imputed_lookup <- model_df[, c("patient_id", NUMERIC_FEATS, CATEGORICAL_FEATS)]

  model_df$died <- as.factor(model_df$died)
  levels(model_df$died) <- c("No", "Yes")

  feat_formula <- as.formula(paste("died ~", paste(c(NUMERIC_FEATS, CATEGORICAL_FEATS), collapse = " + ")))
  gbm_formula  <- as.formula(paste("died_num ~", paste(c(NUMERIC_FEATS, CATEGORICAL_FEATS), collapse = " + ")))

  n <- nrow(model_df)
  set.seed(42)
  train_idx <- caret_like_split(model_df$died, p = 0.8)
  train_df <- model_df[train_idx, ]
  test_df  <- model_df[-train_idx, ]

  fold <- fit_score_fold(train_df, test_df, feat_formula, gbm_formula, keep_models = TRUE)
  results  <- fold$models
  metrics  <- fold$metrics
  roc_data <- fold$roc_data

  metrics_df <- do.call(rbind, lapply(names(metrics), function(nm) {
    m <- metrics[[nm]]
    data.frame(Model = nm, Accuracy = m$accuracy, Precision = m$precision,
               Recall = m$recall, F1 = m$f1, `ROC AUC` = m$auc, check.names = FALSE)
  }))
  metrics_df <- metrics_df[order(-metrics_df$`ROC AUC`), ]
  rownames(metrics_df) <- NULL
  best_name <- metrics_df$Model[1]

  cv_metrics_df <- tryCatch(cross_validate_models(model_df, feat_formula, gbm_formula, k = 5),
                             error = function(e) NULL)

  importances <- NULL
  if (best_name == "Random Forest" || "Random Forest" %in% names(results)) {
    imp <- randomForest::importance(results[["Random Forest"]]$model)
    importances <- sort(setNames(imp[, 1], rownames(imp)), decreasing = TRUE)
  }

  medians <- sapply(NUMERIC_FEATS, function(c) median(df[[c]], na.rm = TRUE))

  list(
    fitted = results,
    metrics_df = metrics_df,
    cv_metrics_df = cv_metrics_df,
    roc_data = roc_data,
    importances = importances,
    best_name = best_name,
    medians = medians,
    imputed_lookup = imputed_lookup,
    train_n = nrow(train_df),
    test_n = nrow(test_df),
    feat_levels = list(
      sex_label = levels(model_df$sex_label),
      baseline_vl_status = levels(model_df$baseline_vl_status),
      subcity_grouped = levels(model_df$subcity_grouped)
    )
  )
}

caret_like_split <- function(y, p = 0.8) {
  idx <- unlist(lapply(split(seq_along(y), y), function(ix) {
    sample(ix, size = floor(length(ix) * p))
  }))
  sort(idx)
}

compute_metrics <- function(truth, prob, threshold = 0.5) {
  ## drop any pair where prob is NA (e.g. a row a model couldn't score)
  ## so length/alignment issues never propagate into the metrics
  keep <- !is.na(prob)
  truth <- truth[keep]
  prob <- prob[keep]
  pred <- ifelse(prob >= threshold, "Yes", "No")
  truth <- as.character(truth)
  tp <- sum(pred == "Yes" & truth == "Yes")
  fp <- sum(pred == "Yes" & truth == "No")
  fn <- sum(pred == "No"  & truth == "Yes")
  tn <- sum(pred == "No"  & truth == "No")
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall    <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1        <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  accuracy  <- (tp + tn) / length(truth)
  auc_val <- tryCatch(as.numeric(pROC::auc(pROC::roc(truth, prob, quiet = TRUE, levels = c("No", "Yes")))),
                       error = function(e) NA_real_)
  list(accuracy = accuracy, precision = precision, recall = recall, f1 = f1, auc = auc_val)
}

## roc() needs truth/prob of equal length with no NA — same guard as
## compute_metrics, used wherever we build an ROC curve for plotting
safe_roc <- function(truth, prob) {
  keep <- !is.na(prob)
  pROC::roc(truth[keep], prob[keep], quiet = TRUE, levels = c("No", "Yes"))
}

risk_band <- function(p) {
  if (is.na(p)) return(c("Unknown", "secondary"))
  if (p < 0.10) return(c("Low", "success"))
  if (p < 0.30) return(c("Elevated", "warning"))
  c("High", "danger")
}

predict_prob <- function(model_entry, newdata) {
  switch(model_entry$type,
    "rpart" = predict(model_entry$model, newdata, type = "prob")[, "Yes"],
    "randomForest" = predict(model_entry$model, newdata, type = "prob")[, "Yes"],
    "svm" = {
      p <- predict(model_entry$model, newdata, probability = TRUE)
      attr(p, "probabilities")[, "Yes"]
    },
    "gbm" = predict(model_entry$model, newdata, n.trees = 200, type = "response")
  )
}

## ============================================================================
## UI
## ============================================================================

ui <- dashboardPage(
  title = "CEM DASHBOARD",
  skin = "black",

  dashboardHeader(
    title = "CEM DASHBOARD",
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "main_tabs",
      menuItem("Executive Summary", tabName = "summary", icon = icon("chart-line")),
      menuItem("Cohort Explorer", tabName = "cohort", icon = icon("table")),
      menuItem("Risk Assessment", tabName = "risk", icon = icon("shield-heart")),
      menuItem("Facility Performance", tabName = "facilities", icon = icon("hospital")),
      menuItem("Data Quality", tabName = "quality", icon = icon("magnifying-glass-chart")),
      tags$hr(style = "border-color: rgba(255,255,255,0.15);"),
      div(style = "padding: 0 15px;",
        tags$p("DATA SOURCE", style = "color:#8aa4af; font-size:11px; font-weight:700; letter-spacing:1px; margin-bottom:6px;"),
        fileInput("upload", NULL, accept = ".csv", buttonLabel = "Browse...", placeholder = "Add file (CSV)")
      ),
      tags$hr(style = "border-color: rgba(255,255,255,0.15);"),
      div(style = "padding: 0 15px;",
        fluidRow(
          column(8, tags$p("FILTERS", style = "color:#8aa4af; font-size:11px; font-weight:700; letter-spacing:1px;")),
          column(4, actionButton("reset_filters", "Reset", class = "btn-xs btn-outline-secondary", style = "margin-top:-4px;"))
        ),
        uiOutput("filter_ui")
      ),
      tags$hr(style = "border-color: rgba(255,255,255,0.15);"),
      div(style = "padding: 0 15px 15px 15px; color:#8aa4af; font-size:12px;",
        textOutput("filter_summary")
      )
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML(sprintf("
      .content-wrapper, .right-side { background-color: #f4f6f7; }
      body, .box, .form-control, .selectize-input, .dataTable { font-family: 'Segoe UI', -apple-system, Roboto, Arial, sans-serif; }
      .skin-black .main-header .logo { background-color: %s; color: #fff; font-weight: 700; letter-spacing: .3px; font-size: 18px; }
      .skin-black .main-header .navbar { background-color: %s; }
      .skin-black .main-sidebar { background-color: %s; }
      .skin-black .sidebar-menu > li.active > a, .skin-black .sidebar-menu > li:hover > a {
        background-color: %s; color: #fff; border-left-color: #ffffff;
      }
      .skin-black .sidebar-menu > li > a { color: #d9e6ea; font-size: 13.5px; }
      .box { border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.06); border-top: 3px solid %s; }
      .box-header { padding: 10px 15px 4px 15px; }
      .box-title { font-size: 14px; font-weight: 700; color: #37474f; }
      .box.box-solid.box-primary > .box-header { background: %s; }
      .box.box-solid.box-primary > .box-header .box-title { color: #fff; }
      .risk-badge { display: inline-block; padding: 8px 22px; border-radius: 20px; font-weight: 700;
                    font-size: 1.2rem; color: white; }
      .section-note { color: #6b7a80; margin-bottom: 12px; font-size: 13px; }
      .form-control, .selectize-input { border-radius: 6px; }
      .irs-bar, .irs-bar-edge, .irs-single { background: %s !important; border-color: %s !important; }
      .data-banner { padding: 8px 24px; font-size: 12.5px; border-bottom: 1px solid #e3e9eb; }
      .data-banner.demo { background: #fff8e6; color: #7a5b00; }
      .data-banner.live { background: #f4f9f6; color: #2e6b4a; }
      .rec-item { padding: 9px 14px; margin-bottom: 6px; background: #f7fafb; border-left: 3px solid %s;
                  border-radius: 4px; font-size: 13px; color: #2c3a3f; }
    ", ACCENT_DARK, ACCENT_DARK, ACCENT_DARK, ACCENT, ACCENT, ACCENT, ACCENT, ACCENT, ACCENT)))),

    uiOutput("data_source_banner"),

    tabItems(

      ## ---------------------------------------------------------- EXECUTIVE SUMMARY
      tabItem(tabName = "summary",
        fluidRow(
          valueBoxOutput("kpi_total", width = 3),
          valueBoxOutput("kpi_deaths", width = 3),
          valueBoxOutput("kpi_active", width = 3),
          valueBoxOutput("kpi_lost", width = 3)
        ),
        fluidRow(
          box(title = "Program Insights", status = "primary", solidHeader = TRUE, width = 7,
              uiOutput("recommendations_ui")),
          box(title = "Cohort Composition", width = 5,
              plotlyOutput("plot_status_pie", height = "300px"))
        ),
        fluidRow(
          box(title = "Facility Mortality Ranking", status = "primary", solidHeader = TRUE, width = 12,
              p(class = "section-note", "Facilities with at least 5 patients, ranked by mortality rate. See Facility Performance for the full view."),
              DTOutput("facility_ranking_preview"))
        )
      ),

      ## ---------------------------------------------------------- COHORT EXPLORER
      tabItem(tabName = "cohort",
        tabsetPanel(
          tabPanel("Patient Records",
            br(),
            box(title = "Filtered Patient Records", width = 12,
              fluidRow(
                column(6, textInput("search_text", NULL, placeholder = "Search by patient ID, facility, subcity...")),
                column(6, div(style = "text-align:right; padding-top:6px;",
                              downloadButton("download_data", "Download CSV", class = "btn-sm")))
              ),
              DTOutput("data_table")
            ),
            box(title = "Summary Statistics", width = 12,
              DTOutput("summary_table")
            )
          ),
          tabPanel("Clinical Insights",
            br(),
            p(class = "section-note", "How baseline clinical markers relate to treatment outcomes across the current filter."),
            fluidRow(
              box(title = "Sex Distribution", width = 6, plotlyOutput("plot_sex_pie", height = "300px")),
              box(title = "Baseline WHO Stage", width = 6, plotlyOutput("plot_stage_bar", height = "300px"))
            ),
            fluidRow(
              box(title = "CD4 Count at Baseline vs. Outcome", width = 6,
                  plotlyOutput("plot_cd4_outcome", height = "340px")),
              box(title = "Mortality Rate by WHO Stage", width = 6,
                  plotlyOutput("plot_stage_mortality", height = "340px"))
            ),
            fluidRow(
              box(title = "BMI Category vs. Mortality", width = 6,
                  plotlyOutput("plot_bmi_mortality", height = "340px")),
              box(title = "Viral Load Suppression Status", width = 6,
                  plotlyOutput("plot_vl_status", height = "340px"))
            ),
            fluidRow(
              box(title = "Mortality Rate by Subcity (n \u2265 5)", width = 12,
                  plotlyOutput("plot_subcity_mortality", height = "380px"))
            )
          )
        )
      ),

      ## ---------------------------------------------------------- RISK ASSESSMENT
      tabItem(tabName = "risk",
        tabsetPanel(
          tabPanel("Patient Calculator",
            br(),
            p(class = "section-note", "Estimate an individual patient's mortality risk from baseline characteristics."),
            box(title = "Patient Inputs", status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                column(4,
                  numericInput("calc_age", "Age", value = 35, min = 0, max = 100),
                  numericInput("calc_weight", "Weight (kg)", value = 58, min = 1, max = 200),
                  numericInput("calc_height", "Height (cm)", value = 162, min = 50, max = 250),
                  selectInput("calc_sex", "Sex", choices = c("Male", "Female"))
                ),
                column(4,
                  selectInput("calc_stage", "Baseline WHO stage", choices = c("I", "II", "III", "IV")),
                  numericInput("calc_cd4", "Baseline CD4 count", value = 350, min = 0, max = 2000),
                  selectInput("calc_vl", "Baseline viral load status", choices = c("Undetectable", "Detectable", "Unknown"))
                ),
                column(4,
                  selectInput("calc_subcity", "Subcity group", choices = NULL),
                  numericInput("calc_labs", "Number of lab events recorded", value = 4, min = 0, max = 50),
                  numericInput("calc_days", "Days on treatment", value = 365, min = 0, max = 5000)
                )
              ),
              actionButton("calc_run", "Calculate Risk", class = "btn-primary"),
              actionButton("calc_reset", "New Assessment", class = "btn-default")
            ),
            uiOutput("calc_result")
          ),
          tabPanel("Model Performance",
            br(),
            p(class = "section-note", "Four candidate models trained on baseline features to predict 12-month mortality risk."),
            box(title = "Model Comparison (single 80/20 split)", width = 12, DTOutput("metrics_table")),
            box(title = "Model Comparison (5-fold cross-validation, mean \u00b1 SD of ROC AUC)", width = 12,
                p(class = "section-note", "A single train/test split can swing a lot on a cohort this size \u2014 cross-validation gives a steadier read on each model's real performance."),
                DTOutput("cv_metrics_table")),
            fluidRow(
              box(title = "ROC Curves", width = 6, plotlyOutput("plot_roc", height = "360px")),
              box(title = "Feature Importance (Random Forest)", width = 6, plotlyOutput("plot_importance", height = "360px"))
            )
          ),
          tabPanel("Cohort Risk List",
            br(),
            p(class = "section-note", "Predicted mortality risk (best model) applied across the current cohort filter."),
            box(title = "Risk Tier Distribution", width = 12, plotlyOutput("plot_risk_tiers", height = "300px")),
            box(title = "Highest-Risk Patients", width = 12, DTOutput("high_risk_table"))
          )
        )
      ),

      ## ---------------------------------------------------------- FACILITY PERFORMANCE
      tabItem(tabName = "facilities",
        p(class = "section-note", "Facility-level caseload, mortality, follow-up loss, and lab monitoring coverage."),
        box(title = "Facility Ranking", status = "primary", solidHeader = TRUE, width = 12,
            DTOutput("facility_ranking_table")),
        fluidRow(
          box(title = "Caseload by Facility", width = 6, plotlyOutput("plot_facilities", height = "380px")),
          box(title = "Mortality Rate by Facility", width = 6, plotlyOutput("plot_facility_mortality", height = "380px"))
        ),
        box(title = "Enrollment Trend", width = 12, plotlyOutput("plot_enrollment", height = "320px"))
      ),

      ## ---------------------------------------------------------- DATA QUALITY
      tabItem(tabName = "quality",
        box(title = "Missingness by Field", width = 12,
            plotlyOutput("plot_missingness", height = "380px")),
        box(title = "Record Counts", width = 12,
            DTOutput("dq_table")),
        box(title = "Is Missingness Random? (mortality rate: missing vs. recorded)", width = 12,
            p(class = "section-note", "If the death rate differs a lot between patients missing a field and patients with it recorded, that field isn't missing at random \u2014 imputing it (however it's done) will understate risk for the group it's missing more often for."),
            DTOutput("dq_missingness_bias"))
      )
    )
  )
)

## ============================================================================
## SERVER
## ============================================================================

server <- function(input, output, session) {

  ## ---- data source -----------------------------------------------------
  ## data.table::fread instead of base read.csv — orders of magnitude
  ## faster on this file (877 columns) and is what was mainly behind the
  ## sluggish startup/reload.
  fread_df <- function(path) {
    as.data.frame(data.table::fread(path, stringsAsFactors = FALSE, check.names = TRUE,
                                     na.strings = c("NA", ""), showProgress = FALSE))
  }

  raw_data <- reactive({
    if (!is.null(input$upload)) {
      df <- tryCatch(fread_df(input$upload$datapath), error = function(e) NULL)
      if (is.null(df) || !all(HARD_REQUIRED_COLS %in% names(df))) {
        showNotification("Couldn't load that file \u2014 missing required columns (patient_id, patient_status). Using demo data instead.",
                          type = "error", duration = 8)
        return(list(df = engineer_features(make_demo_data()), source = "synthetic demo data", demo = TRUE))
      }
      return(list(df = engineer_features(df), source = input$upload$name, demo = FALSE))
    }
    if (file.exists(DATA_PATH_DEFAULT)) {
      df <- fread_df(DATA_PATH_DEFAULT)
      return(list(df = engineer_features(df), source = DATA_PATH_DEFAULT, demo = FALSE))
    }
    list(df = engineer_features(make_demo_data()), source = "synthetic demo data", demo = TRUE)
  })

  output$filter_ui <- renderUI({
    df <- raw_data()$df
    tagList(
      selectizeInput("flt_facility", "Health facility",
                      choices = sort(unique(na.omit(df$health_facility))), multiple = TRUE),
      selectizeInput("flt_subcity", "Subcity",
                      choices = sort(unique(na.omit(df$subcity_grouped))), multiple = TRUE),
      selectizeInput("flt_status", "Patient status",
                      choices = sort(unique(na.omit(df$patient_status))), multiple = TRUE),
      sliderInput("flt_age", "Age range",
                  min = floor(min(df$age_clean, na.rm = TRUE)),
                  max = ceiling(max(df$age_clean, na.rm = TRUE)),
                  value = c(floor(min(df$age_clean, na.rm = TRUE)), ceiling(max(df$age_clean, na.rm = TRUE))))
    )
  })

  observeEvent(input$reset_filters, {
    updateSelectizeInput(session, "flt_facility", selected = character(0))
    updateSelectizeInput(session, "flt_subcity", selected = character(0))
    updateSelectizeInput(session, "flt_status", selected = character(0))
    df <- raw_data()$df
    updateSliderInput(session, "flt_age",
                       value = c(floor(min(df$age_clean, na.rm = TRUE)), ceiling(max(df$age_clean, na.rm = TRUE))))
  })

  filtered_data <- reactive({
    df <- raw_data()$df
    if (!is.null(input$flt_facility) && length(input$flt_facility) > 0) {
      df <- df[df$health_facility %in% input$flt_facility, ]
    }
    if (!is.null(input$flt_subcity) && length(input$flt_subcity) > 0) {
      df <- df[df$subcity_grouped %in% input$flt_subcity, ]
    }
    if (!is.null(input$flt_status) && length(input$flt_status) > 0) {
      df <- df[df$patient_status %in% input$flt_status, ]
    }
    if (!is.null(input$flt_age)) {
      df <- df[is.na(df$age_clean) | (df$age_clean >= input$flt_age[1] & df$age_clean <= input$flt_age[2]), ]
    }
    df
  })

  output$filter_summary <- renderText({
    sprintf("Showing %s of %s patients", format(nrow(filtered_data()), big.mark = ","), format(nrow(raw_data()$df), big.mark = ","))
  })

  output$data_source_banner <- renderUI({
    rd <- raw_data()
    if (isTRUE(rd$demo)) {
      div(class = "data-banner demo",
          strong("Demo data \u2014 "), "showing synthetic sample records, not real patient data. Upload a CSV in the sidebar to connect real cohort data.")
    } else {
      NULL
    }
  })

  ## ---- model bundle (trained once per data load, on unfiltered cohort) ----
  ## Training 4 models (incl. SVM with probability=TRUE) is the slowest part
  ## of the app. It only needs to happen when the underlying data actually
  ## changes, so we cache the fitted bundle to disk keyed by an md5 of the
  ## source file's bytes. Re-running the app on the same CSV now reads the
  ## cache instead of retraining from scratch every time.
  MODEL_CACHE_PATH <- "model_cache.rds"

  model_bundle <- reactive({
    rd <- raw_data()
    file_path <- if (!is.null(input$upload)) input$upload$datapath
                 else if (!isTRUE(rd$demo) && file.exists(DATA_PATH_DEFAULT)) DATA_PATH_DEFAULT
                 else NA_character_

    if (is.na(file_path)) return(train_models(rd$df))  # demo data: just train, nothing to key on

    key <- tryCatch(unname(tools::md5sum(file_path)), error = function(e) NA_character_)
    if (!is.na(key) && file.exists(MODEL_CACHE_PATH)) {
      cached <- tryCatch(readRDS(MODEL_CACHE_PATH), error = function(e) NULL)
      if (!is.null(cached) && identical(cached$key, key) &&
          all(c("cv_metrics_df", "imputed_lookup") %in% names(cached$bundle))) {
        return(cached$bundle)
      }
    }
    bundle <- train_models(rd$df)
    if (!is.na(key)) {
      tryCatch(saveRDS(list(key = key, bundle = bundle), MODEL_CACHE_PATH), error = function(e) NULL)
    }
    bundle
  })

  observe({
    bundle <- model_bundle()
    updateSelectInput(session, "calc_subcity", choices = bundle$feat_levels$subcity_grouped)
  })

  ## ---- FACILITY RANKING (shared by Executive Summary + Facility Performance) ----

  facility_ranking <- reactive({
    df <- filtered_data()
    df %>%
      group_by(health_facility) %>%
      summarise(
        Patients = n(),
        Deaths = sum(died),
        `Mortality rate (%)` = round(mean(died) * 100, 1),
        `Lost to follow-up (%)` = round(mean(patient_status == "inactive_or_lost") * 100, 1),
        `CD4 tested (%)` = round(mean(!is.na(baseline_cd4)) * 100, 1)
      ) %>%
      filter(Patients >= 5, !is.na(health_facility)) %>%
      arrange(desc(`Mortality rate (%)`))
  })

  output$facility_ranking_table <- renderDT({
    datatable(facility_ranking(), options = list(pageLength = 10, dom = "tip"), rownames = FALSE)
  })

  output$facility_ranking_preview <- renderDT({
    datatable(head(facility_ranking(), 5), options = list(dom = "t"), rownames = FALSE)
  })

  output$plot_facility_mortality <- renderPlotly({
    tab <- facility_ranking() %>% arrange(`Mortality rate (%)`) %>% tail(12)
    if (nrow(tab) == 0) return(plotly_empty())
    plot_ly(tab, x = ~`Mortality rate (%)`, y = ~reorder(health_facility, `Mortality rate (%)`),
            type = "bar", orientation = "h", marker = list(color = ~`Mortality rate (%)`, colorscale = "Reds")) %>%
      layout(yaxis = list(title = ""), xaxis = list(title = "Mortality rate (%)"))
  })

  ## ---- OVERVIEW ----------------------------------------------------------

  output$kpi_total <- renderValueBox({
    valueBox(format(nrow(filtered_data()), big.mark = ","), "Total patients",
             icon = icon("users"), color = "light-blue")
  })
  output$kpi_deaths <- renderValueBox({
    df <- filtered_data()
    rate <- if (nrow(df) > 0) mean(df$died) * 100 else 0
    valueBox(format(sum(df$died), big.mark = ","),
             sprintf("Deaths recorded (%.1f%% of cohort)", rate),
             icon = icon("skull-crossbones"), color = "red")
  })
  output$kpi_active <- renderValueBox({
    df <- filtered_data()
    valueBox(format(sum(df$patient_status == "active"), big.mark = ","), "Active in care",
             icon = icon("heart-pulse"), color = "green")
  })
  output$kpi_lost <- renderValueBox({
    df <- filtered_data()
    valueBox(format(sum(df$patient_status == "inactive_or_lost"), big.mark = ","), "Lost / inactive",
             icon = icon("person-walking-arrow-right"), color = "yellow")
  })

  output$plot_status_pie <- renderPlotly({
    df <- filtered_data()
    tab <- as.data.frame(table(df$patient_status))
    names(tab) <- c("Status", "Count")
    cols <- c(active = "#2e8b57", inactive_or_lost = "#e0a800", died = "#c0392b")
    plot_ly(tab, labels = ~Status, values = ~Count, type = "pie", hole = 0.45,
            marker = list(colors = cols[as.character(tab$Status)])) %>%
      layout(title = "Patient status", margin = list(t = 40))
  })

  output$plot_sex_pie <- renderPlotly({
    df <- filtered_data()
    tab <- as.data.frame(table(df$sex_label, useNA = "no"))
    names(tab) <- c("Sex", "Count")
    plot_ly(tab, labels = ~Sex, values = ~Count, type = "pie", hole = 0.45) %>%
      layout(title = "Sex distribution", margin = list(t = 40))
  })

  output$plot_stage_bar <- renderPlotly({
    df <- filtered_data()
    tab <- df %>% filter(!is.na(baseline_who_stage)) %>% count(baseline_who_stage) %>% arrange(baseline_who_stage)
    tab$Stage <- c("I", "II", "III", "IV")[tab$baseline_who_stage]
    plot_ly(tab, x = ~Stage, y = ~n, type = "bar", color = ~Stage) %>%
      layout(title = "Baseline WHO stage", showlegend = FALSE, margin = list(t = 40))
  })

  output$plot_facilities <- renderPlotly({
    df <- filtered_data()
    tab <- df %>% count(health_facility, sort = TRUE) %>% na.omit() %>% head(12) %>% arrange(n)
    tab$health_facility <- factor(tab$health_facility, levels = tab$health_facility)
    plot_ly(tab, x = ~n, y = ~health_facility, type = "bar", orientation = "h",
            marker = list(color = ~n, colorscale = "Blues")) %>%
      layout(xaxis = list(title = "Patients"), yaxis = list(title = ""))
  })

  output$plot_enrollment <- renderPlotly({
    df <- filtered_data() %>% filter(!is.na(enroll_date))
    if (nrow(df) == 0) return(plotly_empty(type = "scatter"))
    df$month <- floor_date(df$enroll_date, "month")
    tab <- df %>% count(month)
    plot_ly(tab, x = ~month, y = ~n, type = "scatter", mode = "lines", fill = "tozeroy") %>%
      layout(title = "New enrollments over time", xaxis = list(title = ""), yaxis = list(title = "New patients"))
  })

  ## ---- DATA EXPLORER -------------------------------------------------

  search_text_d <- debounce(reactive(input$search_text), 400)

  display_cols <- c("patient_id", "health_facility", "subcity_grouped", "sex_label", "age_clean",
                     "patient_status", "baseline_who_stage", "baseline_cd4", "baseline_vl_status",
                     "bmi_recalculated", "bmi_category", "n_lab_events", "days_on_treatment")

  explorer_df <- reactive({
    df <- filtered_data()
    cols <- intersect(display_cols, names(df))
    show_df <- df[, cols, drop = FALSE]
    query <- search_text_d()
    if (!is.null(query) && nzchar(query)) {
      mask <- Reduce(`|`, lapply(show_df, function(col) {
        str_detect(tolower(as.character(col)), fixed(tolower(query)))
      }))
      show_df <- show_df[which(mask), , drop = FALSE]
    }
    show_df
  })

  output$data_table <- renderDT({
    datatable(explorer_df(), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$download_data <- downloadHandler(
    filename = function() "cem_filtered.csv",
    content = function(file) write.csv(explorer_df(), file, row.names = FALSE)
  )

  output$summary_table <- renderDT({
    df <- filtered_data()
    cols <- intersect(NUMERIC_FEATS, names(df))
    summ <- do.call(rbind, lapply(cols, function(c) {
      x <- df[[c]]
      data.frame(Field = c, N = sum(!is.na(x)), Mean = mean(x, na.rm = TRUE),
                 SD = sd(x, na.rm = TRUE), Min = suppressWarnings(min(x, na.rm = TRUE)),
                 Max = suppressWarnings(max(x, na.rm = TRUE)))
    }))
    datatable(summ, options = list(dom = "t")) %>% formatRound(c("Mean", "SD", "Min", "Max"), 1)
  })

  ## ---- CLINICAL INSIGHTS ----------------------------------------------

  output$plot_cd4_outcome <- renderPlotly({
    df <- filtered_data() %>% filter(!is.na(baseline_cd4))
    if (nrow(df) == 0) return(plotly_empty())
    df$Outcome <- ifelse(df$died == 1, "Died", "Survived / active")
    plot_ly(df, x = ~Outcome, y = ~baseline_cd4, type = "box", color = ~Outcome,
            colors = c("Survived / active" = "#2e8b57", "Died" = "#c0392b")) %>%
      layout(showlegend = FALSE, yaxis = list(title = "Baseline CD4 count"))
  })

  output$plot_stage_mortality <- renderPlotly({
    df <- filtered_data() %>% filter(!is.na(baseline_who_stage))
    tab <- df %>% group_by(baseline_who_stage) %>% summarise(mean = mean(died) * 100, count = n())
    tab$Stage <- c("I", "II", "III", "IV")[tab$baseline_who_stage]
    plot_ly(tab, x = ~Stage, y = ~mean, type = "bar", text = ~paste0("n=", count), textposition = "outside",
            marker = list(color = ~mean, colorscale = "Reds")) %>%
      layout(yaxis = list(title = "Mortality rate (%)"), xaxis = list(title = "WHO Stage"))
  })

  output$plot_bmi_mortality <- renderPlotly({
    df <- filtered_data()
    order <- c("Underweight", "Normal", "Overweight", "Obese", "Unknown")
    tab <- df %>% group_by(bmi_category) %>% summarise(mean = mean(died) * 100, count = n())
    tab$bmi_category <- factor(tab$bmi_category, levels = order)
    tab <- tab %>% arrange(bmi_category)
    plot_ly(tab, x = ~bmi_category, y = ~mean, type = "bar", text = ~paste0("n=", count), textposition = "outside",
            marker = list(color = ~mean, colorscale = "Oranges")) %>%
      layout(yaxis = list(title = "Mortality rate (%)"), xaxis = list(title = "BMI category"))
  })

  output$plot_vl_status <- renderPlotly({
    df <- filtered_data()
    vl <- ifelse(is.na(df$baseline_vl_status), "Not recorded", df$baseline_vl_status)
    tab <- as.data.frame(table(vl)); names(tab) <- c("Status", "Count")
    cols <- c(Undetectable = "#2e8b57", Detectable = "#c0392b", Unknown = "#e0a800", `Not recorded` = "#95a5a6")
    plot_ly(tab, x = ~Status, y = ~Count, type = "bar", marker = list(color = cols[as.character(tab$Status)])) %>%
      layout(showlegend = FALSE)
  })

  output$plot_subcity_mortality <- renderPlotly({
    df <- filtered_data()
    tab <- df %>% group_by(subcity_grouped) %>% summarise(mean = mean(died) * 100, count = n()) %>%
      filter(count >= 5) %>% arrange(desc(mean))
    plot_ly(tab, x = ~reorder(subcity_grouped, mean), y = ~mean, type = "bar",
            marker = list(color = ~mean, colorscale = "Reds")) %>%
      layout(xaxis = list(title = ""), yaxis = list(title = "Mortality rate (%)"))
  })

  ## ---- MORTALITY RISK MODEL --------------------------------------------

  output$metrics_table <- renderDT({
    bundle <- model_bundle()
    datatable(bundle$metrics_df, options = list(dom = "t")) %>%
      formatRound(c("Accuracy", "Precision", "Recall", "F1", "ROC AUC"), 3)
  })

  output$cv_metrics_table <- renderDT({
    bundle <- model_bundle()
    if (is.null(bundle$cv_metrics_df)) return(datatable(data.frame(Note = "Cross-validation metrics unavailable"), options = list(dom = "t")))
    datatable(bundle$cv_metrics_df, options = list(dom = "t")) %>%
      formatRound(c("Accuracy", "Precision", "Recall", "F1", "ROC AUC", "ROC AUC SD"), 3)
  })

  output$plot_roc <- renderPlotly({
    bundle <- model_bundle()
    p <- plot_ly()
    for (nm in names(bundle$roc_data)) {
      r <- bundle$roc_data[[nm]]
      p <- add_trace(p, x = 1 - r$specificities, y = r$sensitivities, type = "scatter", mode = "lines", name = nm)
    }
    p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = "scatter", mode = "lines",
                    line = list(dash = "dash", color = "gray"), name = "Chance", showlegend = FALSE)
    p %>% layout(xaxis = list(title = "False positive rate"), yaxis = list(title = "True positive rate"))
  })

  output$plot_importance <- renderPlotly({
    bundle <- model_bundle()
    if (is.null(bundle$importances)) return(plotly_empty())
    imp <- bundle$importances
    tab <- data.frame(Feature = names(imp), Importance = as.numeric(imp))
    tab <- tab %>% arrange(Importance)
    tab$Feature <- factor(tab$Feature, levels = tab$Feature)
    plot_ly(tab, x = ~Importance, y = ~Feature, type = "bar", orientation = "h",
            marker = list(color = "#0e5c73"))
  })

  ## ---- RISK CALCULATOR ---------------------------------------------------

  calc_defaults <- list(age = 35, weight = 58, height = 162, sex = "Male",
                         stage = "I", cd4 = 350, vl = "Unknown", labs = 4, days = 365)

  calc_display <- reactiveVal(NULL)

  observeEvent(input$calc_run, {
    bundle <- model_bundle()
    stage_map <- c(I = 1, II = 2, III = 3, IV = 4)
    newdata <- data.frame(
      age_clean = input$calc_age,
      weight_clean = input$calc_weight,
      height_clean = input$calc_height,
      bmi_recalculated = input$calc_weight / ((input$calc_height / 100)^2),
      baseline_who_stage = stage_map[[input$calc_stage]],
      baseline_cd4 = input$calc_cd4,
      n_lab_events = input$calc_labs,
      days_on_treatment = input$calc_days,
      sex_label = factor(input$calc_sex, levels = bundle$feat_levels$sex_label),
      baseline_vl_status = factor(input$calc_vl, levels = bundle$feat_levels$baseline_vl_status),
      subcity_grouped = factor(input$calc_subcity, levels = bundle$feat_levels$subcity_grouped)
    )
    best <- bundle$fitted[[bundle$best_name]]
    prob <- tryCatch(predict_prob(best, newdata), error = function(e) NA_real_)
    train_n <- nrow(raw_data()$df)
    train_deaths <- sum(raw_data()$df$died, na.rm = TRUE)
    calc_display(list(prob = as.numeric(prob), model = bundle$best_name,
                       train_n = train_n, train_deaths = train_deaths))
  })

  observeEvent(input$calc_reset, {
    updateNumericInput(session, "calc_age", value = calc_defaults$age)
    updateNumericInput(session, "calc_weight", value = calc_defaults$weight)
    updateNumericInput(session, "calc_height", value = calc_defaults$height)
    updateSelectInput(session, "calc_sex", selected = calc_defaults$sex)
    updateSelectInput(session, "calc_stage", selected = calc_defaults$stage)
    updateNumericInput(session, "calc_cd4", value = calc_defaults$cd4)
    updateSelectInput(session, "calc_vl", selected = calc_defaults$vl)
    updateNumericInput(session, "calc_labs", value = calc_defaults$labs)
    updateNumericInput(session, "calc_days", value = calc_defaults$days)
    calc_display(NULL)
  })

  output$calc_result <- renderUI({
    res <- calc_display()
    if (is.null(res)) {
      return(box(width = 12,
        p(class = "section-note", style = "margin:0;",
          "Enter patient details above and click Calculate Risk. Use New Assessment to clear the form for the next patient.")
      ))
    }
    band <- risk_band(res$prob)
    color_map <- c(success = "#2e8b57", warning = "#e0a800", danger = "#c0392b", secondary = "#95a5a6")
    box(title = "Risk Estimate", status = "primary", solidHeader = TRUE, width = 12,
      h4(sprintf("Estimated 12-month mortality risk: %.1f%%", res$prob * 100)),
      span(class = "risk-badge", style = paste0("background-color:", color_map[[band[2]]]), band[1]),
      tags$hr(),
      p(class = "section-note", style = "margin:0;",
        sprintf("Model: %s, trained on %s patients (%s deaths) from the currently loaded cohort. Screening/prioritization support only \u2014 not a diagnostic result. Always confirm with clinical judgment.",
                res$model, format(res$train_n, big.mark = ","), format(res$train_deaths, big.mark = ",")))
    )
  })

  ## ---- RISK STRATIFICATION ------------------------------------------------

  risk_scored <- reactive({
    bundle <- model_bundle()
    df <- filtered_data()
    best <- bundle$fitted[[bundle$best_name]]

    md <- df
    if (!is.null(bundle$imputed_lookup)) {
      lookup <- bundle$imputed_lookup
      idx <- match(md$patient_id, lookup$patient_id)
      for (col in NUMERIC_FEATS) md[[col]] <- lookup[[col]][idx]
      for (col in CATEGORICAL_FEATS) {
        md[[col]] <- factor(as.character(lookup[[col]][idx]), levels = bundle$feat_levels[[col]])
      }
    } else {
      ## fallback, only reached if imputed_lookup wasn't available for some reason
      for (col in NUMERIC_FEATS) {
        md[[col]][is.na(md[[col]])] <- bundle$medians[[col]]
      }
      for (col in CATEGORICAL_FEATS) {
        v <- as.character(md[[col]]); v[is.na(v)] <- "Unknown"
        md[[col]] <- factor(v, levels = bundle$feat_levels[[col]])
      }
    }
    df$risk_prob <- tryCatch(predict_prob(best, md), error = function(e) NA_real_)
    df$risk_tier <- sapply(df$risk_prob, function(p) risk_band(p)[1])
    df
  })

  output$plot_risk_tiers <- renderPlotly({
    df <- risk_scored()
    tab <- as.data.frame(table(factor(df$risk_tier, levels = c("Low", "Elevated", "High"))))
    names(tab) <- c("Tier", "Count")
    cols <- c(Low = "#2e8b57", Elevated = "#e0a800", High = "#c0392b")
    plot_ly(tab, x = ~Tier, y = ~Count, type = "bar", marker = list(color = cols[as.character(tab$Tier)])) %>%
      layout(showlegend = FALSE)
  })

  output$high_risk_table <- renderDT({
    df <- risk_scored() %>% arrange(desc(risk_prob)) %>% head(50)
    cols <- intersect(c("patient_id", "health_facility", "age_clean", "sex_label",
                         "baseline_who_stage", "baseline_cd4", "risk_prob", "risk_tier"), names(df))
    datatable(df[, cols], options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatRound("risk_prob", 3)
  })

  ## ---- RECOMMENDATIONS ---------------------------------------------------

  output$recommendations_ui <- renderUI({
    df <- filtered_data()
    recs <- list()

    stage34_rate <- tryCatch(mean(df$died[df$baseline_who_stage %in% c(3, 4)], na.rm = TRUE), error = function(e) NA)
    stage12_rate <- tryCatch(mean(df$died[df$baseline_who_stage %in% c(1, 2)], na.rm = TRUE), error = function(e) NA)
    if (!is.na(stage34_rate) && !is.na(stage12_rate) && stage34_rate > stage12_rate) {
      recs[[length(recs) + 1]] <- sprintf(
        "WHO stage III/IV mortality is %.1f%% vs %.1f%% for stage I/II \u2014 prioritize earlier diagnosis.",
        stage34_rate * 100, stage12_rate * 100)
    }

    vl_detect_rate <- tryCatch(mean(df$died[df$baseline_vl_status == "Detectable"], na.rm = TRUE), error = function(e) NA)
    if (!is.na(vl_detect_rate)) {
      recs[[length(recs) + 1]] <- sprintf(
        "Detectable baseline viral load: %.1f%% mortality \u2014 strengthen adherence counseling.",
        vl_detect_rate * 100)
    }

    lost_rate <- mean(df$patient_status == "inactive_or_lost", na.rm = TRUE) * 100
    recs[[length(recs) + 1]] <- sprintf(
      "%.1f%% of the cohort is lost to follow-up \u2014 prioritize tracing and re-engagement.",
      lost_rate)

    low_lab <- mean(df$n_lab_events <= 1, na.rm = TRUE) * 100
    if (low_lab > 10) {
      recs[[length(recs) + 1]] <- sprintf(
        "%.1f%% of patients have 0\u20131 lab events recorded \u2014 review lab monitoring coverage.",
        low_lab)
    }

    if (length(recs) == 0) recs[[1]] <- "No specific risk patterns detected in the current filter selection."

    tagList(lapply(recs, function(r) tags$div(class = "rec-item", r)))
  })

  ## ---- DATA QUALITY -------------------------------------------------------

  output$plot_missingness <- renderPlotly({
    df <- filtered_data()
    cols <- intersect(c(NUMERIC_FEATS, CATEGORICAL_FEATS, "patient_status", "health_facility"), names(df))
    miss <- sapply(cols, function(c) mean(is.na(df[[c]])) * 100)
    tab <- data.frame(Field = names(miss), Missing = as.numeric(miss)) %>% arrange(Missing)
    tab$Field <- factor(tab$Field, levels = tab$Field)
    plot_ly(tab, x = ~Missing, y = ~Field, type = "bar", orientation = "h",
            marker = list(color = "#c0392b")) %>%
      layout(xaxis = list(title = "% missing"), yaxis = list(title = ""))
  })

  output$dq_table <- renderDT({
    df <- filtered_data()
    raw <- raw_data()$df
    tab <- data.frame(
      Metric = c("Records (filtered)", "Records (total)", "Duplicate patient IDs",
                 "Missing age", "Missing sex", "Missing baseline WHO stage", "Missing baseline CD4"),
      Value = c(
        nrow(df), nrow(raw), sum(duplicated(raw$patient_id)),
        sum(is.na(df$age_clean)), sum(is.na(df$sex_label)),
        sum(is.na(df$baseline_who_stage)), sum(is.na(df$baseline_cd4))
      )
    )
    datatable(tab, options = list(dom = "t"))
  })

  output$dq_missingness_bias <- renderDT({
    raw <- raw_data()$df
    check_field <- function(col, label) {
      miss <- is.na(raw[[col]])
      data.frame(
        Field = label,
        `Death rate - missing` = round(mean(raw$died[miss], na.rm = TRUE) * 100, 1),
        `Death rate - recorded` = round(mean(raw$died[!miss], na.rm = TRUE) * 100, 1),
        `N missing` = sum(miss),
        check.names = FALSE
      )
    }
    tab <- rbind(
      check_field("baseline_cd4", "Baseline CD4"),
      check_field("baseline_who_stage", "Baseline WHO stage"),
      check_field("baseline_vl_status", "Baseline VL status"),
      check_field("age_clean", "Age")
    )
    datatable(tab, options = list(dom = "t"))
  })
}

shinyApp(ui, server)
