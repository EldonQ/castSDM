#' Evaluate Fitted Models
#'
#' Computes discrimination, threshold-dependent, calibration, and suitability
#' metrics for fitted models on test data. Thresholds are selected on the
#' training reference stored in the fit and then evaluated on the test data.
#'
#' @param fit A [cast_fit] object.
#' @param test_data A `data.frame` with `presence` column and the same
#'   predictor variables used in fitting.
#' @param response Character. Response column name. Default `"presence"`.
#' @param threshold_method Binary threshold rule selected using training
#'   predictions. See [cast_threshold()]. Default `"max_tss"`.
#'
#' @return A `cast_eval` object (S3 class).
#'
#' @details
#' - **AUC**: Area Under the Receiver Operating Characteristic curve,
#'   measuring discrimination ability. Computed via [pROC::roc()] with
#'   `direction = "<"` fixed, so a worse-than-random model correctly reports
#'   AUC < 0.5 instead of being mirrored by `direction = "auto"`.
#' - **TSS**: evaluated at a threshold selected on the training fit, not tuned
#'   on the held-out labels.
#' - **PR-AUC/SEDI**: complementary rare-event discrimination metrics; PR-AUC's
#'   baseline depends on the sampled prevalence.
#' - **Boyce**: moving-window background-based estimate; `cbi_mean` is retained
#'   as a legacy fixed-bin diagnostic.
#'
#' @seealso [cast_fit()], [cast_predict()]
#'
#' @export
cast_evaluate <- function(fit, test_data, response = "presence",
                          threshold_method = "max_tss") {
  check_suggested("pROC", "for AUC computation")

  .cast_check_response(test_data[[response]], response)
  Y_test <- test_data[[response]]
  env_vars <- fit$env_vars

  X_test_raw <- as.data.frame(test_data[, env_vars, drop = FALSE], check.names = FALSE)
  # Same guard as the fit path: reject non-numeric predictors instead of
  # silently coercing factor level codes with as.numeric().
  .cast_check_numeric_predictors(X_test_raw, arg = "test_data")
  for (col in names(X_test_raw)) {
    X_test_raw[[col]] <- as.numeric(X_test_raw[[col]])
  }
  X_test_raw <- .cast_impute(X_test_raw, fit$scaling$impute)

  results <- list()

  for (mdl_name in names(fit$models)) {
    mdl_info <- fit$models[[mdl_name]]

    preds <- tryCatch(
      predict_single_model(mdl_info, X_test_raw),
      error = function(e) rep(NA_real_, nrow(test_data))
    )

    y_train <- fit$scaling$response
    train_ref <- fit$scaling$reference
    threshold <- NULL
    if (!is.null(y_train) && !is.null(train_ref) &&
        length(y_train) == nrow(train_ref)) {
      p_train <- tryCatch(predict_single_model(mdl_info, train_ref),
                          error = function(e) rep(NA_real_, length(y_train)))
      threshold <- tryCatch(cast_threshold(p_train, y_train,
        method = threshold_method), error = function(e) NA_real_)
    } else {
      cli::cli_warn("Fit lacks stored training responses; selecting the threshold on evaluation data is optimistic.")
    }
    ev <- evaluate_model_full(preds, Y_test, threshold = threshold,
                              threshold_method = threshold_method)
    results[[mdl_name]] <- data.frame(
      model       = mdl_name,
      auc_mean    = ev["auc"],
      pr_auc_mean = ev["pr_auc"],
      tss_mean    = ev["tss"],
      sedi_mean   = ev["sedi"],
      brier_mean  = ev["brier"],
      logloss_mean = ev["logloss"],
      boyce_mean  = ev["boyce"],
      cbi_mean    = ev["cbi"],
      tss_threshold = ev["tss_threshold"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  metrics_df <- do.call(rbind, results)
  rownames(metrics_df) <- NULL

  new_cast_eval(metrics = metrics_df, cv_source = FALSE)
}


#' Predict with a Single Model (Internal)
#'
#' @param mdl_info Model info list from cast_fit.
#' @param X_raw Raw (unscaled) test data.
#' @param clamp Logical. Explicitly controls the MaxEnt engine's internal
#'   clamping of predictors to the training range. Default `FALSE`: engine
#'   clamping is off so that every engine (RF, BRT, GAM, MaxEnt) extrapolates
#'   identically outside the training envelope, and clamping is decided once,
#'   at the package level (`.cast_clamp()` on the inputs), never silently
#'   inside one engine. Default `FALSE` keeps predictions comparable across
#'   engines and keeps MESS extrapolation flags meaningful.
#' @return Numeric vector of predictions.
#' @keywords internal
#' @noRd
predict_single_model <- function(mdl_info, X_raw, clamp = FALSE) {
  if (is.null(mdl_info$model)) return(rep(NA_real_, nrow(X_raw)))

  nm <- mdl_info$name
  # Engine namespaces are Suggests: fitting loads them, but a fit object
  # reloaded in a fresh session (e.g. readRDS) predicts into bare
  # stats::predict() with no S3 method registered. Loading the namespace
  # here registers the methods; without the package there is nothing to
  # predict with, so NA preserves the per-model robustness contract.
  engine_pkg <- switch(nm, rf = "ranger", maxent = "maxnet", brt = "gbm",
                       gam = "mgcv", NULL)
  if (!is.null(engine_pkg) && !requireNamespace(engine_pkg, quietly = TRUE)) {
    return(rep(NA_real_, nrow(X_raw)))
  }
  if (nm == "rf") {
    return(stats::predict(mdl_info$model, data = X_raw)$predictions[, "1"])
  } else if (nm == "maxent") {
    # clamp is passed explicitly: maxnet::predict.maxnet defaults to
    # clamp = TRUE, which would silently freeze MaxEnt to the training
    # range while RF/BRT/GAM extrapolate - an engine inconsistency the
    # package never asked for.
    return(as.numeric(stats::predict(mdl_info$model, X_raw, type = "logistic",
                                     clamp = isTRUE(clamp))))
  } else if (nm == "brt") {
    bt <- mdl_info$best_trees %||% 500L
    return(stats::predict(mdl_info$model, X_raw, n.trees = bt, type = "response"))
  } else if (nm == "gam") {
    return(as.numeric(stats::predict(mdl_info$model, newdata = X_raw,
                                     type = "response")))
  }
  rep(NA_real_, nrow(X_raw))
}
