#' Evaluate Fitted Models
#'
#' Computes AUC (Area Under ROC Curve), TSS (True Skill Statistic), and
#' CBI (Continuous Boyce Index) for fitted models on test data.
#'
#' @param fit A [cast_fit] object.
#' @param test_data A `data.frame` with `presence` column and the same
#'   predictor variables used in fitting.
#' @param response Character. Response column name. Default `"presence"`.
#'
#' @return A `cast_eval` object (S3 class).
#'
#' @details
#' - **AUC**: Area Under the Receiver Operating Characteristic curve,
#'   measuring discrimination ability. Computed via [pROC::roc()].
#' - **TSS**: True Skill Statistic = Sensitivity + Specificity - 1.
#' - **CBI**: Continuous Boyce Index, measuring predicted–expected ratio.
#'
#' @seealso [cast_fit()], [cast_predict()]
#'
#' @export
cast_evaluate <- function(fit, test_data, response = "presence") {
  check_suggested("pROC", "for AUC computation")

  Y_test <- test_data[[response]]
  env_vars <- fit$env_vars

  X_test_raw <- as.data.frame(test_data[, env_vars, drop = FALSE])
  for (col in names(X_test_raw)) {
    X_test_raw[[col]] <- as.numeric(X_test_raw[[col]])
  }
  X_test_raw[is.na(X_test_raw)] <- 0

  results <- list()

  for (mdl_name in names(fit$models)) {
    mdl_info <- fit$models[[mdl_name]]

    preds <- tryCatch(
      predict_single_model(mdl_info, X_test_raw),
      error = function(e) rep(NA_real_, nrow(test_data))
    )

    ev <- evaluate_model_full(preds, Y_test)
    results[[mdl_name]] <- data.frame(
      model       = mdl_name,
      auc_mean    = ev["auc"],
      tss_mean    = ev["tss"],
      cbi_mean    = ev["cbi"],
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
#' @return Numeric vector of predictions.
#' @keywords internal
#' @noRd
predict_single_model <- function(mdl_info, X_raw) {
  if (is.null(mdl_info$model)) return(rep(NA_real_, nrow(X_raw)))

  nm <- mdl_info$name
  if (nm == "rf") {
    return(stats::predict(mdl_info$model, data = X_raw)$predictions[, "1"])
  } else if (nm == "maxent") {
    return(as.numeric(stats::predict(mdl_info$model, X_raw, type = "logistic")))
  } else if (nm == "brt") {
    bt <- mdl_info$best_trees %||% 500L
    return(stats::predict(mdl_info$model, X_raw, n.trees = bt, type = "response"))
  } else if (nm == "gam") {
    return(as.numeric(stats::predict(mdl_info$model, newdata = X_raw,
                                     type = "response")))
  } else if (nm == "esm") {
    return(predict_cast_esm(mdl_info$model, X_raw))
  }
  rep(NA_real_, nrow(X_raw))
}
