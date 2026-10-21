#' Evaluate Fitted Models
#'
#' Computes discrimination, threshold-dependent, calibration, and suitability
#' metrics for fitted models on test data. Threshold-dependent metrics (TSS,
#' SEDI, kappa) are evaluated at thresholds supplied by the caller -- by
#' convention the pooled out-of-fold thresholds from [cast_cv()], the same
#' operating point used by [cast_ensemble()] and downstream products. Without
#' supplied thresholds the function falls back to selecting them on this
#' fit's in-sample training predictions, which is optimistic for flexible
#' learners; the fallback is recorded in the result.
#'
#' @param fit A [cast_fit] object.
#' @param test_data A `data.frame` with `presence` column and the same
#'   predictor variables used in fitting.
#' @param response Character. Response column name. Default `"presence"`.
#' @param threshold_method Binary threshold rule when `thresholds` is not
#'   supplied (or lacks an entry for a model). See [cast_threshold()].
#'   Default `"max_tss"`.
#' @param thresholds Optional named numeric vector of per-model binary
#'   thresholds, applied to the held-out predictions. Pass
#'   `cv$thresholds` from [cast_cv()]: those are selected on pooled
#'   out-of-fold predictions, so the held-out labels never influence the
#'   operating point. `NA`/missing entries per model fall back to the
#'   in-sample selection for that model.
#'
#' @return A `cast_eval` object (S3 class) whose `threshold_source` records
#'   whether thresholds were caller-supplied (`"out_of_fold_cv"`) or selected
#'   on in-sample training predictions (`"in_sample_training"`).
#'
#' @details
#' - **AUC**: Area Under the Receiver Operating Characteristic curve,
#'   measuring discrimination ability. Computed via [pROC::roc()] with
#'   `direction = "<"` fixed, so a worse-than-random model correctly reports
#'   AUC < 0.5 instead of being mirrored by `direction = "auto"`.
#' - **TSS**: evaluated at the supplied (out-of-fold) threshold, never tuned
#'   on the held-out labels.
#' - **PR-AUC/SEDI**: complementary rare-event discrimination metrics; PR-AUC's
#'   baseline depends on the sampled prevalence.
#' - **Boyce**: moving-window background-based estimate; `cbi_mean` is retained
#'   as a legacy fixed-bin diagnostic.
#'
#' @seealso [cast_fit()], [cast_predict()], [cast_cv()]
#'
#' @export
cast_evaluate <- function(fit, test_data, response = "presence",
                          threshold_method = "max_tss", thresholds = NULL) {
  check_suggested("pROC", "for AUC computation")

  # `NULL[["presence"]]` is NULL, `sum(is.na(NULL))` is 0 and
  # `all(logical(0) %in% c(0, 1))` is TRUE, so the response check below
  # silently PASSED for a NULL test set and the function returned all-NA
  # metrics as if the evaluation had merely been inconclusive.
  if (!is.data.frame(test_data)) {
    cli::cli_abort(
      "{.arg test_data} must be a data.frame, not {.obj_type_friendly {test_data}}."
    )
  }
  if (!response %in% names(test_data)) {
    cli::cli_abort("{.arg test_data} is missing the response column {.val {response}}.")
  }
  if (!is.null(thresholds) &&
      (!is.numeric(thresholds) || is.null(names(thresholds)))) {
    cli::cli_abort("{.arg thresholds} must be NULL or a named numeric vector.")
  }
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
  fallback_used <- FALSE
  used_thresholds <- stats::setNames(rep(NA_real_, length(fit$models)),
                                     names(fit$models))

  for (mdl_name in names(fit$models)) {
    mdl_info <- fit$models[[mdl_name]]

    preds <- tryCatch(
      predict_single_model(mdl_info, X_test_raw),
      error = function(e) rep(NA_real_, nrow(test_data))
    )

    threshold <- if (!is.null(thresholds)) {
      unname(thresholds[[mdl_name]]) %||% NA_real_
    } else {
      NA_real_
    }
    if (!is.finite(threshold)) {
      # Fallback (or missing entry): select on this fit's training
      # predictions. Optimistic for flexible learners; recorded via
      # threshold_source when it is the only source available.
      fallback_used <- TRUE
      y_train <- fit$scaling$response
      train_ref <- fit$scaling$reference
      if (!is.null(y_train) && !is.null(train_ref) &&
          length(y_train) == nrow(train_ref)) {
        p_train <- tryCatch(predict_single_model(mdl_info, train_ref),
                            error = function(e) rep(NA_real_, length(y_train)))
        threshold <- tryCatch(cast_threshold(p_train, y_train,
          method = threshold_method), error = function(e) NA_real_)
      } else {
        cli::cli_warn("Fit lacks stored training responses; selecting the threshold on evaluation data is optimistic.")
      }
    }
    used_thresholds[[mdl_name]] <- if (is.finite(threshold)) threshold else NA_real_
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
      kappa_mean  = ev["kappa"],
      omission_5_mean  = ev["omission_5"],
      omission_10_mean = ev["omission_10"],
      mpa_mean    = ev["mpa"],
      tss_threshold = ev["tss_threshold"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  metrics_df <- do.call(rbind, results)
  rownames(metrics_df) <- NULL

  threshold_source <- if (is.null(thresholds)) {
    "in_sample_training"
  } else if (fallback_used) {
    "mixed"
  } else {
    "out_of_fold_cv"
  }
  new_cast_eval(metrics = metrics_df, cv_source = FALSE,
                thresholds = used_thresholds,
                threshold_source = threshold_source)
}


#' Predict with a Single Model (Internal)
#'
#' @param mdl_info Model info list from cast_fit.
#' @param X_raw Raw (unscaled) test data.
#' @param clamp Logical. Explicitly controls the MaxEnt engine's internal
#'   clamping of predictors to the training range. Default `FALSE`: engine
#'   clamping is off so that every engine (RF, GAM, MaxEnt) extrapolates
#'   identically outside the training envelope, and clamping is decided once,
#'   at the package level (`.cast_clamp()` on the inputs), never silently
#'   inside one engine. Default `FALSE` keeps predictions comparable across
#'   engines and keeps MESS extrapolation flags meaningful.
#' @return Numeric vector of predictions.
#' @keywords internal
#' @noRd
predict_single_model <- function(mdl_info, X_raw, clamp = FALSE) {
  offset <- mdl_info$offset %||% 0
  if (identical(mdl_info$type, "external")) {
    p <- mdl_info$predict_function(mdl_info$model, X_raw)
    if (!is.numeric(p) || !is.null(dim(p)) || length(p) != nrow(X_raw) ||
        any(is.infinite(p)) || any(p[is.finite(p)] < 0 | p[is.finite(p)] > 1)) {
      cli::cli_abort("External model {.val {mdl_info$name}} must return one suitability in [0, 1] per row.")
    }
    return(.cast_apply_offset(as.numeric(p), offset))
  }
  if (is.null(mdl_info$model)) return(rep(NA_real_, nrow(X_raw)))

  nm <- mdl_info$name
  # Engine namespaces are Suggests: fitting loads them, but a fit object
  # reloaded in a fresh session (e.g. readRDS) predicts into bare
  # stats::predict() with no S3 method registered. Loading the namespace
  # here registers the methods; without the package there is nothing to
  # predict with, so NA preserves the per-model robustness contract.
  engine_pkg <- .cast_engine_package(nm)
  if (!is.null(engine_pkg) && !requireNamespace(engine_pkg, quietly = TRUE)) {
    return(rep(NA_real_, nrow(X_raw)))
  }
  p <- if (nm == "rf") {
    # predict.ranger() draws a random seed from R's RNG when `seed` is not
    # supplied, which would advance the caller's .Random.seed on every
    # prediction call. Pin it: predictions stay deterministic and the
    # caller's RNG state untouched.
    stats::predict(mdl_info$model, data = X_raw, seed = 42L)$predictions[, "1"]
  } else if (nm == "maxent") {
    # clamp is passed explicitly: maxnet::predict.maxnet defaults to
    # clamp = TRUE, which would silently freeze MaxEnt to the training
    # range while RF/GAM extrapolate - an engine inconsistency the
    # package never asked for.
    as.numeric(stats::predict(mdl_info$model, X_raw, type = "logistic",
                              clamp = isTRUE(clamp)))
  } else if (nm == "gam") {
    as.numeric(stats::predict(mdl_info$model, newdata = X_raw,
                              type = "response"))
  } else {
    rep(NA_real_, nrow(X_raw))
  }
  .cast_apply_offset(p, offset)
}
