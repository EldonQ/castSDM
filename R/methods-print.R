# Print Methods -----------------------------------------------------------

#' @export
print.cast_select <- function(x, ...) {
  n_selected <- length(x$selected)
  cli::cli_h1("castSDM Variable Selection")
  cli::cli_ul(c(
    "Method: {x$method %||% 'unknown'}",
    "Selected variables: {n_selected}"
  ))
  d <- x$diagnostics
  if (identical(d$algorithm, "elastic_net")) {
    cli::cli_text("Spatial elastic net; {d$status}. No prescribed predictor count.")
  } else if (!is.null(d$engine)) {
    cli::cli_text("{d$engine}")
  }
  cli::cli_text("Variables: {.val {x$selected}}")
  invisible(x)
}

#' @export
print.cast_importance <- function(x, ...) {
  eff <- x$effects
  cli::cli_h1("castSDM Saved Selection Description")
  cli::cli_text("Measure: {x$metric %||% 'selection status'}; retained: {sum(eff$selected)} / {nrow(eff)}")
  print(utils::head(eff, 10L), row.names = FALSE)
  cli::cli_text("Descriptive fitted-screen signal is not causal importance.")
  invisible(x)
}

#' @export
print.cast_fit <- function(x, ...) {
  model_names <- names(x$models)
  cli::cli_h1("castSDM Model Fit")
  cli::cli_ul(c(
    "Models: {.val {model_names}}",
    "Variables: {length(x$cast_vars)}"
  ))
  invisible(x)
}

#' @export
print.cast_eval <- function(x, ...) {
  src <- if (isTRUE(x$cv_source)) "Spatial CV" else "Hold-out test set"
  cli::cli_h1("castSDM Model Evaluation ({src})")
  metrics <- x$metrics
  scores <- x$predictions
  metrics$average_precision <- vapply(metrics$model, function(model) {
    column <- paste0("HSS_", model)
    if (!is.data.frame(scores) || !all(c("obs", column) %in% names(scores))) return(NA_real_)
    compute_average_precision(scores$obs, scores[[column]])
  }, numeric(1))
  fields <- intersect(c("model", "auc_mean", "average_precision", "tss_mean",
    "brier_mean", "logloss_mean", "omission_at_threshold_mean"), names(metrics))
  print(metrics[, fields, drop = FALSE], row.names = FALSE)
  cli::cli_text("AP is non-interpolated average precision; its constant-score baseline is the sampled positive fraction.")
  cli::cli_text("Trapezoidal PR-AUC and legacy top-fraction omissions remain in x$metrics; they are different diagnostics.")
  invisible(x)
}

#' @export
print.cast_cv <- function(x, ...) {
  tm <- x$threshold_method %||% "max TSS"
  cli::cli_h1("castSDM Spatial Cross-Validation")
  cli::cli_ul(c(
    "Folds (k): {x$k}",
    "Block method: {x$block_method}",
    "Models: {.val {x$metrics$model}}"
  ))
  cli::cli_h2("Aggregated metrics (mean +/- SD)")
  m <- x$metrics
  for (i in seq_len(nrow(m))) {
    cli::cli_li(paste0(
      m$model[i],
      " | AUC=", round(m$auc_mean[i], 3), " (", round(m$auc_sd[i], 3), ")",
      " TSS=", round(m$tss_mean[i], 3),
      " CBI=", round(m$cbi_mean[i], 3)
    ))
  }
  cli::cli_h2("Optimal thresholds ({tm})")
  for (nm in names(x$thresholds)) {
    cli::cli_li("{nm}: {round(x$thresholds[nm], 3)}")
  }
  invisible(x)
}

#' @export
print.cast_predict <- function(x, ...) {
  n_sites <- nrow(x$predictions)
  cli::cli_h1("castSDM Spatial Predictions")
  cli::cli_ul(c("Sites: {n_sites}", "Models: {.val {x$models}}"))
  invisible(x)
}

#' @export
print.cast_ensemble <- function(x, ...) {
  cli::cli_h1("castSDM Ensemble Prediction")
  cli::cli_ul(c(
    "Method: {x$method}",
    "Threshold: {round(x$threshold, 3)}",
    "Sites: {nrow(x$predictions)}"
  ))
  if (!is.null(x$weights) && length(x$weights) > 0) {
    cli::cli_text("Weights: {paste(names(x$weights), round(x$weights, 3), sep='=', collapse=', ')}")
  }
  if (isTRUE(x$weights_fallback)) {
    cli::cli_warn("Weight fallback used: {x$fallback_reason}")
  }
  if (identical(x$method, "committee")) {
    cli::cli_text("Ensemble surface is the fraction of models voting suitable; threshold = 0.5.")
  }
  invisible(x)
}

#' @export
print.cast_ensemble_grouped <- function(x, ...) {
  cli::cli_h1("castSDM Grouped Ensemble")
  cli::cli_ul(c(
    "Groups: {length(x$groups)} ({.val {names(x$groups)}})",
    "Sites: {nrow(x$predictions)}"
  ))
  for (nm in names(x$ensembles)) {
    ens <- x$ensembles[[nm]]
    w <- ens$weights[ens$weights > 0]
    cli::cli_li(
      "{.field {nm}}: method {.val {ens$method}}, threshold {round(ens$threshold, 3)}, weights {paste(names(w), round(w, 3), sep='=', collapse=', ')}"
    )
  }
  invisible(x)
}

#' @export
print.cast_ensemble_importance <- function(x, ...) {
  cli::cli_h1("castSDM Ensemble Permutation Importance")
  cli::cli_ul(c(
    "Permutations: {x$n_perm} per predictor",
    "Correlation: {x$metric} (importance = 1 - cor)"
  ))
  print(x$importance, row.names = FALSE)
  invisible(x)
}

#' @export
print.cast_project <- function(x, ...) {
  n_scenarios <- length(x$future)
  cli::cli_h1("castSDM Future Projection")
  cli::cli_ul(c(
    "Scenarios: {n_scenarios} ({.val {names(x$future)}})",
    "Current range cells: {sum(x$current$predictions$binary_ensemble == 1, na.rm = TRUE)}"
  ))
  if (!is.null(x$stats) && nrow(x$stats) > 0) {
    cli::cli_h2("Range change summary")
    for (i in seq_len(nrow(x$stats))) {
      s <- x$stats[i, ]
      cli::cli_li("{s$scenario}: gain={s$n_gain} loss={s$n_loss} stable={s$n_stable_present} shift={round(s$centroid_shift_km, 1)}km")
    }
  }
  invisible(x)
}

#' @export
print.cast_result <- function(x, ...) {
  cli::cli_h1("castSDM Pipeline Result")
  cli::cli_ul(c(
    "Selected: {length(x$screen$selected)} variables",
    "Models: {.val {names(x$fit$models)}}",
    "Full-data refit: {if (!is.null(x$fit_full)) 'Yes' else 'No'}",
    "Predictions: {if (!is.null(x$predict)) 'Yes' else 'No'}",
    "Ensemble: {if (!is.null(x$ensemble)) 'Yes' else 'No'}"
  ))
  invisible(x)
}
