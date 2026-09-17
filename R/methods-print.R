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
  if (!is.null(d$engine)) {
    cli::cli_text("{d$engine}")
  }
  cli::cli_text("Variables: {.val {x$selected}}")
  invisible(x)
}

#' @export
print.cast_importance <- function(x, ...) {
  eff <- x$effects
  n_sig <- sum(eff$selected, na.rm = TRUE)
  cli::cli_h1("castSDM Predictor Attribution")
  bullets <- c(
    "Interventional effect (shift 1 SD, other predictors fixed)",
    "Calibrated against a permuted-response null",
    "Above the null (p < {x$alpha}): {n_sig} / {nrow(eff)}"
  )
  ag <- x$diagnostics$importance_agreement
  if (!is.null(ag) && is.finite(ag)) {
    bullets <- c(bullets,
                 "Spearman agreement with permutation importance: {round(ag, 3)}")
  }
  cli::cli_ul(bullets)
  show <- utils::head(eff, 10L)
  disp <- data.frame(
    variable = show$variable,
    effect = signif(show$interventional_effect, 4),
    null = signif(show$null_threshold, 4),
    perm_imp = signif(show$perm_importance, 4),
    p_value = signif(show$p_value, 3),
    sig = ifelse(show$selected, "*", ""),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  cli::cli_text("Each predictor is compared with its own permuted-response null ({.field null}).")
  invisible(x)
}

#' @export
print.cast_sensitivity <- function(x, ...) {
  s <- x$summary
  cli::cli_h1("castSDM Sensitivity What-If")
  cli::cli_ul(c(
    "Intervention: {x$variable} + {x$shift} ({x$shift_type})",
    "Models averaged: {.val {x$models}}",
    "Cells with suitability gain: {round(100 * s$frac_positive, 1)}%"
  ))
  cli::cli_text(
    "Delta HSS: mean = {round(s$mean_delta, 4)}, range = [{round(s$max_loss, 3)}, {round(s$max_gain, 3)}]"
  )
  invisible(x)
}

#' @export
print.cast_dose_response <- function(x, ...) {
  cli::cli_h1("castSDM Dose-Response")
  cli::cli_ul(c(
    "Intervention: {x$variable} shift in {x$unit}",
    "Models averaged: {.val {x$models}}"
  ))
  keep <- x$curve[x$curve$estimable, , drop = FALSE]
  if (nrow(keep)) {
    imax <- which.max(keep$mean_abs_delta)
    cli::cli_text(
      "Largest mean |delta| = {round(keep$mean_abs_delta[imax], 4)} at shift {round(keep$shift[imax], 2)} ({x$unit}); support {round(keep$support[imax], 3)}"
    )
  } else {
    cli::cli_text("No shift in the requested range had adequate support.")
  }
  cli::cli_text("{round(100 * mean(x$curve$estimable), 0)}% of the shifted range is on support (>= {attr(x$curve, 'min_support')}).")
  invisible(x)
}


#' @export
print.cast_support <- function(x, ...) {
  cli::cli_h1("castSDM Shift Support (positivity)")
  cli::cli_text("Fraction of observed predictor vectors still inside the training support.")
  tab <- x$support
  tab$support <- round(tab$support, 3)
  tab$shift_raw <- signif(tab$shift_raw, 4)
  print(tab, row.names = FALSE)
  cli::cli_text("Low support means the effect is answered by extrapolation; report it or reduce the shift.")
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
  print(x$metrics)
  invisible(x)
}

#' @export
print.cast_cv <- function(x, ...) {
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
  cli::cli_h2("Optimal thresholds (max TSS)")
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
