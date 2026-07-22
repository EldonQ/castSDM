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
    if (identical(x$method, "dml")) {
      cli::cli_text(
        "{d$engine} | FDR = {d$fdr_method} at alpha = {d$alpha} | {d$n_folds}-fold cross-fitting"
      )
    } else {
      cli::cli_text("{d$engine}")
    }
  }
  cli::cli_text("Variables: {.val {x$selected}}")
  invisible(x)
}

#' @export
print.cast_effect <- function(x, ...) {
  eff <- x$effects
  n_sig <- sum(eff$selected, na.rm = TRUE)
  is_cpi <- identical(x$diagnostics$measure, "cpi")
  if (is_cpi) {
    cli::cli_h1("castSDM Conditional Predictive Impact (CPI)")
    cli::cli_ul(c(
      "Conditional importance (log-loss knockoff), {round(100 * x$conf_level)}% CI",
      "Significant (FDR < {x$alpha}): {n_sig} / {nrow(eff)}"
    ))
  } else {
    cli::cli_h1("castSDM Causal Effects (DML)")
    cli::cli_ul(c(
      "Partial-linear effect per +1 SD, {round(100 * x$conf_level)}% CI",
      "Significant (FDR < {x$alpha}): {n_sig} / {nrow(eff)}"
    ))
  }
  show <- utils::head(eff, 10L)
  disp <- data.frame(
    variable = show$variable,
    estimate = round(show$estimate, 4),
    ci = sprintf("[%.3f, %.3f]", show$conf_low, show$conf_high),
    p_adjusted = signif(show$p_adjusted, 3),
    sig = ifelse(show$selected, "*", ""),
    stringsAsFactors = FALSE
  )
  if (is_cpi) names(disp)[names(disp) == "estimate"] <- "cpi"
  print(disp, row.names = FALSE)
  invisible(x)
}

#' @export
print.cast_counterfactual <- function(x, ...) {
  s <- x$summary
  cli::cli_h1("castSDM Counterfactual What-If")
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
    "Predictions: {if (!is.null(x$predict)) 'Yes' else 'No'}",
    "Ensemble: {if (!is.null(x$ensemble)) 'Yes' else 'No'}"
  ))
  invisible(x)
}

#' @export
print.cast_batch <- function(x, ...) {
  cli::cli_h1("castSDM Batch Results")
  cli::cli_ul(c(
    "Species: {.val {x$species}}",
    "Models: {.val {x$models}}"
  ))
  if (!is.null(x$output_dir)) {
    cli::cli_text("  Output: {x$output_dir}")
  }
  if (!is.null(x$species_metrics) && nrow(x$species_metrics) > 0) {
    cli::cli_h2("Summary metrics")
    print(x$species_metrics[, intersect(
      c("species", "model", "auc_mean", "tss_mean", "cbi_mean"),
      names(x$species_metrics)
    ), drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}
