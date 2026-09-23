# Quantile-box coverage diagnoses range extrapolation, not joint positivity.

.cast_check_fit <- function(fit) {
  if (!inherits(fit, "cast_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls cast_fit} object.")
  }
  if (!length(fit$models)) {
    cli::cli_abort("{.arg fit} carries no fitted models.")
  }
  invisible(fit)
}

# A marginal-quantile box detects extrapolation, not holes in joint support.
.cast_check_support_probs <- function(probs) {
  if (!is.numeric(probs) || length(probs) != 2L || any(!is.finite(probs)) ||
      probs[1] < 0 || probs[2] > 1 || probs[1] >= probs[2]) {
    cli::cli_abort("{.arg support_probs} must be two increasing probabilities in [0, 1].")
  }
  invisible(probs)
}

.cast_support_bounds <- function(ref, probs) {
  .cast_check_support_probs(probs)
  lapply(ref, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(NA_real_, NA_real_))
    unname(stats::quantile(v, probs = probs, names = FALSE))
  })
}

.cast_support_rows <- function(X, driver, steps, bounds) {
  baseline <- rep(TRUE, nrow(X))
  for (v in names(bounds)) {
    b <- bounds[[v]]
    if (any(!is.finite(b))) return(matrix(NA, nrow(X), length(steps)))
    baseline <- baseline & is.finite(X[[v]]) & X[[v]] >= b[1] & X[[v]] <= b[2]
  }
  b <- bounds[[driver]]
  out <- vapply(steps, function(s) {
    z <- X[[driver]] + s
    baseline & is.finite(z) & z >= b[1] & z <= b[2]
  }, logical(nrow(X)))
  matrix(out, nrow = nrow(X), ncol = length(steps))
}

.cast_support_fraction <- function(ref, driver, steps,
                                   probs = c(0.01, 0.99), newdata = ref) {
  colMeans(.cast_support_rows(newdata, driver, steps,
                              .cast_support_bounds(ref, probs)))
}

# Refuse a shift that leaves every observed value of the driver outside the
# training range: the contrast would be pure extrapolation everywhere.
.cast_check_shift_reachable <- function(ref, driver, steps) {
  v <- ref[[driver]]
  v <- v[is.finite(v)]
  if (!length(v)) return(invisible(TRUE))
  lim <- range(v)
  frac <- vapply(steps, function(s)
    mean(v + s >= lim[1] & v + s <= lim[2]), numeric(1))
  if (all(frac <= 0)) {
    cli::cli_abort(c(
      "Every {.val {driver}} value shifted by {.val {steps}} leaves the training range.",
      "i" = "These contrasts would require complete range extrapolation; reduce the shifts."))
  }
  invisible(TRUE)
}

.cast_reference <- function(fit, env_vars) {
  ref <- fit$scaling$reference
  if (is.null(ref)) {
    cli::cli_abort(c(
      "{.arg fit} carries no training predictor reference.",
      "i" = "Refit with {.fn cast_fit}; support cannot be judged without it."))
  }
  ref <- as.data.frame(ref)
  miss <- setdiff(env_vars, names(ref))
  if (length(miss)) {
    cli::cli_abort("{.arg fit} reference lacks predictor{?s}: {.val {miss}}.")
  }
  ref[, env_vars, drop = FALSE]
}

#' Dose-Response (removed in 0.12.0)
#'
#' Removed: use [cast_effect_table()] / [cast_effect_map()] with a single
#' raw-unit shift and hard range-masking.
#' @param ... Ignored.
#' @return Never returns; always aborts.
#' @export
cast_dose_response <- function(...) {
  cli::cli_abort("`cast_dose_response()` was removed in 0.12.0; use `cast_effect_table()` / `cast_effect_map()`.")
}


#' Support Fraction (removed in 0.12.0)
#'
#' Removed: range-masking is now hard-wired into [cast_effect_table()] /
#' [cast_effect_map()] (`support < 0.5` masked).
#' @param ... Ignored.
#' @return Never returns; always aborts.
#' @export
cast_effect_support <- function(...) {
  cli::cli_abort("`cast_effect_support()` was removed in 0.12.0; masking now lives in `cast_effect_table()` / `cast_effect_map()`.")
}

# ---- internals ------------------------------------------------------------

.cast_effect_assumptions <- function() paste(
  "Model-based interventional contrast (g-computation / standardization).",
  "Causal interpretation requires consistency, a justified adjustment set,",
  "no uncontrolled confounding, joint positivity, and adequate response and",
  "observation models. Quantile-box coverage only diagnoses range",
  "extrapolation; it does not establish joint positivity or identification.",
  "On presence/background data the scale is relative suitability, not",
  "occurrence probability.")

# ---- end of effect-grid helpers --------------------------------------------
