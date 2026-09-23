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

# One predictor frame, all shifts, one prediction batch. Keeping the shifts in
# a single stack means one predict() pass per engine rather than one per shift.
# Results are read back out as an n-by-length(steps) matrix.
.cast_shift_stack <- function(X, driver, steps) {
  n <- nrow(X)
  stacked <- X[rep(seq_len(n), times = length(steps)), , drop = FALSE]
  stacked[[driver]] <- stacked[[driver]] +
    rep(as.numeric(steps), each = n)
  list(stacked = stacked, n = n, n_steps = length(steps))
}

.cast_shift_effects <- function(fit, X, driver, steps) {
  base <- .cast_predict_matrix(fit, X, names(fit$models))
  st <- .cast_shift_stack(X, driver, steps)
  cf <- .cast_predict_matrix(fit, st$stacked, names(fit$models))
  paired_base <- base[rep(seq_len(st$n), st$n_steps), , drop = FALSE]
  delta <- cf - paired_base
  delta[!is.finite(cf) | !is.finite(paired_base)] <- NA_real_
  delta <- rowMeans(delta, na.rm = TRUE)
  delta[!is.finite(delta)] <- NA_real_
  matrix(delta, nrow = st$n, ncol = st$n_steps)
}

# A marginal-quantile box detects extrapolation, not holes in joint support.
.cast_support_bounds <- function(ref, probs) {
  if (!is.numeric(probs) || length(probs) != 2L || any(!is.finite(probs)) ||
      probs[1] < 0 || probs[2] > 1 || probs[1] >= probs[2]) {
    cli::cli_abort("{.arg support_probs} must be two increasing probabilities in [0, 1].")
  }
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

.cast_driver_sd <- function(fit, driver) {
  sds <- fit$scaling$sds
  v <- if (!is.null(sds)) sds[[driver]] else NULL
  if (is.null(v) || !is.finite(v) || v <= 0) {
    ref <- fit$scaling$reference
    v <- if (!is.null(ref)) stats::sd(ref[[driver]], na.rm = TRUE) else NA_real_
  }
  if (!is.finite(v) || v <= 0) {
    cli::cli_abort(c("No usable scale for {.val {driver}}.",
                     "i" = "The stored training SD is missing or zero."))
  }
  v
}

#' Dose-Response of Predicted Probability to a Predictor Shift
#'
#' Sweeps the size of an additive intervention on one predictor, holding every
#' other predictor at its observed value, and reports the mean change in
#' predicted probability at each shift size. This is the response curve a
#' single summary number cannot show: a bounded or saturating response is
#' invisible once the curve is collapsed to one signed value.
#'
#' @section What the numbers mean:
#' `mean_delta` is the pointwise g-computation contrast
#' \eqn{E[Y \mid A = a + \delta, C = c] - E[Y \mid A = a, C = c]} averaged over
#' the rows supplied in `newdata`; `mean_abs_delta` is its absolute value,
#' which does not cancel when a driver raises suitability in some places and
#' lowers it in others. Both are contrasts of the \emph{fitted} model on the
#' scale it predicts. On presence/background data that scale is relative
#' suitability, not occurrence probability.
#'
#' `support` is the fraction of evaluated rows whose baseline and shifted
#' vectors both lie inside every predictor's training quantile bounds. This
#' box can contain unobserved combinations, so high coverage does not establish
#' conditional positivity. `range_supported` marks coverage of at least 0.5:
#' a descriptive range screen, not causal estimability. Contrasts are computed
#' over all complete evaluated rows, including those outside the box.
#'
#' @param fit A `cast_fit` object.
#' @param variable Single predictor to intervene on.
#' @param shift Numeric vector of shift sizes, in training standard deviations
#'   unless `shift_type = "raw"`. A zero shift is not an intervention and is
#'   dropped from the scan. Default scans -3 to 3 SD.
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param newdata Optional data frame defining the reference population.
#'   Defaults to training predictors stored in `fit`. Incomplete or non-finite
#'   rows are excluded from both effects and coverage.
#' @param support_probs Length-2 numeric. Increasing training quantiles in
#'   `[0, 1]` defining a box across all predictors. Default `c(0.01, 0.99)`.
#' @param max_rows Positive integer limiting rows sampled evenly from the
#'   complete reference population, or `Inf` for all rows. Default 5000.
#'
#' @return A `cast_dose_response` object with `curve` (data frame: `shift`,
#'   `shift_raw`, `mean_delta`, `mean_abs_delta`, `support`, `range_supported`),
#'   `shift_type`, `variable`, `unit`, `models`, and `assumptions`.
#' @seealso [cast_effect_support()], [cast_effect_table()]
#' @export
cast_dose_response <- function(fit, variable, shift = seq(-3, 3, by = 0.25),
                               shift_type = c("sd", "raw"), newdata = NULL,
                               support_probs = c(0.01, 0.99),
                               max_rows = 5000L) {
  .cast_check_fit(fit)
  shift_type <- match.arg(shift_type)
  env_vars <- fit$env_vars
  if (length(variable) != 1L || !variable %in% env_vars) {
    cli::cli_abort(c("{.arg variable} must be one fitted predictor.",
                     "i" = "Available: {.val {env_vars}}."))
  }
  shift <- as.numeric(shift)
  shift <- shift[is.finite(shift) & shift != 0]
  if (!length(shift)) {
    cli::cli_abort("{.arg shift} must supply at least one non-zero finite value.")
  }
  min_support <- 0.5
  .cast_check_shift_reachable_ref(fit, variable, shift, shift_type)

  ref <- .cast_reference(fit, env_vars)
  X <- .cast_reference_rows(fit, newdata, env_vars, max_rows)
  shift_raw <- if (identical(shift_type, "raw")) shift else
    shift * .cast_driver_sd(fit, variable)

  sup <- .cast_support_fraction(ref, variable, shift_raw,
                                probs = support_probs, newdata = X)
  eff <- .cast_shift_effects(fit, X, variable, shift_raw)

  curve <- data.frame(
    shift = shift,
    shift_raw = shift_raw,
    mean_delta = colMeans(eff, na.rm = TRUE),
    mean_abs_delta = colMeans(abs(eff), na.rm = TRUE),
    support = as.numeric(sup),
    stringsAsFactors = FALSE)
  curve$range_supported <- is.finite(curve$support) & curve$support >= min_support
  attr(curve, "shift_type") <- shift_type
  attr(curve, "support_probs") <- support_probs
  attr(curve, "min_support") <- min_support

  new_cast_dose_response(
    curve = curve, variable = variable, shift_type = shift_type,
    unit = if (identical(shift_type, "raw")) "raw predictor units" else
      "training SD",
    models = names(fit$models),
    assumptions = .cast_effect_assumptions())
}


#' Support Fraction of an Interventional Shift
#'
#' Reports, per driver and shift, the fraction of complete evaluated rows whose
#' baseline and shifted vectors both lie inside every predictor's training
#' quantile bounds. Low coverage flags range extrapolation. Unobserved predictor
#' combinations inside this box remain undetected, so even full coverage does
#' not establish conditional positivity or causal identification.
#'
#' @inheritParams cast_dose_response
#' @param variables Predictors to assess. Default: every fitted predictor.
#' @param shift Numeric vector of shift sizes, in training standard deviations
#'   unless `shift_type = "raw"`.
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param max_rows Positive integer limiting complete rows sampled evenly, or
#'   `Inf` for all rows. Default 20000.
#'
#' @return A `cast_support` object with `support`: a data frame of
#'   `driver`, `shift`, `shift_raw`, and `support`.
#' @seealso [cast_effect_map()], [cast_dose_response()]
#' @export
cast_effect_support <- function(fit, variables = NULL,
                                shift = c(1, 2, 3), shift_type = c("sd", "raw"),
                                newdata = NULL, support_probs = c(0.01, 0.99),
                                max_rows = 20000L) {
  .cast_check_fit(fit)
  shift_type <- match.arg(shift_type)
  env_vars <- fit$env_vars
  variables <- variables %||% env_vars
  if (!length(variables) || any(!variables %in% env_vars)) {
    cli::cli_abort("Unknown {.arg variables}: {.val {setdiff(variables, env_vars)}}.")
  }
  shift <- as.numeric(shift)
  if (!length(shift) || any(!is.finite(shift)) || any(shift == 0)) {
    cli::cli_abort("{.arg shift} must be non-zero finite numbers.")
  }
  ref <- .cast_reference(fit, env_vars)
  X <- .cast_reference_rows(fit, newdata, env_vars, max_rows)

  pieces <- lapply(variables, function(v) {
    sr <- if (identical(shift_type, "raw")) shift else shift * .cast_driver_sd(fit, v)
    frac <- .cast_support_fraction(ref, v, sr, probs = support_probs, newdata = X)
    data.frame(driver = v, shift = shift, shift_raw = sr,
               support = as.numeric(frac),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, pieces)
  rownames(tab) <- NULL
  new_cast_support(support = tab, variables = variables, shift = shift,
                   shift_type = shift_type, support_probs = support_probs,
                   models = names(fit$models),
                   assumptions = .cast_effect_assumptions())
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

.cast_reference_rows <- function(fit, newdata, env_vars, max_rows) {
  if (!is.numeric(max_rows) || length(max_rows) != 1L || is.na(max_rows) ||
      max_rows < 1 || (is.finite(max_rows) && max_rows != floor(max_rows))) {
    cli::cli_abort("{.arg max_rows} must be a positive integer or Inf.")
  }
  X <- if (is.null(newdata)) {
    .cast_reference(fit, env_vars)
  } else {
    miss <- setdiff(env_vars, names(newdata))
    if (length(miss)) {
      cli::cli_abort("{.arg newdata} lacks fitted predictor{?s}: {.val {miss}}.")
    }
    as.data.frame(newdata[, env_vars, drop = FALSE])
  }
  .cast_check_numeric_predictors(X, arg = "newdata")
  ok <- rowSums(!is.finite(as.matrix(X))) == 0L
  X <- X[ok, , drop = FALSE]
  if (!nrow(X)) cli::cli_abort("No complete finite predictor rows to evaluate.")
  if (is.finite(max_rows) && nrow(X) > max_rows) {
    X <- X[unique(round(seq(1, nrow(X), length.out = max_rows))), , drop = FALSE]
  }
  rownames(X) <- NULL
  X
}

.cast_check_shift_reachable_ref <- function(fit, driver, shifts, shift_type) {
  ref <- .cast_reference(fit, fit$env_vars)
  sr <- if (identical(shift_type, "raw")) shifts else
    shifts * .cast_driver_sd(fit, driver)
  .cast_check_shift_reachable(ref, driver, sr)
}
