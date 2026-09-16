# ==========================================================================
# Effect grids: interrogating an interventional effect over the space it is
# claimed to hold, instead of collapsing it to one number.
#
#   cast_dose_response()   - how the effect scales with the size of the shift
#   cast_effect_heatmap()  - where in (driver, effect modifier) space the
#                            effect is large, and which cells are even
#                            answerable given the observed data support
#
# The support rule is the positivity assumption of causal inference made
# operational: a shift is on support when the shifted predictor vector stays
# inside the observed data. Cells outside the support are reported, never
# silently coloured in. See Petersen et al. (2012) for the diagnostic
# prescription this follows.
# ==========================================================================

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
.cast_shift_stack <- function(fit, X, driver, steps) {
  n <- nrow(X)
  stacked <- X[rep(seq_len(n), times = length(steps)), , drop = FALSE]
  stacked[[driver]] <- stacked[[driver]] +
    rep(as.numeric(steps), each = n)
  list(stacked = stacked, n = n, n_steps = length(steps))
}

.cast_shift_effects <- function(fit, X, driver, steps) {
  base <- .cast_predict_matrix(fit, X, names(fit$models))
  base_p <- rowMeans(base, na.rm = TRUE)
  st <- .cast_shift_stack(fit, X, driver, steps)
  cf <- .cast_predict_matrix(fit, st$stacked, names(fit$models))
  cf_p <- rowMeans(cf, na.rm = TRUE)
  delta <- cf_p - base_p
  delta[!is.finite(delta)] <- NA_real_
  matrix(delta, nrow = st$n, ncol = st$n_steps)
}

# Joint-support rule. `quantile` (default) treats the training 1st-99th
# percentile box in every predictor as the supported region; that is a
# marginal-quantile (hyper-rectangle) approximation to the joint support, so
# it flags gross off-support shifts and not subtler holes inside the box.
.cast_support_bounds <- function(ref, probs) {
  lapply(ref, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(NA_real_, NA_real_))
    unname(stats::quantile(v, probs = probs, names = FALSE))
  })
}

.cast_support_fraction <- function(ref, driver, steps,
                                   probs = c(0.01, 0.99)) {
  bounds <- .cast_support_bounds(ref, probs)
  out <- lapply(names(ref), function(v) {
    b <- bounds[[v]]
    if (any(!is.finite(b)) || b[1] >= b[2]) {
      return(stats::setNames(rep(NA_real_, length(steps)), as.character(steps)))
    }
    x <- ref[[v]]
    if (identical(v, driver)) {
      vals <- vapply(steps, function(s) {
        z <- x + s
        mean(z >= b[1] & z <= b[2])
      }, numeric(1))
    } else {
      vals <- rep(mean(x >= b[1] & x <= b[2]), length(steps))
    }
    stats::setNames(vals, as.character(steps))
  })
  names(out) <- names(ref)
  do.call(rbind, out)
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
      "i" = "The effect is not estimable anywhere for this shift; reduce it."))
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
#' `support` is the fraction of rows whose shifted predictor vector remains
#' inside the observed multivariable support. A shift with low `support` is
#' answered mostly by extrapolation; report it as such or reduce it.
#'
#' @param fit A `cast_fit` object.
#' @param variable Single predictor to intervene on.
#' @param shift Numeric vector of shift sizes, in training standard deviations
#'   unless `shift_type = "raw"`. A zero shift is not an intervention and is
#'   dropped from the scan. Default scans -3 to 3 SD.
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param newdata Optional data frame of observed predictor values used as the
#'   reference population. Defaults to the training data stored in `fit`.
#' @param support_probs Length-2 numeric. Training quantiles treated as the
#'   supported range of each predictor. Default `c(0.01, 0.99)`.
#' @param max_rows Rows sampled (evenly) from `newdata`. Default 5000.
#'
#' @return A `cast_dose_response` object with `curve` (data frame: `shift`,
#'   `shift_raw`, `mean_delta`, `mean_abs_delta`, `support`, `estimable`),
#'   `shift_type`, `variable`, `unit`, `models`, and `assumptions`.
#' @seealso [cast_effect_heatmap()], [cast_effect_table()]
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

  eff <- .cast_shift_effects(fit, X, variable, shift_raw)
  sup <- .cast_support_fraction(ref, variable, shift_raw, probs = support_probs)

  curve <- data.frame(
    shift = shift,
    shift_raw = shift_raw,
    mean_delta = colMeans(eff, na.rm = TRUE),
    mean_abs_delta = colMeans(abs(eff), na.rm = TRUE),
    support = as.numeric(sup[variable, ]),
    stringsAsFactors = FALSE)
  est <- is.finite(curve$support) & curve$support >= min_support
  curve$estimable <- est
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

#' Effect Heatmap over Driver and Effect-Modifier Space
#'
#' Bins the observed data by the intervened predictor and one effect modifier,
#' and reports the mean interventional effect and the data support in every
#' bin. The result answers two questions at once: where the shift matters, and
#' where it is answerable at all.
#'
#' @section Reading the grid:
#' `effect` is the mean shift in predicted probability inside the bin.
#' `support` is the fraction of rows in the bin whose shifted predictor vector
#' stays inside the observed multivariable support. Bins flagged
#' `supported = FALSE` are answered mostly by extrapolation and should be left
#' blank in a figure rather than coloured in.
#'
#' The grid is retrospective: bins are built from observed rows, so it shows
#' how the fitted effect varies across the conditions actually sampled. It is
#' not a model-free response surface and it cannot certify effect modification
#' that the data never sampled.
#'
#' @param fit A `cast_fit` object.
#' @param variable Single predictor to intervene on.
#' @param modifier Single predictor defining the second axis.
#' @param shift Shift size, in training standard deviations unless
#'   `shift_type = "raw"`. Default `1`.
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param n_bins Number of quantile bins per axis. Default 8.
#' @param min_n Minimum rows per bin to report. Default 20.
#' @param newdata Optional observed predictor frame. Defaults to training data.
#' @param support_probs Length-2 numeric. Default `c(0.01, 0.99)`.
#' @param max_rows Rows sampled (evenly). Default 20000.
#'
#' @return A `cast_effect_heatmap` object with `grid` (data frame: `x_mid`,
#'   `y_mid`, `x_lo`, `x_hi`, `y_lo`, `y_hi`, `n`, `effect`, `abs_effect`,
#'   `support`, `supported`), plus `variable`, `modifier`, `shift`,
#'   `shift_type`, `unit`, `min_support`, and `assumptions`.
#' @seealso [cast_dose_response()], [cast_effect_map()]
#' @export
cast_effect_heatmap <- function(fit, variable, modifier, shift = 1,
                                shift_type = c("sd", "raw"), n_bins = 8L,
                                min_n = 20L, newdata = NULL,
                                support_probs = c(0.01, 0.99),
                                max_rows = 20000L) {
  .cast_check_fit(fit)
  shift_type <- match.arg(shift_type)
  env_vars <- fit$env_vars
  for (nm in c("variable", "modifier")) {
    v <- get(nm)
    if (length(v) != 1L || !v %in% env_vars) {
      cli::cli_abort(c("{.arg {nm}} must be one fitted predictor.",
                       "i" = "Available: {.val {env_vars}}."))
    }
  }
  if (identical(variable, modifier)) {
    cli::cli_abort("{.arg modifier} must differ from {.arg variable}.")
  }
  shift <- as.numeric(shift)
  if (length(shift) != 1L || !is.finite(shift) || shift == 0) {
    cli::cli_abort("{.arg shift} must be one non-zero finite number.")
  }
  n_bins <- as.integer(n_bins)
  if (is.na(n_bins) || n_bins < 2L) cli::cli_abort("{.arg n_bins} must be >= 2.")
  min_n <- as.integer(min_n)
  if (is.na(min_n) || min_n < 1L) cli::cli_abort("{.arg min_n} must be >= 1.")
  min_support <- 0.5

  shift_raw <- if (identical(shift_type, "raw")) shift else
    shift * .cast_driver_sd(fit, variable)
  .cast_check_shift_reachable_ref(fit, variable, shift_raw, "raw")

  ref <- .cast_reference(fit, env_vars)
  X <- .cast_reference_rows(fit, newdata, env_vars, max_rows)
  eff <- .cast_shift_effects(fit, X, variable, shift_raw)[, 1L]
  sup <- .cast_support_fraction(ref, variable, shift_raw,
                                probs = support_probs)[variable, 1L]

  bx <- .cast_bin_index(X[[variable]], n_bins)
  by <- .cast_bin_index(X[[modifier]], n_bins)
  key <- (bx - 1L) * n_bins + by
  n_bin <- tabulate(key, nbins = n_bins * n_bins)

  rows <- lapply(seq_len(n_bins * n_bins), function(k) {
    sel <- which(key == k)
    if (!length(sel)) return(NULL)
    xb <- bx[sel[1L]]; yb <- by[sel[1L]]
    data.frame(
      x_bin = xb, y_bin = yb,
      x_mid = stats::median(X[[variable]][sel]),
      y_mid = stats::median(X[[modifier]][sel]),
      x_lo = min(X[[variable]][sel]), x_hi = max(X[[variable]][sel]),
      y_lo = min(X[[modifier]][sel]), y_hi = max(X[[modifier]][sel]),
      n = length(sel),
      effect = mean(eff[sel], na.rm = TRUE),
      abs_effect = mean(abs(eff[sel]), na.rm = TRUE),
      support = sup,
      stringsAsFactors = FALSE)
  })
  grid <- do.call(rbind, rows)
  if (is.null(grid)) cli::cli_abort("No usable bins in {.arg newdata}.")
  grid <- grid[grid$n >= min_n, , drop = FALSE]
  if (!nrow(grid)) {
    cli::cli_abort(c("Every bin holds fewer than {.arg min_n} rows.",
                     "i" = "Lower {.arg min_n} or raise {.arg n_bins} span."))
  }
  grid$supported <- is.finite(grid$support) & grid$support >= min_support
  rownames(grid) <- NULL

  new_cast_effect_heatmap(
    grid = grid, variable = variable, modifier = modifier,
    shift = shift, shift_type = shift_type,
    unit = if (identical(shift_type, "raw")) "raw predictor units" else
      "training SD",
    min_support = min_support, n_bins = n_bins, min_n = min_n,
    support_probs = support_probs, models = names(fit$models),
    assumptions = .cast_effect_assumptions())
}

#' Support Fraction of an Interventional Shift
#'
#' Reports, per driver and shift size, the fraction of observed predictor
#' vectors that remain inside the training support after the shift. This is the
#' positivity diagnostic that belongs beside any interventional effect
#' product: where it is low, the effect is answered by extrapolation.
#'
#' @inheritParams cast_dose_response
#' @param variables Predictors to assess. Default: every fitted predictor.
#' @param shift Numeric vector of shift sizes, in training standard deviations
#'   unless `shift_type = "raw"`.
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param newdata Optional observed predictor frame. Defaults to training data.
#' @param support_probs Length-2 numeric. Default `c(0.01, 0.99)`.
#' @param max_rows Rows sampled (evenly). Default 20000.
#'
#' @return A `cast_support` object with `support`: a data frame of
#'   `driver`, `shift`, `shift_raw`, `support`, and `per_variable` (a named
#'   list of per-predictor support fractions behind each joint number).
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
  .cast_reference_rows(fit, newdata, env_vars, max_rows)

  pieces <- lapply(variables, function(v) {
    sr <- if (identical(shift_type, "raw")) shift else shift * .cast_driver_sd(fit, v)
    frac <- .cast_support_fraction(ref, v, sr, probs = support_probs)
    data.frame(driver = v, shift = shift, shift_raw = sr,
               support = as.numeric(frac[v, ]),
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
  "Read as a causal effect only under consistency, no unobserved",
  "confounding given the adjustment set, and positivity; the support column",
  "reports the third. On presence/background data the scale is relative",
  "suitability, not occurrence probability.")

.cast_reference_rows <- function(fit, newdata, env_vars, max_rows) {
  X <- if (is.null(newdata)) {
    as.data.frame(fit$scaling$reference)
  } else {
    miss <- setdiff(env_vars, names(newdata))
    if (length(miss)) {
      cli::cli_abort("{.arg newdata} lacks fitted predictor{?s}: {.val {miss}}.")
    }
    as.data.frame(newdata[, env_vars, drop = FALSE])
  }
  X <- X[, env_vars, drop = FALSE]
  for (col in names(X)) X[[col]] <- as.numeric(X[[col]])
  X <- .cast_impute(X, fit$scaling$impute)
  ok <- rowSums(!is.finite(as.matrix(X))) == 0L
  X <- X[ok, , drop = FALSE]
  if (!nrow(X)) cli::cli_abort("No complete finite predictor rows to evaluate.")
  max_rows <- as.integer(max_rows)
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

.cast_bin_index <- function(x, n_bins) {
  qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n_bins + 1L),
                               na.rm = TRUE, names = FALSE))
  if (length(qs) < 2L) return(rep(1L, length(x)))
  idx <- findInterval(x, qs, all.inside = TRUE)
  as.integer(idx)
}
