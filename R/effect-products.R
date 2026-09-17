# ==========================================================================
# Interventional effect products: per-driver contrast sizes and rasters.
#   cast_effect_table(): per-driver interventional effect sizes on a table.
#   cast_effect_map():   per-driver delta-suitability rasters over a stack.
# Both answer one question: "if we intervene on driver X (shift it by each
# s in a symmetric shift set, all other drivers held fixed), what happens
# to predicted suitability?" Magnitude and direction are reported apart.
# Model-based interventional contrast (g-computation / standardization);
# requires consistency, no unobserved confounding given the adjustment set,
# and positivity. cast_effect_support() reports the third of those;
# cast_dose_response() reports the response shape; cast_effect_support()
# reports how much of the shifted range the observed data can answer.
# Not doubly robust and not TMLE; not a proof of a manipulable mechanism.
# Read alongside cast_necessity(): response without necessity does not
# identify a driver.
# ==========================================================================

.effect_predictions <- function(fit, X, driver, steps, base = NULL) {
  models <- names(fit$models)
  if (!length(models)) cli::cli_abort("No fitted models available.")
  if (is.null(base)) base <- .cast_predict_matrix(fit, X, models)
  aligned <- vapply(steps, function(st) {
    shifted <- X
    shifted[[driver]] <- shifted[[driver]] + st
    cf <- .cast_predict_matrix(fit, shifted, models)
    # Average paired model changes BEFORE taking the absolute value. Taking
    # abs per engine instead measures disagreement, not the ensemble effect.
    delta <- cf - base
    delta[!is.finite(cf) | !is.finite(base)] <- NA_real_
    rowMeans(delta, na.rm = TRUE) * sign(st)
  }, numeric(nrow(X)))
  matrix(aligned, nrow = nrow(X), ncol = length(steps))
}

.steps_from_fit <- function(fit, drivers, shifts, shift_type) {
  shifts <- as.numeric(shifts)
  shift_type <- match.arg(shift_type, c("sd", "raw"))
  if (!length(shifts) || any(!is.finite(shifts)) || any(shifts == 0)) {
    cli::cli_abort("{.arg shifts} must be non-zero finite numbers.")
  }
  if (identical(shift_type, "raw")) {
    return(stats::setNames(rep(list(shifts), length(drivers)), drivers))
  }
  sds <- fit$scaling$sds
  if (is.null(sds)) {
    cli::cli_abort(c(
      "The fit does not store predictor SDs, so SD-based steps cannot be derived.",
      "i" = "Provide {.arg steps} in raw units instead."))
  }
  miss <- setdiff(drivers, names(sds))
  if (length(miss)) cli::cli_abort("No stored SD for: {.val {miss}}.")
  stats::setNames(lapply(drivers, function(v) as.numeric(sds[[v]]) * shifts), drivers)
}

.check_steps <- function(steps, drivers) {
  if (is.numeric(steps)) {
    steps <- if (is.null(names(steps))) {
      stats::setNames(rep(list(as.numeric(steps)), length(drivers)), drivers)
    } else {
      lapply(as.list(steps), as.numeric)
    }
  }
  if (!is.list(steps)) {
    cli::cli_abort("{.arg steps} must be a named list (or numeric vector) of raw-unit shift sets.")
  }
  miss <- setdiff(drivers, names(steps))
  if (length(miss)) cli::cli_abort("{.arg steps} lacks entries for: {.val {miss}}.")
  steps <- steps[drivers]
  bad <- vapply(steps, function(s) !is.numeric(s) || !length(s) ||
                  any(!is.finite(s)) || any(s == 0), logical(1))
  if (any(bad)) cli::cli_abort("{.arg steps} entries must be non-zero finite numbers: {.val {drivers[bad]}}.")
  lapply(steps, as.numeric)
}

# V2 estimator: deltas are aligned to the direction of their own shift
# (delta * sign(step)) so a symmetric shift set does not cancel out.
.effect_summarise <- function(aligned, n_shifts, abs_values = NULL) {
  keep <- is.finite(aligned)
  if (!is.null(abs_values)) keep <- keep & is.finite(abs_values)
  if (!any(keep)) {
    return(data.frame(mean_abs_dHSS = NA_real_, mean_signed_dHSS = NA_real_,
                      median_signed_dHSS = NA_real_, pct_gain = NA_real_,
                      pct_loss = NA_real_, max_gain = NA_real_,
                      max_loss = NA_real_, n = 0L, n_shifts = n_shifts))
  }
  a <- aligned[keep]
  mag <- if (is.null(abs_values)) abs(a) else abs_values[keep]
  data.frame(mean_abs_dHSS = mean(mag), mean_signed_dHSS = mean(a),
             median_signed_dHSS = stats::median(a),
             pct_gain = mean(a > 0.05), pct_loss = mean(a < -0.05),
             max_gain = max(a), max_loss = min(a),
             n = length(a), n_shifts = n_shifts)
}

#' Model Response Table for Species Distribution Models
#'
#' Intervenes on each driver in turn (do-operator: shift the driver, hold
#' every other driver fixed), and summarises the resulting change in
#' predicted suitability on `newdata`: the model-based interventional
#' response of the ensemble to each driver.
#'
#' The effect is averaged over a symmetric shift set (default
#' \eqn{\pm 1, \pm 2} SD) rather than read off one arbitrary step, and is
#' reported as a magnitude (`mean_abs_dHSS`) separately from a direction
#' (`mean_signed_dHSS`). A single signed step at one magnitude is not a
#' stable driver ranking.
#'
#' @section Read with the necessity diagnostic:
#' Pair this table with [cast_necessity()] to describe predictive behavior.
#' The latter uses a separate RF estimator, whereas this table averages the
#' fitted engines equally. Neither agreement nor disagreement identifies a
#' causal effect. A small knockout cost can reflect redundancy, limited power,
#' estimator choice or AUC's insensitivity to calibration changes.
#' Causal interpretation additionally needs a justified adjustment set,
#' consistency, joint support and an adequate response model. With presence/
#' background data the output is relative suitability, not occurrence risk.
#'
#' @param fit A `cast_fit` object (from [cast_fit()] or [cast()]).
#' @param newdata Data frame with the fitted predictors; rows with missing
#'   predictors are dropped consistently with the raster API. Coordinates
#'   are not required.
#' @param drivers Character vector of drivers to intervene on. Default:
#'   every fitted predictor.
#' @param shifts Numeric vector of intervention sizes, in SD units unless
#'   `shift_type = "raw"`. Default `c(-2, -1, 1, 2)`: the effect is averaged
#'   over a symmetric shift set so it does not depend on the sign or the
#'   magnitude of one arbitrary step.
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param steps Optional raw-unit shift sets, overriding
#'   `shifts`/`shift_type`. Either a numeric vector (used for every driver)
#'   or a named list of numeric vectors, one per driver.
#' @param support_probs Length-2 numeric. Training quantiles treated as the
#'   supported range of each predictor. Default `c(0.01, 0.99)`. The `support`
#'   column reports the worst (minimum) such fraction over the driver's shift
#'   set; drivers with `support < 0.5` are answered mostly by extrapolation.
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return A `cast_effect_table` data.frame, one row per driver.
#'   `mean_abs_dHSS` is the effect magnitude (mean absolute delta over the
#'   shift set); `mean_signed_dHSS` is the direction (deltas aligned to the
#'   sign of their own shift, so positive means raising the driver raises
#'   suitability). Rank drivers by `mean_abs_dHSS`; read
#'   `mean_signed_dHSS` for the sign.
#'   Summaries first average over shifts per row; `n` counts complete rows.
#'   `outside_range_fraction` is the fraction of row/shift pairs outside the
#'   predictor's training range. It does not test multivariate support.
#'   `support` is the minimum 1-99% support fraction over the shift set, the
#'   same positivity diagnostic [cast_effect_support()] reports per shift.
#' @export
cast_effect_table <- function(fit, newdata, drivers = NULL,
                              shifts = c(-2, -1, 1, 2), shift_type = "sd",
                              steps = NULL, support_probs = c(0.01, 0.99),
                              verbose = TRUE) {
  if (!inherits(fit, "cast_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls cast_fit} object.")
  }
  env_vars <- fit$env_vars
  drivers <- drivers %||% env_vars
  if (!length(drivers) || anyDuplicated(drivers)) {
    cli::cli_abort("{.arg drivers} must be non-empty and unique.")
  }
  if (!all(drivers %in% env_vars)) {
    cli::cli_abort("Unknown drivers: {.val {setdiff(drivers, env_vars)}}.")
  }
  steps_given <- !is.null(steps)
  steps <- if (steps_given) {
    .check_steps(steps, drivers)
  } else {
    .steps_from_fit(fit, drivers, shifts, shift_type)
  }
  missing_vars <- setdiff(env_vars, names(newdata))
  if (length(missing_vars)) cli::cli_abort("Missing fitted predictors: {.val {missing_vars}}.")
  X <- as.data.frame(newdata[, env_vars, drop = FALSE])
  .cast_check_numeric_predictors(X, arg = "newdata")
  ok <- rowSums(!is.finite(as.matrix(X))) == 0L
  X <- X[ok, , drop = FALSE]
  if (!nrow(X)) cli::cli_abort("No complete finite predictor rows in {.arg newdata}.")
  base <- .cast_predict_matrix(fit, X, names(fit$models))
  ref <- tryCatch(.cast_reference(fit, env_vars), error = function(e) NULL)
  out <- vector("list", length(drivers))
  failed <- character(0)
  for (k in seq_along(drivers)) {
    v <- drivers[k]
    sv <- steps[[v]]
    if (verbose) cli::cli_inform("Intervening on {.val {v}} ({.val {sv}})...")
    aligned <- .effect_predictions(fit, X, v, sv, base)
    row <- .effect_summarise(rowMeans(aligned, na.rm = TRUE), length(sv),
                             rowMeans(abs(aligned), na.rm = TRUE))
    reference <- fit$scaling$reference[[v]]
    row$outside_range_fraction <- if (length(reference)) {
      lim <- range(reference, finite = TRUE)
      mean(vapply(sv, function(st) mean(X[[v]] + st < lim[1] |
                                       X[[v]] + st > lim[2]), numeric(1)))
    } else NA_real_
    # Positivity linkage: worst 1-99% support fraction over this driver's
    # own shift set, on the same scale cast_effect_support() reports.
    row$support <- if (!is.null(ref)) {
      frac <- suppressWarnings(tryCatch(
        .cast_support_fraction(ref, v, sv, probs = support_probs)[v, ],
        error = function(e) NA_real_))
      suppressWarnings(min(frac, na.rm = TRUE))
    } else NA_real_
    if (!is.finite(row$support)) row$support <- NA_real_
    if (row$n == 0L) failed <- c(failed, v)
    out[[k]] <- cbind(data.frame(driver = v, stringsAsFactors = FALSE), row)
  }
  if (length(failed)) {
    cli::cli_warn(c(
      "No finite interventional effect for {.val {failed}}.",
      "i" = "These drivers are reported as {.code NA}, not as zero effect."))
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "shift_type") <- if (steps_given) "raw" else shift_type
  attr(res, "steps") <- steps
  attr(res, "rows_input") <- nrow(newdata)
  attr(res, "rows_complete") <- nrow(X)
  attr(res, "intervention") <- sprintf(
    "do(driver += s) for s in {%s} %s, all other drivers fixed; deltas aligned to sign(s)",
    paste(if (steps_given) unique(unlist(steps)) else shifts, collapse = ", "),
    if (steps_given || identical(shift_type, "raw")) "raw units" else "SD")
  attr(res, "assumptions") <- paste(
    "Model-based interventional effect (g-computation). Assumes no",
    "unobserved confounding; not a proof of a manipulable mechanism.")
  class(res) <- c("cast_effect_table", "data.frame")
  res
}

#' @export
print.cast_effect_table <- function(x, ...) {
  cli::cli_text("{.strong Model response table} (predictor shifts, others fixed)")
  iv <- attr(x, "intervention")
  if (!is.null(iv)) cli::cli_text("{.emph {iv}}")
  print(as.data.frame(x))
  cli::cli_text("Rank by {.field mean_abs_dHSS}; read {.field mean_signed_dHSS} for direction; {.field support} < 0.5 means the shift is answered mostly by extrapolation.")
  cli::cli_text("Pair with {.fn cast_necessity} for predictive diagnostics; neither product establishes causal identification.")
  invisible(x)
}

#' Model Response Maps for Species Distribution Models
#'
#' Intervenes on each driver across every valid cell of `current_stack`
#' (all other drivers held fixed) and returns the spatially explicit
#' equally weighted ensemble-mean change in suitability. Spatial
#' counterpart of [cast_effect_table()], using the same symmetric shift set.
#'
#' @param fit A `cast_fit` object whose predictors are all layers of
#'   `current_stack`.
#' @param current_stack A `SpatRaster` containing every fitted predictor.
#' @param drivers Drivers to map. Default: every fitted predictor.
#' @param shifts Numeric vector of intervention sizes, in SD units unless
#'   `shift_type = "raw"`. Default `c(-2, -1, 1, 2)`.
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param steps Optional raw-unit shift sets, overriding
#'   `shifts`/`shift_type`. Either a numeric vector (used for every driver)
#'   or a named list of numeric vectors, one per driver.
#' @param block_rows Integer rows per processing block. Default 512.
#' @param filename Optional GeoTIFF path for the delta stack (LZW).
#' @param overwrite Logical. Overwrite `filename`. Default `FALSE`.
#' @param verbose Print block progress. Default `TRUE`.
#'
#' @return A `SpatRaster` with two layers per driver (NA where predictors
#'   are missing): `dHSS_<driver>`, the direction (deltas aligned to the
#'   sign of their own shift), and `absdHSS_<driver>`, the magnitude. The
#'   per-driver summary table is attached as attribute `effect_table`.
#' @export
cast_effect_map <- function(fit, current_stack, drivers = NULL,
                            shifts = c(-2, -1, 1, 2), shift_type = "sd",
                            steps = NULL, block_rows = 512L, filename = NULL,
                            overwrite = FALSE, verbose = TRUE) {
  if (!inherits(fit, "cast_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls cast_fit} object.")
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg terra} is required for effect maps.")
  }
  env_vars <- fit$env_vars
  if (!all(env_vars %in% names(current_stack))) {
    cli::cli_abort("{.arg current_stack} lacks fitted predictors: {.val {setdiff(env_vars, names(current_stack))}}.")
  }
  drivers <- drivers %||% env_vars
  if (!length(drivers) || anyDuplicated(drivers)) {
    cli::cli_abort("{.arg drivers} must be non-empty and unique.")
  }
  if (length(block_rows) != 1L || !is.finite(block_rows) ||
      block_rows < 1 || block_rows != as.integer(block_rows)) {
    cli::cli_abort("{.arg block_rows} must be a positive integer.")
  }
  if (!all(drivers %in% env_vars)) {
    cli::cli_abort("Unknown drivers: {.val {setdiff(drivers, env_vars)}}.")
  }
  steps_given <- !is.null(steps)
  steps <- if (steps_given) {
    .check_steps(steps, drivers)
  } else {
    .steps_from_fit(fit, drivers, shifts, shift_type)
  }

  n_cells <- terra::ncell(current_stack)
  nr <- terra::nrow(current_stack)
  ncl <- terra::ncol(current_stack)
  blank <- stats::setNames(rep(list(rep(NA_real_, n_cells)), length(drivers)), drivers)
  acc_signed <- blank
  acc_abs <- blank
  starts <- seq(1L, nr, by = block_rows)
  for (b in seq_along(starts)) {
    r0 <- starts[b]; r1 <- min(nr, r0 + block_rows - 1L)
    cells <- ((r0 - 1L) * ncl + 1L):(r1 * ncl)
    # Read only this block: reading the whole stack per block turned an
    # O(n) pass into O(n * n_blocks) on national grids.
    X <- terra::values(current_stack, mat = TRUE, row = r0,
                       nrows = r1 - r0 + 1L)[, env_vars, drop = FALSE]
    ok <- rowSums(!is.finite(X)) == 0L
    if (!any(ok)) next
    # predict.gbm() and mgcv's predict() reject a matrix, so hand every engine a
    # data.frame. Shifting one column by name also avoids copying the whole
    # block on each intervention step.
    Xi <- as.data.frame(X[ok, , drop = FALSE], check.names = FALSE)
    base <- .cast_predict_matrix(fit, Xi, names(fit$models))
    for (v in drivers) {
      aligned <- .effect_predictions(fit, Xi, v, steps[[v]], base)
      s <- rep(NA_real_, length(cells)); s[ok] <- rowMeans(aligned, na.rm = TRUE)
      a <- rep(NA_real_, length(cells)); a[ok] <- rowMeans(abs(aligned), na.rm = TRUE)
      acc_signed[[v]][cells] <- s
      acc_abs[[v]][cells] <- a
    }
    if (verbose) cli::cli_inform("block {b}/{length(starts)} done")
  }
  to_rast <- function(d) {
    r <- terra::rast(current_stack[[1]])
    terra::values(r) <- d
    r
  }
  out <- c(lapply(acc_signed, to_rast), lapply(acc_abs, to_rast))
  names(out) <- c(sprintf("dHSS_%s", drivers), sprintf("absdHSS_%s", drivers))
  res <- terra::rast(out)
  tab <- do.call(rbind, lapply(drivers, function(v) {
    cbind(data.frame(driver = v, stringsAsFactors = FALSE),
          .effect_summarise(acc_signed[[v]], length(steps[[v]]),
                            abs_values = acc_abs[[v]]))
  }))
  rownames(tab) <- NULL
  if (!is.null(filename)) {
    terra::writeRaster(res, filename, overwrite = overwrite,
                       gdal = c("COMPRESS=LZW", "TILED=YES"))
  }
  attr(res, "effect_table") <- tab
  attr(res, "steps") <- steps
  attr(res, "intervention") <- sprintf(
    "do(driver += s) for s in {%s} %s, all other drivers fixed; deltas aligned to sign(s)",
    paste(if (steps_given) unique(unlist(steps)) else shifts, collapse = ", "),
    if (steps_given || identical(shift_type, "raw")) "raw units" else "SD")
  attr(res, "assumptions") <- paste(
    "Model-based interventional effect (g-computation). Assumes no",
    "unobserved confounding; not a proof of a manipulable mechanism.")
  res
}
