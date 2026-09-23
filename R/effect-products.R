# Fitted response contrasts require additional assumptions for causal interpretation.

# Single-shift stub: one raw-unit shift per driver, signed delta, paired
# model averaging BEFORE any absolute value (per-engine abs would measure
# disagreement, not the ensemble effect).
.effect_predictions <- function(fit, X, driver, step, base = NULL) {
  models <- names(fit$models)
  if (!length(models)) cli::cli_abort("No fitted models available.")
  if (length(step) != 1L || !is.finite(step) || step == 0) {
    cli::cli_abort("{.arg step} must be one non-zero finite raw-unit shift.")
  }
  if (is.null(base)) base <- .cast_predict_matrix(fit, X, models)
  shifted <- X
  shifted[[driver]] <- shifted[[driver]] + step
  cf <- .cast_predict_matrix(fit, shifted, models)
  delta <- cf - base
  delta[!is.finite(cf) | !is.finite(base)] <- NA_real_
  rowMeans(delta, na.rm = TRUE)
}

# Resolve one raw-unit shift per driver. Accepts a single unnamed numeric
# (applied to every driver; only sensible when drivers share units), a named
# numeric vector, or a named list of length-1 numerics. SD units are converted
# through the fit's stored training SDs.
.shift_from_fit <- function(fit, drivers, shift, shift_type) {
  shift_type <- match.arg(shift_type, c("raw", "sd"))
  if (is.list(shift)) {
    miss <- setdiff(drivers, names(shift))
    if (length(miss)) cli::cli_abort("{.arg shift} lacks entries for: {.val {miss}}.")
    out <- lapply(drivers, function(v) {
      s <- shift[[v]]
      if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s == 0) {
        cli::cli_abort("{.arg shift} entries must be one non-zero finite number per driver: {.val {v}}.")
      }
      as.numeric(s)
    })
    names(out) <- drivers
  } else {
    nm <- names(shift)
    shift <- as.numeric(shift)
    if (!length(shift) || any(!is.finite(shift)) || any(shift == 0)) {
      cli::cli_abort("{.arg shift} must be non-zero finite numbers.")
    }
    if (!is.null(nm)) {
      miss <- setdiff(drivers, nm)
      if (length(miss)) cli::cli_abort("{.arg shift} lacks entries for: {.val {miss}}.")
      out <- lapply(drivers, function(v) shift[match(v, nm)])
      names(out) <- drivers
    } else {
      if (length(shift) != 1L) {
        cli::cli_abort(c(
          "Unnamed {.arg shift} must be length 1 (applied to every driver).",
          "i" = "Use a named vector or list to give each driver its own raw-unit shift."))
      }
      out <- stats::setNames(rep(list(as.numeric(shift)), length(drivers)), drivers)
    }
  }
  if (identical(shift_type, "sd")) {
    sds <- fit$scaling$sds
    if (is.null(sds)) {
      cli::cli_abort(c(
        "The fit does not store predictor SDs, so SD-based shifts cannot be derived.",
        "i" = "Provide {.arg shift} in raw units instead."))
    }
    miss <- setdiff(drivers, names(sds))
    if (length(miss)) cli::cli_abort("No stored SD for: {.val {miss}}.")
    out <- lapply(drivers, function(v) as.numeric(sds[[v]]) * out[[v]])
    names(out) <- drivers
  }
  out
}

# Masked mean over supported rows only. Returns NA estimates when nothing is
# supported or the driver is masked by the coverage rule.
.effect_masked_mean <- function(delta, supported) {
  d <- delta[is.finite(delta) & supported]
  if (!length(d)) {
    return(list(mean_signed = NA_real_, mean_abs = NA_real_, n_supported = 0L))
  }
  list(mean_signed = mean(d), mean_abs = mean(abs(d)), n_supported = length(d))
}

#' Shift Effect Table for Species Distribution Models
#'
#' Intervenes on each driver in turn with a single raw-unit shift (do-operator:
#' shift the driver, hold every other driver fixed), and summarises the
#' resulting change in predicted suitability on `newdata`: the model-based
#' shift effect (shift estimand) of the ensemble for each driver.
#'
#' Rows outside the training quantile box - at baseline or after shifting -
#' are masked: estimates average over supported rows only, and drivers whose
#' box coverage falls below `min_support` report `NA` estimates instead of an
#' extrapolated ranking. A zero shift is not an intervention and is refused,
#' as is a shift that leaves the training range everywhere.
#'
#' @section Interpretation (read before citing):
#' The table is a *model-based* shift contrast (g-computation /
#' standardization), reported on the relative suitability scale with
#' presence/background data. A causal reading additionally needs consistency
#' of the stated shift, a justified adjustment set, no uncontrolled
#' confounding, joint support and adequate response and observation models.
#' Quantile-box masking only diagnoses range extrapolation; it does not
#' establish joint positivity or causal identification.
#'
#' @param fit A `cast_fit` object (from [cast_fit()] or [cast()]).
#' @param newdata Data frame with the fitted predictors; rows with missing
#'   predictors are dropped. Coordinates are not required.
#' @param drivers Character vector of drivers to intervene on. Default:
#'   every fitted predictor.
#' @param shift Single raw-unit shift per driver (e.g. `list(bio01 = 2)` for
#'   +2C). Either one unnamed numeric value applied to every driver, or a
#'   named numeric vector / named list with one value per driver. Default `1`.
#' @param shift_type `"raw"` (default, predictor units) or `"sd"` (training
#'   standard deviations, converted through the fit).
#' @param support_probs Length-2 numeric. Increasing training quantiles in
#'   `[0, 1]` defining a box across all predictors. Default `c(0.01, 0.99)`.
#'   This diagnostic cannot detect holes inside the box and does not establish
#'   joint positivity.
#' @param min_support Numeric in `[0, 1]`. Drivers with box coverage below
#'   this threshold are masked (`NA` estimates). Default `0.5`.
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return A `cast_effect_table` data.frame, one row per driver:
#'   `shift_raw` (applied raw-unit shift), `mean_abs_dHSS` (magnitude: mean
#'   absolute delta over supported rows), `mean_signed_dHSS` (direction:
#'   mean signed delta), `n` (complete rows), `n_supported` (rows inside the
#'   box before and after shifting), `support` (`n_supported / n`), and
#'   `masked` (box coverage below `min_support`). Rank supported drivers by
#'   `mean_abs_dHSS`; read `mean_signed_dHSS` for the sign.
#' @export
cast_effect_table <- function(fit, newdata, drivers = NULL,
                              shift = 1, shift_type = "raw",
                              support_probs = c(0.01, 0.99),
                              min_support = 0.5,
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
  if (!is.numeric(min_support) || length(min_support) != 1L ||
      !is.finite(min_support) || min_support < 0 || min_support > 1) {
    cli::cli_abort("{.arg min_support} must be one finite number in [0, 1].")
  }
  .cast_check_support_probs(support_probs)
  steps <- .shift_from_fit(fit, drivers, shift, shift_type)
  missing_vars <- setdiff(env_vars, names(newdata))
  if (length(missing_vars)) cli::cli_abort("Missing fitted predictors: {.val {missing_vars}}.")
  X <- as.data.frame(newdata[, env_vars, drop = FALSE])
  .cast_check_numeric_predictors(X, arg = "newdata")
  ok <- rowSums(!is.finite(as.matrix(X))) == 0L
  X <- X[ok, , drop = FALSE]
  if (!nrow(X)) cli::cli_abort("No complete finite predictor rows in {.arg newdata}.")
  base <- .cast_predict_matrix(fit, X, names(fit$models))
  ref <- tryCatch(.cast_reference(fit, env_vars), error = function(e) NULL)
  bounds <- if (!is.null(ref)) .cast_support_bounds(ref, support_probs) else NULL
  out <- vector("list", length(drivers))
  for (k in seq_along(drivers)) {
    v <- drivers[k]
    st <- steps[[v]]
    .cast_check_shift_reachable(ref, v, st)
    if (verbose) cli::cli_inform("Intervening on {.val {v}} ({.val {st}} raw units)...")
    delta <- .effect_predictions(fit, X, v, st, base)
    covered <- if (!is.null(bounds)) {
      as.logical(.cast_support_rows(X, v, st, bounds)[, 1L])
    } else {
      rep(NA, nrow(X))
    }
    covered[!is.finite(delta)] <- FALSE
    est <- .effect_masked_mean(delta, is.finite(covered) & covered)
    support <- if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE)
    masked <- is.na(support) || support < min_support
    out[[k]] <- data.frame(
      driver = v, shift_raw = st,
      mean_abs_dHSS = if (masked) NA_real_ else est$mean_abs,
      mean_signed_dHSS = if (masked) NA_real_ else est$mean_signed,
      n = nrow(X), n_supported = est$n_supported,
      support = support, masked = masked,
      stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "shift_type") <- shift_type
  attr(res, "steps") <- steps
  attr(res, "min_support") <- min_support
  attr(res, "support_probs") <- support_probs
  attr(res, "rows_input") <- nrow(newdata)
  attr(res, "rows_complete") <- nrow(X)
  attr(res, "intervention") <- "do(driver += shift_raw) in raw units, all other drivers fixed"
  attr(res, "assumptions") <- .cast_effect_assumptions()
  class(res) <- c("cast_effect_table", "data.frame")
  res
}

#' @export
print.cast_effect_table <- function(x, ...) {
  cli::cli_text("{.strong Shift effect table} (single raw-unit shift, others fixed)")
  iv <- attr(x, "intervention")
  if (!is.null(iv)) cli::cli_text("{.emph {iv}}")
  print(as.data.frame(x))
  cli::cli_text("Rank supported drivers by {.field mean_abs_dHSS}; read {.field mean_signed_dHSS} for direction. Masked drivers ({.field masked}) carry {.code NA} estimates: the shift is not answerable inside the training box.")
  invisible(x)
}

#' Shift Effect Maps for Species Distribution Models
#'
#' Intervenes on each driver with a single raw-unit shift across every valid
#' cell of `current_stack` (all other drivers held fixed) and returns the
#' spatially explicit equally weighted ensemble-mean change in suitability.
#' Spatial counterpart of [cast_effect_table()] with the same hard
#' range-masking: cells outside the training quantile box - at baseline or
#' after shifting - carry `NA` effects.
#'
#' @param fit A `cast_fit` object whose predictors are all layers of
#'   `current_stack`.
#' @param current_stack A `SpatRaster` containing every fitted predictor.
#' @param drivers Drivers to map. Default: every fitted predictor.
#' @param shift Single raw-unit shift per driver. Either one unnamed numeric
#'   value applied to every driver, or a named numeric vector / named list
#'   with one value per driver. Default `1`.
#' @param shift_type `"raw"` (default) or `"sd"`.
#' @param block_rows Integer rows per processing block. Default 512.
#' @param filename Optional GeoTIFF path for the effect and coverage stack (LZW).
#' @param overwrite Logical. Overwrite `filename`. Default `FALSE`.
#' @param verbose Print block progress. Default `TRUE`.
#' @param support_probs Length-2 numeric. Increasing training quantiles in
#'   `[0, 1]` defining a box across all predictors. Default `c(0.01, 0.99)`.
#' @param min_support Numeric in `[0, 1]`. Drivers with grid coverage below
#'   this threshold report `NA` in the attached table. Default `0.5`.
#'
#' @return A `SpatRaster` with three layers per driver (NA where any fitted
#'   predictor is missing or non-finite, or where the cell is masked):
#'   `dHSS_<driver>` (signed change), `absdHSS_<driver>` (magnitude), and
#'   `support_<driver>` (1 where the cell's baseline and shifted vectors both
#'   lie inside the training quantile box, 0 otherwise).
#'
#'   The per-driver summary table is attached as attribute `effect_table`
#'   with the same columns as [cast_effect_table()] evaluated on the grid.
#'   Masking is hard: masked cells never contribute to the table means.
#'   Coverage is a range-extrapolation diagnostic, not a test of conditional
#'   positivity or causal identification; holes inside the box are not
#'   detected.
#' @export
cast_effect_map <- function(fit, current_stack, drivers = NULL,
                            shift = 1, shift_type = "raw",
                            block_rows = 512L, filename = NULL,
                            overwrite = FALSE, verbose = TRUE,
                            support_probs = c(0.01, 0.99),
                            min_support = 0.5) {
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
  if (!is.numeric(min_support) || length(min_support) != 1L ||
      !is.finite(min_support) || min_support < 0 || min_support > 1) {
    cli::cli_abort("{.arg min_support} must be one finite number in [0, 1].")
  }
  .cast_check_support_probs(support_probs)
  steps <- .shift_from_fit(fit, drivers, shift, shift_type)

  ref <- tryCatch(.cast_reference(fit, env_vars), error = function(e) NULL)
  bounds <- if (!is.null(ref)) .cast_support_bounds(ref, support_probs) else NULL
  # Degenerate bounds (no finite training values) leave coverage unknown:
  # nothing is counted as supported and every driver masks.
  bounds_ok <- !is.null(bounds) &&
    all(vapply(bounds, function(b) all(is.finite(b)), logical(1)))
  for (v in drivers) {
    if (!is.null(ref)) .cast_check_shift_reachable(ref, v, steps[[v]])
  }
  n_cells <- terra::ncell(current_stack)
  nr <- terra::nrow(current_stack)
  ncl <- terra::ncol(current_stack)
  blank <- stats::setNames(rep(list(rep(NA_real_, n_cells)), length(drivers)), drivers)
  acc_signed <- blank
  acc_abs <- blank
  acc_support <- blank
  n_supported <- stats::setNames(rep(0, length(drivers)), drivers)
  n_complete <- 0
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
    n_complete <- n_complete + nrow(Xi)
    base <- .cast_predict_matrix(fit, Xi, names(fit$models))
    for (v in drivers) {
      delta <- .effect_predictions(fit, Xi, v, steps[[v]], base)
      covered <- if (bounds_ok) {
        as.logical(.cast_support_rows(Xi, v, steps[[v]], bounds)[, 1L])
      } else {
        # Unknown coverage: never counted as supported.
        rep(NA, nrow(Xi))
      }
      covered[!is.finite(delta)] <- FALSE
      s <- rep(NA_real_, length(cells))
      a <- rep(NA_real_, length(cells))
      cvg <- rep(NA_real_, length(cells))
      s[ok] <- ifelse(covered, delta, NA_real_)
      a[ok] <- ifelse(covered, abs(delta), NA_real_)
      cvg[ok] <- as.numeric(covered)
      acc_signed[[v]][cells] <- s
      acc_abs[[v]][cells] <- a
      acc_support[[v]][cells] <- cvg
      n_supported[[v]] <- n_supported[[v]] + sum(covered, na.rm = TRUE)
    }
    if (verbose) cli::cli_inform("block {b}/{length(starts)} done")
  }
  to_rast <- function(d) {
    r <- terra::rast(current_stack[[1]])
    terra::values(r) <- d
    r
  }
  out <- c(lapply(acc_signed, to_rast), lapply(acc_abs, to_rast),
           lapply(acc_support, to_rast))
  names(out) <- c(sprintf("dHSS_%s", drivers), sprintf("absdHSS_%s", drivers),
                  sprintf("support_%s", drivers))
  res <- terra::rast(out)
  tab <- do.call(rbind, lapply(drivers, function(v) {
    svec <- acc_signed[[v]]
    avec <- acc_abs[[v]]
    keep <- is.finite(svec)
    support <- if (!bounds_ok || n_complete <= 0) {
      NA_real_
    } else {
      n_supported[[v]] / n_complete
    }
    masked <- is.na(support) || support < min_support
    data.frame(
      driver = v, shift_raw = steps[[v]],
      mean_abs_dHSS = if (masked || !any(keep)) NA_real_ else mean(avec[keep], na.rm = TRUE),
      mean_signed_dHSS = if (masked || !any(keep)) NA_real_ else mean(svec[keep], na.rm = TRUE),
      n = n_complete, n_supported = n_supported[[v]],
      support = support, masked = masked,
      stringsAsFactors = FALSE)
  }))
  rownames(tab) <- NULL
  if (!is.null(filename)) {
    terra::writeRaster(res, filename, overwrite = overwrite,
                       gdal = c("COMPRESS=LZW", "TILED=YES"))
  }
  attr(res, "effect_table") <- tab
  attr(res, "steps") <- steps
  attr(res, "support_probs") <- support_probs
  attr(res, "min_support") <- min_support
  attr(res, "intervention") <- "do(driver += shift_raw) in raw units, all other drivers fixed; masked cells carry NA"
  attr(res, "assumptions") <- .cast_effect_assumptions()
  res
}
