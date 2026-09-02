# ==========================================================================
# Causal effect products: the interpretive core of castSDM.
#   cast_effect_table(): per-driver interventional effect sizes on a table.
#   cast_effect_map():   per-driver delta-suitability rasters over a stack.
# Both answer one question: "if we intervene on driver X (+shift SD, all
# other drivers held fixed), what happens to predicted suitability?"
# Model-based interventional effect (g-computation); requires the standard
# no-unobserved-confounding assumption; not a proof of a manipulable
# causal mechanism.
# ==========================================================================

.pred_num_engine <- function(engine, model, X) {
  p <- switch(engine,
    rf = {
      pp <- stats::predict(model, X)$predictions
      if (is.matrix(pp)) pp <- pp[, "1"]
      as.numeric(pp)
    },
    brt = as.numeric(stats::predict(model, X, n.trees = model$n.trees,
                                    type = "response")),
    gam = as.numeric(stats::predict(model, X, type = "response")),
    maxent = as.numeric(stats::predict(model, X, type = "cloglog",
                                       clamp = FALSE)),
    cli::cli_abort("Unsupported engine {.val {engine}} for effect products.")
  )
  as.numeric(p)
}

.steps_from_fit <- function(fit, drivers, shift, shift_type) {
  if (identical(shift_type, "raw")) return(stats::setNames(rep(shift, length(drivers)), drivers))
  sds <- fit$scaling$sds
  if (is.null(sds)) {
    cli::cli_abort(c(
      "The fit does not store predictor SDs, so SD-based steps cannot be derived.",
      "i" = "Provide {.arg steps} in raw units instead."))
  }
  miss <- setdiff(drivers, names(sds))
  if (length(miss)) cli::cli_abort("No stored SD for: {.val {miss}}.")
  stats::setNames(unname(sds[drivers]) * shift, drivers)
}

#' Causal Effect Table for Species Distribution Models
#'
#' Intervenes on each driver in turn (do-operator: shift the driver, hold
#' every other driver fixed), and summarises the resulting change in
#' predicted suitability on `newdata`: the model-based interventional
#' response of the ensemble to each driver.
#'
#' @param fit A `cast_fit` object (from [cast_fit()] or [cast()]).
#' @param newdata Data frame with the fitted predictors; rows with missing
#'   predictors are dropped per model.
#' @param drivers Character vector of drivers to intervene on. Default:
#'   every fitted predictor.
#' @param shift Numeric. Intervention size. Default 1 (in SD units unless
#'   `steps` is given).
#' @param shift_type `"sd"` (default) or `"raw"`.
#' @param steps Optional named numeric vector of raw-unit intervention
#'   sizes (overrides `shift`/`shift_type`).
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return A `cast_effect_table` data.frame: one row per driver with mean,
#'   median, max, min delta suitability and the share of rows gaining
#'   (> +0.05) or losing (< -0.05).
#' @export
cast_effect_table <- function(fit, newdata, drivers = NULL,
                              shift = 1, shift_type = "sd", steps = NULL,
                              verbose = TRUE) {
  if (!inherits(fit, "cast_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls cast_fit} object.")
  }
  env_vars <- fit$env_vars
  drivers <- drivers %||% env_vars
  if (!all(drivers %in% env_vars)) {
    cli::cli_abort("Unknown drivers: {.val {setdiff(drivers, env_vars)}}.")
  }
  steps <- steps %||% .steps_from_fit(fit, drivers, shift, shift_type)
  out <- vector("list", length(drivers))
  for (k in seq_along(drivers)) {
    v <- drivers[k]
    if (verbose) cli::cli_inform("Intervening on {.val {v}} ({steps[[k]]})...")
    s <- tryCatch(
      cast_sensitivity(fit, newdata = newdata, variable = v,
                       shift = steps[[k]], shift_type = "raw"),
      error = function(e) NULL
    )
    d <- if (is.null(s)) NA_real_ else s$predictions$delta_hss
    d <- d[is.finite(d)]
    out[[k]] <- if (!length(d)) {
      data.frame(driver = v, mean_dHSS = NA_real_, median_dHSS = NA_real_,
                 pct_gain = NA_real_, pct_loss = NA_real_, max_gain = NA_real_,
                 max_loss = NA_real_, n = 0L)
    } else {
      data.frame(driver = v, mean_dHSS = mean(d), median_dHSS = stats::median(d),
                 pct_gain = mean(d > 0.05), pct_loss = mean(d < -0.05),
                 max_gain = max(d), max_loss = min(d), n = length(d))
    }
  }
  res <- do.call(rbind, out)
  attr(res, "intervention") <- "do(driver = +1 step), all other drivers fixed"
  attr(res, "assumptions") <- paste(
    "Model-based interventional effect (g-computation). Assumes no",
    "unobserved confounding; not a proof of a manipulable mechanism.")
  class(res) <- c("cast_effect_table", "data.frame")
  res
}

#' @export
print.cast_effect_table <- function(x, ...) {
  cli::cli_text("{.strong Causal effect table} (do-intervention, others fixed)")
  print(as.data.frame(x))
  invisible(x)
}

#' Causal Effect Maps for Species Distribution Models
#'
#' Intervenes on each driver across every valid cell of `current_stack`
#' (+`steps` raw units, all other drivers held fixed) and returns the
#' spatially explicit ensemble-mean change in suitability: the causal
#' effect map. Spatial counterpart of [cast_effect_table()].
#'
#' @param fit A `cast_fit` object whose predictors are all layers of
#'   `current_stack`.
#' @param current_stack A `SpatRaster` containing every fitted predictor.
#' @param drivers Drivers to map. Default: every fitted predictor.
#' @param shift,shift_type Intervention size and unit (SD-based unless
#'   `steps` is given).
#' @param steps Optional named numeric vector of raw-unit intervention
#'   sizes per driver (overrides `shift`/`shift_type`).
#' @param block_rows Integer rows per processing block. Default 512.
#' @param filename Optional GeoTIFF path for the delta stack (LZW).
#' @param overwrite Logical. Overwrite `filename`. Default `FALSE`.
#' @param verbose Print block progress. Default `TRUE`.
#'
#' @return A `SpatRaster` with one `dHSS_<driver>` layer per driver (NA
#'   where predictors are missing). The per-driver summary table is
#'   attached as attribute `effect_table`.
#' @export
cast_effect_map <- function(fit, current_stack, drivers = NULL,
                            shift = 1, shift_type = "sd", steps = NULL,
                            block_rows = 512L, filename = NULL,
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
  if (!all(drivers %in% env_vars)) {
    cli::cli_abort("Unknown drivers: {.val {setdiff(drivers, env_vars)}}.")
  }
  steps <- steps %||% .steps_from_fit(fit, drivers, shift, shift_type)
  eng <- vapply(fit$models, function(m) m$name, character(1))
  mdl <- lapply(fit$models, function(m) m$model)

  n_cells <- terra::ncell(current_stack)
  nr <- terra::nrow(current_stack)
  acc <- setNames(rep(list(rep(NA_real_, n_cells)), length(drivers)), drivers)
  starts <- seq(1L, nr, by = block_rows)
  for (b in seq_along(starts)) {
    r0 <- starts[b]; r1 <- min(nr, r0 + block_rows - 1L)
    cells <- ((r0 - 1L) * terra::ncol(current_stack) + 1L):(r1 * terra::ncol(current_stack))
    X <- terra::values(current_stack, mat = TRUE)[cells, env_vars, drop = FALSE]
    ok <- stats::complete.cases(X)
    if (!any(ok)) next
    Xi <- X[ok, , drop = FALSE]
    dsum <- matrix(0, nrow = nrow(Xi), ncol = length(drivers),
                   dimnames = list(NULL, drivers))
    for (m in seq_along(eng)) {
      p0 <- .pred_num_engine(eng[[m]], mdl[[m]], Xi)
      for (v in drivers) {
        Xs <- Xi; Xs[, v] <- Xs[, v] + steps[[v]]
        dsum[, v] <- dsum[, v] + (.pred_num_engine(eng[[m]], mdl[[m]], Xs) - p0)
      }
    }
    dsum <- dsum / length(eng)
    for (v in drivers) {
      tmp <- rep(NA_real_, length(cells)); tmp[ok] <- dsum[, v]
      acc[[v]][cells] <- tmp
    }
    if (verbose) cli::cli_inform("block {b}/{length(starts)} done")
  }
  out <- lapply(acc, function(d) {
    r <- terra::rast(current_stack[[1]])
    terra::values(r) <- d
    r
  })
  names(out) <- sprintf("dHSS_%s", drivers)
  res <- terra::rast(out)
  tab <- do.call(rbind, lapply(drivers, function(v) {
    d <- acc[[v]]; d <- d[is.finite(d)]
    data.frame(driver = v, mean_dHSS = mean(d),
               pct_gain = mean(d > 0.05), pct_loss = mean(d < -0.05))
  }))
  if (!is.null(filename)) {
    terra::writeRaster(res, filename, overwrite = overwrite,
                       gdal = c("COMPRESS=LZW", "TILED=YES"))
  }
  attr(res, "effect_table") <- tab
  attr(res, "assumptions") <- paste(
    "Model-based interventional effect (g-computation). Assumes no",
    "unobserved confounding; not a proof of a manipulable mechanism.")
  res
}
