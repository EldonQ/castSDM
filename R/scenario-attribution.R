# Paired environmental scenarios share one fitted model and frozen weights.

.cast_attribution_weights <- function(fit, cv, models, weights, method,
                                      min_score, decay, fallback,
                                      min_metric_folds) {
  models <- if (is.null(models)) names(fit$models) else models
  if (!length(models) || anyDuplicated(models) ||
      !all(models %in% names(fit$models))) {
    cli::cli_abort("{.arg models} must name unique fitted models.")
  }
  if (!is.null(weights)) {
    if (!is.numeric(weights) || is.null(names(weights)) ||
        anyDuplicated(names(weights)) || !all(models %in% names(weights)) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      cli::cli_abort("{.arg weights} must be finite non-negative named model weights.")
    }
    w <- weights[models]
  } else if (is.null(cv)) {
    if (method != "equal") {
      cli::cli_abort("Supply {.arg cv} or {.arg weights} for a weighted/best attribution.")
    }
    w <- stats::setNames(rep(1, length(models)), models)
  } else {
    if (!all(models %in% cv$metrics$model)) {
      cli::cli_abort("{.arg cv} lacks evaluation metrics for the requested models.")
    }
    scores <- .cast_ensemble_scores(cv$metrics, models, min_metric_folds)
    pw <- .cast_power_scores(scores, decay, method, min_score)
    w <- .cast_ensemble_weights(pw$scores, method, pw$min_score, fallback)
  }
  if (sum(w) <= 0) cli::cli_abort("At least one model must have positive weight.")
  w <- w[w > 0]
  w / sum(w)
}

.cast_attribution_inside <- function(X, bounds) {
  inside <- rep(TRUE, nrow(X))
  for (v in names(bounds)) {
    b <- bounds[[v]]
    inside <- inside & is.finite(X[[v]]) & X[[v]] >= b[1] & X[[v]] <= b[2]
  }
  inside
}

.cast_attribution_block <- function(fit, current, future, drivers, weights,
                                    bounds, reference = NULL, clamp = FALSE) {
  if (isTRUE(clamp) && !is.null(reference)) {
    current <- .cast_clamp(current, reference)
    future <- .cast_clamp(future, reference)
  }
  predict <- function(X) {
    p <- vapply(names(weights), function(m) {
      z <- tryCatch(predict_single_model(fit$models[[m]], X), error = function(e) {
        cli::cli_abort("Attribution prediction failed for {.val {m}}: {conditionMessage(e)}")
      })
      if (!is.numeric(z) || length(z) != nrow(X)) {
        cli::cli_abort("Model {.val {m}} returned an invalid prediction vector.")
      }
      as.numeric(z)
    }, numeric(nrow(X)))
    matrix(p, nrow(X), length(weights), dimnames = list(NULL, names(weights)))
  }
  # Requiring the same finite active models preserves a paired comparison.
  average <- function(p) {
    z <- as.numeric(p %*% weights)
    z[rowSums(!is.finite(p)) > 0] <- NA_real_
    z
  }
  spread <- function(p) {
    z <- .ensemble_rowsd(p)
    z[rowSums(!is.finite(p)) > 0] <- NA_real_
    z
  }
  base <- predict(current)
  fut <- predict(future)
  total <- fut - base
  bc <- .cast_attribution_inside(current, bounds)
  fc <- .cast_attribution_inside(future, bounds)
  out <- cbind(current = average(base), future = average(fut),
               total_change = average(total),
               support_total = as.numeric(bc & fc),
               total_sd = spread(total))
  effects <- matrix(NA_real_, nrow(current), length(drivers))
  changed <- names(current)[!vapply(names(current), function(v)
    identical(current[[v]], future[[v]]), logical(1))]
  for (j in seq_along(drivers)) {
    v <- drivers[j]
    only <- current; only[[v]] <- future[[v]]
    held <- future; held[[v]] <- current[[v]]
    # Unchanged drivers avoid two redundant model evaluations.
    if (identical(current[[v]], future[[v]])) {
      po <- base; ph <- fut
    } else if (identical(changed, v)) {
      po <- fut; ph <- base
    } else {
      po <- predict(only); ph <- predict(held)
    }
    effect <- po - base
    retained <- fut - ph
    effects[, j] <- average(effect)
    z <- cbind(average(effect), average(retained),
               as.numeric(bc & .cast_attribution_inside(only, bounds)),
               as.numeric(fc & .cast_attribution_inside(held, bounds)),
               spread(effect), spread(retained))
    colnames(z) <- paste0(c("effect_", "retained_", "support_effect_",
                           "support_retained_", "effect_sd_", "retained_sd_"), v)
    out <- cbind(out, z)
  }
  cbind(out, remainder = out[, "total_change"] - rowSums(effects))
}

#' Attribute Habitat Change to Environmental Scenarios
#'
#' Compares aligned current and future environments with one fitted SDM and
#' frozen ensemble weights. For each driver, computes its future-value effect
#' on the current environment and its contribution relative to holding that
#' driver at its current values in the otherwise future environment.
#'
#' @param fit A [cast_fit] object. Use [cast_as_fit()] for models fitted by
#'   other tools with an explicit prediction callback and training reference.
#' @param current,future Aligned predictor data frames or `terra::SpatRaster`
#'   stacks. Tables with `lon`/`lat` are matched by coordinates; otherwise rows
#'   must already be paired. Rasters must have identical geometry and CRS.
#' @param drivers Unique fitted predictor names to compare. Default all.
#' @param cv Optional [cast_cv] object for ensemble weights.
#' @param weights Optional named non-negative model weights, overriding `cv`.
#' @param models Optional subset of fitted model names.
#' @param method Ensemble strategy: `"equal"`, `"weighted"` or `"best"`.
#' @param min_score,decay,fallback,min_metric_folds Weight configuration as in
#'   [cast_ensemble()]. Defaults: `0.5`, `1`, `"best"`, `2`.
#' @param support_probs Training quantiles used for range flags. Default
#'   `c(0.01, 0.99)`. Raw contrasts are retained and supported summaries are
#'   reported separately; range flags do not establish ecological causality.
#' @param clamp Logical. Clamp current and future predictors to the training
#'   range before prediction, matching [cast_project_raster()]. Default
#'   `FALSE` (raw scenario values). Use `TRUE` when attribution maps must share
#'   the same extrapolation rule as projected HSS maps.
#' @param block_size Maximum prediction rows per block. Default `50000`.
#' @param filename GeoTIFF output path for raster inputs. When omitted a
#'   temporary file is used. Ignored for tables.
#' @param mask Optional raster mask; non-finite or zero mask cells are excluded.
#' @param overwrite Whether to replace an existing raster output.
#' @param verbose Whether to report progress.
#'
#' @return A `cast_attribution` list with `summary`, frozen `weights`,
#'   `drivers`, `data` (tables) or `raster_file` (rasters), and `layers`.
#'   `effect_<driver>` is driver-only future minus current; `retained_<driver>`
#'   is full future minus future with that driver held current. `remainder`
#'   is total change minus the sum of driver-only effects: it includes joint
#'   response and changes in any omitted drivers. `*_sd` layers describe
#'   cross-model disagreement, not sampling uncertainty or confidence intervals.
#'
#' @details Comparisons are on the fitted suitability scale, using raw scenario
#'   values without clamping or independent rescaling. Single-driver contrasts
#'   need not add to total change. Missing predictor pairs remain missing.
#'   Active models are fixed; a non-finite active-model prediction makes that
#'   comparison missing rather than switching weights. Scenario attribution is
#'   a predictive explanation, not identification of ecological causal effects,
#'   realized occupancy, colonisation, or dispersal. Presence/background scores
#'   remain relative suitability. All summary means weight rows/cells equally;
#'   use an equal-area grid when interpreting raster summaries spatially.
#' @seealso [cast_effect_map()], [cast_project()], [plot.cast_attribution()]
#' @export
cast_attribution <- function(fit, current, future, drivers = NULL, cv = NULL,
                             weights = NULL, models = NULL,
                             method = c("equal", "weighted", "best"),
                             min_score = 0.5, decay = 1, fallback = "best",
                             min_metric_folds = 2L,
                             support_probs = c(0.01, 0.99), clamp = FALSE,
                             block_size = 50000L,
                             filename = NULL, mask = NULL, overwrite = FALSE,
                             verbose = TRUE) {
  .cast_check_fit(fit)
  method <- match.arg(method)
  decay <- .cast_check_decay(decay)
  .cast_check_scalar(block_size, "block_size", lower = 1, integer = TRUE)
  vars <- fit$env_vars
  if (is.null(drivers)) drivers <- vars
  if (!is.character(drivers) || !length(drivers) || anyNA(drivers) ||
      anyDuplicated(drivers) || !all(drivers %in% vars)) {
    cli::cli_abort("{.arg drivers} must name unique fitted predictors.")
  }
  weights <- .cast_attribution_weights(fit, cv, models, weights, method,
    min_score, decay, fallback, min_metric_folds)
  reference <- .cast_reference(fit, vars)
  bounds <- .cast_support_bounds(reference, support_probs)
  if (any(!vapply(bounds, function(b) all(is.finite(b)), logical(1)))) {
    cli::cli_abort("Training reference contains no finite values for a predictor.")
  }
  layer_names <- c("current", "future", "total_change", "support_total", "total_sd",
    unlist(lapply(drivers, function(v) paste0(c("effect_", "retained_",
      "support_effect_", "support_retained_", "effect_sd_", "retained_sd_"), v))),
    "remainder")
  raster <- inherits(current, "SpatRaster")
  if (raster != inherits(future, "SpatRaster")) {
    cli::cli_abort("{.arg current} and {.arg future} must have the same input type.")
  }
  if (raster) {
    check_suggested("terra", "for raster attribution")
    if (!all(vars %in% names(current)) || !all(vars %in% names(future))) {
      cli::cli_abort("Both scenario rasters must contain every fitted predictor.")
    }
    terra::compareGeom(current, future, lyrs = FALSE, stopOnError = TRUE)
    current <- current[[vars]]; future <- future[[vars]]
    mask <- .cast_raster_mask(mask, current)
    if (is.null(filename)) filename <- tempfile("cast_attribution_", fileext = ".tif")
    if (file.exists(filename) && !overwrite) cli::cli_abort("{.arg filename} already exists.")
    writer <- .cast_raster_writer(current, filename, layer_names,
      datatype = "FLT8S", overwrite = overwrite)
    on.exit(.cast_writer_abort(writer), add = TRUE)
    n <- terra::ncell(current)
    rows <- max(1L, as.integer(block_size / terra::ncol(current)))
    starts <- seq.int(1L, terra::nrow(current), by = rows)
  } else {
    if (!is.data.frame(current) || !is.data.frame(future) || !nrow(current) ||
        !all(vars %in% names(current)) || !all(vars %in% names(future))) {
      cli::cli_abort("Scenario tables must be non-empty and contain every fitted predictor.")
    }
    if (all(c("lon", "lat") %in% names(current)) &&
        all(c("lon", "lat") %in% names(future))) {
      future <- .cast_align_future_env(current, future, "attribution")
    } else if (nrow(current) != nrow(future)) {
      cli::cli_abort("Scenario tables must have equal row counts.")
    }
    .cast_check_numeric_predictors(current[, vars, drop = FALSE])
    .cast_check_numeric_predictors(future[, vars, drop = FALSE])
    if (!is.null(mask) || !is.null(filename)) cli::cli_abort("{.arg mask}/{.arg filename} require raster inputs.")
    n <- nrow(current)
    starts <- seq.int(1L, n, by = block_size)
    out <- matrix(NA_real_, n, length(layer_names), dimnames = list(NULL, layer_names))
  }
  # Online moments avoid retaining national effect surfaces in RAM.
  sums <- matrix(0, 1 + 2 * length(drivers), 8,
    dimnames = list(NULL, c("n", "ns", "sum", "abs", "sums", "abss", "sd", "nsd")))
  accumulate <- function(z) {
    cols <- c("total_change", paste0("effect_", drivers), paste0("retained_", drivers))
    for (j in seq_along(cols)) {
      col <- cols[j]; d <- z[, col]
      sc <- if (j == 1) "support_total" else paste0("support_", col)
      dc <- if (j == 1) "total_sd" else sub("^(effect|retained)_", "\\1_sd_", col)
      k <- is.finite(d); s <- k & !is.na(z[, sc]) & z[, sc] == 1
      ds <- z[, dc]; ks <- k & is.finite(ds)
      sums[j, ] <<- sums[j, ] + c(sum(k), sum(s), sum(d[k]), sum(abs(d[k])),
        sum(d[s]), sum(abs(d[s])), sum(ds[ks]), sum(ks))
    }
  }
  for (i in seq_along(starts)) {
    first <- starts[i]
    if (raster) {
      count <- min(rows, terra::nrow(current) - first + 1L)
      a <- .cast_raster_block(current, first, count)
      b <- .cast_raster_block(future, first, count)
    } else {
      idx <- seq.int(first, min(n, first + block_size - 1L))
      a <- current[idx, vars, drop = FALSE]; b <- future[idx, vars, drop = FALSE]
    }
    keep <- .cast_finite_rows(a) & .cast_finite_rows(b)
    if (raster) keep <- keep & .cast_mask_rows(mask, first, count, nrow(a))
    z <- matrix(NA_real_, nrow(a), length(layer_names), dimnames = list(NULL, layer_names))
    if (any(keep)) z[keep, ] <- .cast_attribution_block(fit, a[keep, , drop = FALSE],
      b[keep, , drop = FALSE], drivers, weights, bounds,
      reference = reference, clamp = clamp)
    accumulate(z)
    if (raster) .cast_writer_block(writer, z, first, count) else out[idx, ] <- z
    if (verbose) cli::cli_inform("Scenario attribution: block {i}/{length(starts)}")
  }
  divide <- function(x, n) ifelse(n > 0, x / n, NA_real_)
  summary <- data.frame(driver = c(NA_character_, drivers, drivers),
    contrast = c("total", rep("only", length(drivers)), rep("retained", length(drivers))),
    n = sums[, "n"], n_supported = sums[, "ns"],
    support = divide(sums[, "ns"], sums[, "n"]),
    mean_signed = divide(sums[, "sum"], sums[, "n"]),
    mean_abs = divide(sums[, "abs"], sums[, "n"]),
    mean_signed_supported = divide(sums[, "sums"], sums[, "ns"]),
    mean_abs_supported = divide(sums[, "abss"], sums[, "ns"]),
    mean_model_sd = divide(sums[, "sd"], sums[, "nsd"]), row.names = NULL)
  data <- NULL
  if (raster) {
    .cast_writer_finish(writer)
    filename <- normalizePath(filename, winslash = "/", mustWork = TRUE)
  } else {
    xy <- intersect(c("lon", "lat"), names(current))
    data <- cbind(current[, xy, drop = FALSE], as.data.frame(out))
  }
  structure(list(summary = summary, data = data,
    raster_file = if (raster) filename else NULL, layers = layer_names,
    weights = weights, drivers = drivers, support_probs = support_probs,
    interpretation = "Fitted-model environmental scenario comparisons on the suitability scale"),
    class = c("cast_attribution", "list"))
}

#' @rdname cast_attribution
#' @param x A `cast_attribution` result for printing.
#' @param ... Ignored.
#' @export
print.cast_attribution <- function(x, ...) {
  cli::cli_h2("Environmental attribution of habitat change")
  print(x$summary, row.names = FALSE)
  cli::cli_text("Driver-only and retained-driver contrasts are non-additive scenario comparisons.")
  invisible(x)
}
