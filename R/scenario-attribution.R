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
                                    bounds, reference = NULL, clamp = FALSE,
                                    diagnostics = FALSE) {
  context <- .cast_prediction_context(fit, weights, current, future, clamp)
  predict <- context$predict
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
  base <- predict()
  fut <- predict(names(future))
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
      changed_rows <- current[[v]] != future[[v]]
      po <- predict(v, base = base, changed_rows = changed_rows)
      ph <- predict(setdiff(names(future), v), base = fut, changed_rows = changed_rows)
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
    if (diagnostics) {
      # Existing per-model predictions suffice: no extra fitted-model calls.
      oe <- average(effect); re <- average(retained)
      agree <- function(p, mean) {
        value <- rowMeans(.cast_response_sign(p) == .cast_response_sign(mean))
        value[rowSums(!is.finite(p)) > 0L] <- NA_real_
        value
      }
      direction <- .cast_response_direction(oe, re)
      detail <- cbind(agree(effect, oe), agree(retained, re), direction, re - oe)
      colnames(detail) <- paste0(c("effect_agreement_", "retained_agreement_",
        "context_direction_", "context_gap_"), v)
      out <- cbind(out, detail)
    }
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
#'   Use `character()` for only the paired total response, without additional
#'   driver-only and retained-driver predictions.
#' @param shifts Optional named numeric vector of finite raw-unit changes,
#'   added to the corresponding future predictors inside processing blocks.
#'   With `future = current` this computes finite-change responses without
#'   materialising another environmental raster stack. These are explicitly
#'   hypothetical predictor changes, not observed data or causal effects.
#' @param diagnostics Also retain per-driver model sign agreement, background
#'   direction classes and background response gaps. Default `FALSE`. Uses
#'   predictions already calculated, without further fitted-model calls.
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
#' @param overwrite Permit replacement of an existing raster when its receipt
#'   does not match. A verified receipt is reused when `resume = TRUE`.
#' @param verbose Whether to report progress.
#' @param resources Optional global resource budget; `cores` controls RF
#'   prediction threads in this single-process operation. See [cast_compare()].
#' @param resume Reuse a verified result receipt for an existing filename,
#'   including when `overwrite = TRUE`. Set `FALSE` to force recomputation.
#'
#' @return A `cast_attribution` list with `summary`, frozen `weights`,
#'   optional `consistency` counts/fractions on both-background supported cells,
#'   `drivers`, `data` (tables) or `raster_file` (rasters), and `layers`.
#'   `effect_<driver>` is driver-only future minus current; `retained_<driver>`
#'   is full future minus future with that driver held current. `remainder`
#'   is total change minus the sum of driver-only effects: it includes joint
#'   response and changes in any omitted drivers. `*_sd` layers describe
#'   cross-model disagreement, not sampling uncertainty or confidence intervals.
#'   With `diagnostics = TRUE`, `*_agreement_<driver>` is the fraction of
#'   active models whose response sign equals the frozen ensemble's sign;
#'   models count equally for this diagnostic. `context_direction_<driver>`
#'   is -1 (both negative), 0 (both zero), 1 (both positive), 2 (opposite
#'   signs), or 3 (zero in one background). `context_gap_<driver>` is retained
#'   minus driver-only response. These diagnostics do not apply a minimum
#'   ecological effect threshold and do not identify ecological causal effects.
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
                             verbose = TRUE, resources = list(), resume = TRUE,
                             shifts = NULL, diagnostics = FALSE) {
  .cast_check_fit(fit)
  budget <- .cast_resources(resources)
  old_options <- options(castSDM.predict_threads = budget$cores)
  on.exit(options(old_options), add = TRUE)
  method <- match.arg(method)
  decay <- .cast_check_decay(decay)
  .cast_check_scalar(block_size, "block_size", lower = 1, integer = TRUE)
  vars <- fit$env_vars
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics))
    cli::cli_abort("{.arg diagnostics} must be TRUE or FALSE.")
  if (is.null(drivers)) drivers <- vars
  if (!is.character(drivers) || anyNA(drivers) ||
      anyDuplicated(drivers) || !all(drivers %in% vars)) {
    cli::cli_abort("{.arg drivers} must name unique fitted predictors.")
  }
  if (!is.null(shifts) && (!is.numeric(shifts) || !length(shifts) ||
      is.null(names(shifts)) || anyNA(names(shifts)) || anyDuplicated(names(shifts)) ||
      !all(names(shifts) %in% vars) || any(!is.finite(shifts))))
    cli::cli_abort("{.arg shifts} must be a finite numeric vector named by fitted predictors.")
  weights <- .cast_attribution_weights(fit, cv, models, weights, method,
    min_score, decay, fallback, min_metric_folds)
  reference <- .cast_reference(fit, vars)
  bounds <- .cast_support_bounds(reference, support_probs)
  if (any(!vapply(bounds, function(b) all(is.finite(b)), logical(1)))) {
    cli::cli_abort("Training reference contains no finite values for a predictor.")
  }
  layer_names <- c("current", "future", "total_change", "support_total", "total_sd",
    unlist(lapply(drivers, function(v) paste0(c("effect_", "retained_",
      "support_effect_", "support_retained_", "effect_sd_", "retained_sd_",
      if (diagnostics) c("effect_agreement_", "retained_agreement_",
        "context_direction_", "context_gap_")), v))),
    "remainder")
  # The memory budget bounds block buffers; fitted learner objects are resident.
  cells_budget <- max(1, floor(budget$memory_mb * 1e6 /
    (8 * (3 * length(vars) + 2 * length(weights) * (length(drivers) + 1) + length(layer_names)) * 4)))
  block_size <- min(block_size, cells_budget)
  raster <- inherits(current, "SpatRaster")
  if (raster != inherits(future, "SpatRaster")) {
    cli::cli_abort("{.arg current} and {.arg future} must have the same input type.")
  }
  checkpoint_dir <- NULL
  if (raster) {
    check_suggested("terra", "for raster attribution")
    old_temp <- terra::terraOptions(print = FALSE)$tempdir
    terra::terraOptions(tempdir = budget$scratch_dir)
    on.exit(terra::terraOptions(tempdir = old_temp), add = TRUE)
    if (!all(vars %in% names(current)) || !all(vars %in% names(future))) {
      cli::cli_abort("Both scenario rasters must contain every fitted predictor.")
    }
    terra::compareGeom(current, future, lyrs = FALSE, stopOnError = TRUE)
    current <- current[[vars]]; future <- future[[vars]]
    mask <- .cast_raster_mask(mask, current)
    if (is.null(filename)) filename <- tempfile("cast_attribution_", fileext = ".tif")
    lock <- .cast_lock(paste0(filename, ".lock"))
    on.exit(filelock::unlock(lock), add = TRUE)
    signature <- list(schema = "cast_attribution/2", fit = .cast_fit_identity(fit), weights = weights,
      current = .cast_data_identity(current), future = .cast_data_identity(future),
      mask = if (!is.null(mask)) .cast_data_identity(mask), drivers = drivers,
      support_probs = support_probs, clamp = clamp,
      core = .cast_code_identity(.cast_attribution_block))
    if (!is.null(shifts)) signature$shifts <- shifts
    if (diagnostics) {
      signature$diagnostics <- TRUE
      signature$sign_tolerance <- 1e-12
      signature$direction <- .cast_code_identity(.cast_response_direction)
    }
    key <- .cast_digest(signature)
    receipt <- paste0(filename, ".receipt.rds")
    cached <- if (resume) .cast_cache_read(receipt, key) else NULL
    if (!is.null(cached)) {
      cached$raster_file <- normalizePath(filename, winslash = "/")
      return(cached)
    }
    if (file.exists(filename) && !overwrite) cli::cli_abort("{.arg filename} already exists.")
    checkpoint_dir <- file.path(dirname(filename), paste0(".", basename(filename), ".blocks_", substr(key, 1L, 12L)))
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
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
  sums <- matrix(0, 2 + 2 * length(drivers), 8,
    dimnames = list(NULL, c("n", "ns", "sum", "abs", "sums", "abss", "sd", "nsd")))
  consistency <- matrix(0, length(drivers), 7L, dimnames = list(drivers,
    c("n_finite", "n_supported", "negative", "zero", "positive", "reversal", "one_zero")))
  accumulate <- function(z) {
    cols <- c("total_change", if (length(drivers)) paste0("effect_", drivers),
      if (length(drivers)) paste0("retained_", drivers), "remainder")
    for (j in seq_along(cols)) {
      col <- cols[j]; d <- z[, col]
      sc <- if (col %in% c("total_change", "remainder")) "support_total" else
        paste0("support_", col)
      dc <- if (col %in% c("total_change", "remainder")) "total_sd" else
        sub("^(effect|retained)_", "\\1_sd_", col)
      k <- is.finite(d); s <- k & !is.na(z[, sc]) & z[, sc] == 1
      ds <- z[, dc]; ks <- k & is.finite(ds)
      sums[j, ] <<- sums[j, ] + c(sum(k), sum(s), sum(d[k]), sum(abs(d[k])),
        sum(d[s]), sum(abs(d[s])), sum(ds[ks]), sum(ks))
    }
    if (diagnostics) for (v in drivers) {
      direction <- z[, paste0("context_direction_", v)]
      finite <- is.finite(direction)
      ok <- finite & z[, paste0("support_effect_", v)] == 1 &
        z[, paste0("support_retained_", v)] == 1
      ok[is.na(ok)] <- FALSE
      consistency[v, ] <<- consistency[v, ] + c(sum(finite), sum(ok),
        vapply(c(-1, 0, 1, 2, 3), function(value) sum(direction[ok] == value), numeric(1)))
    }
  }
  started <- proc.time()[["elapsed"]]
  if (verbose) cli::cli_inform("Environmental response maps: {length(drivers)} predictors, up to {2L + 2L*length(drivers)} environmental combinations; completed raster batches can resume.")
  for (i in seq_along(starts)) {
    first <- starts[i]
    if (raster) {
      count <- min(rows, terra::nrow(current) - first + 1L)
      active <- .cast_mask_rows(mask, first, count, count * terra::ncol(current))
      if (!any(active)) {
        .cast_writer_block(writer, matrix(NA_real_, length(active), length(layer_names)), first, count)
        next
      }
      checkpoint <- file.path(checkpoint_dir, paste0("block_", first, "_", count, ".rds"))
      checkpoint_key <- .cast_digest(list(key, first, count, layer_names))
      saved <- if (resume) .cast_cache_read(checkpoint, checkpoint_key) else NULL
      if (!is.null(saved)) {
        z <- matrix(NA_real_, length(active), length(layer_names), dimnames = list(NULL, layer_names))
        z[saved$indices, ] <- saved$values
        accumulate(z); .cast_writer_block(writer, z, first, count)
        if (verbose) cli::cli_inform("Environmental response maps: raster batch {i}/{length(starts)} reused ({length(saved$indices)} finite cells).")
        next
      }
      a <- .cast_raster_block(current, first, count)
      b <- .cast_raster_block(future, first, count)
    } else {
      idx <- seq.int(first, min(n, first + block_size - 1L))
      a <- current[idx, vars, drop = FALSE]; b <- future[idx, vars, drop = FALSE]
    }
    for (v in names(shifts)) b[[v]] <- b[[v]] + shifts[[v]]
    keep <- .cast_finite_rows(a) & .cast_finite_rows(b)
    if (raster) keep <- keep & .cast_mask_rows(mask, first, count, nrow(a))
    z <- matrix(NA_real_, nrow(a), length(layer_names), dimnames = list(NULL, layer_names))
    if (any(keep)) z[keep, ] <- .cast_attribution_block(fit, a[keep, , drop = FALSE],
      b[keep, , drop = FALSE], drivers, weights, bounds,
      reference = reference, clamp = clamp, diagnostics = diagnostics)
    if (raster) .cast_cache_write(list(indices = which(keep), values = z[keep, , drop = FALSE]),
      checkpoint, checkpoint_key)
    accumulate(z)
    if (raster) .cast_writer_block(writer, z, first, count) else out[idx, ] <- z
    if (verbose) cli::cli_inform("Environmental response maps: raster batch {i}/{length(starts)}; {sum(keep)} finite cells; elapsed {round(proc.time()[['elapsed']] - started)} s.")
  }
  divide <- function(x, n) ifelse(n > 0, x / n, NA_real_)
  n_driver_rows <- 1L + 2L * length(drivers)
  summary <- data.frame(
    driver = c(NA_character_, drivers, drivers, NA_character_),
    contrast = c("total", rep("only", length(drivers)),
      rep("retained", length(drivers)), "remainder"),
    n = sums[, "n"], n_supported = sums[, "ns"],
    support = divide(sums[, "ns"], sums[, "n"]),
    mean_signed = divide(sums[, "sum"], sums[, "n"]),
    mean_abs = divide(sums[, "abs"], sums[, "n"]),
    mean_signed_supported = divide(sums[, "sums"], sums[, "ns"]),
    mean_abs_supported = divide(sums[, "abss"], sums[, "ns"]),
    mean_model_sd = c(divide(sums[seq_len(n_driver_rows), "sd"],
      sums[seq_len(n_driver_rows), "nsd"]), NA_real_),
    row.names = NULL)
  data <- NULL
  if (raster) {
    .cast_writer_finish(writer)
    filename <- normalizePath(filename, winslash = "/", mustWork = TRUE)
  } else {
    xy <- intersect(c("lon", "lat"), names(current))
    data <- cbind(current[, xy, drop = FALSE], as.data.frame(out))
  }
  evidence <- if (diagnostics && length(drivers)) data.frame(driver = drivers, consistency,
    row.names = NULL) else NULL
  if (!is.null(evidence)) for (name in c("negative", "zero", "positive", "reversal", "one_zero"))
    evidence[[paste0(name, "_fraction")]] <- ifelse(evidence$n_supported > 0,
      evidence[[name]]/evidence$n_supported, NA_real_)
  result <- structure(list(summary = summary, data = data, consistency = evidence,
    raster_file = if (raster) filename else NULL, layers = layer_names,
    weights = weights, drivers = drivers, support_probs = support_probs,
    shifts = shifts, diagnostics = diagnostics, sign_tolerance = if (diagnostics) 1e-12 else NULL,
    interpretation = "Fitted-model environmental scenario comparisons on the suitability scale"),
    class = c("cast_attribution", "list"))
  if (raster) {
    .cast_cache_write(result, receipt, key, filename)
    # Remove only owned numerical checkpoints after successful publication.
    parent <- normalizePath(dirname(filename), winslash = "/", mustWork = TRUE)
    expected <- file.path(parent, paste0(".", basename(filename), ".blocks_", substr(key, 1L, 12L)))
    resolved <- normalizePath(checkpoint_dir, winslash = "/", mustWork = TRUE)
    if (identical(resolved, expected)) {
      unlink(list.files(resolved, pattern = "^block_[0-9]+_[0-9]+[.]rds$", full.names = TRUE))
      if (!length(list.files(resolved, all.files = TRUE, no.. = TRUE))) unlink(resolved, recursive = TRUE)
    }
  }
  result
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
