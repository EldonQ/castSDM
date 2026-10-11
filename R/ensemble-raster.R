# Ensemble Raster Engine ------------------------------------------------------
#
# The block-wise raster half of the ensemble layer: it reads the covariate
# stack in blocks, predicts every contributing model per block, combines
# them with the same weighting convention as cast_ensemble() (see
# ensemble-core.R), and writes the output layers atomically.
#
# Kept separate from ensemble-core.R because it is the only module that
# performs low-level terra I/O; nothing here is reused by the pointwise
# API, and the two halves change for different reasons.

#' Raster-Based Ensemble Prediction
#'
#' Generates ensemble habitat suitability (HSS) and binary suitability
#' rasters from a fitted model and cross-validation results. This is the
#' raster-native equivalent of [cast_ensemble()]: it reads a
#' `SpatRaster` stack, runs all models in memory, computes the weighted
#' ensemble, and writes GeoTIFF outputs. Extrapolation control is built
#' in: an optional per-cell MESS layer flags out-of-envelope cells, and
#' `clamp` can cap predictors at their training range.
#'
#' @param fit A [cast_fit] object.
#' @param cv A [cast_cv] object providing per-model evaluation metrics.
#' @param raster_stack A `terra::SpatRaster` whose layer names cover the
#'   environmental variables in `fit$env_vars`.
#' @param output_dir Character. Directory for output rasters. Created if
#'   it does not exist.
#' @param method Character. Ensemble strategy: `"weighted"` (default),
#'   `"best"`, `"equal"`, `"median"`, or `"committee"`. See [cast_ensemble()].
#' @param models Character vector or `NULL`. Models to use. Default all.
#' @param min_score Minimum composite score for weighted inclusion. Default
#'   `0.5`.
#' @param fallback Weighted-score fallback: `"best"`, `"equal"`, or `"error"`.
#'   Default `"best"`.
#' @param min_metric_folds Minimum finite spatial-CV folds for a metric to
#'   contribute to the composite score. Default `2`.
#' @param decay Numeric >= 0. Score-power exponent for the `weighted` method,
#'   as in [cast_ensemble()]. Default `1`.
#' @param aoa_cv Optional `cast_cv` object from [cast_cv()] with `aoa = TRUE`. When
#'   supplied, each valid cell is scored against the **area of applicability**
#'   (Meyer & Pebesma 2021): the DI of the cell's standardized predictors to
#'   the nearest training row is compared with the calibrated threshold, and
#'   `<prefix>_aoa.tif` (flag, 1 = inside the AOA) and `<prefix>_aoa_di.tif`
#'   (raw dissimilarity index) are written. The DI is computed on the
#'   unclamped block input, so clamping never softens the flag. Requires the
#'   \pkg{FNN} package.
#' @param mask A `terra::SpatRaster` or `NULL`. If provided, prediction
#'   is restricted to cells where mask is non-NA. Only the first layer is
#'   used; a mask whose geometry cannot be matched to `raster_stack` is
#'   ignored with a warning.
#' @param clamp Logical. Clamp predictors to the training range before
#'   prediction. Default `FALSE`. The MESS layer (see `extrapolation`) is
#'   always computed on the unclamped input, so clamping never hides
#'   extrapolation.
#' @param extrapolation Logical. Compute a per-cell MESS layer (Elith et
#'   al. 2010) block by block and write `<prefix>_mess.tif`. Default
#'   `TRUE` (skipped, with `mess_path = NULL`, if the fit lacks a stored
#'   training reference).
#' @param max_memory_mb Numeric. Approximate per-block memory budget in MB
#'   used to size the row blocks. Default `200`.
#' @param prefix Character. Filename prefix. Default `""`.
#' @param overwrite Logical. Overwrite existing output. Default `FALSE`.
#' @param compression Character. GeoTIFF compression. Default `"LZW"`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{hss_path}{File path to the HSS raster.}
#'   \item{hss_sd_path}{File path to the cross-model uncertainty (SD) raster.}
#'   \item{binary_path}{File path to the binary raster.}
#'   \item{mess_path}{File path to the MESS raster, or `NULL` when
#'     `extrapolation = FALSE` or no training reference is available.}
#'   \item{aoa_path, aoa_di_path}{File paths to the AOA flag and DI rasters,
#'     or `NULL` when `aoa_cv` is not supplied.}
#'   \item{weights}{Named numeric vector of per-model weights.}
#'   \item{threshold}{Binary classification threshold.}
#'   \item{n_valid_cells}{Number of cells predicted.}
#' }
#'
#' @details
#' Models are combined cell by cell: a model that is non-finite at some cells
#' of a block still contributes to the rest of that block, and the positive
#' weights are renormalised per cell over the models that produced a value
#' there - the same convention as [cast_ensemble()]. The cross-model SD layer
#' likewise uses only the cell-level contributors (`NA` where fewer than two
#' models contribute). Cells with `NA` covariates (or masked out) stay `NA` in
#' every output layer. If no model produces a finite prediction for a block,
#' those cells are `NA` and a warning is issued once. `n_valid_cells` counts
#' the cells that received an ensemble value.
#'
#' @seealso [cast_ensemble()], [cast_predict()], [cast_project_raster()]
#'
#' @export
cast_ensemble_raster <- function(fit, cv, raster_stack,
                                 output_dir,
                                 method = c("weighted", "best", "equal", "median", "committee"),
                                 models = NULL,
                                 min_score = 0.5,
                                 decay = 1,
                                 fallback = c("best", "equal", "error"),
                                 min_metric_folds = 2L,
                                 mask = NULL,
                                 clamp = FALSE,
                                 extrapolation = TRUE,
                                 aoa_cv = NULL,
                                 max_memory_mb = 200,
                                 prefix = "",
                                 overwrite = FALSE,
                                 compression = "LZW",
                                 verbose = TRUE) {
  check_suggested("terra", "for raster prediction")
  method <- match.arg(method)
  fallback <- match.arg(fallback)
  decay <- .cast_check_decay(decay)

  if (!inherits(raster_stack, "SpatRaster")) {
    if (is.character(raster_stack)) {
      raster_stack <- terra::rast(raster_stack)
    } else {
      cli::cli_abort("{.arg raster_stack} must be a SpatRaster or file path.")
    }
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  env_vars <- fit$env_vars
  missing_layers <- setdiff(env_vars, names(raster_stack))
  if (length(missing_layers) > 0) {
    cli::cli_abort("Raster missing required layer{?s}: {.val {missing_layers}}.")
  }

  # ---- Output paths -----------------------------------------------------------
  # `prefix` reaches the filesystem through six paste0() calls below, so a
  # value like "../../tmp/evil" would write outside `output_dir`. Sanitise it
  # once; an empty prefix is legitimate and stays empty.
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
    cli::cli_abort("{.arg prefix} must be a single character string.")
  }
  if (nzchar(prefix)) {
    safe_prefix <- .cast_safe_name(prefix, arg = "prefix")
    if (!identical(safe_prefix, prefix)) {
      cli::cli_warn(c(
        "{.arg prefix} was changed to a safe filename stem:",
        "{.val {prefix}} -> {.val {safe_prefix}}."
      ))
    }
    prefix <- safe_prefix
  }

  hss_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "hss_ensemble.tif")
  )
  hss_sd_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "hss_sd.tif")
  )
  bin_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "binary_ensemble.tif")
  )
  mess_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "mess.tif")
  )
  aoa_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "aoa.tif")
  )
  aoa_di_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "aoa_di.tif")
  )

  reference <- fit$scaling$reference
  want_mess <- isTRUE(extrapolation) && !is.null(reference)

  # AOA setup: validate the calibration and project the training reference
  # into the standardized DI space once, reused block by block.
  want_aoa <- FALSE
  aoa_threshold <- NA_real_
  aoa_params <- NULL
  if (!is.null(aoa_cv)) {
    if (!isTRUE(aoa_cv$aoa$enabled) ||
        !is.finite(aoa_cv$aoa$threshold %||% NA_real_)) {
      cli::cli_abort(c(
        "{.arg aoa_cv} must come from {.code cast_cv(aoa = TRUE)} with a finite threshold.",
        "i" = "The supplied object carries no usable AOA calibration."
      ))
    }
    if (is.null(reference)) {
      cli::cli_abort("The fit carries no stored training reference; AOA cannot be computed.")
    }
    check_suggested("FNN", "for the area of applicability (aoa_cv)")
    ref_m <- as.data.frame(reference, check.names = FALSE)
    ref_m <- ref_m[, intersect(fit$env_vars, names(ref_m)), drop = FALSE]
    aoa_params <- .cast_aoa_params(ref_m)
    aoa_threshold <- aoa_cv$aoa$threshold
    want_aoa <- TRUE
  }

  # ---- Compute ensemble weights and threshold ---------------------------------
  # Computed before the skip check so a cached run still reports the
  # configuration its rasters were written with.
  mdl_names <- if (is.null(models)) names(fit$models) else models
  mdl_names <- intersect(mdl_names, names(fit$models))
  if (length(mdl_names) == 0) {
    cli::cli_abort("No matching models found in {.arg fit}.")
  }

  cv_metrics <- cv$metrics
  mdl_names <- intersect(mdl_names, cv_metrics$model)
  if (length(mdl_names) == 0) {
    cli::cli_abort("No models found in both {.arg fit} and {.arg cv}.")
  }
  cv_sub <- cv_metrics[cv_metrics$model %in% mdl_names, , drop = FALSE]

  scores <- .cast_ensemble_scores(cv_sub, mdl_names,
                                  min_metric_folds = min_metric_folds)
  pw <- .cast_power_scores(scores, decay, method, min_score)
  weights <- .cast_ensemble_weights(pw$scores, method, pw$min_score, fallback)
  weights_fallback <- isTRUE(attr(weights, "fallback_used"))
  fallback_reason <- attr(weights, "fallback_reason")
  score_components <- attr(scores, "components")
  threshold <- if (identical(method, "committee")) 0.5 else
    .ensemble_threshold(cv, mdl_names, weights, method)

  outputs_exist <- file.exists(hss_path) && file.exists(hss_sd_path) &&
    file.exists(bin_path) && (!want_mess || file.exists(mess_path)) &&
    (!want_aoa || (file.exists(aoa_path) && file.exists(aoa_di_path)))
  if (!overwrite && outputs_exist) {
    if (verbose) cli::cli_inform("Ensemble rasters exist; skipping (overwrite = FALSE).")
    return(invisible(list(
      hss_path = hss_path, hss_sd_path = hss_sd_path,
       binary_path = bin_path,
       mess_path = if (want_mess) mess_path else NULL,
       aoa_path = if (want_aoa) aoa_path else NULL,
       aoa_di_path = if (want_aoa) aoa_di_path else NULL,
       weights = weights, threshold = threshold, n_valid_cells = NA_integer_,
       model_scores = scores, weights_fallback = weights_fallback,
       fallback_reason = fallback_reason, score_components = score_components,
       decay = decay
     )))
  }

  if (verbose) {
    w_str <- paste0(mdl_names, "=", round(weights, 3), collapse = ", ")
    cli::cli_inform(c(
      "Ensemble config: {.val {method}} method",
      " " = "Weights: {w_str}",
      " " = "Threshold: {round(threshold, 4)}"
    ))
  }

  # ---- Block-based raster prediction (memory-safe) ----------------------------
  # Uses terra::crop() for block reading instead of readStart/readValues,
  # which avoids both the full-grid segfault and the readStart NULL-pointer
  # issue on some terra installations.
  r <- raster_stack

  nr    <- terra::nrow(r)
  nc    <- terra::ncol(r)
  nl    <- terra::nlyr(r)
  res_y <- terra::yres(r)

  bytes_per_row <- as.double(nc) * nl * 8 * 5
  rows_per_block <- max(10L, as.integer(max_memory_mb * 1e6 / bytes_per_row))
  rows_per_block <- min(rows_per_block, nr)
  n_blocks <- ceiling(nr / rows_per_block)

  if (verbose) {
    cli::cli_inform(c(
      "Grid: {nr} x {nc}  |  {nl} vars  |  {n_blocks} block{?s} of ~{rows_per_block} rows"
    ))
  }

  # Validate the mask geometry once; values are read per block below so a
  # national-scale mask never has to be held in memory as a whole.
  mask_ok <- FALSE
  if (!is.null(mask)) {
    geom_ok <- tryCatch(
      terra::compareGeom(mask, r, ext = TRUE, rowcol = TRUE, res = TRUE,
                         crs = TRUE, stopOnError = FALSE),
      error = function(e) FALSE
    )
    if (!isTRUE(geom_ok)) {
      # A mask in a different CRS (or grid) is reprojected onto the stack
      # geometry before giving up and predicting the full extent.
      mask_proj <- tryCatch(terra::project(mask, r), error = function(e) NULL)
      if (!is.null(mask_proj)) {
        geom_ok <- tryCatch(
          terra::compareGeom(mask_proj, r, ext = TRUE, rowcol = TRUE,
                             res = TRUE, crs = TRUE, stopOnError = FALSE),
          error = function(e) FALSE
        )
        if (isTRUE(geom_ok)) mask <- mask_proj
      }
    }
    if (!isTRUE(geom_ok)) {
      cli::cli_warn(
        "Mask geometry does not match the raster stack; predicting the full extent."
      )
    } else {
      if (terra::nlyr(mask) > 1L) {
        cli::cli_warn(
          "{.arg mask} has {terra::nlyr(mask)} layers; using the first one."
        )
        mask <- mask[[1L]]
      }
      mask_ok <- TRUE
    }
  }

  # Pre-allocate output vectors
  n_cells_total <- as.double(nr) * as.double(nc)
  hss_vec  <- rep(NA_real_,    n_cells_total)
  hss_sd_vec <- rep(NA_real_,  n_cells_total)
  bin_vec  <- rep(NA_integer_, n_cells_total)
  mess_vec <- if (want_mess) rep(NA_real_, n_cells_total) else NULL
  aoa_vec <- if (want_aoa) rep(NA_integer_, n_cells_total) else NULL
  aoa_di_vec <- if (want_aoa) rep(NA_real_, n_cells_total) else NULL

  n_valid <- 0L
  warned <- stats::setNames(rep(FALSE, length(mdl_names)), mdl_names)
  warned_empty_block <- FALSE
  warned_mask <- FALSE

  for (bi in seq_len(n_blocks)) {
    row_start <- (bi - 1L) * rows_per_block + 1L
    row_count <- min(rows_per_block, nr - row_start + 1L)

    y_top <- terra::yFromRow(r, row_start) + res_y / 2
    y_bot <- terra::yFromRow(r, row_start + row_count - 1L) - res_y / 2
    block_ext <- terra::ext(terra::xmin(r), terra::xmax(r), y_bot, y_top)

    block_r <- terra::crop(r, block_ext)
    v <- terra::as.data.frame(block_r, na.rm = FALSE)
    rm(block_r)

    if (!all(env_vars %in% names(v))) {
      cli::cli_abort(c(
        "Block values lost expected layer name{?s}: {.val {setdiff(env_vars, names(v))}}.",
        "i" = "Refusing to match predictors by position; check that raster layer names are stable."
      ))
    }
    v <- v[, env_vars, drop = FALSE]

    cell_start <- (row_start - 1L) * nc + 1L

    if (mask_ok) {
      m <- tryCatch(
        terra::values(mask, mat = FALSE, row = row_start, nrows = row_count),
        error = function(e) NULL
      )
      if (is.null(m) || length(m) != nrow(v)) {
        if (!warned_mask) {
          cli::cli_warn(
            "Mask read failed for at least one block; those cells are predicted unmasked."
          )
          warned_mask <- TRUE
        }
      } else {
        v[is.na(m), ] <- NA
      }
      rm(m)
    }

    ok <- stats::complete.cases(v)
    n_ok <- sum(ok)

    if (n_ok > 0) {
      X_raw <- as.data.frame(v[ok, , drop = FALSE])
      for (col in names(X_raw)) X_raw[[col]] <- as.numeric(X_raw[[col]])

      valid_idx <- cell_start - 1L + which(ok)

      # MESS on the (unclamped) block input; clamping happens afterwards and
      # therefore never masks extrapolation.
      if (want_mess) {
        mess_vec[valid_idx] <- tryCatch(
          .cast_mess(reference, X_raw),
          error = function(e) rep(NA_real_, n_ok)
        )
      }
      # AOA likewise on the unclamped input (same convention as MESS).
      if (want_aoa) {
        aoa_di_vec[valid_idx] <- tryCatch({
          new_s <- .cast_aoa_apply(aoa_params, X_raw[, env_vars, drop = FALSE])
          .cast_aoa_di(aoa_params$train, new_s)
        }, error = function(e) rep(NA_real_, n_ok))
        di_blk <- aoa_di_vec[valid_idx]
        aoa_vec[valid_idx] <- as.integer(!is.na(di_blk) & di_blk <= aoa_threshold)
      }
      if (isTRUE(clamp) && !is.null(reference)) {
        X_raw <- .cast_clamp(X_raw, reference)
      }

      pred_mat <- matrix(NA_real_, nrow = n_ok, ncol = length(mdl_names),
                         dimnames = list(NULL, mdl_names))
      contrib <- stats::setNames(rep(FALSE, length(mdl_names)), mdl_names)
      for (mdl_name in mdl_names) {
        if (weights[mdl_name] == 0) next
        mdl_info <- fit$models[[mdl_name]]
        preds <- tryCatch(
          predict_single_model(mdl_info, X_raw),
          error = function(e) {
            cli::cli_warn("Prediction failed for {.val {mdl_name}} in block {bi}: {e$message}")
            rep(NA_real_, n_ok)
          }
        )
        finite <- is.finite(preds)
        if (!any(finite)) {
          if (!warned[mdl_name]) {
            cli::cli_warn("{.val {mdl_name}} is non-finite across at least one whole block; the other models cover those cells.")
            warned[mdl_name] <- TRUE
          }
          next
        }
        if (!all(finite) && !warned[mdl_name]) {
          cli::cli_warn(c(
            "{.val {mdl_name}} is non-finite at some cells.",
            i = "Those cells are averaged over the remaining models; the model still contributes elsewhere."
          ))
          warned[mdl_name] <- TRUE
        }
        preds[!finite] <- NA_real_
        pred_mat[, mdl_name] <- preds
        contrib[mdl_name] <- TRUE
      }

      if (any(contrib)) {
        # Weights renormalise per cell over the models that produced a value
        # there (same convention as cast_ensemble()); the cross-model SD uses
        # exactly those cell-level contributors.
        w_blk <- weights[contrib]
        w_blk <- w_blk / sum(w_blk)
        sub <- pred_mat[, names(w_blk), drop = FALSE]
        if (identical(method, "median")) {
          ens_hss <- .ensemble_rowmedian(sub)
        } else if (identical(method, "committee")) {
          thr <- cv$thresholds[names(w_blk)]
          if (length(thr) != ncol(sub) || any(!is.finite(thr))) {
            cli::cli_abort("{.arg method = 'committee'} requires finite outer-CV thresholds for contributing models.")
          }
          votes <- sweep(sub, 2L, thr, FUN = ">=")
          votes[!is.finite(sub)] <- NA
          ens_hss <- rowMeans(votes, na.rm = TRUE)
          ens_hss[!is.finite(ens_hss)] <- NA_real_
        } else {
          ens_hss <- .ensemble_rowmean(sub, w_blk)
        }
        ens_sd <- .ensemble_rowsd(sub)
        hss_vec[valid_idx] <- ens_hss
        hss_sd_vec[valid_idx] <- ens_sd
        bin_vec[valid_idx] <- as.integer(ens_hss >= threshold)
        n_valid <- n_valid + sum(is.finite(ens_hss))
        rm(ens_hss, ens_sd, w_blk, sub)
      } else if (!warned_empty_block) {
        cli::cli_warn(
          "No model produced finite predictions for at least one block; those cells are NA."
        )
        warned_empty_block <- TRUE
      }
      rm(X_raw, pred_mat, contrib, valid_idx)
    }

    rm(v)
    if (bi %% 5 == 0) invisible(gc())
  }

  if (n_valid == 0L) {
    cli::cli_abort("No valid (non-NA) cells found in raster stack.")
  }

  # Write output rasters
  template <- terra::rast(
    nrows = nr, ncols = nc,
    xmin = terra::xmin(r), xmax = terra::xmax(r),
    ymin = terra::ymin(r), ymax = terra::ymax(r),
    crs = terra::crs(r)
  )

  hss_out <- terra::setValues(template, hss_vec)
  names(hss_out) <- "hss_ensemble"
  .cast_write_raster_atomic(hss_out, hss_path, compression = compression)
  rm(hss_out, hss_vec)
  invisible(gc())

  hss_sd_out <- terra::setValues(template, hss_sd_vec)
  names(hss_sd_out) <- "hss_sd"
  .cast_write_raster_atomic(hss_sd_out, hss_sd_path, compression = compression)
  rm(hss_sd_out, hss_sd_vec)
  invisible(gc())

  if (want_mess) {
    mess_out <- terra::setValues(terra::rast(template), mess_vec)
    names(mess_out) <- "mess"
    .cast_write_raster_atomic(mess_out, mess_path, compression = compression)
    rm(mess_out)
  }
  rm(mess_vec)

  if (want_aoa) {
    aoa_di_out <- terra::setValues(terra::rast(template), aoa_di_vec)
    names(aoa_di_out) <- "aoa_di"
    .cast_write_raster_atomic(aoa_di_out, aoa_di_path, compression = compression)
    rm(aoa_di_out)

    aoa_out <- terra::setValues(terra::rast(template), aoa_vec)
    names(aoa_out) <- "aoa"
    .cast_write_raster_atomic(aoa_out, aoa_path, compression = compression,
                              datatype = "INT1U")
    rm(aoa_out)
  }
  rm(aoa_di_vec, aoa_vec)

  bin_out <- terra::setValues(terra::rast(template), bin_vec)
  names(bin_out) <- "binary_ensemble"
  .cast_write_raster_atomic(bin_out, bin_path, compression = compression,
                            datatype = "INT1U")
  rm(bin_out, bin_vec)

  if (verbose) {
    msg <- c(
      "v" = "Ensemble rasters saved ({format(n_valid, big.mark = ',')} valid cells):",
      " " = "HSS: {.path {hss_path}}",
      " " = "HSS SD: {.path {hss_sd_path}}",
      " " = "Binary: {.path {bin_path}}"
    )
    if (want_mess) msg <- c(msg, " " = "MESS: {.path {mess_path}}")
    if (want_aoa) msg <- c(msg,
      " " = "AOA: {.path {aoa_path}}",
      " " = "AOA DI: {.path {aoa_di_path}}"
    )
    cli::cli_inform(msg)
  }

  invisible(list(
    hss_path     = hss_path,
    hss_sd_path  = hss_sd_path,
    binary_path  = bin_path,
    mess_path    = if (want_mess) mess_path else NULL,
    aoa_path     = if (want_aoa) aoa_path else NULL,
    aoa_di_path  = if (want_aoa) aoa_di_path else NULL,
    weights      = weights,
    threshold    = threshold,
    model_scores = scores,
    weights_fallback = weights_fallback,
    fallback_reason = fallback_reason,
    score_components = score_components,
    decay = decay,
    n_valid_cells = n_valid
  ))
}


