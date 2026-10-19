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
#'   is restricted to finite non-zero cells. The mask must have one layer
#'   and exactly match the prediction grid. Geometry or read failures abort;
#'   align categorical masks explicitly with nearest-neighbour resampling.
#' @param clamp Logical. Clamp predictors to the training range before
#'   prediction. Default `FALSE`. The MESS layer (see `extrapolation`) is
#'   always computed on the unclamped input, so clamping never hides
#'   extrapolation.
#' @param extrapolation Logical. Compute a per-cell MESS layer (Elith et
#'   al. 2010) block by block and write `<prefix>_mess.tif`. Default
#'   `TRUE` (skipped, with `mess_path = NULL`, if the fit lacks a stored
#'   training reference).
#' @param max_memory_mb Numeric. Approximate per-block memory budget in MB
#'   used to size the row blocks. Default `200`. This is not a whole-process
#'   limit: fitted models, prediction workspaces and GDAL add memory costs.
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
#' models contribute). Cells with non-finite covariates (or masked out) stay `NA` in
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
  .cast_check_fit(fit)
  max_memory_mb <- .cast_check_scalar(max_memory_mb, "max_memory_mb", lower = 0, min_inclusive = FALSE)

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

  raster_stack <- raster_stack[[env_vars]]
  mask <- .cast_raster_mask(mask, raster_stack)

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

  # Only block-sized output buffers are retained; each layer streams to disk.
  r <- raster_stack
  nr <- terra::nrow(r); nc <- terra::ncol(r); nl <- terra::nlyr(r)
  bytes_per_row <- as.double(nc) * (nl + length(mdl_names) + 6L) * 8 * 5
  rows_per_block <- max(1L, floor(min(nr, max_memory_mb * 1e6 / bytes_per_row)))
  n_blocks <- ceiling(nr / rows_per_block)
  if (verbose) cli::cli_inform("Grid: {nr} x {nc} | {n_blocks} blocks of ~{rows_per_block} rows")
  paths <- c(hss_ensemble = hss_path, hss_sd = hss_sd_path, binary_ensemble = bin_path,
    if (want_mess) c(mess = mess_path),
    if (want_aoa) c(aoa = aoa_path, aoa_di = aoa_di_path))
  writers <- list()
  on.exit(for (writer in writers) .cast_writer_abort(writer), add = TRUE)
  for (name in names(paths)) writers[[name]] <- .cast_raster_writer(r, paths[[name]], name,
    datatype = if (name %in% c("binary_ensemble", "aoa")) "INT1U" else "FLT4S",
    compression = compression, overwrite = TRUE)

  n_valid <- 0L
  warned <- stats::setNames(rep(FALSE, length(mdl_names)), mdl_names)
  warned_empty_block <- FALSE

  for (bi in seq_len(n_blocks)) {
    row_start <- (bi - 1L) * rows_per_block + 1L
    row_count <- min(rows_per_block, nr - row_start + 1L)

    v <- .cast_raster_block(r, row_start, row_count)
    if (!identical(names(v), env_vars)) cli::cli_abort("Raster block lost fitted predictor names.")
    ok <- .cast_finite_rows(v) & .cast_mask_rows(mask, row_start, row_count, nrow(v))
    block <- matrix(NA_real_, nrow(v), length(writers), dimnames = list(NULL, names(writers)))
    n_ok <- sum(ok)

    if (n_ok > 0) {
      X_raw <- as.data.frame(v[ok, , drop = FALSE])
      for (col in names(X_raw)) X_raw[[col]] <- as.numeric(X_raw[[col]])

      valid_idx <- which(ok)

      # MESS on the (unclamped) block input; clamping happens afterwards and
      # therefore never masks extrapolation.
      if (want_mess) {
        block[valid_idx, "mess"] <- tryCatch(
          .cast_mess(reference, X_raw),
          error = function(e) rep(NA_real_, n_ok)
        )
      }
      # AOA likewise on the unclamped input (same convention as MESS).
      if (want_aoa) {
        block[valid_idx, "aoa_di"] <- tryCatch({
          new_s <- .cast_aoa_apply(aoa_params, X_raw[, env_vars, drop = FALSE])
          .cast_aoa_di(aoa_params$train, new_s)
        }, error = function(e) rep(NA_real_, n_ok))
        di_blk <- block[valid_idx, "aoa_di"]
        block[valid_idx, "aoa"] <- ifelse(is.finite(di_blk), as.integer(di_blk <= aoa_threshold), NA_integer_)
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
        if (!is.numeric(preds) || length(preds) != n_ok) {
          cli::cli_abort("Model {.val {mdl_name}} returned an invalid prediction vector in block {bi}.")
        }
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
        block[valid_idx, "hss_ensemble"] <- ens_hss
        block[valid_idx, "hss_sd"] <- ens_sd
        block[valid_idx, "binary_ensemble"] <- as.integer(ens_hss >= threshold)
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

    for (name in names(writers)) .cast_writer_block(writers[[name]], block[, name], row_start, row_count)
    rm(v, block)
    if (bi %% 5 == 0) invisible(gc())
  }

  if (n_valid == 0L) {
    cli::cli_abort("No valid (non-NA) cells found in raster stack.")
  }

  for (writer in writers) .cast_writer_finish(writer)

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


