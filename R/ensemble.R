#' Performance-Weighted Ensemble Prediction
#'
#' Combines predictions from multiple SDM algorithms into a single
#' ensemble habitat suitability map. Model weights are derived from
#' cross-validation performance scores.
#'
#' @param fit A [cast_fit] object.
#' @param cv A [cast_cv] object providing per-model evaluation metrics.
#' @param new_data A `data.frame` with `lon`, `lat`, and environmental
#'   variables matching the training data.
#' @param method Character. Ensemble strategy:
#'   - `"weighted"` (default): weight = Score / sum(Score), zero out
#'     models with Score < 0.5.
#'   - `"best"`: use the single highest-scoring model.
#'   - `"equal"`: simple average of all models.
#' @param threshold_method Character. Method for binary thresholding.
#'   `"maxTSS"` (default) uses the threshold that maximises TSS from CV.
#' @param models Character vector. Subset of models to include. Default
#'   `NULL` (all fitted models).
#'
#' @return A `cast_ensemble` object with components:
#' \describe{
#'   \item{predictions}{A `data.frame` with `lon`, `lat`, `hss_ensemble`,
#'     and `binary_ensemble` columns.}
#'   \item{weights}{Named numeric vector of per-model weights.}
#'   \item{method}{The ensemble method used.}
#'   \item{threshold}{Binary classification threshold.}
#'   \item{model_scores}{Named numeric vector of per-model composite scores.}
#' }
#'
#' @details
#' The composite score for each model is:
#'
#' \deqn{Score = \frac{1}{3}(2 \times AUC - 1 + maxTSS + CBI)}
#'
#' following the N-SDM (Nfiles Species Distribution Modeling) framework.
#' Models with Score < 0.5 are excluded from the weighted ensemble.
#'
#' @seealso [cast_cv()], [cast_predict()], [cast_project()]
#'
#' @export
cast_ensemble <- function(fit, cv, new_data,
                          method = c("weighted", "best", "equal"),
                          threshold_method = "maxTSS",
                          models = NULL) {
  method <- match.arg(method)

  # ---- Compute per-model composite scores from CV -------------------------
  cv_metrics <- cv$metrics
  mdl_names <- models %||% names(fit$models)
  mdl_names <- intersect(mdl_names, names(fit$models))
  mdl_names <- intersect(mdl_names, cv_metrics$model)
  if (length(mdl_names) == 0) {
    cli::cli_abort("No models found in both {.arg fit} and {.arg cv}.")
  }

  cv_sub <- cv_metrics[cv_metrics$model %in% mdl_names, , drop = FALSE]
  scores <- vapply(mdl_names, function(m) {
    row <- cv_sub[cv_sub$model == m, , drop = FALSE]
    if (nrow(row) == 0) return(NA_real_)
    auc_val <- row$auc_mean[1]
    tss_val <- row$tss_mean[1]
    cbi_val <- if ("cbi_mean" %in% names(row)) row$cbi_mean[1] else 0
    mean(c(2 * auc_val - 1, tss_val, cbi_val), na.rm = TRUE)
  }, numeric(1))
  names(scores) <- mdl_names

  # ---- Determine weights --------------------------------------------------
  weights <- switch(method,
    weighted = {
      w <- scores
      w[is.na(w) | w < 0.5] <- 0
      total <- sum(w)
      if (total > 0) w / total else rep(1 / length(w), length(w))
    },
    best = {
      w <- rep(0, length(mdl_names))
      names(w) <- mdl_names
      best_idx <- which.max(scores)
      if (length(best_idx) > 0) w[best_idx] <- 1
      w
    },
    equal = {
      rep(1 / length(mdl_names), length(mdl_names))
    }
  )
  names(weights) <- mdl_names

  # ---- Generate per-model predictions -------------------------------------
  pred_obj <- cast_predict(fit, new_data, models = mdl_names)
  pred_df <- pred_obj$predictions

  # ---- Combine into ensemble HSS ------------------------------------------
  hss_cols <- paste0("HSS_", mdl_names)
  ensemble_hss <- rep(0, nrow(pred_df))
  for (i in seq_along(mdl_names)) {
    col <- hss_cols[i]
    if (col %in% names(pred_df) && weights[i] > 0) {
      vals <- pred_df[[col]]
      vals[is.na(vals)] <- 0
      ensemble_hss <- ensemble_hss + weights[i] * vals
    }
  }

  # ---- Binary threshold ---------------------------------------------------
  threshold <- .ensemble_threshold(cv, mdl_names, weights, method)

  # ---- Build output -------------------------------------------------------
  has_coords <- all(c("lon", "lat") %in% names(pred_df))
  out_df <- if (has_coords) {
    data.frame(lon = pred_df$lon, lat = pred_df$lat)
  } else {
    data.frame(site = seq_len(nrow(pred_df)))
  }
  out_df$hss_ensemble <- ensemble_hss
  out_df$binary_ensemble <- as.integer(ensemble_hss >= threshold)

  new_cast_ensemble(
    predictions  = out_df,
    weights      = weights,
    method       = method,
    threshold    = threshold,
    model_scores = scores
  )
}


#' Compute Ensemble Binary Threshold
#' @keywords internal
#' @noRd
.ensemble_threshold <- function(cv, mdl_names, weights, method) {
  # Use weighted average of per-model optimal thresholds from CV
  thresholds <- cv$thresholds
  if (is.null(thresholds) || length(thresholds) == 0) return(0.5)

  avail <- intersect(names(thresholds), mdl_names)
  if (length(avail) == 0) return(0.5)

  if (method == "best") {
    best_mdl <- names(which.max(weights))
    if (best_mdl %in% names(thresholds)) return(thresholds[best_mdl])
    return(0.5)
  }

  w <- weights[avail]
  t <- thresholds[avail]
  if (sum(w) > 0) {
    sum(w * t) / sum(w)
  } else {
    mean(t, na.rm = TRUE)
  }
}


#' Raster-Based Ensemble Prediction
#'
#' Generates ensemble habitat suitability (HSS) and binary suitability
#' rasters from a fitted model and cross-validation results. This is the
#' raster-native equivalent of [cast_ensemble()]: it reads a
#' `SpatRaster` stack, runs all models in memory, computes the weighted
#' ensemble, and writes GeoTIFF outputs.
#'
#' @param fit A [cast_fit] object.
#' @param cv A [cast_cv] object providing per-model evaluation metrics.
#' @param raster_stack A `terra::SpatRaster` whose layer names cover the
#'   environmental variables in `fit$env_vars`.
#' @param output_dir Character. Directory for output rasters. Created if
#'   it does not exist.
#' @param method Character. Ensemble strategy: `"weighted"` (default),
#'   `"best"`, `"equal"`. See [cast_ensemble()].
#' @param threshold_method Character. Binary threshold method. Default
#'   `"maxTSS"`.
#' @param models Character vector or `NULL`. Models to use. Default all.
#' @param mask A `terra::SpatRaster` or `NULL`. If provided, prediction
#'   is restricted to cells where mask is non-NA.
#' @param prefix Character. Filename prefix. Default `""`.
#' @param overwrite Logical. Overwrite existing output. Default `FALSE`.
#' @param compression Character. GeoTIFF compression. Default `"LZW"`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{hss_path}{File path to the HSS raster.}
#'   \item{binary_path}{File path to the binary raster.}
#'   \item{weights}{Named numeric vector of per-model weights.}
#'   \item{threshold}{Binary classification threshold.}
#'   \item{n_valid_cells}{Number of cells predicted.}
#' }
#'
#' @seealso [cast_ensemble()], [cast_predict()], [cast_project_raster()]
#'
#' @export
cast_ensemble_raster <- function(fit, cv, raster_stack,
                                 output_dir,
                                 method = c("weighted", "best", "equal"),
                                 threshold_method = "maxTSS",
                                 models = NULL,
                                 mask = NULL,
                                 prefix = "",
                                 overwrite = FALSE,
                                 compression = "LZW",
                                 verbose = TRUE) {
  check_suggested("terra", "for raster prediction")
  method <- match.arg(method)

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
  hss_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "hss_ensemble.tif")
  )
  bin_path <- file.path(
    output_dir,
    paste0(prefix, if (nzchar(prefix)) "_" else "", "binary_ensemble.tif")
  )

  if (!overwrite && file.exists(hss_path) && file.exists(bin_path)) {
    if (verbose) cli::cli_inform("Ensemble rasters exist; skipping (overwrite = FALSE).")
    return(invisible(list(
      hss_path = hss_path, binary_path = bin_path,
      weights = NULL, threshold = NULL, n_valid_cells = NA_integer_
    )))
  }

  # ---- Compute ensemble weights and threshold ---------------------------------
  mdl_names <- models %||% names(fit$models)
  mdl_names <- intersect(mdl_names, names(fit$models))
  if (length(mdl_names) == 0) {
    cli::cli_abort("No matching models found in {.arg fit}.")
  }

  cv_metrics <- cv$metrics
  mdl_names <- intersect(mdl_names, cv_metrics$model)
  cv_sub <- cv_metrics[cv_metrics$model %in% mdl_names, , drop = FALSE]

  scores <- vapply(mdl_names, function(m) {
    row <- cv_sub[cv_sub$model == m, , drop = FALSE]
    if (nrow(row) == 0) return(NA_real_)
    auc_val <- row$auc_mean[1]
    tss_val <- row$tss_mean[1]
    cbi_val <- if ("cbi_mean" %in% names(row)) row$cbi_mean[1] else 0
    mean(c(2 * auc_val - 1, tss_val, cbi_val), na.rm = TRUE)
  }, numeric(1))
  names(scores) <- mdl_names

  weights <- switch(method,
    weighted = {
      w <- scores
      w[is.na(w) | w < 0.5] <- 0
      total <- sum(w)
      if (total > 0) w / total else rep(1 / length(w), length(w))
    },
    best = {
      w <- rep(0, length(mdl_names))
      names(w) <- mdl_names
      best_idx <- which.max(scores)
      if (length(best_idx) > 0) w[best_idx] <- 1
      w
    },
    equal = {
      w <- rep(1 / length(mdl_names), length(mdl_names))
      names(w) <- mdl_names
      w
    }
  )
  names(weights) <- mdl_names

  threshold <- .ensemble_threshold(cv, mdl_names, weights, method)

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
  rows_per_block <- max(10L, as.integer(200e6 / bytes_per_row))
  rows_per_block <- min(rows_per_block, nr)
  n_blocks <- ceiling(nr / rows_per_block)

  if (verbose) {
    cli::cli_inform(c(
      "Grid: {nr} x {nc}  |  {nl} vars  |  {n_blocks} block{?s} of ~{rows_per_block} rows"
    ))
  }

  # Pre-load mask values (single layer, safe for memory)
  mask_vals <- NULL
  if (!is.null(mask)) {
    mask_vals <- tryCatch(
      terra::values(mask, mat = FALSE),
      error = function(e) {
        cli::cli_warn("Mask read failed ({e$message}); predicting full extent")
        NULL
      }
    )
  }

  # Pre-allocate output vectors
  n_cells_total <- as.double(nr) * as.double(nc)
  hss_vec <- rep(NA_real_,    n_cells_total)
  bin_vec <- rep(NA_integer_, n_cells_total)

  n_valid <- 0L

  for (bi in seq_len(n_blocks)) {
    row_start <- (bi - 1L) * rows_per_block + 1L
    row_count <- min(rows_per_block, nr - row_start + 1L)

    y_top <- terra::yFromRow(r, row_start) + res_y / 2
    y_bot <- terra::yFromRow(r, row_start + row_count - 1L) - res_y / 2
    block_ext <- terra::ext(terra::xmin(r), terra::xmax(r), y_bot, y_top)

    block_r <- terra::crop(r, block_ext)
    v <- terra::as.data.frame(block_r, na.rm = FALSE)
    rm(block_r)

    if (!all(env_vars %in% names(v))) names(v)[seq_along(env_vars)] <- env_vars
    v <- v[, env_vars, drop = FALSE]

    cell_start <- (row_start - 1L) * nc + 1L
    cell_end   <- (row_start + row_count - 1L) * nc

    if (!is.null(mask_vals)) {
      m <- mask_vals[cell_start:cell_end]
      v[is.na(m), ] <- NA
      rm(m)
    }

    ok <- stats::complete.cases(v)
    n_ok <- sum(ok)

    if (n_ok > 0) {
      X_raw <- as.data.frame(v[ok, , drop = FALSE])
      for (col in names(X_raw)) X_raw[[col]] <- as.numeric(X_raw[[col]])

      ens_hss <- rep(0, n_ok)
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
        preds[is.na(preds)] <- 0
        ens_hss <- ens_hss + weights[mdl_name] * preds
      }

      valid_idx <- cell_start - 1L + which(ok)
      hss_vec[valid_idx] <- ens_hss
      bin_vec[valid_idx] <- as.integer(ens_hss >= threshold)
      n_valid <- n_valid + n_ok
      rm(X_raw, ens_hss, valid_idx)
    }

    rm(v)
    if (bi %% 5 == 0) invisible(gc())
  }

  rm(mask_vals)

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
  terra::writeRaster(hss_out, hss_path, overwrite = TRUE,
    gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
    wopt = list(datatype = "FLT4S"))
  rm(hss_out, hss_vec)
  invisible(gc())

  bin_out <- terra::setValues(terra::rast(template), bin_vec)
  names(bin_out) <- "binary_ensemble"
  terra::writeRaster(bin_out, bin_path, overwrite = TRUE,
    gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
    wopt = list(datatype = "INT1U"))
  rm(bin_out, bin_vec)

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Ensemble rasters saved ({format(n_valid, big.mark = ',')} valid cells):",
      " " = "HSS: {.path {hss_path}}",
      " " = "Binary: {.path {bin_path}}"
    ))
  }

  invisible(list(
    hss_path     = hss_path,
    binary_path  = bin_path,
    weights      = weights,
    threshold    = threshold,
    n_valid_cells = n_valid
  ))
}
