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
#' @param models Character vector. Subset of models to include. Default
#'   `NULL` (all fitted models).
#'
#' @return A `cast_ensemble` object with components:
#' \describe{
#'   \item{predictions}{A `data.frame` with `lon`, `lat`, `hss_ensemble`,
#'     `hss_sd` (cross-model standard deviation over the models contributing
#'     to that cell, `NA` when fewer than two contribute), and
#'     `binary_ensemble` columns. When the fit carries a training
#'     reference, the MESS columns (`mess`, `extrapolating`) from
#'     [cast_predict()] are retained.}
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
#' following the N-SDM nested-modelling framework (Adde et al. 2020). All
#' three components are required: a model whose CV metrics are incomplete
#' scores `NA` and is dropped with a warning, because averaging over
#' whichever components happen to be present would rescale the score and
#' make models incomparable. Models with Score < 0.5 are excluded from the
#' weighted ensemble (a warning is issued when this excludes every model and
#' equal weights are used as a fallback).
#'
#' Models are combined cell by cell. A model that is non-finite at some cells
#' still contributes everywhere else, and the weights are renormalised per
#' cell over the models that produced a value there; only a model that is
#' non-finite everywhere is dropped outright. Cells where no model produced a
#' value are `NA`.
#'
#' The binary threshold is the TSS-maximising cut of the ensemble built from
#' the out-of-fold predictions stored by [cast_cv()], i.e. of the surface that
#' is actually thresholded. When `cv` carries no out-of-fold predictions the
#' function falls back to a weighted mean of the per-model thresholds and
#' warns that this is only an approximation.
#'
#' @references
#' Adde, A., Rey, C., Brun, P., et al. (2020). N-SDM: a high-performance
#' computing pipeline for Nested Species Distribution Modelling.
#' *Ecography*, 43(2), 331-334.
#'
#' @seealso [cast_cv()], [cast_predict()], [cast_project()]
#'
#' @export
cast_ensemble <- function(fit, cv, new_data,
                          method = c("weighted", "best", "equal"),
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
  scores <- .cast_ensemble_scores(cv_sub, mdl_names)

  # ---- Determine weights --------------------------------------------------
  weights <- .cast_ensemble_weights(scores, method)

  # ---- Generate per-model predictions -------------------------------------
  pred_obj <- cast_predict(fit, new_data, models = mdl_names)
  pred_df <- pred_obj$predictions

  # ---- Combine into ensemble HSS + cross-model uncertainty --------------
  # Only models with positive weight AND finite predictions take part, in
  # both the weighted mean (renormalised) and the cross-model SD; the
  # raster path (cast_ensemble_raster) uses the same convention.
  hss_cols <- paste0("HSS_", mdl_names)
  include <- rep(TRUE, length(mdl_names))
  for (i in seq_along(mdl_names)) {
    col <- hss_cols[i]
    if (!col %in% names(pred_df)) {
      cli::cli_warn(
        "No {.field {col}} column in the predictions; {.val {mdl_names[i]}} is excluded from the ensemble."
      )
      include[i] <- FALSE
      next
    }
    vals <- pred_df[[col]]
    n_bad <- sum(!is.finite(vals))
    if (n_bad == length(vals)) {
      cli::cli_warn(
        "{.val {mdl_names[i]}} is non-finite everywhere; excluded from the ensemble."
      )
      include[i] <- FALSE
      next
    }
    if (n_bad > 0L) {
      cli::cli_warn(c(
        "{.val {mdl_names[i]}} is non-finite at {n_bad} cell{?s}.",
        i = "Those cells are averaged over the remaining models; the model still contributes elsewhere."
      ))
    }
    # Zero-weight models neither average in nor count towards the
    # cross-model SD (same convention as cast_ensemble_raster()).
    if (weights[mdl_names[i]] <= 0) include[i] <- FALSE
  }
  if (!any(include)) {
    cli::cli_abort("No model produced finite predictions for the ensemble.")
  }
  w <- weights[include]
  if (sum(w) <= 0) w <- rep(1, sum(include))
  w <- w / sum(w)

  pred_mat <- as.matrix(pred_df[, hss_cols[include], drop = FALSE])
  ensemble_hss <- .ensemble_rowmean(pred_mat, w)
  hss_sd <- .ensemble_rowsd(pred_mat)

  # ---- Binary threshold ---------------------------------------------------
  # Use the same model set that actually contributes to the ensemble.
  threshold <- .ensemble_threshold(cv, mdl_names[include], weights, method)

  # ---- Build output -------------------------------------------------------
  has_coords <- all(c("lon", "lat") %in% names(pred_df))
  out_df <- if (has_coords) {
    data.frame(lon = pred_df$lon, lat = pred_df$lat)
  } else {
    data.frame(site = seq_len(nrow(pred_df)))
  }
  out_df$hss_ensemble <- ensemble_hss
  out_df$hss_sd <- hss_sd
  out_df$binary_ensemble <- as.integer(ensemble_hss >= threshold)
  # Keep the MESS extrapolation flags computed by cast_predict() instead of
  # silently dropping them.
  if (all(c("mess", "extrapolating") %in% names(pred_df))) {
    out_df$mess <- pred_df$mess
    out_df$extrapolating <- pred_df$extrapolating
  }

  new_cast_ensemble(
    predictions  = out_df,
    weights      = weights,
    method       = method,
    threshold    = threshold,
    model_scores = scores
  )
}


#' Row-Wise Weighted Mean Over Model Columns
#'
#' Renormalises the weights per row over the models that produced a finite
#' value there, so one model failing at one cell neither voids the cell nor
#' costs that model its contribution to every other cell.
#'
#' @param pred_mat Numeric matrix, one column per model.
#' @param w Numeric weights, one per column.
#' @return Numeric vector of length `nrow(pred_mat)`; `NA` where no column
#'   carries a finite value.
#' @keywords internal
#' @noRd
.ensemble_rowmean <- function(pred_mat, w) {
  if (!nrow(pred_mat)) return(numeric(0))
  finite <- is.finite(pred_mat)
  wt <- finite * rep(as.numeric(w), each = nrow(pred_mat))
  denom <- rowSums(wt)
  filled <- pred_mat
  filled[!finite] <- 0
  out <- rowSums(filled * wt) / denom
  out[denom <= 0] <- NA_real_
  out
}

#' Row-Wise Cross-Model Standard Deviation
#'
#' @param pred_mat Numeric matrix, one column per model.
#' @return Numeric vector; `NA` where fewer than two models contribute.
#' @keywords internal
#' @noRd
.ensemble_rowsd <- function(pred_mat) {
  if (!nrow(pred_mat)) return(numeric(0))
  if (ncol(pred_mat) < 2L) return(rep(NA_real_, nrow(pred_mat)))
  apply(pred_mat, 1L, function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2L) NA_real_ else stats::sd(z)
  })
}


#' Composite per-model scores from CV metrics (N-SDM score)
#' @keywords internal
#' @noRd
.cast_ensemble_scores <- function(cv_sub, mdl_names) {
  need <- c("auc_mean", "tss_mean", "cbi_mean")
  miss <- setdiff(need, names(cv_sub))
  if (length(miss)) {
    cli::cli_abort(c(
      "{.arg cv} metrics lack the column{?s} {.val {miss}}.",
      i = "The N-SDM score is (2*AUC - 1 + maxTSS + CBI) / 3; no component can be dropped."
    ))
  }
  scores <- vapply(mdl_names, function(m) {
    row <- cv_sub[cv_sub$model == m, , drop = FALSE]
    if (nrow(row) == 0) return(NA_real_)
    parts <- c(2 * row$auc_mean[1] - 1, row$tss_mean[1], row$cbi_mean[1])
    # No na.rm: averaging over whichever components survived would change the
    # divisor and make models with different missing metrics incomparable.
    if (anyNA(parts)) return(NA_real_)
    sum(parts) / 3
  }, numeric(1))
  stats::setNames(scores, mdl_names)
}

#' Ensemble weights from composite scores
#' @keywords internal
#' @noRd
.cast_ensemble_weights <- function(scores, method) {
  mdl_names <- names(scores)
  na_mdl <- mdl_names[is.na(scores)]
  if (length(na_mdl) == length(mdl_names)) {
    cli::cli_abort(c(
      "No model has a complete composite score (AUC, maxTSS and CBI).",
      i = "Check {.code cv$metrics} for missing values."
    ))
  }
  if (length(na_mdl)) {
    cli::cli_warn(
      "Dropping {.val {na_mdl}} from the ensemble: incomplete CV metrics."
    )
  }
  weights <- switch(method,
    weighted = {
      w <- scores
      w[is.na(w) | w < 0.5] <- 0
      total <- sum(w)
      if (total > 0) {
        w / total
      } else {
        cli::cli_warn(
          "All model scores are below 0.5; falling back to equal ensemble weights."
        )
        w <- as.numeric(!is.na(scores))
        w / sum(w)
      }
    },
    best = {
      w <- rep(0, length(mdl_names))
      names(w) <- mdl_names
      best_idx <- which.max(scores)
      if (length(best_idx) > 0) w[best_idx] <- 1
      w
    },
    equal = {
      w <- as.numeric(!is.na(scores))
      w / sum(w)
    }
  )
  stats::setNames(weights, mdl_names)
}


#' Compute Ensemble Binary Threshold
#'
#' Optimises TSS on the ensemble of out-of-fold CV predictions, i.e. on the
#' surface that is actually thresholded. The mean of per-model thresholds is
#' not the optimal cut of the averaged surface, so it is used only as a
#' fallback when `cv` carries no out-of-fold predictions.
#'
#' @keywords internal
#' @noRd
.ensemble_threshold <- function(cv, mdl_names, weights, method) {
  if (!length(mdl_names)) return(0.5)
  w <- weights[mdl_names]
  w[!is.finite(w)] <- 0
  if (identical(method, "best")) {
    keep <- which.max(w)
    if (length(keep)) {
      mdl_names <- mdl_names[keep]
      w <- w[keep]
    }
  }

  oof <- cv$oof
  cols <- paste0("HSS_", mdl_names)
  if (is.data.frame(oof) && "obs" %in% names(oof) &&
      all(cols %in% names(oof)) && sum(w) > 0) {
    ens <- .ensemble_rowmean(as.matrix(oof[, cols, drop = FALSE]), w / sum(w))
    ok <- is.finite(ens) & !is.na(oof$obs)
    if (sum(ok) >= 10L && length(unique(oof$obs[ok])) > 1L) {
      return(find_tss_threshold(ens[ok], oof$obs[ok]))
    }
  }

  thresholds <- cv$thresholds
  avail <- intersect(names(thresholds %||% numeric(0)), mdl_names)
  if (!length(avail)) return(0.5)
  cli::cli_warn(c(
    "{.arg cv} carries no usable out-of-fold predictions for these models.",
    i = "Thresholding on the weighted mean of per-model thresholds, which only approximates the ensemble-surface optimum."
  ))
  wa <- w[avail]
  ta <- thresholds[avail]
  out <- if (sum(wa) > 0) sum(wa * ta) / sum(wa) else mean(ta, na.rm = TRUE)
  if (!is.finite(out)) 0.5 else unname(out)
}


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
#'   `"best"`, `"equal"`. See [cast_ensemble()].
#' @param models Character vector or `NULL`. Models to use. Default all.
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
                                 method = c("weighted", "best", "equal"),
                                 models = NULL,
                                 mask = NULL,
                                 clamp = FALSE,
                                 extrapolation = TRUE,
                                 max_memory_mb = 200,
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

  reference <- fit$scaling$reference
  want_mess <- isTRUE(extrapolation) && !is.null(reference)

  # ---- Compute ensemble weights and threshold ---------------------------------
  # Computed before the skip check so a cached run still reports the
  # configuration its rasters were written with.
  mdl_names <- models %||% names(fit$models)
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

  scores <- .cast_ensemble_scores(cv_sub, mdl_names)
  weights <- .cast_ensemble_weights(scores, method)
  threshold <- .ensemble_threshold(cv, mdl_names, weights, method)

  outputs_exist <- file.exists(hss_path) && file.exists(hss_sd_path) &&
    file.exists(bin_path) && (!want_mess || file.exists(mess_path))
  if (!overwrite && outputs_exist) {
    if (verbose) cli::cli_inform("Ensemble rasters exist; skipping (overwrite = FALSE).")
    return(invisible(list(
      hss_path = hss_path, hss_sd_path = hss_sd_path,
      binary_path = bin_path,
      mess_path = if (want_mess) mess_path else NULL,
      weights = weights, threshold = threshold, n_valid_cells = NA_integer_
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
        ens_hss <- .ensemble_rowmean(sub, w_blk)
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
  terra::writeRaster(hss_out, hss_path, overwrite = TRUE,
    gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
    wopt = list(datatype = "FLT4S"))
  rm(hss_out, hss_vec)
  invisible(gc())

  hss_sd_out <- terra::setValues(template, hss_sd_vec)
  names(hss_sd_out) <- "hss_sd"
  terra::writeRaster(hss_sd_out, hss_sd_path, overwrite = TRUE,
    gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
    wopt = list(datatype = "FLT4S"))
  rm(hss_sd_out, hss_sd_vec)
  invisible(gc())

  if (want_mess) {
    mess_out <- terra::setValues(terra::rast(template), mess_vec)
    names(mess_out) <- "mess"
    terra::writeRaster(mess_out, mess_path, overwrite = TRUE,
      gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
      wopt = list(datatype = "FLT4S"))
    rm(mess_out)
  }
  rm(mess_vec)

  bin_out <- terra::setValues(terra::rast(template), bin_vec)
  names(bin_out) <- "binary_ensemble"
  terra::writeRaster(bin_out, bin_path, overwrite = TRUE,
    gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
    wopt = list(datatype = "INT1U"))
  rm(bin_out, bin_vec)

  if (verbose) {
    msg <- c(
      "v" = "Ensemble rasters saved ({format(n_valid, big.mark = ',')} valid cells):",
      " " = "HSS: {.path {hss_path}}",
      " " = "HSS SD: {.path {hss_sd_path}}",
      " " = "Binary: {.path {bin_path}}"
    )
    if (want_mess) msg <- c(msg, " " = "MESS: {.path {mess_path}}")
    cli::cli_inform(msg)
  }

  invisible(list(
    hss_path     = hss_path,
    hss_sd_path  = hss_sd_path,
    binary_path  = bin_path,
    mess_path    = if (want_mess) mess_path else NULL,
    weights      = weights,
    threshold    = threshold,
    n_valid_cells = n_valid
  ))
}
