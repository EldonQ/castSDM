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
#'     models with Score < `min_score`.
#'   - `"best"`: use the single highest-scoring model.
#'   - `"equal"`: simple average of all models.
#'   - `"median"`: cell-wise median of contributing suitability scores.
#'   - `"committee"`: fraction of contributing models voting suitable at
#'     their outer-CV thresholds; the ensemble threshold is 0.5.
#' @param models Character vector. Subset of models to include. Default
#'   `NULL` (all fitted models).
#' @param min_score Minimum composite score for weighted inclusion. Default
#'   `0.5`.
#' @param decay Numeric >= 0. Exponent applied to the composite scores before
#'   weighting (`w = score^decay`, negative scores clamped at 0), then
#'   renormalised. Only used by `method = "weighted"`: `1` (default)
#'   reproduces linear score weighting, `0` degenerates to equal weights, and
#'   values > 1 sharpen the ensemble towards the best-scoring models. The
#'   raw (unpowered) scores are still reported in `model_scores`, and
#'   `min_score` filters on the raw scores.
#' @param fallback When no model reaches `min_score`: `"best"` retains the
#'   highest-scoring finite model, `"equal"` uses equal weights, or `"error"`
#'   aborts. The choice and reason are stored in the result. Default `"best"`.
#' @param min_metric_folds Minimum finite spatial-CV folds required for a
#'   metric to contribute to the composite score. Default `2`.
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
#'   \item{weights_fallback}{Whether weighted selection needed a fallback.}
#'   \item{fallback_reason}{Reason for the fallback, or `NULL`.}
#'   \item{decay}{The decay exponent used for weighting.}
#' }
#'
#' @details
#' The composite score for each model is the mean of its finite score
#' components among `2 x AUC - 1` (AUC skill), `maxTSS`, and Boyce/CBI:
#'
#' \deqn{Score = mean(2 \times AUC - 1, maxTSS, Boyce)}
#'
#' averaged over whichever components are finite, following the N-SDM
#' nested-modelling framework (Adde et al. 2023), whose reference
#' implementation takes the `na.omit` mean of the available components
#' rather than a fixed divisor (the docstring formula \eqn{\frac{1}{3}}
#' applies only when all three components are present). When the CV metrics
#' carry Boyce it substitutes CBI as the ordination-aware component. A
#' model with fewer than two finite score components scores `NA` and is
#' dropped with a warning, because averaging over
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
#' Adde, A., Rey, P.-L., Brun, P., Külling, N., Fopp, F., Altermatt, F.,
#' Broennimann, O., Lehmann, A., Petitpierre, B., Zimmermann, N. E.,
#' Pellissier, L., Guisan, A. (2023). N-SDM: a high-performance computing
#' pipeline for Nested Species Distribution Modelling.
#' *Ecography*, 2023(6), e06540. \doi{10.1111/ecog.06540}
#'
#' @seealso [cast_cv()], [cast_predict()], [cast_project()]
#'
#' @export
cast_ensemble <- function(fit, cv, new_data,
                          method = c("weighted", "best", "equal", "median", "committee"),
                          models = NULL, min_score = 0.5,
                          decay = 1,
                          fallback = c("best", "equal", "error"),
                          min_metric_folds = 2L) {
  method <- match.arg(method)
  fallback <- match.arg(fallback)
  decay <- .cast_check_decay(decay)

  # ---- Compute per-model composite scores from CV -------------------------
  cv_metrics <- cv$metrics
  mdl_names <- models %||% names(fit$models)
  mdl_names <- intersect(mdl_names, names(fit$models))
  mdl_names <- intersect(mdl_names, cv_metrics$model)
  if (length(mdl_names) == 0) {
    cli::cli_abort("No models found in both {.arg fit} and {.arg cv}.")
  }

  cv_sub <- cv_metrics[cv_metrics$model %in% mdl_names, , drop = FALSE]
  scores <- .cast_ensemble_scores(cv_sub, mdl_names,
                                  min_metric_folds = min_metric_folds)

  # ---- Determine weights --------------------------------------------------
  pw <- .cast_power_scores(scores, decay, method, min_score)
  weights <- .cast_ensemble_weights(pw$scores, method, pw$min_score, fallback)
  weights_fallback <- isTRUE(attr(weights, "fallback_used"))
  fallback_reason <- attr(weights, "fallback_reason")
  score_components <- attr(scores, "components")

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
  if (identical(method, "median")) {
    ensemble_hss <- apply(pred_mat, 1L, function(z) {
      z <- z[is.finite(z)]
      if (length(z)) stats::median(z) else NA_real_
    })
  } else if (identical(method, "committee")) {
    thr <- cv$thresholds[mdl_names[include]]
    if (length(thr) != ncol(pred_mat) || any(!is.finite(thr))) {
      cli::cli_abort("{.arg method = 'committee'} requires finite outer-CV thresholds for every contributing model.")
    }
    votes <- sweep(pred_mat, 2L, thr, FUN = ">=")
    votes[!is.finite(pred_mat)] <- NA
    ensemble_hss <- rowMeans(votes, na.rm = TRUE)
    ensemble_hss[!is.finite(ensemble_hss)] <- NA_real_
  } else {
    ensemble_hss <- .ensemble_rowmean(pred_mat, w)
  }
  hss_sd <- .ensemble_rowsd(pred_mat)

  # ---- Binary threshold ---------------------------------------------------
  # Use the same model set that actually contributes to the ensemble.
  threshold <- if (identical(method, "committee")) 0.5 else
    .ensemble_threshold(cv, mdl_names[include], weights, method)

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
    model_scores = scores,
    weights_fallback = weights_fallback,
    fallback_reason = fallback_reason,
    score_components = score_components,
    decay = decay
  )
}

#' Validate the Decay Exponent
#' @keywords internal
#' @noRd
.cast_check_decay <- function(decay) {
  if (!is.numeric(decay) || length(decay) != 1L || !is.finite(decay) ||
      decay < 0) {
    cli::cli_abort("{.arg decay} must be a single finite non-negative number.")
  }
  decay
}

#' Power-Weight the Composite Scores (decay weighting)
#'
#' For `method = "weighted"` and `decay != 1`: `w = max(score, 0)^decay`,
#' non-finite scores stay non-finite. The `min_score` filter keeps working on
#' the **raw** scores (documented semantics), so models below the threshold
#' are zeroed here and the caller passes `min_score = 0` to the weight
#' builder, which then only excludes negative powered values. Any other
#' method/decay combination returns the raw scores unchanged (the exponent
#' is monotone for `best`, irrelevant for `equal`/`median`/`committee`).
#'
#' @return List with `scores` (possibly powered and pre-filtered) and the
#'   `min_score` to use in [.cast_ensemble_weights()].
#' @keywords internal
#' @noRd
.cast_power_scores <- function(scores, decay, method, min_score = 0.5) {
  if (!identical(method, "weighted") || decay == 1) {
    return(list(scores = scores, min_score = min_score))
  }
  pw <- ifelse(is.finite(scores), pmax(scores, 0), NA_real_)^decay
  pw[is.finite(scores) & scores < min_score] <- 0
  list(scores = pw, min_score = 0)
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
.cast_ensemble_scores <- function(cv_sub, mdl_names, min_metric_folds = 2L) {
  need <- c("auc_mean", "tss_mean")
  miss <- setdiff(need, names(cv_sub))
  if (length(miss)) {
    cli::cli_abort(c(
      "{.arg cv} metrics lack the column{?s} {.val {miss}}.",
      i = "Composite scoring needs AUC and TSS; Boyce contributes only when sufficiently evaluable."
    ))
  }
  min_metric_folds <- as.integer(min_metric_folds)
  if (is.na(min_metric_folds) || min_metric_folds < 1L) {
    cli::cli_abort("{.arg min_metric_folds} must be a positive integer.")
  }
  components <- stats::setNames(vector("list", length(mdl_names)), mdl_names)
  get_metric <- function(row, value_name, count_name = NULL) {
    if (!value_name %in% names(row) || !is.finite(row[[value_name]][1])) return(NA_real_)
    if (!is.null(count_name) && count_name %in% names(row) &&
        (!is.finite(row[[count_name]][1]) || row[[count_name]][1] < min_metric_folds)) {
      return(NA_real_)
    }
    unname(row[[value_name]][1])
  }
  scores <- vapply(mdl_names, function(m) {
    row <- cv_sub[cv_sub$model == m, , drop = FALSE]
    if (nrow(row) == 0) return(NA_real_)
    auc <- get_metric(row, "auc_mean", "auc_n_folds")
    tss <- get_metric(row, "tss_mean", "tss_n_folds")
    if ("boyce_mean" %in% names(row)) {
      boyce <- get_metric(row, "boyce_mean", "boyce_n_folds")
    } else {
      boyce <- get_metric(row, "cbi_mean", "cbi_n_folds")
    }
    parts <- c(AUC_S = if (is.finite(auc)) 2 * auc - 1 else NA_real_,
               TSS = tss, Boyce = boyce)
    parts <- parts[is.finite(parts)]
    components[[m]] <<- names(parts)
    if (length(parts) < 2L) return(NA_real_)
    mean(parts)
  }, numeric(1))
  scores <- stats::setNames(scores, mdl_names)
  attr(scores, "components") <- components
  scores
}

#' Ensemble weights from composite scores
#' @keywords internal
#' @noRd
.cast_ensemble_weights <- function(scores, method, min_score = 0.5,
                                   fallback = c("best", "equal", "error")) {
  fallback <- match.arg(fallback)
  if (!is.numeric(min_score) || length(min_score) != 1L || !is.finite(min_score)) {
    cli::cli_abort("{.arg min_score} must be one finite number.")
  }
  mdl_names <- names(scores)
  finite <- is.finite(scores)
  na_mdl <- mdl_names[!finite]
  if (length(na_mdl)) {
    cli::cli_warn("Composite score unavailable for {.val {na_mdl}}; inspect per-metric fold counts.")
  }
  fallback_used <- FALSE
  fallback_reason <- NULL
  weights <- switch(method,
    weighted = {
      w <- scores
      w[!finite | w < min_score] <- 0
      total <- sum(w)
      if (total > 0) {
        w / total
      } else {
        fallback_used <<- TRUE
        fallback_reason <<- sprintf("No finite model score reached min_score = %g", min_score)
        if (fallback == "error") cli::cli_abort("{fallback_reason}.")
        w[] <- 0
        candidates <- which(finite)
        if (fallback == "best") {
          if (!length(candidates)) cli::cli_abort("No finite composite score is available for fallback = 'best'.")
          w[candidates[which.max(scores[candidates])]] <- 1
          cli::cli_warn(c(
            "{fallback_reason}; using the highest-scoring finite model only.",
            i = "The fallback is recorded in the cast_ensemble object."
          ))
        } else {
          if (!length(candidates)) candidates <- seq_along(scores)
          w[candidates] <- 1 / length(candidates)
          cli::cli_warn(c(
            "{fallback_reason}; using equal weights.",
            i = "The fallback is recorded in the cast_ensemble object."
          ))
        }
        w
      }
    },
    best = {
      if (!any(finite)) cli::cli_abort("No finite composite score is available for {.arg method = 'best'}.")
      w <- rep(0, length(mdl_names)); names(w) <- mdl_names
      w[which.max(ifelse(finite, scores, -Inf))] <- 1
      w
    },
    equal = rep(1 / length(mdl_names), length(mdl_names)),
    median = rep(1 / length(mdl_names), length(mdl_names)),
    committee = rep(1 / length(mdl_names), length(mdl_names))
  )
  weights <- stats::setNames(as.numeric(weights), mdl_names)
  attr(weights, "fallback_used") <- fallback_used
  attr(weights, "fallback_reason") <- fallback_reason
  weights
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
      return(cast_threshold(ens[ok], oof$obs[ok],
                            method = cv$threshold_method %||% "max_tss"))
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
          ens_hss <- apply(sub, 1L, function(z) {
            z <- z[is.finite(z)]
            if (length(z)) stats::median(z) else NA_real_
          })
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

  if (want_aoa) {
    aoa_di_out <- terra::setValues(terra::rast(template), aoa_di_vec)
    names(aoa_di_out) <- "aoa_di"
    terra::writeRaster(aoa_di_out, aoa_di_path, overwrite = TRUE,
      gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
      wopt = list(datatype = "FLT4S"))
    rm(aoa_di_out)

    aoa_out <- terra::setValues(terra::rast(template), aoa_vec)
    names(aoa_out) <- "aoa"
    terra::writeRaster(aoa_out, aoa_path, overwrite = TRUE,
      gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
      wopt = list(datatype = "INT1U"))
    rm(aoa_out)
  }
  rm(aoa_di_vec, aoa_vec)

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


#' Grouped Ensemble Prediction
#'
#' Builds one ensemble per caller-specified group of models: the castSDM
#' analogue of biomod2's `em.by` grouping, made explicit by the caller
#' instead of inferred from model attributes. Each group is combined with
#' [cast_ensemble()] independently, so every group carries its own weights,
#' threshold, cross-model SD and binary map. Groups may overlap (a model can
#' belong to several groups, e.g. a family ensemble plus a full ensemble).
#'
#' @param fit A [cast_fit] object.
#' @param cv A [cast_cv] object providing per-model evaluation metrics.
#' @param new_data A `data.frame` with `lon`, `lat`, and environmental
#'   variables matching the training data.
#' @param groups A **named** list of character vectors of model names; each
#'   name becomes an output column suffix and must be unique and non-empty.
#' @param ... Further arguments passed to [cast_ensemble()] (`method`,
#'   `min_score`, `decay`, `fallback`, `min_metric_folds`).
#'
#' @return A `cast_ensemble_grouped` object with components:
#' \describe{
#'   \item{groups}{The groups as supplied (model names resolved).}
#'   \item{ensembles}{Named list of `cast_ensemble` objects, one per group.}
#'   \item{predictions}{A `data.frame` with the coordinate columns and, per
#'     group, `hss_<group>`, `hss_sd_<group>` and `binary_<group>` columns.}
#' }
#'
#' @seealso [cast_ensemble()], [cast_ensemble_raster()]
#'
#' @examples
#' \dontrun{
#' # One ensemble per algorithm family plus a consensus view
#' cast_ensemble_by(fit, cv, preds,
#'   groups = list(
#'     trees   = c("rf", "brt"),
#'     static  = c("glm", "maxent"),
#'     full    = c("rf", "brt", "glm", "maxent")
#'   ))
#' }
#'
#' @export
cast_ensemble_by <- function(fit, cv, new_data, groups, ...) {
  if (!is.list(groups) || !length(groups) || is.null(names(groups)) ||
      any(!nzchar(names(groups))) || anyDuplicated(names(groups))) {
    cli::cli_abort(
      "{.arg groups} must be a non-empty, uniquely named list of character vectors."
    )
  }
  avail <- intersect(names(fit$models), cv$metrics$model)
  resolved <- lapply(groups, function(g) {
    if (!is.character(g) || !length(g)) {
      cli::cli_abort("Every {.arg groups} element must be a non-empty character vector.")
    }
    hit <- intersect(g, avail)
    if (!length(hit)) {
      cli::cli_abort(
        "Group {.val {g}} has no models in both {.arg fit} and {.arg cv}; available: {.val {avail}}."
      )
    }
    hit
  })

  ensembles <- vector("list", length(resolved))
  names(ensembles) <- names(resolved)
  pieces <- vector("list", length(resolved))
  names(pieces) <- names(resolved)
  for (nm in names(resolved)) {
    ens <- cast_ensemble(fit, cv, new_data, models = resolved[[nm]], ...)
    ensembles[[nm]] <- ens
    coords <- ens$predictions[, intersect(
      c("lon", "lat", "site"), names(ens$predictions)), drop = FALSE]
    pieces[[nm]] <- data.frame(
      coords,
      stats::setNames(
        data.frame(
          ens$predictions$hss_ensemble,
          ens$predictions$hss_sd,
          as.integer(ens$predictions$binary_ensemble)
        ),
        c(paste0("hss_", nm), paste0("hss_sd_", nm), paste0("binary_", nm))
      ),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  # Align pieces row by row: all cast_ensemble calls receive the same
  # new_data, so the coordinate blocks are identical; bind the model columns
  # onto the first piece.
  out <- pieces[[1L]]
  if (length(pieces) > 1L) {
    for (nm in names(pieces)[-1L]) {
      extra <- pieces[[nm]][, setdiff(names(pieces[[nm]]), names(out)),
                            drop = FALSE]
      out <- cbind(out, extra)
    }
  }

  structure(
    list(groups = resolved, ensembles = ensembles, predictions = out),
    class = "cast_ensemble_grouped"
  )
}


#' Ensemble-Layer Permutation Importance
#'
#' Variable importance measured on the **ensemble** output rather than on
#' individual models: each predictor is permuted `n_perm` times, the
#' ensemble prediction is recomputed with the frozen per-model weights, and
#' the importance of the predictor is `1 - cor(original, permuted)` averaged
#' over permutations (the biomod2-style ensemble convention). A predictor
#' the ensemble ignores yields values near 0; permuting an influential
#' predictor decorrelates the surface and yields large values.
#'
#' The per-model weights are computed exactly as in [cast_ensemble()]
#' (same composite scores, `min_score`, `decay` and fallback policy) and then
#' frozen, so the importance reflects the ensemble that would actually be
#' deployed. Per-model predictions are evaluated directly on the imputed
#' input, skipping the MESS/clamp machinery of [cast_predict()] - this is a
#' pure re-scoring loop, not a re-prediction of the full predict stack.
#'
#' @param fit A [cast_fit] object.
#' @param cv A [cast_cv] object providing per-model evaluation metrics.
#' @param new_data A `data.frame` with the environmental variables used in
#'   fitting (coordinates optional). Use a manageable sample - background
#'   points, a raster sub-sample or held-out records - not a full prediction
#'   grid: every predictor permutation re-evaluates every model.
#' @param models Character vector. Subset of models to include. Default
#'   `NULL` (all models present in both `fit` and `cv`).
#' @param method,min_score,decay,fallback,min_metric_folds Ensemble
#'   weighting arguments, passed to the same internal logic as
#'   [cast_ensemble()]. Default `method = "weighted"`, `min_score = 0.5`,
#'   `decay = 1`.
#' @param n_perm Integer >= 1. Permutations per predictor. Default `25`.
#' @param metric Correlation used in `1 - cor`: `"pearson"` (default) or
#'   `"spearman"`.
#' @param seed Random seed for the permutations.
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return A `cast_ensemble_importance` object with components:
#' \describe{
#'   \item{importance}{A `data.frame` with `variable`, `importance_mean`
#'     and `importance_sd` (across permutations), sorted descending by mean.}
#'   \item{n_perm, metric, method, decay}{As supplied/resolved.}
#'   \item{weights}{The frozen per-model weights used.}
#'   \item{model_scores}{The raw composite scores.}
#'   \item{baseline}{The unpermuted ensemble prediction.}
#' }
#'
#' @seealso [cast_ensemble()], [cast_ensemble_by()]
#'
#' @export
cast_ensemble_importance <- function(fit, cv, new_data,
                                     models = NULL, method = "weighted",
                                     min_score = 0.5, decay = 1,
                                     fallback = "best",
                                     min_metric_folds = 2L,
                                     n_perm = 25L,
                                     metric = c("pearson", "spearman"),
                                     seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  fallback <- match.arg(fallback)
  metric <- match.arg(metric)
  decay <- .cast_check_decay(decay)
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) {
    cli::cli_abort("{.arg n_perm} must be a positive integer.")
  }
  if (!is.data.frame(new_data) || nrow(new_data) < 3L) {
    cli::cli_abort("{.arg new_data} needs at least 3 rows to correlate surfaces.")
  }

  # ---- Frozen weights, identical to cast_ensemble() ----------------------
  cv_metrics <- cv$metrics
  mdl_names <- models %||% names(fit$models)
  mdl_names <- intersect(mdl_names, names(fit$models))
  mdl_names <- intersect(mdl_names, cv_metrics$model)
  if (!length(mdl_names)) {
    cli::cli_abort("No models found in both {.arg fit} and {.arg cv}.")
  }
  cv_sub <- cv_metrics[cv_metrics$model %in% mdl_names, , drop = FALSE]
  scores <- .cast_ensemble_scores(cv_sub, mdl_names,
                                  min_metric_folds = min_metric_folds)
  pw <- .cast_power_scores(scores, decay, method, min_score)
  weights <- .cast_ensemble_weights(pw$scores, method, pw$min_score, fallback)
  w <- weights[mdl_names]
  w[!is.finite(w)] <- 0
  if (sum(w) <= 0) w[] <- 1
  w <- w / sum(w)

  # ---- Baseline per-model predictions on the imputed input ---------------
  env_vars <- fit$env_vars
  missing <- setdiff(env_vars, names(new_data))
  if (length(missing)) {
    cli::cli_abort("{.arg new_data} is missing fitted predictor{?s}: {.val {missing}}.")
  }
  X_base <- as.data.frame(new_data[, env_vars, drop = FALSE], check.names = FALSE)
  .cast_check_numeric_predictors(X_base)
  for (col in names(X_base)) X_base[[col]] <- as.numeric(X_base[[col]])
  X_base <- .cast_impute(X_base, fit$scaling$impute)

  pred_one <- function(X) {
    out <- matrix(NA_real_, nrow(X), length(mdl_names),
                  dimnames = list(NULL, mdl_names))
    for (m in mdl_names) {
      out[, m] <- tryCatch(
        predict_single_model(fit$models[[m]], X),
        error = function(e) rep(NA_real_, nrow(X))
      )
    }
    out
  }

  ens0 <- .ensemble_rowmean(pred_one(X_base), w)
  ok0 <- is.finite(ens0)
  if (sum(ok0) < 3L ||
      stats::sd(ens0[ok0], na.rm = TRUE) %||% 0 == 0 ||
      !is.finite(stats::sd(ens0[ok0]))) {
    cli::cli_abort(
      "The baseline ensemble prediction is constant or too sparse; permutation importance is undefined."
    )
  }

  corr_fun <- if (identical(metric, "pearson")) {
    function(a, b) stats::cor(a, b, method = "pearson")
  } else {
    function(a, b) stats::cor(a, b, method = "spearman")
  }

  if (!is.null(seed)) set.seed(seed)
  imp <- matrix(NA_real_, nrow = length(env_vars), ncol = n_perm,
                dimnames = list(env_vars, NULL))
  if (verbose) {
    cli::cli_inform(
      "Ensemble permutation importance: {length(env_vars)} predictor{?s} x {n_perm} permutation{?s}."
    )
  }
  for (vi in seq_along(env_vars)) {
    v <- env_vars[vi]
    vals <- X_base[[v]]
    for (p in seq_len(n_perm)) {
      X_p <- X_base
      X_p[[v]] <- sample(vals, size = length(vals), replace = FALSE)
      ens_p <- .ensemble_rowmean(pred_one(X_p), w)
      ok <- ok0 & is.finite(ens_p)
      imp[vi, p] <- if (sum(ok) >= 3L) {
        r <- suppressWarnings(corr_fun(ens0[ok], ens_p[ok]))
        if (is.finite(r)) 1 - r else NA_real_
      } else NA_real_
    }
  }

  importance <- data.frame(
    variable = env_vars,
    importance_mean = vapply(seq_along(env_vars), function(i) {
      z <- imp[i, ][is.finite(imp[i, ])]
      if (length(z)) mean(z) else NA_real_
    }, numeric(1)),
    importance_sd = vapply(seq_along(env_vars), function(i) {
      z <- imp[i, ][is.finite(imp[i, ])]
      if (length(z) >= 2L) stats::sd(z) else NA_real_
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  importance <- importance[order(
    -xtfrm(importance$importance_mean), importance$variable), ,
    drop = FALSE]
  rownames(importance) <- NULL

  structure(
    list(
      importance = importance,
      n_perm = n_perm,
      metric = metric,
      method = method,
      decay = decay,
      weights = weights,
      model_scores = scores,
      baseline = ens0
    ),
    class = "cast_ensemble_importance"
  )
}
