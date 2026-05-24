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
