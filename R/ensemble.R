# Ensemble Core -----------------------------------------------------------------
#
# Pointwise ensemble construction: composite model scoring from nested CV
# metrics, weight construction (score-power / min_score / fallback), the
# per-cell aggregation helpers shared with the raster path, and the
# grouped and permutation-importance entry points.
#
# The raster half lives in ensemble-raster.R.

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
#' averaged over whichever components are finite. The mean uses the
#' available components rather than a fixed divisor (the docstring formula
#' \eqn{\frac{1}{3}} applies only when all three components are present).
#' When the CV metrics
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
  mdl_names <- if (is.null(models)) names(fit$models) else models
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
    ensemble_hss <- .ensemble_rowmedian(pred_mat)
  } else if (identical(method, "committee")) {
    thr <- cv$thresholds[mdl_names[include]]
    if (length(thr) != ncol(pred_mat) || any(!is.finite(thr))) {
      cli::cli_abort("{.arg method = 'committee'} requires finite outer-CV thresholds for every contributing model.")
    }
    votes <- sweep(pred_mat, 2L, thr, FUN = ">=")
    votes[!is.finite(pred_mat)] <- NA
    ensemble_hss <- rowMeans(votes, na.rm = TRUE)
    # A committee of one is not a committee. With a single valid vote the
    # fraction collapses to 0/1, and the 0.5 binarisation then turns one
    # model's opinion into a hard decision. Require two valid votes.
    n_votes <- rowSums(!is.na(votes))
    ensemble_hss[n_votes < 2L] <- NA_real_
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
  .row_sd_finite(pred_mat)
}


#' Row-Wise Median Over Model Columns, Ignoring Non-Finite Values
#'
#' Mirrors the `.ensemble_rowmean()` convention: a model that fails at one
#' cell neither voids the cell nor costs it its contribution elsewhere.
#'
#' @keywords internal
#' @noRd
.ensemble_rowmedian <- function(pred_mat) {
  if (!nrow(pred_mat)) return(numeric(0))
  if (ncol(pred_mat) == 1L) {
    v <- pred_mat[, 1L]
    v[!is.finite(v)] <- NA_real_
    return(v)
  }
  m <- as.matrix(pred_mat)
  storage.mode(m) <- "double"
  # Mask non-finite BEFORE delegating: matrixStats' na.rm drops NA/NaN but
  # still counts +/-Inf as data, whereas the reference implementation this
  # replaces filtered on is.finite(). Without this an Inf would silently
  # shift a cell's median.
  bad <- !is.finite(m)
  if (any(bad)) m[bad] <- NA_real_
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    # C-level; rowMedians uses the same "average the two middle values" rule
    # as stats::median for an even number of finite values.
    out <- matrixStats::rowMedians(m, na.rm = TRUE)
  } else {
    out <- apply(m, 1L, function(z) {
      z <- z[is.finite(z)]
      if (length(z)) stats::median(z) else NA_real_
    })
  }
  out[!is.finite(out)] <- NA_real_
  out
}


#' Row-Wise Sample SD Ignoring Non-Finite Values
#'
#' Vectorised replacement for an `apply(pred_mat, 1, ...)` closure: at raster
#' product scale (~3e7 cells) the closure form spends minutes in R-level
#' interpretation, dwarfing the model predictions themselves.
#'
#' Prefers `matrixStats::rowSds()` (C-level) and otherwise uses the
#' two-pass-free `sqrt((sum(x^2) - n*mean^2)/(n-1))` identity computed from
#' vectorised matrix products.
#'
#' @param pred_mat Numeric matrix, one column per model.
#' @return Numeric vector; `NA` where fewer than two finite values exist.
#' @keywords internal
#' @noRd
.row_sd_finite <- function(pred_mat) {
  m <- as.matrix(pred_mat)
  storage.mode(m) <- "double"
  bad <- !is.finite(m)
  if (any(bad)) m[bad] <- NA_real_
  n <- rowSums(!is.na(m))
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    out <- matrixStats::rowSds(m, na.rm = TRUE)
  } else {
    s1 <- rowSums(m, na.rm = TRUE)
    s2 <- rowSums(m^2, na.rm = TRUE)
    out <- ifelse(n > 1L, sqrt(pmax((s2 - (s1^2) / n) / (n - 1), 0)), NA_real_)
  }
  out[!is.finite(out) | n < 2L] <- NA_real_
  out
}


#' Composite per-model scores from CV metrics
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
        # Plain assignment: switch alternatives evaluate in this function's
        # own frame, so `<<-` would SKIP these locals and store the flags in
        # the namespace/global environment instead (the caller reads them
        # from here, and would silently see FALSE/NULL after a fallback).
        fallback_used <- TRUE
        fallback_reason <- sprintf("No finite model score reached min_score = %g", min_score)
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
    # Mirror the fold-coverage guard in cast_cv(): an OOF surface missing
    # whole folds can still clear the total-count test, and thresholding an
    # ensemble mean computed from a partial surface is not a valid operating
    # point.
    cov_ok <- .cast_oof_fold_coverage(cv, mdl_names)
    if (all(cov_ok$ok)) {
      if (sum(ok) >= 10L && length(unique(oof$obs[ok])) > 1L) {
        return(cast_threshold(ens[ok], oof$obs[ok],
                              method = cv$threshold_method %||% "max_tss"))
      }
    } else {
      cli::cli_warn(c(
        "Out-of-fold coverage is incomplete for {sum(!cov_ok$ok)} of {length(mdl_names)} contributing model{?s}.",
        i = "The ensemble threshold would be derived from a partial surface;",
        i = "falling back to the weighted mean of per-model thresholds."
      ))
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
#'     trees   = c("rf"),
#'     other   = c("maxent", "gam"),
#'     full    = c("rf", "maxent", "gam")
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
  mdl_names <- if (is.null(models)) names(fit$models) else models
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
        error = function(e) {
          # Silence here would let .ensemble_rowmean() silently reallocate
          # this model's weight, so a "4-model ensemble" would actually be
          # scored on 3 models while still being reported as 4.
          cli::cli_warn("Permutation pass: prediction failed for {.val {m}}: {conditionMessage(e)}")
          rep(NA_real_, nrow(X))
        }
      )
    }
    out
  }

  # A model that fails on the baseline would fail on every permutation too,
  # and the importance score would then describe a smaller ensemble than the
  # one reported. Refuse rather than report a silently narrowed average.
  base_mat <- pred_one(X_base)
  base_ok <- stats::setNames(
    vapply(mdl_names, function(m) any(is.finite(base_mat[, m])), logical(1)),
    mdl_names
  )
  dead <- mdl_names[!base_ok]
  if (length(dead)) {
    cli::cli_abort(c(
      "{.field cast_ensemble_importance}: {length(dead)} of {length(mdl_names)} model{?s} failed on the unpermuted baseline ({.val {dead}}).",
      i = "Permutation importance would measure a narrower ensemble than the",
      i = "one reported, so its denominator would be wrong. Fix or drop the",
      i = "failing model{?s} first."
    ))
  }

  ens0 <- .ensemble_rowmean(base_mat, w)
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

  imp <- .cast_with_seed(seed, {
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
        perm_mat <- pred_one(X_p)
        # Keep the denominator honest: score each permutation against the SAME
        # model set the baseline used, so a model that fails only under
        # permutation cannot shrink the ensemble mid-comparison.
        perm_keep <- base_ok & vapply(mdl_names, function(m) any(is.finite(perm_mat[, m])), logical(1))
        if (!all(perm_keep)) {
          perm_keep <- base_ok
        }
        ens_p <- .ensemble_rowmean(perm_mat[, perm_keep, drop = FALSE], w[perm_keep])
        ok <- ok0 & is.finite(ens_p)
        imp[vi, p] <- if (sum(ok) >= 3L) {
          r <- suppressWarnings(corr_fun(ens0[ok], ens_p[ok]))
          if (is.finite(r)) 1 - r else NA_real_
        } else NA_real_
      }
    }
    imp
  })

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
