#' Spatial K-Fold Cross-Validation for SDMs
#'
#' Evaluates fitted models using spatially blocked k-fold cross-validation.
#' Unlike random splitting, spatial blocks ensure test folds are geographically
#' separated from training folds, providing honest estimates of model
#' transferability. The resulting metrics (AUC, TSS, CBI, SEDI, Kappa, PRAUC)
#' feed directly into [cast_predict()] for threshold-calibrated HSS maps.
#'
#' @param data A `data.frame` with `lon`, `lat`, `presence` columns and
#'   environmental variables. Typically the full dataset before any split,
#'   or `split$train` from [cast_prepare()].
#' @param screen A [cast_screen] object, or `NULL`.
#' @param dag A [cast_dag] object, or `NULL`.
#' @param ate A [cast_ate] object, or `NULL`.
#' @param k Integer. Number of spatial folds. Default `5`. Increase to `10`
#'   for large datasets; reduce to `3` for small (<200 presence points).
#' @param models Character vector. Models to cross-validate. Same options as
#'   [cast_fit()]: `"cast"`, `"rf"`, `"brt"`, `"maxent"`. Default `c("rf")`.
#' @param block_method Character. How to create spatial blocks:
#'   \describe{
#'     \item{`"grid"`}{(Default) Divide bounding box into a k-by-k grid, then
#'       assign each point to the nearest grid cell, merging until k non-empty
#'       folds are obtained. No external packages required.}
#'     \item{`"cluster"`}{k-means clustering on lon/lat. Fast but ignores
#'       spatial extent balance.}
#'   }
#' @param response Character. Response column. Default `"presence"`.
#' @param n_epochs Integer. Epochs for CI-MLP per fold. Default `100`.
#' @param n_runs Integer. NN ensemble runs per fold. Default `2`.
#' @param rf_ntree Integer. RF trees per fold. Default `300`.
#' @param brt_n_trees Integer. BRT iterations per fold. Default `500`.
#' @param seed Integer or `NULL`. Base random seed.
#' @param verbose Logical. Print fold progress. Default `TRUE`.
#'
#' @return A `cast_cv` object with components:
#' \describe{
#'   \item{`metrics`}{`data.frame` -- per-model mean +/- SD for all metrics.}
#'   \item{`fold_metrics`}{`data.frame` -- per-fold per-model raw metrics.}
#'   \item{`folds`}{Integer vector -- fold assignment for each row of `data`.}
#'   \item{`k`}{Number of folds used.}
#'   \item{`block_method`}{Blocking method used.}
#'   \item{`thresholds`}{Named numeric -- optimal threshold per model (from
#'     pooled OOF predictions, used by [cast_predict()] for binary maps).}
#' }
#'
#' @details
#' ## Why spatial CV matters for HSS maps
#' Random CV recycles nearby points across train/test splits, causing AUC
#' inflation due to spatial autocorrelation (Roberts et al. 2017 MEE). Spatial
#' CV ensures each test fold lies in a geographic region unseen during training,
#' so metrics approximate true extrapolation ability -- the same challenge faced
#' when predicting the full environmental grid for HSS maps.
#'
#' ## Optimal threshold
#' The function pools out-of-fold (OOF) predictions across all k folds and
#' finds the threshold maximising TSS on the full training set. This threshold
#' is returned in `$thresholds` and used by [cast_predict()] when
#' `threshold_from_cv = TRUE`.
#'
#' @references
#' Roberts, D.R. et al. (2017). Cross-validation strategies for data with
#' temporal, spatial, hierarchical, or phylogenetic structure.
#' *Ecography*, 40(8), 913-929.
#'
#' @seealso [cast_fit()], [cast_evaluate()], [cast_predict()]
#'
#' @export
cast_cv <- function(data,
                    screen  = NULL,
                    dag     = NULL,
                    ate     = NULL,
                    k       = 5L,
                    models  = c("rf"),
                    block_method = c("grid", "cluster"),
                    response     = "presence",
                    n_epochs     = 100L,
                    n_runs       = 2L,
                    rf_ntree     = 300L,
                    brt_n_trees  = 500L,
                    seed         = NULL,
                    verbose      = TRUE) {

  block_method <- match.arg(block_method)
  k <- as.integer(k)
  stopifnot(k >= 2L)

  validate_species_data(data)

  # -- 1. Build spatial folds ------------------------------------------------
  folds <- make_spatial_folds(data$lon, data$lat, k = k,
                               method = block_method, seed = seed)
  fold_sizes <- table(folds)
  if (verbose) {
    cli::cli_inform(c(
      "v" = "Spatial {block_method} blocking: {k} folds",
      "i" = "Fold sizes: {paste(as.integer(fold_sizes), collapse = ' | ')}"
    ))
  }

  # -- 2. Env vars -----------------------------------------------------------
  env_vars <- if (!is.null(dag)) dag$nodes else get_env_vars(data, response)

  # -- 3. Cross-validate -----------------------------------------------------
  all_fold_rows <- list()
  # OOF container: for threshold calibration
  oof_preds  <- lapply(models, function(m) rep(NA_real_, nrow(data)))
  names(oof_preds) <- models
  oof_obs    <- data[[response]]

  for (fold_i in seq_len(k)) {
    test_idx  <- which(folds == fold_i)
    train_idx <- which(folds != fold_i)

    if (length(test_idx) == 0 || sum(data[[response]][train_idx]) < 5) {
      if (verbose) cli::cli_warn("Fold {fold_i}: skipping (too few presences).")
      next
    }

    if (verbose) {
      n_pres_test <- sum(data[[response]][test_idx])
      cli::cli_inform(
        "Fold {fold_i}/{k}: train n={length(train_idx)}, test n={length(test_idx)} (pres={n_pres_test})"
      )
    }

    train_fold <- data[train_idx, , drop = FALSE]
    test_fold  <- data[test_idx,  , drop = FALSE]

    # Fit on this fold
    fold_fit <- tryCatch(
      cast_fit(
        train_fold,
        screen   = screen,
        dag      = dag,
        ate      = ate,
        models   = models,
        response = response,
        n_epochs = n_epochs,
        n_runs   = n_runs,
        rf_ntree = rf_ntree,
        brt_n_trees = brt_n_trees,
        seed     = if (!is.null(seed)) seed + fold_i else NULL,
        verbose  = FALSE
      ),
      error = function(e) {
        cli::cli_warn("Fold {fold_i}: cast_fit failed -- {e$message}")
        NULL
      }
    )
    if (is.null(fold_fit)) next

    # Evaluate each model
    fold_rows <- list()
    for (mdl in models) {
      if (!mdl %in% names(fold_fit$models)) next
      mdl_info <- fold_fit$models[[mdl]]

      X_test_raw <- as.data.frame(test_fold[, env_vars, drop = FALSE])
      X_test_raw[is.na(X_test_raw)] <- 0
      X_test_sc <- as.data.frame(
        scale(X_test_raw,
              center = fold_fit$scaling$means,
              scale  = fold_fit$scaling$sds)
      )
      X_test_sc[is.na(X_test_sc)] <- 0

      preds <- tryCatch(
        predict_single_model(
          mdl_info, X_test_raw, X_test_sc,
          fold_fit$screen, fold_fit$dag, fold_fit$ate
        ),
        error = function(e) rep(NA_real_, nrow(test_fold))
      )

      # Store OOF
      oof_preds[[mdl]][test_idx] <- preds

      m <- evaluate_model_full(preds, test_fold[[response]])
      fold_rows[[mdl]] <- data.frame(
        fold  = fold_i,
        model = mdl,
        auc   = m["auc"],
        tss   = m["tss"],
        cbi   = m["cbi"],
        sedi  = m["sedi"],
        kappa = m["kappa"],
        prauc = m["prauc"],
        row.names = NULL
      )
    }
    all_fold_rows <- c(all_fold_rows, fold_rows)
  }

  # -- 4. Aggregate ----------------------------------------------------------
  fold_df <- if (length(all_fold_rows) > 0) {
    do.call(rbind, all_fold_rows)
  } else {
    data.frame(fold=integer(), model=character(),
               auc=numeric(), tss=numeric(), cbi=numeric(),
               sedi=numeric(), kappa=numeric(), prauc=numeric())
  }

  metric_cols <- c("auc", "tss", "cbi", "sedi", "kappa", "prauc")
  agg_rows <- list()
  for (mdl in models) {
    sub <- fold_df[fold_df$model == mdl, metric_cols, drop = FALSE]
    if (nrow(sub) == 0) next
    means <- colMeans(sub, na.rm = TRUE)
    sds   <- apply(sub, 2, stats::sd, na.rm = TRUE)
    r <- data.frame(
      model       = mdl,
      auc_mean    = means["auc"],   auc_sd    = sds["auc"],
      tss_mean    = means["tss"],   tss_sd    = sds["tss"],
      cbi_mean    = means["cbi"],   cbi_sd    = sds["cbi"],
      sedi_mean   = means["sedi"],  sedi_sd   = sds["sedi"],
      kappa_mean  = means["kappa"], kappa_sd  = sds["kappa"],
      prauc_mean  = means["prauc"], prauc_sd  = sds["prauc"],
      n_folds     = nrow(sub),
      row.names   = NULL
    )
    agg_rows[[mdl]] <- r
  }
  metrics_df <- if (length(agg_rows) > 0) {
    do.call(rbind, agg_rows)
  } else {
    data.frame()
  }
  rownames(metrics_df) <- NULL

  # -- 5. OOF optimal thresholds ----------------------------------------------
  thresholds <- vapply(models, function(mdl) {
    p <- oof_preds[[mdl]]
    o <- oof_obs
    valid <- !is.na(p) & !is.na(o)
    if (sum(valid) < 10 || length(unique(o[valid])) < 2) return(0.5)
    find_tss_threshold(p[valid], o[valid])
  }, numeric(1))

  if (verbose) {
    cli::cli_inform("Spatial CV complete.")
    for (mdl in models) {
      if (mdl %in% metrics_df$model) {
        r <- metrics_df[metrics_df$model == mdl, ]
        cli::cli_inform(
          "  {mdl}: AUC={round(r$auc_mean,3)} TSS={round(r$tss_mean,3)} CBI={round(r$cbi_mean,3)} threshold={round(thresholds[mdl],3)}"
        )
      }
    }
  }

  new_cast_cv(
    metrics      = metrics_df,
    fold_metrics = fold_df,
    folds        = folds,
    k            = k,
    block_method = block_method,
    thresholds   = thresholds
  )
}


# =============================================================================
# Internal: spatial fold creation
# =============================================================================

#' Create Spatial Folds
#' @keywords internal
#' @noRd
make_spatial_folds <- function(lon, lat, k, method = "grid", seed = NULL) {
  n <- length(lon)
  if (!is.null(seed)) set.seed(seed)

  if (method == "cluster") {
    # k-means on scaled coordinates
    coords <- scale(cbind(lon, lat))
    km <- stats::kmeans(coords, centers = k, nstart = 10, iter.max = 100)
    folds <- km$cluster
    # Recode to 1..k sequentially
    folds <- as.integer(factor(folds))
    return(folds)
  }

  # method == "grid": divide into grid, merge small cells
  lon_breaks <- seq(min(lon), max(lon), length.out = k + 1L)
  lat_breaks <- seq(min(lat), max(lat), length.out = k + 1L)
  lon_cell <- findInterval(lon, lon_breaks, rightmost.closed = TRUE)
  lat_cell <- findInterval(lat, lat_breaks, rightmost.closed = TRUE)
  grid_id  <- paste(lon_cell, lat_cell, sep = "_")

  # Assign grid cells to k folds (approx equal size)
  unique_cells <- unique(grid_id)
  cell_counts  <- tabulate(match(grid_id, unique_cells))
  # Sort cells by size descending, round-robin assign to folds
  ord   <- order(cell_counts, decreasing = TRUE)
  cells_sorted <- unique_cells[ord]
  fold_assign  <- integer(length(cells_sorted))
  fold_totals  <- integer(k)
  for (ci in seq_along(cells_sorted)) {
    f <- which.min(fold_totals)
    fold_assign[ci] <- f
    fold_totals[f]  <- fold_totals[f] + cell_counts[ord[ci]]
  }
  cell_to_fold <- stats::setNames(fold_assign, cells_sorted)
  folds <- as.integer(cell_to_fold[grid_id])

  # Ensure no NA
  folds[is.na(folds)] <- 1L
  folds
}


#' Find TSS-Maximising Threshold
#' @keywords internal
#' @noRd
find_tss_threshold <- function(pred, obs) {
  thresholds <- seq(0.01, 0.99, by = 0.01)
  tss_vals <- vapply(thresholds, function(thr) {
    pred_bin <- as.integer(pred >= thr)
    tp <- sum(pred_bin == 1 & obs == 1)
    tn <- sum(pred_bin == 0 & obs == 0)
    fp <- sum(pred_bin == 1 & obs == 0)
    fn <- sum(pred_bin == 0 & obs == 1)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    sens + spec - 1
  }, numeric(1))
  thresholds[which.max(tss_vals)]
}
