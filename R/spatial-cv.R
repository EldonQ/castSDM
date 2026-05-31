#' Spatial K-Fold Cross-Validation for SDMs
#'
#' Evaluates fitted models using spatially blocked k-fold cross-validation.
#' Unlike random splitting, spatial blocks ensure test folds are geographically
#' separated from training folds, providing honest estimates of model
#' transferability.
#'
#' @param data A `data.frame` with `lon`, `lat`, `presence` columns and
#'   environmental variables.
#' @param screen A `cast_select` object from [cast_select()], or `NULL`.
#' @param dag A [cast_dag] object, or `NULL`.
#' @param k Integer. Number of spatial folds. Default `5`.
#' @param models Character vector. Models to cross-validate: `"rf"`, `"brt"`,
#'   `"maxent"`, `"gam"`. Default `c("rf")`.
#' @param block_method Character. How to create spatial blocks:
#'   \describe{
#'     \item{`"grid"`}{(Default) Divide bounding box into a grid and
#'       round-robin assign cells to folds.}
#'     \item{`"cluster"`}{k-means clustering on lon/lat.}
#'   }
#' @param response Character. Response column. Default `"presence"`.
#' @param rf_ntree Integer. RF trees per fold. Default `300`.
#' @param brt_n_trees Integer. BRT iterations per fold. Default `500`.
#' @param parallel Logical. Run folds in parallel via \pkg{future.apply}.
#'   Default `FALSE`.
#' @param seed Integer or `NULL`. Base random seed.
#' @param verbose Logical. Print fold progress. Default `TRUE`.
#'
#' @return A `cast_cv` object with components:
#' \describe{
#'   \item{`metrics`}{`data.frame` -- per-model mean +/- SD.}
#'   \item{`fold_metrics`}{`data.frame` -- per-fold per-model raw metrics.}
#'   \item{`folds`}{Integer vector -- fold assignment for each row.}
#'   \item{`k`}{Number of folds used.}
#'   \item{`block_method`}{Blocking method used.}
#'   \item{`thresholds`}{Named numeric -- optimal TSS threshold per model.}
#' }
#'
#' @references
#' Roberts, D.R. et al. (2017). Cross-validation strategies for data with
#' temporal, spatial, hierarchical, or phylogenetic structure.
#' *Ecography*, 40(8), 913-929.
#'
#' @seealso [cast_fit()], [cast_evaluate()], [cast_ensemble()]
#'
#' @export
cast_cv <- function(data,
                    screen  = NULL,
                    dag     = NULL,
                    k       = 5L,
                    models  = c("rf"),
                    block_method = c("grid", "cluster"),
                    response     = "presence",
                    rf_ntree     = 300L,
                    brt_n_trees  = 500L,
                    parallel     = FALSE,
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
  env_vars <- if (!is.null(screen)) {
    screen$selected
  } else if (!is.null(dag)) {
    setdiff(dag$nodes, c(response, "lon", "lat"))
  } else {
    get_env_vars(data, response)
  }

  # -- 3. Cross-validate -----------------------------------------------------
  oof_obs <- data[[response]]

  cv_one_fold <- function(fold_i, data, folds, screen, dag, models,
                          response, env_vars, rf_ntree,
                          brt_n_trees, seed) {
    test_idx  <- which(folds == fold_i)
    train_idx <- which(folds != fold_i)

    if (length(test_idx) == 0 || sum(data[[response]][train_idx]) < 5) {
      return(NULL)
    }

    train_fold <- data[train_idx, , drop = FALSE]
    test_fold  <- data[test_idx,  , drop = FALSE]

    fold_fit <- tryCatch(
      cast_fit(
        train_fold,
        screen   = screen,
        dag      = dag,
        models   = models,
        response = response,
        rf_ntree = rf_ntree,
        brt_n_trees = brt_n_trees,
        seed     = if (!is.null(seed)) seed + fold_i else NULL,
        verbose  = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fold_fit)) return(NULL)

    fold_rows <- list()
    oof_updates <- list()
    for (mdl in models) {
      if (!mdl %in% names(fold_fit$models)) next
      mdl_info <- fold_fit$models[[mdl]]

      X_test_raw <- as.data.frame(test_fold[, env_vars, drop = FALSE])
      for (col in names(X_test_raw)) {
        X_test_raw[[col]] <- as.numeric(X_test_raw[[col]])
      }
      X_test_raw[is.na(X_test_raw)] <- 0

      preds <- tryCatch(
        predict_single_model(mdl_info, X_test_raw),
        error = function(e) rep(NA_real_, nrow(test_fold))
      )

      oof_updates[[mdl]] <- list(idx = test_idx, preds = preds)

      m <- evaluate_model_full(preds, test_fold[[response]])
      fold_rows[[mdl]] <- data.frame(
        fold  = fold_i,
        model = mdl,
        auc   = m["auc"],
        tss   = m["tss"],
        cbi   = m["cbi"],
        row.names = NULL
      )
    }
    list(fold_rows = fold_rows, oof_updates = oof_updates)
  }

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    if (verbose) cli::cli_inform("Running {k}-fold CV in parallel...")
    fold_results <- future.apply::future_lapply(
      seq_len(k), cv_one_fold,
      data = data, folds = folds, screen = screen, dag = dag,
      models = models, response = response, env_vars = env_vars,
      rf_ntree = rf_ntree, brt_n_trees = brt_n_trees, seed = seed,
      future.seed = TRUE
    )
  } else {
    fold_results <- vector("list", k)
    for (fold_i in seq_len(k)) {
      if (verbose) {
        test_idx_v  <- which(folds == fold_i)
        train_idx_v <- which(folds != fold_i)
        n_pres_test <- sum(data[[response]][test_idx_v])
        cli::cli_inform(
          "Fold {fold_i}/{k}: train n={length(train_idx_v)}, test n={length(test_idx_v)} (pres={n_pres_test})"
        )
      }
      fold_results[[fold_i]] <- cv_one_fold(
        fold_i, data, folds, screen, dag, models, response,
        env_vars, rf_ntree, brt_n_trees, seed
      )
    }
  }

  # Collect results
  all_fold_rows <- list()
  oof_preds <- lapply(models, function(m) rep(NA_real_, nrow(data)))
  names(oof_preds) <- models

  for (fr in fold_results) {
    if (is.null(fr)) next
    all_fold_rows <- c(all_fold_rows, fr$fold_rows)
    for (mdl in names(fr$oof_updates)) {
      upd <- fr$oof_updates[[mdl]]
      oof_preds[[mdl]][upd$idx] <- upd$preds
    }
  }

  # -- 4. Aggregate ----------------------------------------------------------
  fold_df <- if (length(all_fold_rows) > 0) {
    do.call(rbind, all_fold_rows)
  } else {
    data.frame(fold = integer(), model = character(),
               auc = numeric(), tss = numeric(), cbi = numeric())
  }

  metric_cols <- c("auc", "tss", "cbi")
  agg_rows <- list()
  for (mdl in models) {
    sub <- fold_df[fold_df$model == mdl, metric_cols, drop = FALSE]
    if (nrow(sub) == 0) next
    means <- colMeans(sub, na.rm = TRUE)
    sds   <- apply(sub, 2, stats::sd, na.rm = TRUE)
    agg_rows[[mdl]] <- data.frame(
      model     = mdl,
      auc_mean  = means["auc"],   auc_sd  = sds["auc"],
      tss_mean  = means["tss"],   tss_sd  = sds["tss"],
      cbi_mean  = means["cbi"],   cbi_sd  = sds["cbi"],
      n_folds   = nrow(sub),
      row.names = NULL
    )
  }
  metrics_df <- if (length(agg_rows) > 0) {
    do.call(rbind, agg_rows)
  } else {
    data.frame()
  }
  rownames(metrics_df) <- NULL

  # -- 5. OOF optimal thresholds ---------------------------------------------
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

#' @keywords internal
#' @noRd
make_spatial_folds <- function(lon, lat, k, method = "grid", seed = NULL) {
  n <- length(lon)
  if (!is.null(seed)) set.seed(seed)

  if (method == "cluster") {
    coords <- scale(cbind(lon, lat))
    km <- stats::kmeans(coords, centers = k, nstart = 10, iter.max = 100)
    folds <- as.integer(factor(km$cluster))
    return(folds)
  }

  # method == "grid"
  lon_breaks <- seq(min(lon), max(lon), length.out = k + 1L)
  lat_breaks <- seq(min(lat), max(lat), length.out = k + 1L)
  lon_cell <- findInterval(lon, lon_breaks, rightmost.closed = TRUE)
  lat_cell <- findInterval(lat, lat_breaks, rightmost.closed = TRUE)
  grid_id  <- paste(lon_cell, lat_cell, sep = "_")

  unique_cells <- unique(grid_id)
  cell_counts  <- tabulate(match(grid_id, unique_cells))
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
  folds[is.na(folds)] <- 1L
  folds
}


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
