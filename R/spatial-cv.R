#' Nested Spatial Cross-Validation for SDMs
#'
#' Variable selection is re-fitted inside every outer training fold. This
#' prevents the held-out fold from influencing predictor choice.
#'
#' @param data Data frame with `lon`, `lat`, binary response, and predictors.
#' @param screen Optional final-data screen. It is used only when
#'   `select_method = NULL`; it is never re-used as a nested screen. Because
#'   such a screen was fitted on the full data set, evaluating it here leaks
#'   selection into the CV metrics (a warning is issued); treat those metrics
#'   as optimistic.
#' @param select_method Selection method passed to [cast_select()]. Default
#'   `"cpi"`. Set to `NULL` only to evaluate a fixed supplied screen.
#' @param select_args Named list of additional [cast_select()] arguments.
#' @param k Number of outer spatial folds.
#' @param models Models passed to [cast_fit()].
#' @param block_method `"grid"` (default; grid cells grouped into spatially
#'   contiguous folds via k-means on cell centroids), `"grid_random"`
#'   (legacy behaviour: count-balanced greedy packing that ignores cell
#'   position and produces interleaved, non-contiguous folds), or `"cluster"`
#'   (k-means on point coordinates).
#' @param buffer Non-negative number. When `> 0`, training rows closer than
#'   `buffer` (in coordinate units, e.g. degrees for lon/lat) to any test row
#'   are dropped for that fold, thinning spatial-autocorrelation leakage at
#'   the fold boundary. Default `0` (no exclusion band).
#' @param response Binary response column.
#' @param rf_ntree RF trees per fold.
#' @param brt_n_trees BRT iterations per fold.
#' @param parallel Run folds with `future.apply`; requires the user to set a
#'   [future::plan()] first (a warning is issued when none is set).
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_cv` object including fold-level selections.
#' @export
cast_cv <- function(data,
                    screen = NULL,
                    select_method = "cpi",
                    select_args = list(),
                    k = 5L,
                    models = c("rf"),
                    block_method = c("grid", "grid_random", "cluster"),
                    buffer = 0,
                    response = "presence",
                    rf_ntree = 300L,
                    brt_n_trees = 500L,
                    parallel = FALSE,
                    seed = NULL,
                    verbose = TRUE) {
  block_method <- match.arg(block_method)
  check_suggested("pROC", "for fold evaluation metrics")
  k <- as.integer(k)
  if (k < 2L) cli::cli_abort("{.arg k} must be at least 2.")
  if (!is.numeric(buffer) || length(buffer) != 1L || buffer < 0) {
    cli::cli_abort("{.arg buffer} must be a single non-negative number.")
  }
  validate_species_data(data, required_cols = c("lon", "lat", response),
                        response = response)
  if (is.null(select_method) && !is.null(screen)) {
    cli::cli_warn(c(
      "{.code select_method = NULL} reuses the supplied {.arg screen} in every fold.",
      i = "That screen was selected on the full data set, so selection leaks into these CV metrics; treat them as optimistic."
    ))
  }
  folds <- make_spatial_folds(data$lon, data$lat, k, block_method, seed)

  cv_one <- function(fold_i) {
    test_idx <- which(folds == fold_i)
    train_idx <- .cast_buffer_train_idx(data$lon, data$lat, test_idx, buffer)
    train <- data[train_idx, , drop = FALSE]
    test <- data[test_idx, , drop = FALSE]
    if (length(unique(train[[response]])) < 2L ||
        length(unique(test[[response]])) < 2L) {
      return(list(skipped_single_class = TRUE))
    }

    fold_screen <- if (!is.null(select_method)) {
      args <- utils::modifyList(
        list(
          data = train, response = response, method = select_method,
          seed = if (is.null(seed)) NULL else seed + fold_i,
          verbose = FALSE
        ),
        select_args
      )
      tryCatch(do.call(cast_select, args), error = function(e) {
        warning(sprintf("Selection failed in fold %d: %s", fold_i, e$message))
        NULL
      })
    } else {
      screen
    }
    if (is.null(fold_screen) || !length(fold_screen$selected)) return(NULL)

    fit <- tryCatch(
      cast_fit(
        train, screen = fold_screen, models = models, response = response,
        rf_ntree = rf_ntree, brt_n_trees = brt_n_trees,
        seed = if (is.null(seed)) NULL else seed + 100L + fold_i,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    rows <- list()
    updates <- list()
    for (mdl in models) {
      info <- fit$models[[mdl]]
      if (is.null(info) || is.null(info$model)) next
      x_test <- as.data.frame(test[, fit$env_vars, drop = FALSE], check.names = FALSE)
      .cast_check_numeric_predictors(x_test, arg = "data")
      for (nm in names(x_test)) x_test[[nm]] <- as.numeric(x_test[[nm]])
      x_test <- .cast_impute(x_test, fit$scaling$impute)
      pred <- tryCatch(
        predict_single_model(info, x_test),
        error = function(e) rep(NA_real_, nrow(test))
      )
      met <- tryCatch(
        evaluate_model_full(pred, test[[response]]),
        error = function(e) c(auc = NA_real_, tss = NA_real_, cbi = NA_real_)
      )
      rows[[mdl]] <- data.frame(
        fold = fold_i, model = mdl, auc = met["auc"], tss = met["tss"],
        cbi = met["cbi"], n_selected = length(fold_screen$selected),
        stringsAsFactors = FALSE
      )
      updates[[mdl]] <- list(idx = test_idx, pred = pred)
    }
    list(rows = rows, updates = updates, selected = fold_screen$selected,
         screen = fold_screen)
  }

  if (verbose) {
    cli::cli_inform(
      "Nested spatial CV: {k} folds; selector={select_method %||% 'fixed'}."
    )
  }
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    if (requireNamespace("future", quietly = TRUE) &&
        inherits(future::plan(), "sequential")) {
      cli::cli_warn(c(
        "{.code parallel = TRUE} but no {.pkg future} plan is set; folds will run sequentially.",
        i = "Set one first, e.g. {.code future::plan(future::multisession)}."
      ))
    }
    results <- future.apply::future_lapply(
      seq_len(k), function(i) cv_one(i), future.seed = TRUE
    )
  } else {
    results <- lapply(seq_len(k), cv_one)
  }

  row_list <- list()
  oof <- stats::setNames(lapply(models, function(x) rep(NA_real_, nrow(data))), models)
  selections <- vector("list", k)
  screens <- vector("list", k)
  skipped_single <- integer(0)
  for (i in seq_along(results)) {
    res <- results[[i]]
    if (is.null(res)) next
    if (isTRUE(res$skipped_single_class)) {
      skipped_single <- c(skipped_single, i)
      next
    }
    row_list <- c(row_list, res$rows)
    selections[[i]] <- res$selected
    screens[[i]] <- res$screen
    for (mdl in names(res$updates)) {
      upd <- res$updates[[mdl]]
      oof[[mdl]][upd$idx] <- upd$pred
    }
  }

  # Fold-level selection frequency: the fraction of outer folds in which
  # each predictor was retained. The basis of the consensus selector
  # (cast_consensus()) and of the spatial-stability diagnostic. The
  # denominator is ALL k folds: an empty (or failed) fold contributes zero
  # rather than silently shrinking the denominator and inflating frequency.
  all_vars <- unique(unlist(selections))
  selection_freq <- data.frame(variable = character(0), freq = numeric(0),
                               stringsAsFactors = FALSE)
  if (length(all_vars)) {
    freq <- vapply(all_vars, function(v) {
      mean(vapply(selections, function(s) v %in% (s %||% character(0)),
                  logical(1)))
    }, numeric(1))
    selection_freq <- data.frame(
      variable = all_vars, freq = unname(freq), stringsAsFactors = FALSE)
    selection_freq <- selection_freq[order(-selection_freq$freq), , drop = FALSE]
    rownames(selection_freq) <- NULL
  }
  fold_df <- if (length(row_list)) do.call(rbind, row_list) else data.frame()
  if (length(skipped_single)) {
    cli::cli_warn(
      "Skipped {length(skipped_single)}/{k} fold{?s} ({skipped_single}): a single response class in the train or test split."
    )
  }
  if (!nrow(fold_df)) cli::cli_abort("All spatial CV folds failed.")

  agg <- lapply(models, function(mdl) {
    z <- fold_df[fold_df$model == mdl, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    data.frame(
      model = mdl,
      auc_mean = mean(z$auc, na.rm = TRUE), auc_sd = stats::sd(z$auc, na.rm = TRUE),
      tss_mean = mean(z$tss, na.rm = TRUE), tss_sd = stats::sd(z$tss, na.rm = TRUE),
      cbi_mean = mean(z$cbi, na.rm = TRUE), cbi_sd = stats::sd(z$cbi, na.rm = TRUE),
      n_folds = nrow(z), n_selected_mean = mean(z$n_selected, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, Filter(Negate(is.null), agg))
  thresholds <- vapply(models, function(mdl) {
    pred <- oof[[mdl]]
    ok <- is.finite(pred)
    if (sum(ok) < 10L || length(unique(data[[response]][ok])) < 2L) return(0.5)
    find_tss_threshold(pred[ok], data[[response]][ok])
  }, numeric(1))

  new_cast_cv(
    metrics = metrics, fold_metrics = fold_df, folds = folds, k = k,
    block_method = block_method, thresholds = thresholds,
    selections = selections, screens = screens,
    selection_freq = selection_freq
  )
}

#' Assign Spatial Folds
#'
#' @param lon,lat Numeric coordinate vectors.
#' @param k Number of folds.
#' @param method `"grid"` (default): bin coordinates into grid cells and group
#'   cells into `k` spatially contiguous folds via k-means on cell centroids,
#'   so every fold is a connected region. `"grid_random"`: legacy behaviour
#'   (greedy count-balanced packing that ignores cell position; kept for
#'   backwards comparability). `"cluster"`: k-means on the (scaled) point
#'   coordinates.
#' @param seed Random seed.
#' @return Integer fold ids in `1:k` (degenerate inputs collapse to `1L`).
#' @keywords internal
#' @noRd
make_spatial_folds <- function(lon, lat, k,
                               method = c("grid", "grid_random", "cluster"),
                               seed = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  n <- length(lon)
  k_use <- min(k, n)
  if (k_use < 2L) return(rep(1L, n))
  if (method == "cluster") {
    return(as.integer(factor(stats::kmeans(
      scale(cbind(lon, lat)), centers = k_use, nstart = 10L
    )$cluster)))
  }
  # grid blocking; degenerate coordinates (few distinct values) collapse to
  # a coarser grid rather than erroring
  side <- ceiling(sqrt(k * 2))
  xb <- .cast_bin(lon, side)
  yb <- .cast_bin(lat, side)
  cell <- interaction(xb, yb, drop = TRUE)
  counts <- table(cell)
  if (method == "grid_random") {
    # Legacy: greedy count-balanced packing, ignores cell position entirely.
    counts <- sort(counts, decreasing = TRUE)
    totals <- integer(k)
    assignment <- integer(length(counts))
    names(assignment) <- names(counts)
    for (nm in names(counts)) {
      f <- which.min(totals)
      assignment[nm] <- f
      totals[f] <- totals[f] + counts[nm]
    }
    return(as.integer(assignment[as.character(cell)]))
  }
  # "grid": spatially contiguous folds. Cells are the atoms; their centroids
  # are clustered into k spatial groups and each group becomes one fold.
  cell_lon <- as.numeric(tapply(lon, cell, mean))
  cell_lat <- as.numeric(tapply(lat, cell, mean))
  n_cells <- length(cell_lon)
  if (n_cells <= k) {
    return(as.integer(factor(cell, levels = names(counts))))
  }
  km <- stats::kmeans(cbind(cell_lon, cell_lat), centers = k, nstart = 10L)
  fold_of_cell <- stats::setNames(km$cluster, names(counts))
  as.integer(fold_of_cell[as.character(cell)])
}


#' Training Index with a Buffer Exclusion Band
#'
#' Returns the training-row indices for one CV fold: every row not in
#' `test_idx` whose distance to the nearest test row is at least `buffer`
#' (in coordinate units). Rows inside the band are dropped, thinning
#' spatial-autocorrelation leakage at the fold boundary.
#'
#' @param lon,lat Numeric coordinate vectors.
#' @param test_idx Integer indices of the held-out fold.
#' @param buffer Non-negative exclusion distance; `0` keeps all non-test rows.
#' @keywords internal
#' @noRd
.cast_buffer_train_idx <- function(lon, lat, test_idx, buffer = 0) {
  all_idx <- seq_along(lon)
  if (buffer <= 0) return(setdiff(all_idx, test_idx))
  d <- as.matrix(stats::dist(cbind(lon, lat)))
  near <- apply(d[test_idx, , drop = FALSE] < buffer, 2L, any)
  all_idx[!near & !(all_idx %in% test_idx)]
}

#' Bin a coordinate into `side` quantile bins (degenerate-safe)
#' @keywords internal
#' @noRd
.cast_bin <- function(x, side) {
  x <- as.numeric(x)
  u <- sort(unique(x[!is.na(x)]))
  if (length(u) < 2L) return(rep(1L, length(x)))
  breaks <- unique(stats::quantile(x, seq(0, 1, length.out = side + 1L),
                                    na.rm = TRUE, names = FALSE))
  if (length(breaks) < 2L) breaks <- c(u[1], u[length(u)])
  b <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  b[is.na(b)] <- 1L
  b
}

#' @keywords internal
#' @noRd
find_tss_threshold <- function(pred, obs) {
  thresholds <- seq(0.01, 0.99, by = 0.01)
  values <- vapply(thresholds, function(thr) {
    cls <- pred >= thr
    sens <- sum(cls & obs == 1) / max(1, sum(obs == 1))
    spec <- sum(!cls & obs == 0) / max(1, sum(obs == 0))
    sens + spec - 1
  }, numeric(1))
  thresholds[which.max(values)]
}
