#' Nested Spatial Cross-Validation for SDMs
#'
#' Variable selection is re-fitted inside every outer training fold. This
#' prevents the held-out fold from influencing predictor choice.
#'
#' @param data Data frame with `lon`, `lat`, binary response, and predictors.
#' @param screen Optional final-data screen. It is used only when
#'   `select_method = NULL`; it is never re-used as a nested screen.
#' @param select_method Selection method passed to [cast_select()]. Default
#'   `"dml"`. Set to `NULL` only to evaluate a fixed supplied screen.
#' @param select_args Named list of additional [cast_select()] arguments.
#' @param k Number of outer spatial folds.
#' @param models Models passed to [cast_fit()].
#' @param block_method `"grid"` or `"cluster"`.
#' @param response Binary response column.
#' @param rf_ntree RF trees per fold.
#' @param brt_n_trees BRT iterations per fold.
#' @param parallel Run folds with `future.apply`.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_cv` object including fold-level selections.
#' @export
cast_cv <- function(data,
                    screen = NULL,
                    select_method = "dml",
                    select_args = list(),
                    k = 5L,
                    models = c("rf"),
                    block_method = c("grid", "cluster"),
                    response = "presence",
                    rf_ntree = 300L,
                    brt_n_trees = 500L,
                    parallel = FALSE,
                    seed = NULL,
                    verbose = TRUE) {
  block_method <- match.arg(block_method)
  k <- as.integer(k)
  if (k < 2L) cli::cli_abort("{.arg k} must be at least 2.")
  validate_species_data(data, required_cols = c("lon", "lat", response))
  folds <- make_spatial_folds(data$lon, data$lat, k, block_method, seed)

  cv_one <- function(fold_i) {
    test_idx <- which(folds == fold_i)
    train_idx <- which(folds != fold_i)
    train <- data[train_idx, , drop = FALSE]
    test <- data[test_idx, , drop = FALSE]
    if (length(unique(train[[response]])) < 2L ||
        length(unique(test[[response]])) < 2L) return(NULL)

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
      x_test <- as.data.frame(test[, fit$env_vars, drop = FALSE])
      for (nm in names(x_test)) x_test[[nm]] <- as.numeric(x_test[[nm]])
      x_test[is.na(x_test)] <- 0
      pred <- tryCatch(
        predict_single_model(info, x_test),
        error = function(e) rep(NA_real_, nrow(test))
      )
      met <- evaluate_model_full(pred, test[[response]])
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
  for (i in seq_along(results)) {
    res <- results[[i]]
    if (is.null(res)) next
    row_list <- c(row_list, res$rows)
    selections[[i]] <- res$selected
    screens[[i]] <- res$screen
    for (mdl in names(res$updates)) {
      upd <- res$updates[[mdl]]
      oof[[mdl]][upd$idx] <- upd$pred
    }
  }
  fold_df <- if (length(row_list)) do.call(rbind, row_list) else data.frame()
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
    selections = selections, screens = screens
  )
}

#' @keywords internal
#' @noRd
make_spatial_folds <- function(lon, lat, k, method = "grid", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (method == "cluster") {
    return(as.integer(factor(stats::kmeans(
      scale(cbind(lon, lat)), centers = k, nstart = 10L
    )$cluster)))
  }
  side <- ceiling(sqrt(k * 2))
  xb <- cut(lon, breaks = unique(stats::quantile(
    lon, seq(0, 1, length.out = side + 1L), na.rm = TRUE
  )), include.lowest = TRUE, labels = FALSE)
  yb <- cut(lat, breaks = unique(stats::quantile(
    lat, seq(0, 1, length.out = side + 1L), na.rm = TRUE
  )), include.lowest = TRUE, labels = FALSE)
  cell <- interaction(xb, yb, drop = TRUE)
  counts <- sort(table(cell), decreasing = TRUE)
  totals <- integer(k)
  assignment <- integer(length(counts))
  names(assignment) <- names(counts)
  for (nm in names(counts)) {
    f <- which.min(totals)
    assignment[nm] <- f
    totals[f] <- totals[f] + counts[nm]
  }
  as.integer(assignment[as.character(cell)])
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
