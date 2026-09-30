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
#'   `"two_stage"`. Set to `NULL` only to evaluate a fixed supplied screen.
#' @param select_args Named list of additional [cast_select()] arguments.
#'   Cannot override `data`, `response`, or `method`.
#' @param k Number of outer spatial folds.
#' @param models Models passed to [cast_fit()].
#' @param block_method `"grid"` (default; grid cells grouped into spatially
#'   contiguous folds via k-means on cell centroids), `"grid_random"`
#'   (legacy behaviour: count-balanced greedy packing that ignores cell
#'   position and produces interleaved, non-contiguous folds), `"cluster"`
#'   (k-means on point coordinates), or `"env"` (k-means on the scaled
#'   predictor matrix, so every fold is contiguous in **environmental**
#'   space rather than geographic space; Roberts et al. 2017-style
#'   environmental blocking).
#' @param buffer Non-negative number. When `> 0`, training rows closer than
#'   `buffer` to any test row are dropped for that fold, thinning
#'   spatial-autocorrelation leakage at the fold boundary. Units follow
#'   `buffer_unit`. Default `0` (no exclusion band).
#' @param buffer_unit Unit of `buffer`: `"deg"` (default) planar degrees as
#'   before, or `"km"` great-circle kilometres (Haversine).
#'   `"km"` requires decimal-degree lon/lat coordinates and is refused on
#'   projected inputs.
#' @param response Binary response column.
#' @param rf_ntree RF trees per fold. Default `500`.
#' @param brt_n_trees BRT trees per fold. Default `2000`.
#' @param brt_shrinkage BRT learning rate per fold. Default `0.005`.
#' @param tune Logical. Run the per-engine hyperparameter grid search
#'   ([cast_fit()], argument `tune`) inside every outer training fold, so the
#'   held-out folds never influence the chosen hyperparameters. Default
#'   `FALSE`. Raises compute cost by the grid size per fold.
#' @param tune_folds Integer. Folds for the inner grid-search scoring.
#'   Default `3`.
#' @param prevalence_target Numeric in (0, 1), or `NULL`. Passed to
#'   [cast_fit()] in every outer fold: presence/background rows are weighted
#'   so the fitted prevalence equals the target (weight = target/n per class).
#'   Evaluation metrics remain unweighted. Default `NULL`.
#' @param n_repeat Integer >= 1. Number of independent spatial-blocking
#'   repeats (Barbet-Massin et al. 2012-style replicated evaluation).
#'   Each repeat redraws the fold assignment under its own seed, so the
#'   spread across repeats isolates blocking randomness from sampling noise.
#'   With `n_repeat = 1` (default) `metrics` averages and standardises
#'   across folds as before. With `n_repeat > 1`, `metrics` first averages
#'   within each repeat, then reports the mean across repeats and the
#'   **between-repeat** standard deviation; `fold_metrics` gains a `repeat`
#'   column; `selection_freq` pools all repeats; `oof`, `thresholds` and
#'   `folds` refer to the first repeat.
#' @param threshold_method Rule for selecting the binary threshold inside each
#'   outer training fold. That threshold is then applied to the held-out fold.
#'   See [cast_threshold()]. Default `"max_tss"`.
#' @param aoa Logical. Calibrate the area-of-applicability (AOA)
#'   dissimilarity threshold (Meyer & Pebesma 2021). When `TRUE`, every fold
#'   computes the DI of its held-out rows against the fold's
#'   (buffer-thinned) training rows in the standardized predictor space, and
#'   the stored threshold is the 95th percentile of the pooled held-out DI
#'   values. The calibration is model-independent (uniform predictor
#'   weights). The result is stored in the `aoa` component of the returned
#'   object, together with the per-row held-out DI (`di_oof`), and is
#'   consumed by [cast_predict()] (`aoa = TRUE`) and
#'   [cast_ensemble_raster()] (`aoa_cv = TRUE`). Requires \pkg{FNN}.
#' @param parallel Run folds with `future.apply`; requires the user to set a
#'   [future::plan()] first (a warning is issued when none is set).
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_cv` object including fold-level selections and `fold_status`.
#'   Metrics include `auc_n_folds`, `tss_n_folds`, and `cbi_n_folds` counting
#'   finite values separately; `n_folds` counts model result rows.
#'   Empty or failed folds do not contribute predictive metrics. When no fold
#'   is evaluable, the `cast_cv_no_evaluable_folds` error carries `screens` and
#'   `fold_status` for inspection. With `aoa = TRUE` the `aoa` component
#'   carries the calibrated `threshold` (NA when fewer than 10 finite
#'   held-out DI values were collected).
#' @export
cast_cv <- function(data,
                    screen = NULL,
                    select_method = "two_stage",
                    select_args = list(),
                    k = 5L,
                    models = c("rf"),
                    block_method = c("grid", "grid_random", "cluster", "env"),
                    buffer = 0,
                    buffer_unit = c("deg", "km"),
                    response = "presence",
                    rf_ntree = 500L,
                    brt_n_trees = 2000L,
                    brt_shrinkage = 0.005,
                    tune = FALSE,
                    tune_folds = 3L,
                    prevalence_target = NULL,
                    n_repeat = 1L,
                    threshold_method = "max_tss",
                    aoa = FALSE,
                    parallel = FALSE,
                    seed = NULL,
                    verbose = TRUE) {
  block_method <- match.arg(block_method)
  if (!is.logical(aoa) || length(aoa) != 1L || is.na(aoa)) {
    cli::cli_abort("{.arg aoa} must be TRUE or FALSE.")
  }
  if (aoa) check_suggested("FNN", "for the area of applicability (aoa = TRUE)")
  if (!is.list(select_args) ||
      (length(select_args) && (is.null(names(select_args)) ||
       anyNA(names(select_args)) || any(!nzchar(names(select_args))) ||
       anyDuplicated(names(select_args))))) {
    cli::cli_abort("{.arg select_args} must be a uniquely named list.")
  }
  if (any(c("data", "response", "method") %in% names(select_args))) {
    cli::cli_abort("{.arg select_args} cannot override data, response, or method; selection must use the outer training fold.")
  }
  check_suggested("pROC", "for fold evaluation metrics")
  k <- as.integer(k)
  if (k < 2L) cli::cli_abort("{.arg k} must be at least 2.")
  if (!is.numeric(buffer) || length(buffer) != 1L || buffer < 0) {
    cli::cli_abort("{.arg buffer} must be a single non-negative number.")
  }
  buffer_unit <- match.arg(buffer_unit)
  if (buffer_unit == "km" && buffer > 0) {
    # Haversine is defined on decimal degrees; silently treating projected
    # metres as degrees would inflate the exclusion band by orders of
    # magnitude (same guard as cast_thin()).
    bad_deg <- any(data$lon < -180 | data$lon > 180, na.rm = TRUE) ||
      any(data$lat < -90 | data$lat > 90, na.rm = TRUE)
    if (bad_deg) {
      cli::cli_abort(c(
        "{.code buffer_unit = \"km\"} requires decimal-degree lon/lat coordinates.",
        "x" = "The data hold coordinates outside [-180, 180] / [-90, 90].",
        "i" = "Reproject to lon/lat, or keep {.code buffer_unit = \"deg\"} on projected grids."
      ))
    }
  }
  n_repeat <- as.integer(n_repeat)
  if (is.na(n_repeat) || n_repeat < 1L) {
    cli::cli_abort("{.arg n_repeat} must be an integer >= 1.")
  }
  validate_species_data(data, required_cols = c("lon", "lat", response),
                        response = response)
  if (is.null(select_method) && !is.null(screen)) {
    cli::cli_warn(c(
      "{.code select_method = NULL} reuses the supplied {.arg screen} in every fold.",
      i = "That screen was selected on the full data set, so selection leaks into these CV metrics; treat them as optimistic."
    ))
  }

  # Environmental blocking needs the numeric predictor matrix (coordinates
  # and the response excluded); computed once, reused by every repeat.
  env_block <- if (identical(block_method, "env")) {
    cand <- setdiff(names(data), c("lon", "lat", response))
    cand <- cand[vapply(data[cand], is.numeric, logical(1))]
    if (length(cand) < 2L) {
      cli::cli_abort("{.code block_method = \"env\"} needs at least two numeric predictor columns.")
    }
    as.matrix(data[, cand, drop = FALSE])
  } else NULL

  cv_one <- function(folds, fold_i) {
    test_idx <- which(folds == fold_i)
    train_idx <- .cast_buffer_train_idx(data$lon, data$lat, test_idx, buffer,
                                        buffer_unit = buffer_unit)
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
    if (is.null(fold_screen)) return(list(status = "selection_error"))
    if (!length(fold_screen$selected)) {
      return(list(status = fold_screen$diagnostics$status %||% "empty_selection",
                  selected = character(0), screen = fold_screen))
    }

    fit <- tryCatch(
      cast_fit(
        train, screen = fold_screen, models = models, response = response,
        rf_ntree = rf_ntree, brt_n_trees = brt_n_trees,
        brt_shrinkage = brt_shrinkage,
        tune = tune, tune_folds = tune_folds,
        prevalence_target = prevalence_target,
        seed = if (is.null(seed)) NULL else seed + 100L + fold_i,
        verbose = FALSE
      ),
      error = function(e) {
        warning(sprintf("Model fitting failed in fold %d: %s", fold_i, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit)) {
      return(list(status = "model_error", selected = fold_screen$selected,
                  screen = fold_screen))
    }

    # Held-out AOA calibration (model-independent): DI of every test row to
    # its nearest training-row neighbour in the standardized fold space. The
    # training side uses the same buffer-thinned rows the models saw.
    di_fold <- NULL
    if (aoa) {
      di_fold <- tryCatch({
        x_tr <- as.data.frame(train[, fit$env_vars, drop = FALSE],
                              check.names = FALSE)
        x_te <- as.data.frame(test[, fit$env_vars, drop = FALSE],
                              check.names = FALSE)
        for (nm in fit$env_vars) {
          x_tr[[nm]] <- as.numeric(x_tr[[nm]])
          x_te[[nm]] <- as.numeric(x_te[[nm]])
        }
        prep <- .cast_aoa_prep(x_tr, x_te)
        .cast_aoa_di(prep$train, prep$new)
      }, error = function(e) {
        warning(sprintf("AOA DI failed in fold %d: %s", fold_i, e$message))
        rep(NA_real_, nrow(test))
      })
    }

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
      fold_threshold <- NULL
      if (!is.null(fit$scaling$response) && !is.null(fit$scaling$reference)) {
        train_pred <- tryCatch(
          predict_single_model(info, fit$scaling$reference),
          error = function(e) rep(NA_real_, nrow(train)))
        fold_threshold <- tryCatch(cast_threshold(
          train_pred, fit$scaling$response, method = threshold_method),
          error = function(e) NA_real_)
      }
      met <- tryCatch(
        if (is.null(fold_threshold)) {
          evaluate_model_full(pred, test[[response]])
        } else {
          evaluate_model_full(pred, test[[response]], threshold = fold_threshold,
                              threshold_method = threshold_method)
        },
        error = function(e) c(auc = NA_real_, pr_auc = NA_real_, tss = NA_real_,
          sedi = NA_real_, brier = NA_real_, logloss = NA_real_, boyce = NA_real_,
          cbi = NA_real_, kappa = NA_real_, omission_5 = NA_real_,
          omission_10 = NA_real_, mpa = NA_real_, tss_threshold = NA_real_)
      )
      metric_value <- function(nm) {
        if (nm %in% names(met)) unname(met[[nm]]) else NA_real_
      }
      rows[[mdl]] <- data.frame(
        fold = fold_i, model = mdl, auc = metric_value("auc"),
        pr_auc = metric_value("pr_auc"), tss = metric_value("tss"),
        sedi = metric_value("sedi"), brier = metric_value("brier"),
        logloss = metric_value("logloss"), boyce = metric_value("boyce"),
        cbi = metric_value("cbi"), kappa = metric_value("kappa"),
        omission_5 = metric_value("omission_5"),
        omission_10 = metric_value("omission_10"),
        mpa = metric_value("mpa"),
        tss_threshold = metric_value("tss_threshold"),
        n_selected = length(fold_screen$selected),
        stringsAsFactors = FALSE
      )
      updates[[mdl]] <- list(idx = test_idx, pred = pred)
    }
    list(rows = rows, updates = updates, selected = fold_screen$selected,
         screen = fold_screen, di = di_fold,
         status = if (length(rows)) "evaluated" else "no_predictions")
  }

  run_once <- function(rep_i) {
    # Per-repeat fold assignment: seeds separated by a prime so distinct
    # repeats draw unrelated blockings under the same user seed.
    rep_seed <- if (!is.null(seed)) seed + (rep_i - 1L) * 7919L else NULL
    folds <- make_spatial_folds(data$lon, data$lat, k, block_method, rep_seed,
                                env = env_block)
    if (verbose && n_repeat > 1L) {
      cli::cli_inform("Spatial CV repeat {rep_i}/{n_repeat}.")
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
        seq_len(k), function(i) cv_one(folds, i), future.seed = TRUE
      )
    } else {
      results <- lapply(seq_len(k), function(i) cv_one(folds, i))
    }

    row_list <- list()
    oof <- stats::setNames(lapply(models, function(x) rep(NA_real_, nrow(data))), models)
    oof_di <- rep(NA_real_, nrow(data))
    selections <- vector("list", k)
    screens <- vector("list", k)
    fold_status <- rep("failed", k)
    skipped_single <- integer(0)
    for (i in seq_along(results)) {
      res <- results[[i]]
      if (is.null(res)) next
      if (isTRUE(res$skipped_single_class)) {
        skipped_single <- c(skipped_single, i)
        fold_status[i] <- "single_class"
        next
      }
      row_list <- c(row_list, res$rows)
      selections[i] <- list(res$selected)
      screens[i] <- list(res$screen)
      fold_status[i] <- res$status
      if (aoa && !is.null(res$di)) oof_di[which(folds == i)] <- res$di
      for (mdl in names(res$updates)) {
        upd <- res$updates[[mdl]]
        oof[[mdl]][upd$idx] <- upd$pred
      }
    }

    fold_df <- if (length(row_list)) do.call(rbind, row_list) else data.frame()
    # Keep the out-of-fold surface: it is the only labelled ensemble-scale
    # prediction available downstream, and cast_ensemble() needs it to
    # threshold the ensemble rather than averaging per-model thresholds.
    oof_df <- data.frame(obs = data[[response]])
    for (mdl in models) oof_df[[paste0("HSS_", mdl)]] <- oof[[mdl]]
    thresholds <- vapply(models, function(mdl) {
      pred <- oof[[mdl]]
      ok <- is.finite(pred)
      if (sum(ok) < 10L || length(unique(data[[response]][ok])) < 2L) return(0.5)
      cast_threshold(pred[ok], data[[response]][ok], method = threshold_method)
    }, numeric(1))

    list(folds = folds, fold_df = fold_df, oof_df = oof_df,
         thresholds = thresholds, selections = selections,
         screens = screens, fold_status = fold_status,
         skipped_single = skipped_single, di = oof_di)
  }

  if (verbose) {
    cli::cli_inform(
      "Nested spatial CV: {k} fold{?s} x {n_repeat} repeat{?s}; selector={select_method %||% 'fixed'}."
    )
  }
  reps <- lapply(seq_len(n_repeat), run_once)

  total_skipped <- sum(vapply(reps, function(rp) length(rp$skipped_single),
                              integer(1)))
  if (total_skipped) {
    cli::cli_warn(
      "Skipped {total_skipped} fold run{?s} across {n_repeat} repeat{?s} x {k} folds: a single response class in the train or test split."
    )
  }

  # Fold status per fold position across repeats: identical states pass
  # through; disagreement is surfaced as "mixed" rather than hidden.
  fold_status_final <- if (n_repeat == 1L) reps[[1]]$fold_status else
    vapply(seq_len(k), function(i) {
      st <- vapply(reps, function(rp) rp$fold_status[i], character(1))
      if (all(st == st[1])) st[1] else "mixed"
    }, character(1))

  parts <- lapply(seq_len(n_repeat), function(ri) {
    if (!nrow(reps[[ri]]$fold_df)) return(NULL)
    if (n_repeat == 1L) return(reps[[ri]]$fold_df)
    cbind(data.frame(replicate = ri, stringsAsFactors = FALSE),
          reps[[ri]]$fold_df)
  })
  fold_df_all <- do.call(rbind, Filter(Negate(is.null), parts))
  if (is.null(fold_df_all)) fold_df_all <- data.frame()
  if (n_repeat > 1L && nrow(fold_df_all)) {
    fold_df_all <- fold_df_all[order(fold_df_all$replicate,
                                     fold_df_all$fold), , drop = FALSE]
    rownames(fold_df_all) <- NULL
  }

  # Fold-level selection frequency pooled over ALL repeats: the fraction of
  # fold runs in which each predictor was retained. The basis of the
  # consensus selector (cast_consensus()) and of the spatial-stability
  # diagnostic. The denominator is every fold run (k x n_repeat): an empty
  # (or failed) fold contributes zero rather than silently shrinking the
  # denominator and inflating frequency.
  selections_all <- do.call(c, lapply(reps, `[[`, "selections"))
  all_vars <- unique(unlist(selections_all))
  selection_freq <- data.frame(variable = character(0), freq = numeric(0),
                               stringsAsFactors = FALSE)
  if (length(all_vars)) {
    freq <- vapply(all_vars, function(v) {
      mean(vapply(selections_all, function(s) v %in% (s %||% character(0)),
                  logical(1)))
    }, numeric(1))
    selection_freq <- data.frame(
      variable = all_vars, freq = unname(freq), stringsAsFactors = FALSE)
    selection_freq <- selection_freq[order(-selection_freq$freq), , drop = FALSE]
    rownames(selection_freq) <- NULL
  }

  if (!nrow(fold_df_all)) {
    cli::cli_abort("All spatial CV folds failed.",
                   class = "cast_cv_no_evaluable_folds",
                   screens = reps[[1]]$screens,
                   fold_status = fold_status_final)
  }

  metrics <- if (n_repeat == 1L) {
    .cast_fold_metrics(reps[[1]]$fold_df, models)
  } else {
    per_rep <- lapply(reps, function(rp) .cast_fold_metrics(rp$fold_df, models))
    .cast_repeat_metrics(per_rep, models)
  }

  # AOA calibration from the first repeat (same convention as the oof
  # surface): the threshold is the 95th percentile of the pooled held-out
  # fold DI values. Fewer than 10 finite values leaves the threshold NA -
  # too sparse to be anything but noise.
  aoa_info <- NULL
  if (aoa) {
    di <- reps[[1]]$di
    vals <- di[is.finite(di)]
    if (length(vals) < 10L) {
      cli::cli_warn(
        "Only {length(vals)} finite held-out DI value{?s}; the AOA threshold is NA."
      )
    }
    aoa_info <- list(
      enabled = TRUE,
      weighted = FALSE,
      threshold = if (length(vals) >= 10L) {
        as.numeric(stats::quantile(vals, 0.95, names = FALSE))
      } else NA_real_,
      di_oof = di,
      n_di = length(vals),
      method = "cv_heldout_di_p95"
    )
  }

  new_cast_cv(
    metrics = metrics, fold_metrics = fold_df_all, folds = reps[[1]]$folds,
    k = k, block_method = block_method, thresholds = reps[[1]]$thresholds,
    selections = selections_all, screens = do.call(c, lapply(reps, `[[`, "screens")),
    selection_freq = selection_freq, oof = reps[[1]]$oof_df,
    fold_status = fold_status_final, threshold_method = threshold_method,
    n_repeat = n_repeat, aoa = aoa_info
  )
}

#' Aggregate Fold-Level Metrics Per Model (single blocking)
#'
#' Mean/sd across folds for every metric, plus the finite-value counts and
#' the mean tuned threshold. Shared by the single-repeat and repeated CV
#' aggregation paths.
#' @keywords internal
#' @noRd
.cast_fold_metrics <- function(fold_df, models) {
  agg <- lapply(models, function(mdl) {
    z <- fold_df[fold_df$model == mdl, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    metrics <- c("auc", "pr_auc", "tss", "sedi", "brier", "logloss", "boyce", "cbi",
                 "kappa", "omission_5", "omission_10", "mpa")
    out <- list(model = mdl)
    for (metric in metrics) {
      vals <- z[[metric]][is.finite(z[[metric]])]
      out[[paste0(metric, "_mean")]] <- if (length(vals)) mean(vals) else NA_real_
      out[[paste0(metric, "_sd")]] <- if (length(vals) >= 2L) stats::sd(vals) else NA_real_
      out[[paste0(metric, "_n_folds")]] <- length(vals)
    }
    thr <- z$tss_threshold[is.finite(z$tss_threshold)]
    out$tss_threshold_mean <- if (length(thr)) mean(thr) else NA_real_
    out$n_folds <- nrow(z)
    out$n_selected_mean <- mean(z$n_selected, na.rm = TRUE)
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), agg))
}

#' Aggregate Metrics Across Blocking Repeats
#'
#' Two-level aggregation for `n_repeat > 1`: per-repeat fold means are
#' averaged across repeats (`*_mean`), the **between-repeat** standard
#' deviation of those means is reported (`*_sd`; blocking randomness, not
#' fold noise), and fold counts are summed.
#' @keywords internal
#' @noRd
.cast_repeat_metrics <- function(per_rep, models) {
  agg <- lapply(models, function(mdl) {
    rows <- Filter(Negate(is.null), lapply(per_rep, function(m) {
      z <- m[m$model == mdl, , drop = FALSE]
      if (!nrow(z)) NULL else z
    }))
    if (!length(rows)) return(NULL)
    metric_names <- c("auc", "pr_auc", "tss", "sedi", "brier", "logloss",
                      "boyce", "cbi", "kappa", "omission_5", "omission_10", "mpa")
    out <- list(model = mdl)
    for (metric in metric_names) {
      means <- vapply(rows, function(z) z[[paste0(metric, "_mean")]][1],
                      numeric(1))
      out[[paste0(metric, "_mean")]] <- mean(means)
      out[[paste0(metric, "_sd")]] <- if (length(means) >= 2L) {
        stats::sd(means)
      } else NA_real_
      out[[paste0(metric, "_n_folds")]] <-
        sum(vapply(rows, function(z) z[[paste0(metric, "_n_folds")]][1],
                   numeric(1)))
    }
    out$tss_threshold_mean <- mean(vapply(rows, function(z)
      z$tss_threshold_mean[1], numeric(1)))
    out$n_selected_mean <- mean(vapply(rows, function(z)
      z$n_selected_mean[1], numeric(1)))
    out$n_folds <- sum(vapply(rows, function(z) z$n_folds[1], numeric(1)))
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), agg))
}

#' Assign Spatial or Environmental Folds
#'
#' @param lon,lat Numeric coordinate vectors.
#' @param k Number of folds.
#' @param method `"grid"` (default): bin coordinates into grid cells and group
#'   cells into `k` spatially contiguous folds via k-means on cell centroids,
#'   so every fold is a connected region. `"grid_random"`: legacy behaviour
#'   (greedy count-balanced packing that ignores cell position; kept for
#'   backwards comparability). `"cluster"`: k-means on the (scaled) point
#'   coordinates. `"env"`: k-means on the scaled predictor matrix (`env`),
#'   so folds are contiguous in environmental space; gaps are imputed with
#'   column medians and zero-variance columns dropped before scaling.
#' @param seed Random seed.
#' @param env Optional numeric matrix with one row per record, used only by
#'   `method = "env"`.
#' @return Integer fold ids in `1:k` (degenerate inputs collapse to `1L`).
#' @keywords internal
#' @noRd
make_spatial_folds <- function(lon, lat, k,
                               method = c("grid", "grid_random", "cluster", "env"),
                               seed = NULL, env = NULL) {
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
  if (method == "env") {
    if (is.null(env) || nrow(env) != n) {
      cli::cli_abort("{.code block_method = \"env\"} requires a predictor matrix with one row per record.")
    }
    env <- as.matrix(env)
    storage.mode(env) <- "numeric"
    # Impute gaps with column medians and drop zero-variance columns so
    # scale()/kmeans() stay defined on sparse or constant predictors.
    med <- vapply(seq_len(ncol(env)), function(j) stats::median(env[, j], na.rm = TRUE),
                  numeric(1))
    med[!is.finite(med)] <- 0
    for (j in seq_len(ncol(env))) env[!is.finite(env[, j]), j] <- med[j]
    keep <- vapply(seq_len(ncol(env)), function(j) diff(range(env[, j])) > 0, logical(1))
    if (sum(keep) < 2L) {
      cli::cli_abort("{.code block_method = \"env\"} needs at least two predictors with non-zero variance.")
    }
    return(as.integer(factor(stats::kmeans(
      scale(env[, keep, drop = FALSE]), centers = k_use, nstart = 10L
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
#' `test_idx` whose distance to the nearest test row is at least `buffer`.
#' Rows inside the band are dropped, thinning spatial-autocorrelation
#' leakage at the fold boundary.
#'
#' @param lon,lat Numeric coordinate vectors.
#' @param test_idx Integer indices of the held-out fold.
#' @param buffer Non-negative exclusion distance; `0` keeps all non-test rows.
#' @param buffer_unit `"deg"` planar degrees (default) or `"km"`
#'   great-circle kilometres via Haversine.
#' @keywords internal
#' @noRd
.cast_buffer_train_idx <- function(lon, lat, test_idx, buffer = 0,
                                   buffer_unit = "deg") {
  all_idx <- seq_along(lon)
  if (buffer <= 0) return(setdiff(all_idx, test_idx))
  # Only distances to the test rows matter; the full n x n matrix is quadratic
  # in the number of records and exhausts memory on national data sets.
  near <- rep(FALSE, length(all_idx))
  if (identical(buffer_unit, "km")) {
    for (i in test_idx) {
      d <- .haversine_km(lat, lon, lat[i], lon[i])
      near <- near | (d < buffer)
    }
  } else {
    buf2 <- buffer^2
    for (i in test_idx) {
      near <- near | ((lon - lon[i])^2 + (lat - lat[i])^2 < buf2)
    }
  }
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
