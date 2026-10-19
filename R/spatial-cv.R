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
#'   `"covsel"`. Set to `NULL` only to evaluate a fixed supplied screen.
#' @param select_args Named list of additional [cast_select()] arguments.
#'   Cannot override `data`, `response`, or `method`.
#' @param k Number of outer spatial folds.
#' @param models Models passed to [cast_fit()].
#' @param block_method Spatial blocking strategy:
#'   - `"blockcv"` (default): \pkg{blockCV} spatial blocks
#'     ([blockCV::cv_spatial()]) with square blocks. The block size follows
#'     the range of spatial autocorrelation estimated by
#'     [blockCV::cv_spatial_autocor()] (Roberts et al. 2017), or the
#'     user-supplied `block_size`. Falls back to `"grid"` with a warning when
#'     \pkg{blockCV}/\pkg{sf} are unavailable or the design cannot be built.
#'   - `"grid"`: grid cells grouped into spatially contiguous folds via
#'     k-means on cell centroids (dependency-free fallback).
#'   - `"env"`: k-means on the scaled predictor matrix, so every fold is
#'     contiguous in **environmental** space rather than geographic space
#'     (Roberts et al. 2017-style environmental blocking).
#'   - `"external"`: set internally when `fold_ids` is supplied.
#' @param block_size Optional numeric block edge length in **metres**, passed
#'   to [blockCV::cv_spatial()] (it converts to degrees automatically for
#'   geographic coordinates). `NULL` (default) estimates the size from the
#'   range of spatial autocorrelation.
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
#' @param threshold_method Rule for selecting the binary threshold. Thresholds
#'   are computed once per model from the pooled **out-of-fold** predictions
#'   after all folds finish, then applied to every fold's held-out rows and to
#'   downstream products. Selecting a threshold on a fold's own in-sample
#'   training predictions inflates it for flexible learners and collapses the
#'   held-out TSS; see [cast_threshold()]. Default `"max_tss"`.
#' @param aoa Logical. Calibrate the area-of-applicability (AOA)
#'   dissimilarity threshold (Meyer & Pebesma 2021). When `TRUE`, every fold
#'   computes, in the standardized (optionally weighted) predictor space:
#'   the held-out DI of each test row to its nearest training neighbour, and
#'   the training-side DI of each training row to its nearest *other*
#'   training neighbour. The stored threshold is `Q3 + 1.5 * IQR` of the
#'   pooled training-side DI values (the published rule); the held-out DI
#'   (`di_oof`) and the share of held-out rows above the threshold
#'   (`share_outside_oof`) are validation evidence. The calibration is
#'   model-independent (uniform predictor weights). The result is stored in
#'   the `aoa` component of the returned object and is consumed by
#'   [cast_predict()] (`aoa = TRUE`) and
#'   [cast_ensemble_raster()] (`aoa_cv = TRUE`). Requires \pkg{FNN}.
#' @param parallel Run folds with `future.apply`; requires the user to set a
#'   [future::plan()] first (a warning is issued when none is set).
#' @param seed Random seed.
#' @param verbose Print progress.
#' @param rf_mtry,brt_depth,maxent_classes,maxent_regmult,num_threads Engine
#'   settings as in [cast_fit()], applied in every outer training fold.
#'
#' @return A `cast_cv` object including fold-level selections and `fold_status`.
#'   Metrics include `auc_n_folds`, `tss_n_folds`, and `cbi_n_folds` counting
#'   finite values separately; `n_folds` counts model result rows.
#'   Empty or failed folds do not contribute predictive metrics. When no fold
#'   is evaluable, the `cast_cv_no_evaluable_folds` error carries `screens` and
#'   `fold_status` for inspection. With `aoa = TRUE` the `aoa` component
#'   carries the calibrated `threshold` (NA when fewer than 10 finite
#'   training-side DI values were collected), `di_oof`, `di_train` and
#'   `share_outside_oof`.
#' @export
#' @param fold_ids Optional integer fold assignments from an external design
#'   such as `blockCV`, aligned with data rows. Values must cover `1:k`.
#'   Use `n_repeat = 1` for a fixed external design. This allows physical block
#'   sizes or basin groups to be chosen independently of castSDM's automatic
#'   blocking; supplied assignments are returned unchanged.
cast_cv <- function(data,
                    screen = NULL,
                    select_method = "covsel",
                    select_args = list(),
                    k = 5L,
                    models = c("rf"),
                    block_method = c("blockcv", "grid", "env"),
                    block_size = NULL,
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
                    verbose = TRUE,
                    rf_mtry = NULL, brt_depth = 5L,
                    maxent_classes = NULL, maxent_regmult = NULL,
                    num_threads = 1L, fold_ids = NULL) {
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
  k <- .cast_check_scalar(k, "k", lower = 2, integer = TRUE)
  if (!is.numeric(buffer) || length(buffer) != 1L || buffer < 0) {
    cli::cli_abort("{.arg buffer} must be a single non-negative number.")
  }
  # `Inf` passed the non-negative test above and then excluded every training
  # row, silently emptying the model side of each fold.
  if (!is.finite(buffer)) {
    cli::cli_abort("{.arg buffer} must be finite; {.val Inf} would exclude every training row.")
  }
  if (!is.null(block_size) &&
      (!is.numeric(block_size) || length(block_size) != 1L ||
       !is.finite(block_size) || block_size <= 0)) {
    cli::cli_abort("{.arg block_size} must be NULL or a single positive number of metres.")
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
  n_repeat <- .cast_check_scalar(n_repeat, "n_repeat", lower = 1, integer = TRUE)
  validate_species_data(data, required_cols = c("lon", "lat", response),
                        response = response)
  if (!is.null(fold_ids)) {
    if (!is.numeric(fold_ids) || length(fold_ids) != nrow(data) ||
        any(!is.finite(fold_ids)) || any(fold_ids != as.integer(fold_ids)) ||
        !setequal(unique(fold_ids), seq_len(k))) {
      cli::cli_abort("{.arg fold_ids} must assign every data row to an integer in 1:k, with all folds present.")
    }
    if (n_repeat != 1L) cli::cli_abort("Fixed {.arg fold_ids} require n_repeat = 1; evaluate different external designs separately.")
    fold_ids <- as.integer(fold_ids)
    block_method <- "external"
  }
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
        rf_mtry = rf_mtry, brt_depth = brt_depth, brt_shrinkage = brt_shrinkage,
        maxent_classes = maxent_classes, maxent_regmult = maxent_regmult,
        num_threads = num_threads,
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

    # AOA (model-independent, standardized fold space). Two DI sets per fold:
    #   held  -- every test row to its nearest training-row neighbour;
    #   train -- every training row to its nearest OTHER training neighbour
    #            (k = 2 skips the self-match), the distribution Meyer &
    #            Pebesma (2021) threshold.
    di_fold <- NULL
    di_fold_train <- NULL
    if (aoa) {
      aoa_di <- tryCatch({
        x_tr <- as.data.frame(train[, fit$env_vars, drop = FALSE],
                              check.names = FALSE)
        x_te <- as.data.frame(test[, fit$env_vars, drop = FALSE],
                              check.names = FALSE)
        for (nm in fit$env_vars) {
          x_tr[[nm]] <- as.numeric(x_tr[[nm]])
          x_te[[nm]] <- as.numeric(x_te[[nm]])
        }
        prep <- .cast_aoa_prep(x_tr, x_te)
        list(held = .cast_aoa_di(prep$train, prep$new),
             train = .cast_aoa_di_train(prep$train))
      }, error = function(e) {
        warning(sprintf("AOA DI failed in fold %d: %s", fold_i, e$message))
        NULL
      })
      if (!is.null(aoa_di)) {
        di_fold <- aoa_di$held
        di_fold_train <- aoa_di$train
      } else {
        di_fold <- rep(NA_real_, nrow(test))
        di_fold_train <- rep(NA_real_, nrow(train))
      }
    }

    # Fold predictions only. The threshold-dependent metrics (TSS, SEDI,
    # kappa) are computed in run_once() against the pooled out-of-fold
    # threshold for each model -- the in-sample training-prediction
    # threshold used previously is optimistically high for flexible
    # learners (RF/BRT) and made the held-out TSS collapse.
    updates <- list()
    pred_ok <- list()
    for (mdl in models) {
      info <- fit$models[[mdl]]
      if (is.null(info) || is.null(info$model)) next
      x_test <- as.data.frame(test[, fit$env_vars, drop = FALSE], check.names = FALSE)
      .cast_check_numeric_predictors(x_test, arg = "data")
      for (nm in names(x_test)) x_test[[nm]] <- as.numeric(x_test[[nm]])
      x_test <- .cast_impute(x_test, fit$scaling$impute)
      pred <- tryCatch(
        predict_single_model(info, x_test),
        error = function(e) {
          # A silent NA vector here used to flow straight into the OOF
          # surface and then into the ensemble threshold, so a model that
          # never predicted looked like a model that predicted poorly.
          warning(sprintf(
            "Prediction failed in fold %d for %s: %s",
            fold_i, mdl, conditionMessage(e)
          ))
          rep(NA_real_, nrow(test))
        }
      )
      pred_ok[[mdl]] <- length(pred) == nrow(test) && any(is.finite(pred))
      updates[[mdl]] <- list(idx = test_idx, pred = pred)
    }
    list(updates = updates, selected = fold_screen$selected,
         screen = fold_screen, di = di_fold, di_train = di_fold_train,
         pred_ok = pred_ok,
         status = if (any(unlist(pred_ok))) "evaluated" else "no_predictions")
  }

  run_once <- function(rep_i) {
    # Per-repeat fold assignment: seeds separated by a prime so distinct
    # repeats draw unrelated blockings under the same user seed.
    rep_seed <- if (!is.null(seed)) seed + (rep_i - 1L) * 7919L else NULL
    folds <- if (is.null(fold_ids)) {
      make_spatial_folds(data$lon, data$lat, k, block_method, rep_seed,
                         env = env_block, y = data[[response]],
                         block_size = block_size)
    } else fold_ids
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

    oof <- stats::setNames(lapply(models, function(x) rep(NA_real_, nrow(data))), models)
    oof_di <- rep(NA_real_, nrow(data))
    train_di <- numeric(0)
    selections <- vector("list", k)
    screens <- vector("list", k)
    fold_status <- rep("failed", k)
    skipped_single <- integer(0)
    # How many folds actually produced usable predictions per model. A model
    # can be present in every fold's result list and still have failed to
    # predict in some of them; without this counter the threshold guard below
    # only saw the TOTAL number of finite values, so a surface missing 40% of
    # its folds still produced a "valid" threshold.
    pred_folds <- stats::setNames(integer(length(models)), models)
    for (i in seq_along(results)) {
      res <- results[[i]]
      if (is.null(res)) next
      if (isTRUE(res$skipped_single_class)) {
        skipped_single <- c(skipped_single, i)
        fold_status[i] <- "single_class"
        next
      }
      selections[i] <- list(res$selected)
      screens[i] <- list(res$screen)
      fold_status[i] <- res$status
      if (aoa && !is.null(res$di)) oof_di[which(folds == i)] <- res$di
      if (aoa && !is.null(res$di_train)) train_di <- c(train_di, res$di_train)
      for (mdl in names(res$updates)) {
        upd <- res$updates[[mdl]]
        oof[[mdl]][upd$idx] <- upd$pred
      }
      for (mdl in names(res$pred_ok %||% list())) {
        if (isTRUE(res$pred_ok[[mdl]])) pred_folds[[mdl]] <- pred_folds[[mdl]] + 1L
      }
    }

    # Keep the out-of-fold surface: it is the only labelled ensemble-scale
    # prediction available downstream, and cast_ensemble() needs it to
    # threshold the ensemble rather than averaging per-model thresholds.
    oof_df <- data.frame(obs = data[[response]])
    for (mdl in models) oof_df[[paste0("HSS_", mdl)]] <- oof[[mdl]]
    thresholds <- vapply(models, function(mdl) {
      pred <- oof[[mdl]]
      ok <- is.finite(pred)
      # Fold coverage floor: at least 80% of folds (and never fewer than 2)
      # must have produced predictions. `sum(ok)` alone cannot detect a
      # surface that is complete within the folds that worked but missing
      # every fold where the model failed.
      n_folds_ok <- pred_folds[[mdl]]
      min_folds_ok <- max(2L, ceiling(k * 0.8))
      if (n_folds_ok < min_folds_ok) {
        cli::cli_warn(c(
          "{.val {mdl}} produced predictions in only {n_folds_ok}/{k} fold{?s}",
          i = "Fewer than the required {min_folds_ok} usable fold{?s}; the",
          i = "out-of-fold surface is incomplete and its threshold is not",
          i = "trustworthy. Threshold falls back to 0.5, which changes every",
          i = "threshold-based metric and any binarised product."
        ))
        return(0.5)
      }
      if (sum(ok) < 10L || length(unique(data[[response]][ok])) < 2L) {
        cli::cli_warn(c(
          "Fewer than 10 usable {.val {mdl}} out-of-fold values (or one class",
          "only); threshold falls back to 0.5, which silently changes every",
          "threshold-based metric and any binarised product."
        ))
        return(0.5)
      }
      cast_threshold(pred[ok], data[[response]][ok], method = threshold_method)
    }, numeric(1))

    # Fold metrics are assembled here, after the pooled out-of-fold
    # thresholds exist: AUC/Brier/Boyce and the area metrics do not need a
    # threshold, while TSS/SEDI/kappa are evaluated at the model's pooled
    # out-of-fold operating point (never at an in-sample training cut). Each
    # fold row records that one threshold, so `tss_threshold_mean` is the
    # operating point actually used.
    row_list <- list()
    for (i in seq_along(results)) {
      res <- results[[i]]
      if (is.null(res) || isTRUE(res$skipped_single_class)) next
      for (mdl in names(res$updates %||% list())) {
        upd <- res$updates[[mdl]]
        met <- tryCatch(
          evaluate_model_full(upd$pred, data[[response]][upd$idx],
                              threshold = thresholds[[mdl]],
                              threshold_method = threshold_method),
          error = function(e) .cast_na_full_metrics()
        )
        metric_value <- function(nm) {
          if (nm %in% names(met)) unname(met[[nm]]) else NA_real_
        }
        row_list[[length(row_list) + 1L]] <- data.frame(
          fold = i, model = mdl, auc = metric_value("auc"),
          pr_auc = metric_value("pr_auc"), tss = metric_value("tss"),
          sedi = metric_value("sedi"), brier = metric_value("brier"),
          logloss = metric_value("logloss"), boyce = metric_value("boyce"),
          cbi = metric_value("cbi"), kappa = metric_value("kappa"),
          omission_5 = metric_value("omission_5"),
          omission_10 = metric_value("omission_10"),
          mpa = metric_value("mpa"),
          tss_threshold = metric_value("tss_threshold"),
          n_selected = length(res$selected %||% character(0)),
          stringsAsFactors = FALSE
        )
      }
    }
    fold_df <- if (length(row_list)) do.call(rbind, row_list) else data.frame()

    list(folds = folds, fold_df = fold_df, oof_df = oof_df,
         thresholds = thresholds, selections = selections,
         screens = screens, fold_status = fold_status,
         skipped_single = skipped_single, di = oof_di,
         di_train = train_di, pred_folds = pred_folds)
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
  # surface). Threshold rule follows Meyer & Pebesma (2021): the 75th
  # percentile plus 1.5 x IQR of the pooled TRAINING-side nearest-other-
  # neighbour DI values (the distribution the published method thresholds).
  # The held-out DI is validation evidence, not the threshold source: its
  # share above the threshold is reported alongside. Fewer than 10 finite
  # training DI values leaves the threshold NA - too sparse to be anything
  # but noise.
  aoa_info <- NULL
  if (aoa) {
    di <- reps[[1]]$di
    di_tr <- reps[[1]]$di_train
    vals <- di_tr[is.finite(di_tr)]
    if (length(vals) < 10L) {
      cli::cli_warn(
        "Only {length(vals)} finite training-side DI value{?s}; the AOA threshold is NA."
      )
    }
    thr <- if (length(vals) >= 10L) {
      qs <- as.numeric(stats::quantile(vals, c(0.25, 0.75), names = FALSE))
      qs[2] + 1.5 * (qs[2] - qs[1])
    } else NA_real_
    di_held <- di[is.finite(di)]
    aoa_info <- list(
      enabled = TRUE,
      weighted = FALSE,
      threshold = thr,
      di_oof = di,
      di_train = di_tr,
      n_di = length(vals),
      share_outside_oof = if (is.finite(thr) && length(di_held)) {
        mean(di_held > thr)
      } else NA_real_,
      method = "train_di_q3_plus_1.5iqr"
    )
  }

  new_cast_cv(
    metrics = metrics, fold_metrics = fold_df_all, folds = reps[[1]]$folds,
    k = k, block_method = block_method, thresholds = reps[[1]]$thresholds,
    selections = selections_all, screens = do.call(c, lapply(reps, `[[`, "screens")),
    selection_freq = selection_freq, oof = reps[[1]]$oof_df,
    fold_status = fold_status_final, threshold_method = threshold_method,
    n_repeat = n_repeat, aoa = aoa_info,
    pred_folds = if (length(reps)) reps[[1]]$pred_folds else NULL
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
#' @param method `"blockcv"` (default): [blockCV::cv_spatial()] square blocks
#'   over an sf point layer, with class-balanced fold assignment when `y` is
#'   supplied; block size from `block_size` or the autocorrelation range.
#'   `"grid"`: grid cells grouped into `k` spatially contiguous folds via
#'   k-means on cell centroids (dependency-free fallback). `"env"`: k-means
#'   on the scaled predictor matrix (`env`); gaps imputed with column medians
#'   and zero-variance columns dropped before scaling.
#' @param seed Random seed.
#' @param env Optional numeric matrix with one row per record (method "env").
#' @param y Optional binary response vector; passed to blockCV for class-
#'   balanced fold assignment.
#' @param block_size Optional block edge length in metres (method "blockcv").
#' @return Integer fold ids in `1:k` (degenerate inputs collapse to `1L`).
#' @keywords internal
#' @noRd
make_spatial_folds <- function(lon, lat, k,
                               method = c("blockcv", "grid", "env"),
                               seed = NULL, env = NULL,
                               y = NULL, block_size = NULL) {
  method <- match.arg(method)
  n <- length(lon)
  k_use <- min(k, n)
  if (k_use < 2L) return(rep(1L, n))
  if (method == "blockcv") {
    folds <- .cast_blockcv_folds(lon, lat, k_use, seed, y, block_size)
    if (!is.null(folds)) return(folds)
    cli::cli_warn(c(
      "blockCV blocking is unavailable or could not be built for this data set;",
      i = "falling back to grid blocking."
    ))
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
    return(as.integer(factor(.cast_with_seed(seed, stats::kmeans(
      scale(env[, keep, drop = FALSE]), centers = k_use, nstart = 10L
    )$cluster))))
  }
  # "grid" (also the blockcv fallback): spatially contiguous folds. Cells are
  # the atoms; their centroids are clustered into k spatial groups and each
  # group becomes one fold. Degenerate coordinates (few distinct values)
  # collapse to a coarser grid rather than erroring.
  side <- ceiling(sqrt(k * 2))
  xb <- .cast_bin(lon, side)
  yb <- .cast_bin(lat, side)
  cell <- interaction(xb, yb, drop = TRUE)
  counts <- table(cell)
  cell_lon <- as.numeric(tapply(lon, cell, mean))
  cell_lat <- as.numeric(tapply(lat, cell, mean))
  n_cells <- length(cell_lon)
  if (n_cells <= k) {
    return(as.integer(factor(cell, levels = names(counts))))
  }
  km <- .cast_with_seed(seed,
    stats::kmeans(cbind(cell_lon, cell_lat), centers = k, nstart = 10L))
  fold_of_cell <- stats::setNames(km$cluster, names(counts))
  as.integer(fold_of_cell[as.character(cell)])
}

#' blockCV-Backed Fold Assignment
#'
#' Builds square spatial blocks with [blockCV::cv_spatial()] over an sf point
#' layer. Block size: `block_size` when supplied, else the range of spatial
#' autocorrelation estimated by [blockCV::cv_spatial_autocor()] (Roberts et
#' al. 2017); when the autocorrelation range is unavailable the blockCV
#' default sizing is used. Geographic coordinates (lon/lat within their
#' bounds) are tagged with CRS 4326 so blockCV converts metres to degrees;
#' projected coordinates are passed through unchanged (metres assumed).
#'
#' @return Integer fold ids, or `NULL` when blockCV/sf are unavailable or the
#'   design cannot be built (caller falls back to grid blocking).
#' @keywords internal
#' @noRd
.cast_blockcv_folds <- function(lon, lat, k, seed, y, block_size) {
  if (!requireNamespace("blockCV", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    return(NULL)
  }
  ok <- is.finite(lon) & is.finite(lat)
  if (sum(ok) < 2L * k) return(NULL)
  geographic <- all(lon[ok] >= -180 & lon[ok] <= 180) &&
    all(lat[ok] >= -90 & lat[ok] <= 90)
  df <- data.frame(.lon = lon, .lat = lat)
  if (!is.null(y)) df$.response <- as.integer(y)
  pts <- tryCatch(
    sf::st_as_sf(df, coords = c(".lon", ".lat"),
                 crs = if (geographic) 4326L else NA),
    error = function(e) NULL
  )
  if (is.null(pts) || any(sf::st_is_empty(pts))) return(NULL)
  column <- if (!is.null(y)) ".response" else NULL
  size <- block_size
  if (is.null(size)) {
    size <- tryCatch({
      ac <- suppressWarnings(blockCV::cv_spatial_autocor(
        x = pts, column = column,
        num_sample = min(5000L, nrow(pts)),
        plot = FALSE, progress = FALSE
      ))
      rng <- as.numeric(ac$range %||% NA_real_)
      if (is.finite(rng) && rng > 0) rng else NULL
    }, error = function(e) NULL)
  }
  blk <- tryCatch(
    suppressWarnings(blockCV::cv_spatial(
      x = pts, column = column, k = k, size = size,
      hexagon = FALSE, selection = "random", iteration = 100L,
      seed = seed, progress = FALSE, report = FALSE, plot = FALSE
    )),
    error = function(e) NULL
  )
  if (is.null(blk)) return(NULL)
  ids <- as.integer(blk$folds_ids)
  if (length(ids) != length(lon) || anyNA(ids)) return(NULL)
  # blockCV may return fewer than k non-empty folds on degenerate extents.
  if (!all(seq_len(k) %in% unique(ids))) return(NULL)
  ids
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

#' Bin a coordinate into `side` equal-width bins (degenerate-safe)
#'
#' Equal-WIDTH bins keep neighbouring cells neighbours on the raster, which
#' the k-means step downstream needs to recover fold contiguity. Equal-count
#' (quantile) bins were the legacy default but on a clustered distribution
#' they turn into long thin cells that straddle the whole domain.
#'
#' @keywords internal
#' @noRd
.cast_bin <- function(x, side) {
  x <- as.numeric(x)
  u <- sort(unique(x[!is.na(x)]))
  if (length(u) < 2L) return(rep(1L, length(x)))
  breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                length.out = side + 1L)
  if (length(breaks) < 2L) breaks <- c(u[1], u[length(u)])
  b <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  b[is.na(b)] <- 1L
  b
}
