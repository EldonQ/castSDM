#' Fit Species Distribution Models
#'
#' Trains one or more SDM algorithms on prepared data. Supported models:
#' Random Forest (RF), Boosted Regression Trees (BRT), MaxEnt, and
#' Generalised Additive Models (GAM).
#'
#' Variable selection is driven by the `cast_select` object from
#' [cast_select()]. If no screen is provided, all environmental variables
#' detected in `data` are used.
#'
#' @param data A `data.frame` with `presence` column and predictor variables.
#' @param screen A `cast_select` object from [cast_select()], or `NULL`.
#' @param models Character vector. Models to fit: `"rf"`, `"maxent"`, `"brt"`,
#'   `"gam"`. Default `c("rf", "brt", "maxent", "gam")`.
#' @param response Character. Response column name. Default `"presence"`.
#' @param rf_ntree Integer. Number of RF trees. Default `500`.
#' @param rf_mtry Integer. RF variables considered at each split. Default
#'   `NULL` (the ranger default, `floor(sqrt(p))`).
#' @param brt_n_trees Integer. Number of BRT trees. Default `2000`.
#' @param brt_depth Integer. BRT tree depth. Default `5`. Ignored when
#'   `tune = TRUE` (the depth enters the search grid).
#' @param brt_shrinkage Numeric. BRT learning rate. Default `0.005`
#'   (with 2000 trees, the Elith et al. 2008 recipe for presence-background
#'   data; larger rates underfit unless the tree count grows with them).
#'   Ignored when `tune = TRUE` (the rate enters the search grid).
#' @param maxent_classes Character. MaxEnt feature classes as a single string
#'   of class letters, e.g. `"l"`, `"lq"`, `"lqh"`, `"lqhp"`. Default `NULL`
#'   (the maxnet default feature set). When `tune = TRUE`, a supplied value is
#'   held fixed while the regularisation multiplier is searched.
#' @param maxent_regmult Numeric. MaxEnt regularisation multiplier. Default
#'   `NULL` (i.e. `1`, the maxnet default). When `tune = TRUE`, a supplied
#'   value is held fixed while the feature classes are searched.
#' @param tune Logical. Run a small grid search per engine (RF `mtry`; BRT
#'   `interaction.depth` x `shrinkage`; MaxEnt feature classes x regularisation
#'   multiplier) and fit with the best-scoring combination. Default `FALSE`.
#'   See Details.
#' @param tune_folds Integer. Folds used by the inner grid-search scoring
#'   (BRT internal CV; MaxEnt stratified random folds). Default `3`.
#' @param num_threads Integer. Threads for the Random Forest learner. Default
#'   `1` (safe under fold-parallel cross-validation; raise for a single fit).
#' @param seed Integer or `NULL`. Base random seed.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A `cast_fit` object containing fitted models and metadata. Each
#'   entry of the `models` list additionally carries a `tune` component
#'   (grid, per-combination scores, best combination) when `tune = TRUE`.
#'
#' @details
#' ## Supported Models
#' - **RF**: [ranger::ranger()] with probability output.
#' - **MaxEnt**: [maxnet::maxnet()] with logistic output. The engine never
#'   clamps predictions internally; clamping to the training range is a
#'   package-level decision (see [cast_predict()] and
#'   [cast_ensemble_raster()], argument `clamp`).
#' - **BRT**: [gbm::gbm()] with Bernoulli loss and 5-fold CV.
#' - **GAM**: [mgcv::gam()] with thin-plate splines.
#'
#' ## Hyperparameter grids (`tune = TRUE`)
#' Small default grids, informed by the reference SDM toolkits (biomod2,
#' flexsdm) but trimmed for cost:
#' - **RF** searches `mtry` around `sqrt(p)` and is scored by out-of-bag TSS
#'   (no extra refits beyond the grid itself).
#' - **BRT** searches `interaction.depth` x `shrinkage`; the tree count stays
#'   adaptive via [gbm::gbm.perf()]. Each combination is scored by the TSS of
#'   its internal cross-validated predictions.
#' - **MaxEnt** searches feature classes x regularisation multiplier, scored
#'   by stratified random-fold TSS with the threshold selected on the
#'   training part of each fold.
#'
#' Combinations that fail to fit (e.g. very small folds) score `NA`; if no
#' combination scores, the fixed defaults are kept with a warning.
#' Hyperparameters supplied explicitly are held fixed rather than searched.
#' Selecting hyperparameters on inner data and evaluating on outer folds is
#' the nested-selection principle of Cawley & Talbot (2010): pass
#' `tune = TRUE` to [cast_cv()] to run the grid inside every outer training
#' fold, so the held-out folds never influence the chosen hyperparameters.
#'
#' @references
#' Cawley, G. C. & Talbot, N. L. C. (2010). On over-fitting in model selection
#' and subsequent selection bias in performance evaluation. *Journal of
#' Machine Learning Research*, 11, 2079-2107.
#'
#' Elith, J., Leathwick, J. R. & Hastie, T. (2008). A working guide to boosted
#' regression trees. *Journal of Animal Ecology*, 77(4), 802-813.
#'
#' @seealso [cast_select()], [cast_evaluate()], [cast_predict()]
#'
#' @export
cast_fit <- function(data,
                     screen       = NULL,
                     models       = c("rf", "brt", "maxent", "gam"),
                     response     = "presence",
                     rf_ntree     = 500L,
                     rf_mtry      = NULL,
                     brt_n_trees  = 2000L,
                     brt_depth    = 5L,
                     brt_shrinkage = 0.005,
                     maxent_classes = NULL,
                     maxent_regmult = NULL,
                     tune         = FALSE,
                     tune_folds   = 3L,
                     num_threads  = 1L,
                     seed         = NULL,
                     verbose      = TRUE) {
  models <- tolower(models)
  valid_models <- c("rf", "maxent", "brt", "gam")
  bad <- setdiff(models, valid_models)
  if (length(bad) > 0) {
    cli::cli_abort(
      "Unknown model(s): {.val {bad}}. Use one or more of: {.val {valid_models}}."
    )
  }
  if (!is.null(rf_mtry) &&
      (!is.numeric(rf_mtry) || length(rf_mtry) != 1L ||
       !is.finite(rf_mtry) || rf_mtry < 1)) {
    cli::cli_abort("{.arg rf_mtry} must be a single positive integer or NULL.")
  }
  if (!is.null(maxent_classes) &&
      (!is.character(maxent_classes) || length(maxent_classes) != 1L ||
       !grepl("^[lqpht]+$", maxent_classes))) {
    cli::cli_abort(c(
      "{.arg maxent_classes} must be a single string of feature-class letters.",
      i = "Allowed letters: l (linear), q (quadratic), p (product), h (hinge), t (threshold); e.g. {.val \"lqh\"}."
    ))
  }
  if (!is.null(maxent_regmult) &&
      (!is.numeric(maxent_regmult) || length(maxent_regmult) != 1L ||
       !is.finite(maxent_regmult) || maxent_regmult <= 0)) {
    cli::cli_abort("{.arg maxent_regmult} must be a single positive number.")
  }
  if (isTRUE(tune)) {
    tune_folds <- as.integer(tune_folds)
    if (is.na(tune_folds) || tune_folds < 2L) {
      cli::cli_abort("{.arg tune_folds} must be an integer >= 2.")
    }
  }

  # ---- Determine variables ------------------------------------------------
  env_vars <- if (!is.null(screen)) {
    if (!length(screen$selected)) {
      cli::cli_abort(c(
        "The supplied {.arg screen} has an empty {.field selected} set.",
        i = "Refit {.fun cast_select} with {.code method = \"full\"}, or pass {.code screen = NULL} to use all predictors."
      ))
    }
    screen$selected
  } else {
    get_env_vars(data, response)
  }
  cast_vars <- env_vars

  .cast_check_response(data[[response]], response)
  Y <- as.integer(data[[response]])
  X_raw <- as.data.frame(data[, env_vars, drop = FALSE], check.names = FALSE)
  # Reject non-numeric / factor predictors explicitly: silently coercing a
  # factor with as.numeric() would model its level codes, not its values.
  .cast_check_numeric_predictors(X_raw, arg = "data")
  for (col in names(X_raw)) X_raw[[col]] <- as.numeric(X_raw[[col]])

  # -- Training-set median imputation, reused by evaluate/predict/CV --
  X_impute <- vapply(X_raw, function(v) {
    m <- stats::median(v, na.rm = TRUE)
    if (is.finite(m)) m else 0
  }, numeric(1))
  X_raw <- .cast_impute(X_raw, X_impute)
  # -- Standardize (stored for prediction) --
  X_means <- colMeans(X_raw, na.rm = TRUE)
  X_sds   <- apply(X_raw, 2, stats::sd, na.rm = TRUE)
  X_sds[X_sds < 1e-10] <- 1

  # ---- Fit each model -----------------------------------------------------
  fitted_models <- list()
  for (mdl in models) {
    if (verbose) cli::cli_inform("Training {.val {mdl}}...")

    # Per-engine hyperparameters; explicit overrides win, tuning replaces
    # whichever dimensions were not pinned down.
    p_rf_mtry <- rf_mtry
    p_brt_depth <- brt_depth
    p_brt_shrinkage <- brt_shrinkage
    p_me_classes <- maxent_classes
    p_me_regmult <- maxent_regmult
    if (isTRUE(tune) && mdl %in% c("rf", "brt", "maxent")) {
      tune_info <- tryCatch(
        .cast_tune(
          mdl, X_raw, Y,
          tune_folds = tune_folds,
          rf_ntree = rf_ntree, brt_n_trees = brt_n_trees,
          fixed = list(maxent_classes = maxent_classes,
                       maxent_regmult = maxent_regmult),
          seed = seed, num_threads = num_threads, verbose = verbose
        ),
        error = function(e) {
          cli::cli_warn(
            "Hyperparameter tuning failed for {.val {mdl}} ({e$message}); keeping fixed defaults."
          )
          NULL
        }
      )
      if (!is.null(tune_info)) {
        p_rf_mtry <- tune_info$best$mtry
        p_brt_depth <- tune_info$best$depth %||% p_brt_depth
        p_brt_shrinkage <- tune_info$best$shrinkage %||% p_brt_shrinkage
        p_me_classes <- tune_info$best$classes %||% p_me_classes
        p_me_regmult <- tune_info$best$regmult %||% p_me_regmult
      }
    }

    fitted_models[[mdl]] <- tryCatch(
      fit_traditional(mdl, X_raw, Y, rf_ntree, brt_n_trees,
                      p_brt_depth, p_brt_shrinkage, seed, num_threads,
                      maxent_classes = p_me_classes,
                      maxent_regmult = p_me_regmult,
                      rf_mtry = p_rf_mtry),
      error = function(e) {
        cli::cli_abort(c(
          "Model {.val {mdl}} failed to fit.",
          "x" = "{e$message}",
          "i" = "Refit without this engine or inspect the training data."
        ))
      }
    )
    if (isTRUE(tune) && mdl %in% c("rf", "brt", "maxent") &&
        !is.null(tune_info)) {
      fitted_models[[mdl]]$tune <- tune_info
    }
  }

  new_cast_fit(
    models    = fitted_models,
    cast_vars = cast_vars,
    env_vars  = env_vars,
    scaling   = list(means = X_means, sds = X_sds, impute = X_impute,
                     reference = X_raw, response = Y, response_name = response),
    screen    = screen
  )
}

#' Impute missing predictor values from stored training statistics
#'
#' Fills `NA`s in a predictor `data.frame` column-by-column using the
#' training-set imputation vector stored in `fit$scaling$impute` (median of
#' each predictor at fit time). Columns absent from `impute` fall back to `0`.
#' Centralizing this keeps train, hold-out, spatial-CV, prediction, and
#' sensitivity pathways on the same, statistically defensible fill.
#'
#' @keywords internal
#' @noRd
.cast_impute <- function(X, impute = NULL) {
  X <- as.data.frame(X)
  for (col in names(X)) {
    if (!is.numeric(X[[col]])) X[[col]] <- as.numeric(X[[col]])
    na <- is.na(X[[col]])
    if (any(na)) {
      fill <- if (!is.null(impute) && col %in% names(impute) &&
                  is.finite(impute[[col]])) impute[[col]] else 0
      X[[col]][na] <- fill
    }
  }
  X
}


# ========================================================================
# Internal: Fit Traditional SDM
# ========================================================================

#' @keywords internal
#' @noRd
fit_traditional <- function(name, X, Y, rf_ntree, brt_n_trees,
                             brt_depth, brt_shrinkage, seed,
                             num_threads = 1L,
                             maxent_classes = NULL,
                             maxent_regmult = NULL,
                             rf_mtry = NULL) {
  switch(name,
    "rf" = {
      check_suggested("ranger", "for Random Forest")
      if (!is.null(seed)) set.seed(seed)
      m <- ranger::ranger(
        presence ~ .,
        data = cbind(presence = as.factor(Y), X),
        num.trees = rf_ntree, probability = TRUE, seed = seed %||% 42L,
        mtry = rf_mtry,  # NULL keeps the ranger default floor(sqrt(nvar))
        num.threads = as.integer(num_threads), verbose = FALSE
      )
      list(type = "traditional", model = m, name = "rf")
    },
    "maxent" = {
      m <- .fit_maxent(X, Y, classes = maxent_classes, regmult = maxent_regmult)
      list(type = "traditional", model = m, name = "maxent")
    },
    "brt" = {
      check_suggested("gbm", "for BRT")
      if (!is.null(seed)) set.seed(seed)
      m <- gbm::gbm(
        presence ~ .,
        data = cbind(presence = Y, X),
        distribution = "bernoulli",
        n.trees = brt_n_trees,
        interaction.depth = brt_depth,
        shrinkage = brt_shrinkage,
        cv.folds = 5L,
        # gbm() defaults n.cores to a parallel cluster. On a small CI runner
        # (or inside an already-parallel cross-validation) spawning workers
        # makes the fit fail outright rather than run slower, so the fold
        # parallelism the package controls is the only one used.
        n.cores = 1L,
        verbose = FALSE
      )
      bt <- gbm::gbm.perf(m, method = "cv", plot.it = FALSE)
      list(type = "traditional", model = m, name = "brt", best_trees = bt)
    },
    "gam" = {
      check_suggested("mgcv", "for GAM")
      df <- cbind(presence = Y, X)
      # Backtick-quote names: non-syntactic predictor names (spaces, leading
      # digits) otherwise break formula parsing and the error was swallowed
      # by the caller's tryCatch, silently dropping GAM from every fold.
      qname <- function(nm) sprintf("`%s`", gsub("`", "", nm))
      pred_terms <- vapply(seq_len(ncol(X)), function(i) {
        v <- X[, i]
        nm <- qname(colnames(X)[i])
        if (is.numeric(v) && length(unique(v)) >= 8L)
          sprintf("s(%s, k = 5)", nm) else nm
      }, character(1))
      f <- stats::as.formula(
        paste("presence ~", paste(pred_terms, collapse = " + "))
      )
      m <- tryCatch(
        mgcv::gam(f, data = df, family = stats::binomial(),
                  method = "REML"),
        error = function(e) {
          flin <- stats::as.formula(
            paste("presence ~",
                  paste(vapply(colnames(X), qname, character(1)),
                        collapse = " + "))
          )
          mgcv::gam(flin, data = df, family = stats::binomial())
        }
      )
      list(type = "traditional", model = m, name = "gam")
    }
  )
}

#' Fit a MaxEnt (maxnet) Model
#'
#' Single construction path shared by [cast_fit()] and the MaxEnt grid
#' search, so the single-predictor augmentation trick applies in both.
#'
#' @param X Predictor data.frame (imputed, raw units).
#' @param Y Binary 0/1 response.
#' @param classes Character feature-class string (e.g. `"lqh"`), or `NULL`
#'   for the maxnet default set.
#' @param regmult Numeric regularisation multiplier, or `NULL` for `1`.
#' @keywords internal
#' @noRd
.fit_maxent <- function(X, Y, classes = NULL, regmult = NULL) {
  check_suggested("maxnet", "for MaxEnt")
  f <- if (is.null(classes)) {
    maxnet::maxnet.formula(p = Y, data = X)
  } else {
    maxnet::maxnet.formula(p = Y, data = X, classes = classes)
  }
  single_predictor <- ncol(X) == 1L
  if (single_predictor) {
    # maxnet's background augmentation drops one-column data frames to vectors.
    pres <- X[Y == 1L, , drop = FALSE]
    add <- !pres[[1]] %in% X[[1]][Y == 0L]
    X <- rbind(X, pres[add, , drop = FALSE])
    Y <- c(Y, rep(0L, sum(add)))
    # glmnet requires two columns; its excluded constant leaves the feature set unchanged.
    if (ncol(stats::model.matrix(f, X)) == 1L) f <- stats::update(f, ~ . + 1)
  }
  maxnet::maxnet(
    p = Y, data = X, f = f,
    regmult = as.numeric(regmult %||% 1),
    addsamplestobackground = !single_predictor
  )
}

# ========================================================================
# Internal: Hyperparameter grid search (C1/C2)
# ========================================================================

#' TSS at the max-TSS threshold of the supplied predictions
#' @keywords internal
#' @noRd
.tss_score <- function(pred, obs) {
  ok <- is.finite(pred) & !is.na(obs)
  pred <- as.numeric(pred[ok]); obs <- as.integer(obs[ok])
  if (!length(pred) || !all(c(0L, 1L) %in% unique(obs))) return(NA_real_)
  thr <- tryCatch(cast_threshold(pred, obs, method = "max_tss"),
                  error = function(e) NA_real_)
  if (!is.finite(thr)) return(NA_real_)
  hit <- pred >= thr
  sum(hit & obs == 1L) / sum(obs == 1L) +
    sum(!hit & obs == 0L) / sum(obs == 0L) - 1
}

#' Run the Per-Engine Hyperparameter Grid Search
#'
#' Grids (biomod2/flexsdm-informed, trimmed for cost):
#' - rf: `mtry` around `sqrt(p)`, scored by out-of-bag TSS.
#' - brt: `interaction.depth` x `shrinkage`, scored by internal-CV TSS.
#' - maxent: feature classes x regularisation multiplier (each dimension held
#'   fixed when supplied via `fixed`), scored by stratified random-fold TSS.
#'
#' Returns `NULL` (keeping fixed defaults) when no combination scores.
#'
#' @keywords internal
#' @noRd
.cast_tune <- function(name, X, Y, tune_folds = 3L,
                       rf_ntree = 500L, brt_n_trees = 2000L,
                       fixed = list(), seed = NULL,
                       num_threads = 1L, verbose = TRUE) {
  p <- ncol(X)
  grid <- switch(name,
    rf = data.frame(
      mtry = unique(pmax(1L, pmin(p,
        as.integer(round(sqrt(p) * c(0.5, 1, 2)))))),
      stringsAsFactors = FALSE
    ),
    brt = expand.grid(depth = c(2, 5), shrinkage = c(0.005, 0.01)),
    maxent = {
      cls <- fixed$maxent_classes %||% c("l", "lq", "lqh", "lqhp")
      rms <- fixed$maxent_regmult %||% c(0.5, 1, 2)
      expand.grid(classes = cls, regmult = as.numeric(rms),
                  stringsAsFactors = FALSE)
    },
    cli::cli_abort("No tuning grid for engine {.val {name}}.")
  )

  scores <- switch(name,
    rf   = .tune_rf(X, Y, grid, rf_ntree = rf_ntree, seed = seed,
                    num_threads = num_threads),
    brt  = .tune_brt(X, Y, grid, brt_n_trees = brt_n_trees,
                     tune_folds = tune_folds, seed = seed),
    maxent = .tune_maxent(X, Y, grid, tune_folds = tune_folds, seed = seed)
  )
  names(scores) <- apply(grid, 1L, function(r) paste(names(r), r, sep = "=",
                                                     collapse = ", "))
  ok <- is.finite(scores)
  if (!any(ok)) {
    if (verbose) {
      cli::cli_warn("No hyperparameter combination scored for {.val {name}}; keeping fixed defaults.")
    }
    return(NULL)
  }
  best_i <- which.max(scores)  # which.max ignores NA; any(ok) is TRUE here
  best <- as.list(grid[best_i, , drop = FALSE])
  if (verbose) {
    cli::cli_inform(
      "Tuned {.val {name}}: {.val {names(scores)[best_i]}} (TSS = {round(scores[best_i], 3)})."
    )
  }
  metric <- switch(name, rf = "tss_oob", brt = "tss_cv", maxent = "tss_cv")
  list(
    engine = name, metric = metric, grid = grid, scores = scores,
    best = best, best_score = unname(scores[best_i])
  )
}

#' Grid-score RF by out-of-bag TSS (ranger probability forests carry OOB
#' predictions, so no extra evaluation splits are needed)
#' @keywords internal
#' @noRd
.tune_rf <- function(X, Y, grid, rf_ntree, seed, num_threads) {
  check_suggested("ranger", "for Random Forest tuning")
  vapply(seq_len(nrow(grid)), function(i) {
    tryCatch({
      if (!is.null(seed)) set.seed(seed + grid$mtry[i])
      m <- ranger::ranger(
        presence ~ ., data = cbind(presence = as.factor(Y), X),
        num.trees = rf_ntree, probability = TRUE, mtry = grid$mtry[i],
        seed = seed %||% 42L, num.threads = as.integer(num_threads),
        verbose = FALSE
      )
      .tss_score(m$predictions[, "1"], Y)
    }, error = function(e) NA_real_)
  }, numeric(1))
}

#' Grid-score BRT by internal-CV TSS (each combination is fit once with
#' `cv.folds`; `gbm.perf()` keeps the tree count adaptive)
#' @keywords internal
#' @noRd
.tune_brt <- function(X, Y, grid, brt_n_trees, tune_folds, seed) {
  check_suggested("gbm", "for BRT tuning")
  vapply(seq_len(nrow(grid)), function(i) {
    tryCatch({
      if (!is.null(seed)) set.seed(seed + i)
      m <- gbm::gbm(
        presence ~ ., data = cbind(presence = Y, X),
        distribution = "bernoulli", n.trees = brt_n_trees,
        interaction.depth = grid$depth[i], shrinkage = grid$shrinkage[i],
        cv.folds = tune_folds, n.cores = 1L, verbose = FALSE
      )
      if (is.null(m$cv.fitted)) return(NA_real_)
      .tss_score(as.numeric(m$cv.fitted), Y)
    }, error = function(e) NA_real_)
  }, numeric(1))
}

#' Grid-score MaxEnt by stratified random-fold TSS
#'
#' Presence and background rows are split into `tune_folds` separately, so
#' every fold sees both classes. The threshold is selected on the training
#' part of each fold and applied to the held-out part.
#' @keywords internal
#' @noRd
.tune_maxent <- function(X, Y, grid, tune_folds, seed) {
  if (!is.null(seed)) set.seed(seed)
  idx1 <- which(Y == 1L); idx0 <- which(Y == 0L)
  fold <- integer(length(Y))
  fold[idx1] <- sample(rep(seq_len(tune_folds), length.out = length(idx1)))
  fold[idx0] <- sample(rep(seq_len(tune_folds), length.out = length(idx0)))

  vapply(seq_len(nrow(grid)), function(i) {
    fold_scores <- vapply(seq_len(tune_folds), function(f) {
      tr <- which(fold != f); te <- which(fold == f)
      m <- tryCatch(
        .fit_maxent(X[tr, , drop = FALSE], Y[tr],
                    classes = grid$classes[i], regmult = grid$regmult[i]),
        error = function(e) NULL
      )
      if (is.null(m)) return(NA_real_)
      pred <- function(idx) tryCatch(
        as.numeric(stats::predict(m, X[idx, , drop = FALSE],
                                  type = "logistic", clamp = FALSE)),
        error = function(e) rep(NA_real_, length(idx))
      )
      p_tr <- pred(tr); p_te <- pred(te)
      thr <- tryCatch(cast_threshold(p_tr, Y[tr], method = "max_tss"),
                      error = function(e) NA_real_)
      if (!is.finite(thr)) return(NA_real_)
      hit <- p_te >= thr
      n1 <- sum(Y[te] == 1L); n0 <- sum(Y[te] == 0L)
      if (!n1 || !n0) return(NA_real_)
      sum(hit & Y[te] == 1L) / n1 + sum(!hit & Y[te] == 0L) / n0 - 1
    }, numeric(1))
    if (any(is.finite(fold_scores))) {
      mean(fold_scores[is.finite(fold_scores)])
    } else {
      NA_real_
    }
  }, numeric(1))
}
