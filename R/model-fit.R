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
#' @param rf_ntree Integer. Number of RF trees. Default `300`.
#' @param brt_n_trees Integer. Number of BRT iterations. Default `500`.
#' @param brt_depth Integer. BRT tree depth. Default `5`.
#' @param num_threads Integer. Threads for the Random Forest learner. Default
#'   `1` (safe under fold-parallel cross-validation; raise for a single fit).
#' @param seed Integer or `NULL`. Base random seed.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A `cast_fit` object containing fitted models and metadata.
#'
#' @details
#' ## Supported Models
#' - **RF**: [ranger::ranger()] with probability output.
#' - **MaxEnt**: [maxnet::maxnet()] with logistic output.
#' - **BRT**: [gbm::gbm()] with Bernoulli loss and 5-fold CV.
#' - **GAM**: [mgcv::gam()] with thin-plate splines.
#'
#' @seealso [cast_select()], [cast_evaluate()], [cast_predict()]
#'
#' @export
cast_fit <- function(data,
                     screen       = NULL,
                     models       = c("rf", "brt", "maxent", "gam"),
                     response     = "presence",
                     rf_ntree     = 300L,
                     brt_n_trees  = 500L,
                     brt_depth    = 5L,
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

  # ---- Determine variables ------------------------------------------------
  env_vars <- if (!is.null(screen)) {
    if (!length(screen$selected)) {
      cli::cli_abort(c(
        "The supplied {.arg screen} has an empty {.field selected} set.",
        i = "Refit {.fun cast_select} (e.g. lower {.arg alpha} or raise {.arg min_vars}) or pass {.code screen = NULL} to use all predictors."
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
    fitted_models[[mdl]] <- tryCatch(
      fit_traditional(mdl, X_raw, Y, rf_ntree, brt_n_trees, brt_depth, seed,
                      num_threads),
      error = function(e) {
        cli::cli_warn("{mdl} failed: {e$message}")
        list(type = "traditional", model = NULL, name = mdl)
      }
    )
  }

  new_cast_fit(
    models    = fitted_models,
    cast_vars = cast_vars,
    env_vars  = env_vars,
    scaling   = list(means = X_means, sds = X_sds, impute = X_impute,
                     reference = X_raw),
    screen    = screen
  )
}

#' Impute missing predictor values from stored training statistics
#'
#' Fills `NA`s in a predictor `data.frame` column-by-column using the
#' training-set imputation vector stored in `fit$scaling$impute` (median of
#' each predictor at fit time). Columns absent from `impute` fall back to `0`.
#' Centralizing this keeps train, hold-out, spatial-CV, prediction, and
#' counterfactual pathways on the same, statistically defensible fill.
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
                            brt_depth, seed, num_threads = 1L) {
  switch(name,
    "rf" = {
      check_suggested("ranger", "for Random Forest")
      if (!is.null(seed)) set.seed(seed)
      m <- ranger::ranger(
        presence ~ .,
        data = cbind(presence = as.factor(Y), X),
        num.trees = rf_ntree, probability = TRUE, seed = seed %||% 42L,
        num.threads = as.integer(num_threads), verbose = FALSE
      )
      list(type = "traditional", model = m, name = "rf")
    },
    "maxent" = {
      check_suggested("maxnet", "for MaxEnt")
      m <- maxnet::maxnet(
        p = Y, data = X,
        maxnet::maxnet.formula(p = Y, data = X)
      )
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
        shrinkage = 0.01,
        cv.folds = 5L,
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
