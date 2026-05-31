#' Fit Species Distribution Models
#'
#' Trains one or more SDM algorithms on prepared data. Supported models:
#' Random Forest (RF), Boosted Regression Trees (BRT), MaxEnt, and
#' Generalised Additive Models (GAM).
#'
#' Variable selection is driven by the `cast_select` object from
#' [cast_select()]. If no screen is provided, all environmental variables
#' from the DAG (or data) are used.
#'
#' @param data A `data.frame` with `presence` column and predictor variables.
#' @param screen A `cast_select` object from [cast_select()], or `NULL`.
#' @param dag A [cast_dag] object, or `NULL`.
#' @param models Character vector. Models to fit: `"rf"`, `"maxent"`, `"brt"`,
#'   `"gam"`. Default `c("rf", "brt", "maxent", "gam")`.
#' @param response Character. Response column name. Default `"presence"`.
#' @param rf_ntree Integer. Number of RF trees. Default `300`.
#' @param brt_n_trees Integer. Number of BRT iterations. Default `500`.
#' @param brt_depth Integer. BRT tree depth. Default `5`.
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
                     dag          = NULL,
                     models       = c("rf", "brt", "maxent", "gam"),
                     response     = "presence",
                     rf_ntree     = 300L,
                     brt_n_trees  = 500L,
                     brt_depth    = 5L,
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
    screen$selected
  } else if (!is.null(dag)) {
    setdiff(dag$nodes, c(response, "lon", "lat"))
  } else {
    get_env_vars(data, response)
  }
  cast_vars <- env_vars

  Y <- data[[response]]
  X_raw <- as.data.frame(data[, env_vars, drop = FALSE])
  for (col in names(X_raw)) X_raw[[col]] <- as.numeric(X_raw[[col]])
  X_raw[is.na(X_raw)] <- 0

  # -- Standardize (stored for prediction) --
  X_means <- colMeans(X_raw, na.rm = TRUE)
  X_sds   <- apply(X_raw, 2, stats::sd, na.rm = TRUE)
  X_sds[X_sds < 1e-10] <- 1

  # ---- Fit each model -----------------------------------------------------
  fitted_models <- list()
  for (mdl in models) {
    if (verbose) cli::cli_inform("Training {.val {mdl}}...")
    fitted_models[[mdl]] <- tryCatch(
      fit_traditional(mdl, X_raw, Y, rf_ntree, brt_n_trees, brt_depth, seed),
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
    scaling   = list(means = X_means, sds = X_sds),
    dag       = dag,
    screen    = screen
  )
}


# ========================================================================
# Internal: Fit Traditional SDM
# ========================================================================

#' @keywords internal
#' @noRd
fit_traditional <- function(name, X, Y, rf_ntree, brt_n_trees,
                            brt_depth, seed) {
  switch(name,
    "rf" = {
      check_suggested("ranger", "for Random Forest")
      if (!is.null(seed)) set.seed(seed)
      m <- ranger::ranger(
        presence ~ .,
        data = cbind(presence = as.factor(Y), X),
        num.trees = rf_ntree, probability = TRUE, seed = seed %||% 42L
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
      pred_terms <- vapply(seq_len(ncol(X)), function(i) {
        v <- X[, i]
        nm <- colnames(X)[i]
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
            paste("presence ~", paste(colnames(X), collapse = " + "))
          )
          mgcv::gam(flin, data = df, family = stats::binomial())
        }
      )
      list(type = "traditional", model = m, name = "gam")
    }
  )
}
