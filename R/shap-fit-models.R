#' SHAP explanations for fitted RF models
#'
#' Uses the \pkg{fastshap} package (Monte Carlo Shapley with optional
#' sum-to-prediction adjustment) on the **same** fitted RF object produced by
#' [cast_fit()]. Interaction strengths for the circular plot are a **proxy**
#' derived from cross-products of mean absolute SHAP values.
#'
#' @param fit A [cast_fit()] object that already contains the requested model.
#' @param which Character: `"rf"`. The model must appear in
#'   `names(fit$models)`.
#' @param data Training-like `data.frame` with `response` and columns
#'   `fit$env_vars` (e.g. `split$train` from [cast_prepare()]).
#' @param response Binary response column name. Default `"presence"`.
#' @param test_fraction Fraction of rows (after `na.omit`) held out for rows in
#'   the SHAP matrix. Default `0.2`.
#' @param seed Optional integer seed.
#' @param fastshap_nsim Monte Carlo replications per feature column in
#'   [fastshap::explain()]. Default `60`. Increase for stability (slower).
#' @param max_background_rows Maximum training rows passed as `X` to
#'   `fastshap::explain()` (subsampled). Default `500`.
#' @param max_explain_rows Cap on the number of held-out rows explained.
#'   Default `80`.
#' @param verbose Logical.
#'
#' @return A `cast_shap` object compatible with [plot.cast_shap()]:
#'   `shap`, `bias_shap`, `shap_interaction`, `feature_names`, `base_score`,
#'   `expected_value`, `label`, `train_matrix`, plus `engine`, `method`,
#'   `fitted_model`, `shap_engine_note`.
#'
#' @details
#' Shapley values are on the **probability scale** (positive class),
#' matching `ranger` probability predictions.
#'
#' @seealso [cast_shap_xgb()], [plot.cast_shap()], [cast_fit()]
#'
#' @export
cast_shap_fit <- function(fit,
                          which = "rf",
                          data,
                          response = "presence",
                          test_fraction = 0.2,
                          seed = NULL,
                          fastshap_nsim = 60L,
                          max_background_rows = 500L,
                          max_explain_rows = 80L,
                          verbose = FALSE) {
  which <- match.arg(which, choices = "rf")
  fastshap_nsim <- max(2L, as.integer(fastshap_nsim)[1L])
  check_suggested("fastshap", "for cast_shap_fit()")

  if (!inherits(fit, "cast_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls cast_fit} object.")
  }
  if (!which %in% names(fit$models)) {
    cli::cli_abort(
      "Model {.val {which}} not found in {.code fit$models}: {.val {names(fit$models)}}."
    )
  }
  mdl <- fit$models[[which]]
  if (is.null(mdl$model) || !inherits(mdl$model, "ranger")) {
    cli::cli_abort("{.arg fit} does not contain a fitted ranger RF model.")
  }

  env_vars <- fit$env_vars
  df <- data[, unique(c(response, env_vars)), drop = FALSE]
  df <- stats::na.omit(df)
  if (nrow(df) < 30L) {
    cli::cli_abort("Too few complete rows ({nrow(df)}) for SHAP.")
  }

  Y <- as.numeric(df[[response]])
  X_raw <- as.data.frame(df[, env_vars, drop = FALSE])
  X_raw[is.na(X_raw)] <- 0

  n <- nrow(X_raw)
  if (!is.null(seed)) set.seed(seed)
  n_te <- max(5L, floor(n * test_fraction))
  ii <- sample.int(n, size = n - n_te)
  te <- setdiff(seq_len(n), ii)
  n_te <- min(length(te), max_explain_rows)
  te <- te[seq_len(n_te)]
  y_te <- Y[te]

  n_bg <- min(length(ii), max_background_rows)
  ii_use <- sample(ii, n_bg)

  rf_mod <- mdl$model
  X_bg <- X_raw[ii_use, , drop = FALSE]
  X_exp <- X_raw[te, , drop = FALSE]

  pfun_rf <- function(object, newdata) {
    nd <- as.data.frame(newdata)
    pr <- stats::predict(object, data = nd)
    pm <- as.matrix(pr$predictions)
    if (!("1" %in% colnames(pm))) {
      cli::cli_abort("Expected probability column {.val {'1'}} in ranger predictions.")
    }
    as.numeric(pm[, "1", drop = TRUE])
  }

  fs <- .fastshap_rows(
    object = rf_mod,
    X_bg = X_bg,
    X_exp = X_exp,
    pred_wrapper = pfun_rf,
    fastshap_nsim = fastshap_nsim,
    verbose = verbose
  )
  sh_mat <- fs$shap
  bias_sh <- fs$bias_shap
  exp_val <- mean(pfun_rf(rf_mod, X_exp))
  eng <- "ranger"
  note <- "fastshap (MC Shapley, adjust=TRUE) on ranger probability; Vint = proxy from |SHAP|."
  cap <- paste0(
    "SHAP for the ranger model inside cast_fit(). p=", ncol(sh_mat),
    " raw env columns (= fit$env_vars). fastshap on probability scale. ",
    "Vint edges: proxy from mean |SHAP| cross-products (not exact RF interaction SHAP)."
  )
  feat_sp <- "raw_env"
  sh_sc <- "probability"

  mean_inter <- .shap_interaction_proxy(sh_mat)
  base_score <- mean(bias_sh, na.rm = TRUE)
  n_feat <- ncol(sh_mat)

  structure(
    list(
      shap = sh_mat,
      bias_shap = bias_sh,
      shap_interaction = mean_inter,
      feature_names = colnames(sh_mat),
      base_score = base_score,
      xgb_model = NULL,
      expected_value = exp_val,
      train_matrix = X_exp,
      label = y_te,
      response = response,
      engine = eng,
      method = which,
      fitted_model = mdl,
      shap_engine_note = note,
      feature_space = feat_sp,
      shap_scale = sh_sc,
      n_features = n_feat,
      n_interactions = 0L,
      shap_plot_caption = cap
    ),
    class = "cast_shap"
  )
}


#' @keywords internal
#' @noRd
.shap_interaction_proxy <- function(sh_mat) {
  Ms <- abs(as.matrix(sh_mat))
  (t(Ms) %*% Ms) / pmax(1L, nrow(Ms))
}


#' Run fastshap::explain row-wise; attach bias_shap on matrix
#' @keywords internal
#' @noRd
.fastshap_rows <- function(object,
                           X_bg,
                           X_exp,
                           pred_wrapper,
                           fastshap_nsim,
                           verbose) {
  X_bg <- as.data.frame(X_bg)
  X_exp <- as.data.frame(X_exp)
  feat_names <- colnames(X_bg)
  if (!identical(feat_names, colnames(X_exp))) {
    cli::cli_abort("Background and explain feature columns differ.")
  }

  n_ex <- nrow(X_exp)
  sh_rows <- vector("list", n_ex)
  bias_scalar <- NA_real_

  for (ir in seq_len(n_ex)) {
    if (verbose && (ir == 1L || ir %% 5L == 0L)) {
      cli::cli_inform("fastshap row {ir}/{n_ex}...")
    }
    ex <- fastshap::explain(
      object,
      X = X_bg,
      pred_wrapper = pred_wrapper,
      newdata = X_exp[ir, , drop = FALSE],
      nsim = as.integer(fastshap_nsim)[1L],
      adjust = TRUE
    )
    if (is.na(bias_scalar)) {
      bias_scalar <- as.numeric(attr(ex, "baseline"))[1L]
    }
    sh_rows[[ir]] <- as.numeric(ex[1L, feat_names, drop = TRUE])
  }

  sh_mat <- do.call(rbind, sh_rows)
  colnames(sh_mat) <- feat_names
  list(
    shap = sh_mat,
    bias_shap = rep(bias_scalar, n_ex)
  )
}
