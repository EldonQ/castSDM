# Causal interpretation layer -------------------------------------------------
#
# Two functions turn a fitted castSDM workflow into causal-flavoured evidence:
#   * cast_effect()         - reads the conditional-importance estimates already
#                             produced by the causal screen (CPI impacts or DML
#                             partial-linear effects) and returns a tidy table
#                             with confidence intervals.
#   * cast_counterfactual() - g-computation "what-if" on the current climate:
#                             shift one predictor, hold the rest fixed, and map
#                             the change in predicted habitat suitability.
# Both reuse existing machinery (the screen; the fitted models) so no new
# estimation engine or hand-tuned knob is introduced.

#' @keywords internal
#' @noRd
.cast_extract_screen <- function(object) {
  screen <- if (inherits(object, "cast_select")) {
    object
  } else if (inherits(object, "cast_fit")) {
    object$screen
  } else if (inherits(object, "cast_result")) {
    object$screen
  } else {
    NULL
  }
  if (is.null(screen) || !inherits(screen, "cast_select")) {
    cli::cli_abort(
      "{.arg object} must be a {.cls cast_select}, {.cls cast_fit}, or {.cls cast_result} carrying a screen."
    )
  }
  screen
}

#' Causal Effect / Conditional Importance Table from the Screen
#'
#' Turns a causal screen into a tidy per-predictor table with confidence
#' intervals and FDR-adjusted significance.
#'
#' For `method = "cpi"` screens the table reports each predictor's conditional
#' predictive impact (CPI, Watson & Wright 2021): a non-negative measure of how
#' much predictive accuracy is lost when the predictor is replaced by a knockoff
#' given all other predictors. CPI is a magnitude of conditional dependence and
#' carries no sign (direction), so it should be read together with
#' [cast_counterfactual()] or partial-dependence for the shape of the response.
#'
#' For `method = "dml"` screens the table reports each predictor's
#' Neyman-orthogonal partial-linear effect on occurrence, expressed per one
#' within-sample standard deviation of the predictor, so magnitudes are
#' directly comparable and signed.
#'
#' @section Causal interpretation (read before citing effects):
#' CPI and DML both quantify a predictor's contribution *conditional on the
#' other predictors*. A causal reading is licensed only under the usual
#' assumptions - no unobserved confounding, a correctly specified adjustment
#' set, and no reverse causation - which observational, sampling-biased SDM
#' data rarely satisfy. Absent a defensible causal design, report these as
#' conditional-importance / adjusted-association evidence rather than validated
#' causal effects (Byrnes & Dee 2025).
#'
#' @param object A `cast_select` (from `method = "cpi"` or `"dml"`), or a
#'   `cast_fit` / `cast_result` that carries such a screen.
#' @param conf_level Confidence level for the intervals. Default `0.95`.
#'
#' @return A `cast_effect` object.
#' @references
#' Watson, D. S. & Wright, M. N. (2021). Testing conditional independence in
#' supervised learning algorithms. *Machine Learning*, 110(8), 2107-2129.
#'
#' Byrnes, J. E. K. & Dee, L. E. (2025). Causal inference with observational
#' data and unobserved confounding variables. *Ecology Letters*, 28(1), e70023.
#' @seealso [cast_select()], [cast_counterfactual()]
#' @export
cast_effect <- function(object, conf_level = 0.95) {
  screen <- .cast_extract_screen(object)
  if (!screen$method %in% c("cpi", "dml")) {
    cli::cli_abort(c(
      "{.fn cast_effect} needs a causal screen.",
      i = "Run {.code cast_select(..., method = \"cpi\")} (or {.val dml}) first."
    ))
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort("{.arg conf_level} must be a single number in (0, 1).")
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  sc <- screen$scores

  if (identical(screen$method, "cpi")) {
    sc <- sc[is.finite(sc$cpi) & is.finite(sc$std_error), , drop = FALSE]
    if (!nrow(sc)) {
      cli::cli_abort("The CPI screen holds no finite impact estimates.")
    }
    effects <- data.frame(
      variable    = sc$variable,
      estimate    = sc$cpi,
      std_error   = sc$std_error,
      statistic   = sc$statistic,
      p_value     = sc$p_value,
      p_adjusted  = sc$p_adjusted,
      conf_low    = sc$cpi - z * sc$std_error,
      conf_high   = sc$cpi + z * sc$std_error,
      selected    = sc$selected,
      stringsAsFactors = FALSE
    )
    effects <- effects[order(-effects$estimate), , drop = FALSE]
    measure <- "cpi"
  } else {
    sc <- sc[is.finite(sc$estimate) & is.finite(sc$std_error), , drop = FALSE]
    if (!nrow(sc)) {
      cli::cli_abort("The DML screen holds no finite effect estimates.")
    }
    effects <- data.frame(
      variable    = sc$variable,
      estimate    = sc$estimate,
      std_error   = sc$std_error,
      statistic   = sc$statistic,
      p_value     = sc$p_value,
      p_adjusted  = sc$p_adjusted,
      conf_low    = sc$estimate - z * sc$std_error,
      conf_high   = sc$estimate + z * sc$std_error,
      selected    = sc$selected,
      stringsAsFactors = FALSE
    )
    effects <- effects[order(-abs(effects$estimate)), , drop = FALSE]
    measure <- "dml_plr"
  }
  rownames(effects) <- NULL

  diagnostics <- screen$diagnostics
  diagnostics$measure <- measure
  new_cast_effect(
    effects = effects,
    conf_level = conf_level,
    alpha = screen$diagnostics$alpha %||% 0.05,
    diagnostics = diagnostics
  )
}


#' @keywords internal
#' @noRd
.cast_predict_mean <- function(fit, X_raw, models, chunk_size = 200000L) {
  n <- nrow(X_raw)
  predict_block <- function(Xi) {
    preds <- lapply(models, function(m) {
      tryCatch(predict_single_model(fit$models[[m]], Xi),
               error = function(e) rep(NA_real_, nrow(Xi)))
    })
    rowMeans(do.call(cbind, preds), na.rm = TRUE)
  }
  # Predict in row chunks. A single full-grid predict() (especially mgcv GAM,
  # which allocates an n-by-ncoef basis matrix) can need many GB on national
  # 1km grids (~10M cells) and segfault the process. Chunking bounds peak
  # memory to one block at a time; results are identical (row-independent).
  if (n <= chunk_size) return(predict_block(X_raw))
  out <- numeric(n)
  for (s in seq.int(1L, n, by = chunk_size)) {
    idx <- seq.int(s, min(s + chunk_size - 1L, n))
    out[idx] <- predict_block(X_raw[idx, , drop = FALSE])
  }
  out
}

#' Counterfactual What-If Map on the Current Climate
#'
#' Performs a g-computation intervention: one predictor is shifted while every
#' other predictor is held at its observed value, and the change in predicted
#' habitat suitability is mapped cell by cell. This isolates the modelled effect
#' of a single driver on today's landscape - a clean, purely interpretive
#' "what-if", with no future-climate extrapolation.
#'
#' @section Interpretation:
#' The result is a *model-based* what-if conditional on the fitted models and
#' the observed predictor distribution, reported on the relative habitat
#' suitability scale (not calibrated occurrence probability). It is not a
#' validated causal effect: it inherits the same no-unobserved-confounding and
#' correct-functional-form assumptions as the screen (Byrnes & Dee 2025), and
#' large shifts can push cells outside the training envelope - pair it with the
#' MESS/extrapolation flags from [cast_predict()].
#'
#' @param fit A `cast_fit` object.
#' @param newdata A `data.frame` with the fitted predictors and coordinate
#'   columns (`coords`). Typically the prediction grid on the current climate.
#' @param variable Single predictor name to intervene on (must be a fitted
#'   variable in `fit`).
#' @param shift Size of the intervention. Default `1`.
#' @param shift_type How `shift` is interpreted: `"sd"` (default, in predictor
#'   standard deviations), `"raw"` (predictor units), or `"percent"` (percent of
#'   each cell's value).
#' @param model Model name(s) whose predictions are averaged. Default `NULL`
#'   uses every successfully fitted model (an ensemble mean).
#' @param coords Coordinate column names. Default `c("lon", "lat")`.
#'
#' @return A `cast_counterfactual` object.
#' @seealso [cast_effect()], [cast_predict()]
#' @export
cast_counterfactual <- function(fit, newdata, variable,
                                shift = 1, shift_type = c("sd", "raw", "percent"),
                                model = NULL, coords = c("lon", "lat")) {
  if (!inherits(fit, "cast_fit")) cli::cli_abort("{.arg fit} must be a {.cls cast_fit}.")
  shift_type <- match.arg(shift_type)
  env_vars <- fit$env_vars
  if (length(variable) != 1L || !variable %in% env_vars) {
    cli::cli_abort(c(
      "{.arg variable} must be one fitted predictor.",
      i = "Available: {.val {env_vars}}."
    ))
  }
  if (!all(coords %in% names(newdata))) {
    cli::cli_abort("{.arg newdata} must contain coordinate columns {.val {coords}}.")
  }
  missing_vars <- setdiff(env_vars, names(newdata))
  if (length(missing_vars)) {
    cli::cli_abort("{.arg newdata} is missing fitted predictor{?s}: {.val {missing_vars}}.")
  }

  models <- model %||% names(fit$models)
  models <- models[vapply(models, function(m) {
    !is.null(fit$models[[m]]) && !is.null(fit$models[[m]]$model)
  }, logical(1))]
  if (!length(models)) cli::cli_abort("No usable fitted model available for prediction.")

  X_base <- as.data.frame(newdata[, env_vars, drop = FALSE], check.names = FALSE)
  for (col in names(X_base)) X_base[[col]] <- as.numeric(X_base[[col]])
  X_base <- .cast_impute(X_base, fit$scaling$impute)

  delta <- switch(shift_type,
    sd      = shift * (fit$scaling$sds[[variable]] %||% stats::sd(X_base[[variable]])),
    raw     = shift,
    percent = shift / 100 * X_base[[variable]]
  )
  X_cf <- X_base
  X_cf[[variable]] <- X_base[[variable]] + delta

  base_pred <- .cast_predict_mean(fit, X_base, models)
  cf_pred   <- .cast_predict_mean(fit, X_cf, models)

  predictions <- data.frame(
    lon = newdata[[coords[1]]],
    lat = newdata[[coords[2]]],
    baseline = base_pred,
    counterfactual = cf_pred,
    delta_hss = cf_pred - base_pred,
    stringsAsFactors = FALSE
  )

  d <- predictions$delta_hss[is.finite(predictions$delta_hss)]
  summary <- list(
    mean_delta = mean(d),
    median_delta = stats::median(d),
    frac_positive = mean(d > 0),
    max_gain = if (length(d)) max(d) else NA_real_,
    max_loss = if (length(d)) min(d) else NA_real_
  )

  new_cast_counterfactual(
    predictions = predictions,
    variable = variable,
    shift = shift,
    shift_type = shift_type,
    models = models,
    summary = summary
  )
}
