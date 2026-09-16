# Importance reporting + sensitivity -----------------------------------------
#
# Two functions turn a fitted castSDM workflow into interpretable evidence:
#   * cast_importance() - tidies the screen's two attribution columns: the
#                         interventional effect that selected each predictor
#                         and the permutation importance kept as a diagnostic,
#                         with permutation p-values and null thresholds.
#   * cast_sensitivity() - a "what-if" on the current climate: shift one
#                          predictor, hold the rest fixed, and map the change in
#                          predicted habitat suitability.
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

#' Predictor Importance Table from the Screen
#'
#' Turns a two-stage screen into a tidy per-predictor table carrying both
#' attribution vocabularies: the \strong{interventional effect} that selected
#' the predictor (stage 2 of [cast_select()]), and the random-forest
#' \strong{permutation importance} retained as a diagnostic, each with its
#' permutation-null p-value and BH-adjusted p-value.
#'
#' @section Why both columns are reported:
#' The two columns answer different questions. Permutation importance permutes
#' a predictor, which breaks its correlation with every other predictor and
#' pushes collinear rows outside the observed data support, so the score is
#' governed by the model's extrapolation behaviour (Hooker, Mentch & Zhou
#' 2021). The interventional effect shifts one predictor and holds the rest at
#' their observed values, so it stays inside the support. A predictor with high
#' permutation importance but a low interventional effect is one the forest
#' leans on but barely responds to when changed - the signature of a collinear
#' stand-in. [cast_select()] reports their Spearman agreement as
#' `screen$diagnostics$importance_agreement`.
#'
#' Importance carries no sign, so read it together with
#' [cast_effect_table()] for direction.
#'
#' @section Interpretation (read before citing):
#' Both columns describe the \emph{fitted model}. Neither is a causal effect:
#' a causal reading additionally requires no unobserved confounding, a
#' correctly specified adjustment set, and no reverse causation, which
#' observational, sampling-biased SDM data rarely satisfy (Byrnes & Dee 2025).
#' Pair with [cast_necessity()].
#'
#' @param object A `cast_select` from `method = "two_stage"`, or a
#'   `cast_fit` / `cast_result` that carries such a screen.
#'
#' @return A `cast_importance` object.
#' @references
#' Altmann, A., Toloşi, L., Sander, O. & Lengauer, T. (2010). Permutation
#' importance: a corrected feature importance measure. *Bioinformatics*,
#' 26(10), 1340-1347.
#'
#' Byrnes, J. E. K. & Dee, L. E. (2025). Causal inference with observational
#' data and unobserved confounding variables. *Ecology Letters*, 28(1), e70023.
#'
#' Hooker, G., Mentch, L. & Zhou, S. (2021). Unrestricted permutation forces
#' extrapolation: variable importance requires at least one more model, or
#' there is no free variable importance. *WIREs Data Mining and Knowledge
#' Discovery*, 11, e1421.
#' @seealso [cast_select()], [cast_sensitivity()], [cast_necessity()]
#' @export
cast_importance <- function(object) {
  screen <- .cast_extract_screen(object)
  sc <- screen$scores
  needed <- c("interventional_effect", "perm_importance", "p_value", "p_adjusted")
  if (!all(needed %in% names(sc))) {
    cli::cli_abort(c(
      "{.fn cast_importance} needs a two-stage screen carrying a stage-2 statistic.",
      i = "Run {.code cast_select(..., method = \"two_stage\")} first."))
  }
  sc <- sc[is.finite(sc$interventional_effect) | is.finite(sc$perm_importance), ,
           drop = FALSE]
  if (!nrow(sc)) {
    cli::cli_abort(c(
      "The screen holds no finite importance estimates.",
      i = "{.code method = \"full\"} skips stage 2, so there is nothing to report."))
  }
  effects <- data.frame(
    variable              = sc$variable,
    interventional_effect = sc$interventional_effect,
    perm_importance       = sc$perm_importance,
    null_threshold        = sc$null_threshold,
    p_value               = sc$p_value,
    p_adjusted            = sc$p_adjusted,
    selected              = sc$selected,
    stringsAsFactors = FALSE
  )
  effects <- effects[order(-effects$interventional_effect), , drop = FALSE]
  rownames(effects) <- NULL

  diagnostics <- screen$diagnostics
  diagnostics$measure <- "interventional_effect"
  new_cast_importance(
    effects = effects,
    alpha = screen$diagnostics$alpha %||% 0.05,
    threshold = screen$diagnostics$null_threshold %||% NA_real_,
    diagnostics = diagnostics
  )
}


#' @keywords internal
#' @noRd
.cast_predict_matrix <- function(fit, X_raw, models, chunk_size = 200000L) {
  n <- nrow(X_raw)
  predict_block <- function(Xi) {
    preds <- lapply(models, function(m) {
      tryCatch(predict_single_model(fit$models[[m]], Xi),
               error = function(e) rep(NA_real_, nrow(Xi)))
    })
    do.call(cbind, preds)
  }
  # Predict in row chunks. A single full-grid predict() (especially mgcv GAM,
  # which allocates an n-by-ncoef basis matrix) can need many GB on national
  # 1km grids (~10M cells) and segfault the process. Chunking bounds peak
  # memory to one block at a time; results are identical (row-independent).
  if (n <= chunk_size) return(predict_block(X_raw))
  out <- matrix(NA_real_, nrow = n, ncol = length(models))
  for (s in seq.int(1L, n, by = chunk_size)) {
    idx <- seq.int(s, min(s + chunk_size - 1L, n))
    out[idx, ] <- predict_block(X_raw[idx, , drop = FALSE])
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
#' validated causal effect: causal interpretation requires a justified
#' adjustment set, no unobserved confounding, consistency and joint support; and
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
#' @return A `cast_sensitivity` object. The `predictions` data.frame
#'   carries `baseline`, `counterfactual`, `delta_hss` (ensemble-mean change)
#'   and `delta_sd` (cross-model standard deviation of the change; `NA` for
#'   a single model).
#' @seealso [cast_importance()], [cast_predict()]
#' @export
cast_sensitivity <- function(fit, newdata, variable,
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

  base_mat <- .cast_predict_matrix(fit, X_base, models)
  cf_mat   <- .cast_predict_matrix(fit, X_cf, models)
  base_pred <- rowMeans(base_mat, na.rm = TRUE)
  cf_pred   <- rowMeans(cf_mat, na.rm = TRUE)
  delta_mat <- cf_mat - base_mat
  delta_mat[!is.finite(cf_mat) | !is.finite(base_mat)] <- NA_real_
  delta_sd <- if (ncol(delta_mat) > 1L) {
    apply(delta_mat, 1L, stats::sd, na.rm = TRUE)
  } else {
    rep(NA_real_, nrow(delta_mat))
  }

  predictions <- data.frame(
    lon = newdata[[coords[1]]],
    lat = newdata[[coords[2]]],
    baseline = base_pred,
    counterfactual = cf_pred,
    delta_hss = rowMeans(delta_mat, na.rm = TRUE),
    delta_sd = delta_sd,
    stringsAsFactors = FALSE
  )

  d <- predictions$delta_hss[is.finite(predictions$delta_hss)]
  if (!length(d)) {
    # Every model failed on base or counterfactual: refuse to return a
    # silent NaN summary.
    cli::cli_abort(
      "All model predictions failed for the {.val {variable}} counterfactual; no finite deltas to summarise."
    )
  }
  summary <- list(
    mean_delta = mean(d),
    median_delta = stats::median(d),
    frac_positive = mean(d > 0),
    max_gain = if (length(d)) max(d) else NA_real_,
    max_loss = if (length(d)) min(d) else NA_real_,
    mean_abs_delta_sd = mean(abs(predictions$delta_sd), na.rm = TRUE)
  )

  new_cast_sensitivity(
    predictions = predictions,
    variable = variable,
    shift = shift,
    shift_type = shift_type,
    models = models,
    summary = summary
  )
}
