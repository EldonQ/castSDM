# S3 Class Constructors -----------------------------------------------------

#' Create a cast_select Object
#'
#' @param selected Character vector of selected variable names.
#' @param scores A `data.frame` with per-variable scores. For
#'   `method = "two_stage"`: marginal `assoc`, the stage-1
#'   `collinear_thinned` flag, `interventional_effect` (the stage-2 selection
#'   statistic), the `perm_importance` diagnostic, the permutation `p_value`,
#'   the BH-adjusted `p_adjusted`, and the `selected` flag. For
#'   `method = "full"` the score columns are `NA`.
#' @param method Character screening-method identifier.
#' @param diagnostics Named list of method diagnostics.
#'
#' @return A `cast_select` object.
#' @keywords internal
#' @export
new_cast_select <- function(selected, scores, method = NULL, diagnostics = list()) {
  structure(
    list(
      selected = selected,
      scores = scores,
      method = method,
      diagnostics = diagnostics
    ),
    class = "cast_select"
  )
}

#' Create a cast_importance Object
#'
#' @param effects A `data.frame` of per-predictor interventional effect and
#'   permutation importance with permutation-null `p_value` and BH-adjusted
#'   `p_adjusted`.
#' @param alpha Significance level used to flag predictors.
#' @param threshold Numeric. The stage-2 permutation-null threshold, or `NA`
#'   when unknown (one entry per predictor for a two-stage screen).
#' @param diagnostics Named list carried over from the screen.
#'
#' @return A `cast_importance` object.
#' @keywords internal
#' @export
new_cast_importance <- function(effects, alpha = 0.05,
                                threshold = NA_real_, diagnostics = list()) {
  structure(
    list(
      effects = effects,
      alpha = alpha,
      threshold = threshold,
      diagnostics = diagnostics
    ),
    class = "cast_importance"
  )
}

#' Create a cast_sensitivity Object
#'
#' @param predictions A `data.frame` with `lon`, `lat`, `baseline`,
#'   `counterfactual`, and `delta_hss`.
#' @param variable Character. The intervened predictor.
#' @param shift Numeric. Intervention size.
#' @param shift_type Character. How `shift` is interpreted.
#' @param models Character vector of models averaged for prediction.
#' @param summary Named list of change summaries.
#'
#' @return A `cast_sensitivity` object.
#' @keywords internal
#' @export
new_cast_sensitivity <- function(predictions, variable, shift, shift_type,
                                    models, summary = list()) {
  structure(
    list(
      predictions = predictions,
      variable = variable,
      shift = shift,
      shift_type = shift_type,
      models = models,
      summary = summary
    ),
    class = "cast_sensitivity"
  )
}

#' Create a cast_dose_response Object
#'
#' @param curve A `data.frame` with `shift`, `shift_raw`, `mean_delta`,
#'   `mean_abs_delta`, `support`, and `range_supported`.
#' @param variable Character. The intervened predictor.
#' @param shift_type Character. How `shift` is interpreted.
#' @param unit Character. Human-readable shift unit.
#' @param models Character vector of models averaged.
#' @param assumptions Character. The identification assumptions carried by the
#'   estimate.
#'
#' @return A `cast_dose_response` object.
#' @keywords internal
#' @export
new_cast_dose_response <- function(curve, variable, shift_type, unit, models,
                                   assumptions = NULL) {
  structure(
    list(curve = curve, variable = variable, shift_type = shift_type,
         unit = unit, models = models, assumptions = assumptions),
    class = "cast_dose_response"
  )
}

#' Create a cast_support Object
#'
#' @param support A `data.frame` with `driver`, `shift`, `shift_raw`, and
#'   `support` (the fraction of evaluated rows inside the training quantile
#'   box both before and after the shift).
#' @param variables Character vector of assessed predictors.
#' @param shift Numeric vector of assessed shift sizes.
#' @param shift_type Character. How `shift` is interpreted.
#' @param support_probs Numeric length 2. Quantiles defining the support box.
#' @param models Character vector of models in the fit.
#' @param assumptions Character. The identification assumptions carried.
#'
#' @return A `cast_support` object.
#' @keywords internal
#' @export
new_cast_support <- function(support, variables, shift, shift_type,
                             support_probs, models, assumptions = NULL) {
  structure(
    list(support = support, variables = variables, shift = shift,
         shift_type = shift_type, support_probs = support_probs,
         models = models, assumptions = assumptions),
    class = "cast_support"
  )
}

#' Create a cast_necessity Object
#'
#' @param necessity A `data.frame` with one row per driver: `mean_dAUC`,
#'   `sd_dAUC`, `min_dAUC`, `max_dAUC`, `pct_folds_positive`, `n_folds`,
#'   `necessary`.
#' @param fold_dauc Numeric matrix of per-driver (rows) by per-fold
#'   (columns) held-out AUC loss.
#' @param auc_full Numeric vector of full-model held-out AUC per fold.
#' @param folds Integer vector. Spatial fold assignment for each data row.
#' @param k Integer. Number of folds actually used.
#' @param block_method Character. Spatial blocking strategy used.
#' @param diagnostics Named list of diagnostics.
#'
#' @return A `cast_necessity` object.
#' @keywords internal
#' @export
new_cast_necessity <- function(necessity, fold_dauc, auc_full, folds, k,
                               block_method, diagnostics = list()) {
  structure(
    list(
      necessity = necessity,
      fold_dauc = fold_dauc,
      auc_full = auc_full,
      folds = folds,
      k = k,
      block_method = block_method,
      diagnostics = diagnostics
    ),
    class = "cast_necessity"
  )
}

#' Create a cast_fit Object
#'
#' @param models Named list of fitted model objects.
#' @param cast_vars Character vector of variables used for modeling.
#' @param env_vars Character vector of all environmental variable names.
#' @param scaling List of training-set predictor statistics reused across the
#'   prediction stack: `means` and `sds` (for sensitivity SD-based shifts in
#'   [cast_sensitivity()]), `impute` (per-predictor training median used by
#'   the internal `.cast_impute()` helper), and `reference` (the imputed
#'   training predictor frame for MESS/clamp extrapolation control). Models are
#'   trained on the raw predictors; `means`/`sds` are not applied to fitting
#'   inputs.
#' @param screen A `cast_select` object (or `NULL`).
#'
#' @return A `cast_fit` object.
#' @keywords internal
#' @export
new_cast_fit <- function(models, cast_vars, env_vars, scaling, screen = NULL) {
  structure(
    list(
      models = models,
      cast_vars = cast_vars,
      env_vars = env_vars,
      scaling = scaling,
      screen = screen
    ),
    class = "cast_fit"
  )
}

#' Create a cast_eval Object
#'
#' @param metrics A `data.frame` with per-model evaluation metrics.
#'   Columns: `model`, `auc_mean`, `tss_mean`, `cbi_mean`.
#' @param cv_source Logical. Whether metrics came from spatial CV.
#'   Default `FALSE`.
#'
#' @return A `cast_eval` object.
#' @keywords internal
#' @export
new_cast_eval <- function(metrics, cv_source = FALSE) {
  structure(
    list(metrics = metrics, cv_source = cv_source),
    class = "cast_eval"
  )
}

#' Create a cast_cv Object
#'
#' @param metrics `data.frame`. Aggregated per-model metrics (mean +/- SD).
#' @param fold_metrics `data.frame`. Per-fold per-model raw metrics.
#' @param folds Integer vector. Fold assignment for each data row.
#' @param k Integer. Number of folds.
#' @param block_method Character. Blocking strategy used.
#' @param thresholds Named numeric. TSS-optimal threshold per model.
#' @param selections List of selected variables for each outer fold.
#' @param screens List of fold-specific `cast_select` objects.
#' @param selection_freq A `data.frame` with each predictor's fold-level
#'   selection frequency (`variable`, `freq`), sorted descending.
#' @param oof A `data.frame` of out-of-fold predictions with an `obs` column
#'   and one `HSS_<model>` column per model, or `NULL`. This is the labelled
#'   surface [cast_ensemble()] thresholds on.
#' @param fold_status Character vector of evaluation, empty-selection or failure
#'   statuses, one per scheduled fold.
#'
#' @return A `cast_cv` object.
#' @keywords internal
#' @export
new_cast_cv <- function(metrics, fold_metrics, folds,
                        k, block_method, thresholds,
                        selections = list(), screens = list(),
                        selection_freq = NULL, oof = NULL, fold_status = NULL) {
  structure(
    list(
      metrics      = metrics,
      fold_metrics = fold_metrics,
      folds        = folds,
      k            = k,
      block_method = block_method,
      thresholds   = thresholds,
      selections   = selections,
      screens      = screens,
      selection_freq = selection_freq,
      oof          = oof,
      fold_status  = fold_status
    ),
    class = "cast_cv"
  )
}

#' Create a cast_predict Object
#'
#' @param predictions A `data.frame` with columns: `lon`, `lat`, and one
#'   `HSS_*` column per model.
#' @param models Character vector of model names included.
#'
#' @return A `cast_predict` object.
#' @keywords internal
#' @export
new_cast_predict <- function(predictions, models) {
  structure(
    list(
      predictions = predictions,
      models = models
    ),
    class = "cast_predict"
  )
}

#' Create a cast_result Object
#'
#' Container for the full pipeline output.
#'
#' @param screen A `cast_select` object.
#' @param fit A `cast_fit` object.
#' @param eval A `cast_eval` object (hold-out evaluation).
#' @param cv A `cast_cv` object (spatial CV), or `NULL`.
#' @param predict A `cast_predict` object (or `NULL`).
#' @param ensemble A `cast_ensemble` object (or `NULL`).
#' @param fit_full A `cast_fit` object refitted on the full data set for
#'   final prediction maps (or `NULL`).
#' @param call The original function call.
#'
#' @return A `cast_result` object.
#' @keywords internal
#' @export
new_cast_result <- function(screen, fit, eval,
                            cv = NULL, predict = NULL,
                            ensemble = NULL,
                            fit_full = NULL,
                            call = NULL) {
  structure(
    list(
      screen = screen,
      fit = fit,
      eval = eval,
      cv = cv,
      predict = predict,
      ensemble = ensemble,
      fit_full = fit_full,
      call = call
    ),
    class = "cast_result"
  )
}

#' Create a cast_ensemble Object
#'
#' @param predictions A `data.frame` with columns `lon`, `lat`, `hss_ensemble`
#'   and optionally `binary_ensemble`.
#' @param weights Named numeric vector of per-model weights.
#' @param method Character. Ensemble method used (`"weighted"`, `"best"`,
#'   `"equal"`).
#' @param threshold Numeric. Binary classification threshold.
#' @param model_scores Named numeric vector of per-model composite scores.
#'
#' @return A `cast_ensemble` object.
#' @keywords internal
#' @export
new_cast_ensemble <- function(predictions, weights, method,
                              threshold, model_scores) {
  structure(
    list(
      predictions = predictions,
      weights = weights,
      method = method,
      threshold = threshold,
      model_scores = model_scores
    ),
    class = "cast_ensemble"
  )
}

#' Create a cast_project Object
#'
#' @param current A `cast_ensemble` object for the current climate.
#' @param future A named list of `cast_ensemble` objects for future scenarios.
#' @param changes A named list of `data.frame`s with columns `lon`, `lat`,
#'   `change` (`"gain"`, `"loss"`, `"stable_present"`, `"stable_absent"`).
#' @param stats A `data.frame` with summary statistics per scenario.
#'
#' @return A `cast_project` object.
#' @keywords internal
#' @export
new_cast_project <- function(current, future, changes, stats) {
  structure(
    list(
      current = current,
      future = future,
      changes = changes,
      stats = stats
    ),
    class = "cast_project"
  )
}
