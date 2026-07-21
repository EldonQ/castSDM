# S3 Class Constructors -----------------------------------------------------

#' Create a cast_select Object
#'
#' @param selected Character vector of selected variable names.
#' @param scores A `data.frame` with per-variable scores. For `method = "dml"`:
#'   the orthogonalized effect estimate, standard error, statistic, and
#'   FDR-adjusted p-value; for `method = "rf"`: permutation importance.
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

#' Create a cast_effect Object
#'
#' @param effects A `data.frame` of per-predictor causal effects with
#'   confidence intervals and FDR-adjusted significance.
#' @param conf_level Confidence level used for the intervals.
#' @param alpha FDR level used to flag significance.
#' @param diagnostics Named list carried over from the DML screen.
#'
#' @return A `cast_effect` object.
#' @keywords internal
#' @export
new_cast_effect <- function(effects, conf_level = 0.95, alpha = 0.05,
                            diagnostics = list()) {
  structure(
    list(
      effects = effects,
      conf_level = conf_level,
      alpha = alpha,
      diagnostics = diagnostics
    ),
    class = "cast_effect"
  )
}

#' Create a cast_counterfactual Object
#'
#' @param predictions A `data.frame` with `lon`, `lat`, `baseline`,
#'   `counterfactual`, and `delta_hss`.
#' @param variable Character. The intervened predictor.
#' @param shift Numeric. Intervention size.
#' @param shift_type Character. How `shift` is interpreted.
#' @param models Character vector of models averaged for prediction.
#' @param summary Named list of change summaries.
#'
#' @return A `cast_counterfactual` object.
#' @keywords internal
#' @export
new_cast_counterfactual <- function(predictions, variable, shift, shift_type,
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
    class = "cast_counterfactual"
  )
}

#' Create a cast_fit Object
#'
#' @param models Named list of fitted model objects.
#' @param cast_vars Character vector of variables used for modeling.
#' @param env_vars Character vector of all environmental variable names.
#' @param scaling List with `means` and `sds` used for standardization.
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
#'
#' @return A `cast_cv` object.
#' @keywords internal
#' @export
new_cast_cv <- function(metrics, fold_metrics, folds,
                        k, block_method, thresholds,
                        selections = list(), screens = list()) {
  structure(
    list(
      metrics      = metrics,
      fold_metrics = fold_metrics,
      folds        = folds,
      k            = k,
      block_method = block_method,
      thresholds   = thresholds,
      selections   = selections,
      screens      = screens
    ),
    class = "cast_cv"
  )
}

#' Create a cast_predict Object
#'
#' @param predictions A `data.frame` with columns: `lon`, `lat`, and one
#'   `hss_*` column per model.
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
#' @param call The original function call.
#'
#' @return A `cast_result` object.
#' @keywords internal
#' @export
new_cast_result <- function(screen, fit, eval,
                            cv = NULL, predict = NULL,
                            ensemble = NULL,
                            call = NULL) {
  structure(
    list(
      screen = screen,
      fit = fit,
      eval = eval,
      cv = cv,
      predict = predict,
      ensemble = ensemble,
      call = call
    ),
    class = "cast_result"
  )
}

#' Create a cast_batch Object
#'
#' @param species_metrics A `data.frame` with per-species per-model
#'   evaluation metrics.
#' @param species Character vector of species names.
#' @param models Character vector of model names.
#' @param results Named list of per-species pipeline results (optional).
#' @param output_dir Character. Output directory path.
#'
#' @return A `cast_batch` object.
#' @keywords internal
#' @export
new_cast_batch <- function(species_metrics, species, models,
                           results = NULL, output_dir = NULL) {
  structure(
    list(
      species_metrics = species_metrics,
      species = species,
      models = models,
      results = results,
      output_dir = output_dir
    ),
    class = "cast_batch"
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
