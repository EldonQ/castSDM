# S3 Class Constructors -------------------------------------------------------

#' Create a cast_dag Object
#'
#' @param edges A `data.frame` with columns: `from`, `to`, `strength`,
#'   `direction`.
#' @param nodes Character vector of node (variable) names.
#' @param boot_R Integer. Number of bootstrap replicates used.
#' @param strength_threshold Numeric. Edge strength threshold applied.
#' @param direction_threshold Numeric. Direction consistency threshold applied.
#' @param score Character. Scoring criterion used (e.g., `"bic-cg"`).
#' @param structure_method Character. Structure learning method used.
#' @param response_node Character or `NULL`. Name of the response node if
#'   included in DAG learning.
#'
#' @return A `cast_dag` object.
#' @keywords internal
#' @export
new_cast_dag <- function(edges, nodes, boot_R, strength_threshold,
                         direction_threshold, score = "bic-cg",
                         structure_method = "pc",
                         response_node = NULL,
                         metadata = list()) {
  structure(
    list(
      edges = edges,
      nodes = nodes,
      boot_R = boot_R,
      strength_threshold = strength_threshold,
      direction_threshold = direction_threshold,
      score = score,
      structure_method = structure_method,
      response_node = response_node,
      metadata = metadata %||% list()
    ),
    class = "cast_dag"
  )
}

#' Create a cast_select Object
#'
#' @param selected Character vector of selected variable names.
#' @param scores A `data.frame` with per-variable scores (RF importance and
#'   DAG Markov Blanket membership).
#' @param roles A `data.frame` with columns `variable` and `role` indicating
#'   each variable's causal role: `"parent"`, `"child"`, `"co_parent"`, or
#'   `"predictive"`.
#'
#' @return A `cast_select` object.
#' @keywords internal
#' @export
new_cast_select <- function(selected, scores, roles) {
  structure(
    list(
      selected = selected,
      scores = scores,
      roles = roles
    ),
    class = "cast_select"
  )
}

#' Create a cast_fit Object
#'
#' @param models Named list of fitted model objects.
#' @param cast_vars Character vector of variables used for modeling.
#' @param env_vars Character vector of all environmental variable names.
#' @param scaling List with `means` and `sds` used for standardization.
#' @param dag A `cast_dag` object (or `NULL`).
#' @param screen A `cast_select` object (or `NULL`).
#'
#' @return A `cast_fit` object.
#' @keywords internal
#' @export
new_cast_fit <- function(models, cast_vars, env_vars, scaling,
                         dag = NULL, screen = NULL) {
  structure(
    list(
      models = models,
      cast_vars = cast_vars,
      env_vars = env_vars,
      scaling = scaling,
      dag = dag,
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
#'
#' @return A `cast_cv` object.
#' @keywords internal
#' @export
new_cast_cv <- function(metrics, fold_metrics, folds,
                        k, block_method, thresholds) {
  structure(
    list(
      metrics      = metrics,
      fold_metrics = fold_metrics,
      folds        = folds,
      k            = k,
      block_method = block_method,
      thresholds   = thresholds
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

#' Create a cast_cate Object
#'
#' @param effects A `data.frame` with columns: `lon`, `lat`, `variable`,
#'   `cate`.
#' @param variables Character vector of variables for which CATE was estimated.
#' @param n_trees Integer. Number of causal forest trees used.
#'
#' @return A `cast_cate` object.
#' @keywords internal
#' @export
new_cast_cate <- function(effects, variables, n_trees = 1000L) {
  structure(
    list(
      effects = effects,
      variables = variables,
      n_trees = n_trees
    ),
    class = "cast_cate"
  )
}

#' Create a cast_result Object
#'
#' Container for the full pipeline output.
#'
#' @param dag A `cast_dag` object.
#' @param screen A `cast_select` object.
#' @param fit A `cast_fit` object.
#' @param eval A `cast_eval` object (hold-out evaluation).
#' @param cv A `cast_cv` object (spatial CV), or `NULL`.
#' @param predict A `cast_predict` object (or `NULL`).
#' @param ensemble A `cast_ensemble` object (or `NULL`).
#' @param cate A `cast_cate` object (or `NULL`).
#' @param call The original function call.
#'
#' @return A `cast_result` object.
#' @keywords internal
#' @export
new_cast_result <- function(dag, screen, fit, eval,
                            cv = NULL, predict = NULL,
                            ensemble = NULL,
                            cate = NULL,
                            call = NULL) {
  structure(
    list(
      dag = dag,
      screen = screen,
      fit = fit,
      eval = eval,
      cv = cv,
      predict = predict,
      ensemble = ensemble,
      cate = cate,
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
