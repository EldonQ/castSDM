# S3 Class Constructors -------------------------------------------------------

#' Create a cast_dag Object
#'
#' @param edges A `data.frame` with columns: `from`, `to`, `strength`,
#'   `direction`.
#' @param nodes Character vector of node (variable) names.
#' @param boot_R Integer. Number of bootstrap replicates used.
#' @param strength_threshold Numeric. Edge strength threshold applied.
#' @param direction_threshold Numeric. Direction consistency threshold applied.
#' @param score Character. Scoring criterion used (e.g., `"bic-g"`).
#'
#' @return A `cast_dag` object.
#' @keywords internal
#' @export
new_cast_dag <- function(edges, nodes, boot_R, strength_threshold,
                         direction_threshold, score = "bic-g") {
  structure(
    list(
      edges = edges,
      nodes = nodes,
      boot_R = boot_R,
      strength_threshold = strength_threshold,
      direction_threshold = direction_threshold,
      score = score
    ),
    class = "cast_dag"
  )
}

#' Create a cast_ate Object
#'
#' @param estimates A `data.frame` with columns: `variable`, `coef`, `se`,
#'   `p_value`, `significant`.
#' @param K Integer. Number of cross-fitting folds used.
#' @param alpha Numeric. Significance level used.
#'
#' @return A `cast_ate` object.
#' @keywords internal
#' @export
new_cast_ate <- function(estimates, K, alpha = 0.05) {
  structure(
    list(
      estimates = estimates,
      K = K,
      alpha = alpha
    ),
    class = "cast_ate"
  )
}

#' Create a cast_screen Object
#'
#' @param selected Character vector of selected variable names.
#' @param scores A `data.frame` with per-variable composite scores.
#' @param weights Named numeric vector of adaptive weights (`w_dag`, `w_ate`,
#'   `w_imp`).
#'
#' @return A `cast_screen` object.
#' @keywords internal
#' @export
new_cast_screen <- function(selected, scores, weights) {
  structure(
    list(
      selected = selected,
      scores = scores,
      weights = weights
    ),
    class = "cast_screen"
  )
}

#' Create a cast_roles Object
#'
#' @param roles A `data.frame` with columns: `variable`, `role`, `in_degree`,
#'   `out_degree`, `role_score`.
#'
#' @return A `cast_roles` object.
#' @keywords internal
#' @export
new_cast_roles <- function(roles) {
  structure(
    list(roles = roles),
    class = "cast_roles"
  )
}

#' Create a cast_fit Object
#'
#' @param models Named list of fitted model objects.
#' @param cast_vars Character vector of variables used by CAST model.
#' @param env_vars Character vector of all environmental variable names.
#' @param scaling List with `means` and `sds` used for standardization.
#' @param dag A `cast_dag` object (or `NULL`).
#' @param ate A `cast_ate` object (or `NULL`).
#' @param screen A `cast_screen` object (or `NULL`).
#'
#' @return A `cast_fit` object.
#' @keywords internal
#' @export
new_cast_fit <- function(models, cast_vars, env_vars, scaling,
                         dag = NULL, ate = NULL, screen = NULL) {
  structure(
    list(
      models = models,
      cast_vars = cast_vars,
      env_vars = env_vars,
      scaling = scaling,
      dag = dag,
      ate = ate,
      screen = screen
    ),
    class = "cast_fit"
  )
}

#' Create a cast_eval Object
#'
#' @param metrics A `data.frame` with columns: `model`, `auc_mean`, `auc_sd`,
#'   `tss_mean`, `tss_sd`.
#'
#' @return A `cast_eval` object.
#' @keywords internal
#' @export
new_cast_eval <- function(metrics) {
  structure(
    list(metrics = metrics),
    class = "cast_eval"
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
#' @param ate A `cast_ate` object.
#' @param screen A `cast_screen` object.
#' @param roles A `cast_roles` object.
#' @param fit A `cast_fit` object.
#' @param eval A `cast_eval` object.
#' @param predict A `cast_predict` object (or `NULL`).
#' @param cate A `cast_cate` object (or `NULL`).
#' @param call The original function call.
#'
#' @return A `cast_result` object.
#' @keywords internal
#' @export
new_cast_result <- function(dag, ate, screen, roles, fit, eval,
                            predict = NULL, cate = NULL, call = NULL) {
  structure(
    list(
      dag = dag,
      ate = ate,
      screen = screen,
      roles = roles,
      fit = fit,
      eval = eval,
      predict = predict,
      cate = cate,
      call = call
    ),
    class = "cast_result"
  )
}
