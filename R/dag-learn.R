#' Learn Causal DAG via Bootstrap Hill-Climbing
#'
#' Discovers the causal structure among environmental variables using
#' bootstrap-aggregated Hill-Climbing (HC). Edges are retained based on
#' bootstrap strength and direction consistency thresholds.
#'
#' @param data A `data.frame` containing the `presence` column and
#'   environmental variables. Must be numeric (no factors).
#' @param response Character. Name of the response column. Default
#'   `"presence"`.
#' @param R Integer. Number of bootstrap replicates. Default `100`.
#' @param algorithm Character. Structure learning algorithm. Default `"hc"`
#'   (Hill-Climbing). See [bnlearn::boot.strength()] for options.
#' @param env_vars Character vector or `NULL`. Environmental variable names
#'   to include in DAG learning. When `NULL` (default), all numeric columns
#'   in `data` excluding `response`, `lon`, and `lat` are used. **Always pass
#'   this explicitly** when `data` contains columns beyond the intended
#'   predictors (e.g., after [cast_prepare()], which preserves all original
#'   columns). Typical usage: `env_vars = split$env_vars`.
#' @param score Character. Scoring criterion. Default `"bic-g"` (BIC for
#'   Gaussian data). Note: the DAG is learned on **environmental variables
#'   only** (presence excluded), so `bic-g` is appropriate for the continuous
#'   predictors.
#' @param strength_threshold Numeric. Minimum edge strength (proportion of
#'   bootstrap replicates). Default `0.7`.
#' @param direction_threshold Numeric. Minimum direction consistency. Default
#'   `0.6`.
#' @param max_rows Integer. Maximum rows for DAG learning; subsample if
#'   exceeded. Default `8000`.
#' @param seed Integer or `NULL`. Random seed. Default `NULL`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A [cast_dag] object.
#'
#' @details
#' ## Why presence is excluded from DAG learning
#' The DAG is learned on **environmental variables only**, not on `presence`.
#' The original design included `presence` in the network but then discarded
#' all edges touching it. This was inconsistent: the `bic-g` (Gaussian BIC)
#' score treats `presence` as a continuous variable, corrupting edge weights
#' involving nearby nodes. Excluding `presence` from the start keeps the
#' Gaussian assumption valid for the continuous environmental predictors
#' and avoids spurious confounding of the env-env edge strengths.
#'
#' The bootstrap procedure (Friedman et al. 1999) provides confidence measures:
#' - **Strength**: Proportion of bootstrap replicates containing the edge.
#' - **Direction**: Proportion of times the edge points from A -> B.
#'
#' @references
#' Pearl, J. (2000). *Causality: Models, Reasoning, and Inference*.
#' Cambridge University Press.
#'
#' Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package.
#' *Journal of Statistical Software*, 35(3), 1-22.
#'
#' @seealso [bnlearn::boot.strength()], [cast_ate()], [cast_screen()]
#'
#' @export
cast_dag <- function(data,
                     response            = "presence",
                     env_vars            = NULL,
                     R                   = 100L,
                     algorithm           = "hc",
                     score               = "bic-g",
                     strength_threshold  = 0.7,
                     direction_threshold = 0.6,
                     max_rows            = 8000L,
                     seed                = NULL,
                     verbose             = TRUE) {
  check_suggested("bnlearn", "for DAG structure learning")

  env_vars <- env_vars %||% get_env_vars(data, response = response)
  if (length(env_vars) < 3) {
    cli::cli_abort("Need at least 3 environmental variables for DAG learning.")
  }

  # Learn DAG on environmental variables ONLY.
  # Including the binary `presence` column violates the bic-g Gaussian
  # assumption and distorts edge strengths among neighbouring nodes.
  dag_df <- as.data.frame(data[, env_vars, drop = FALSE])
  for (col in names(dag_df)) {
    dag_df[[col]] <- as.numeric(dag_df[[col]])
  }
  dag_df <- stats::na.omit(dag_df)

  if (nrow(dag_df) < 10) {
    cli::cli_abort(
      "Fewer than 10 complete cases for DAG learning ({nrow(dag_df)} available)."
    )
  }

  # Subsample if too large
  if (nrow(dag_df) > max_rows) {
    if (!is.null(seed)) set.seed(seed)
    dag_df <- dag_df[sample.int(nrow(dag_df), max_rows), ]
  }

  if (verbose) {
    cli::cli_inform(
      "Learning DAG: {length(env_vars)} env vars (presence excluded), {nrow(dag_df)} obs, R={R}..."
    )
  }

  if (!is.null(seed)) set.seed(seed)
  boot_str <- bnlearn::boot.strength(
    dag_df, R = R, algorithm = algorithm,
    algorithm.args = list(score = score)
  )

  # Filter by thresholds
  strong <- boot_str[
    boot_str$strength  >= strength_threshold &
    boot_str$direction >= direction_threshold, ,
    drop = FALSE
  ]
  env_edges <- as.data.frame(strong)
  rownames(env_edges) <- NULL

  if (verbose) {
    cli::cli_inform(
      "DAG: {nrow(env_edges)} edges retained (strength >= {strength_threshold})."
    )
  }

  new_cast_dag(
    edges               = env_edges,
    nodes               = env_vars,
    boot_R              = R,
    strength_threshold  = strength_threshold,
    direction_threshold = direction_threshold,
    score               = score
  )
}
