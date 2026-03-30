#' Learn Causal DAG via Bootstrap Hill-Climbing
#'
#' Discovers species-specific causal structure among environmental variables
#' and species presence using bootstrap-aggregated Hill-Climbing with the
#' BIC-Gaussian score. Edges are retained based on strength and direction
#' consistency thresholds.
#'
#' @param data A `data.frame` containing the `presence` column and
#'   environmental variables. Must be numeric (no factors).
#' @param response Character. Name of the response column. Default
#'   `"presence"`.
#' @param R Integer. Number of bootstrap replicates. Default `100`.
#' @param algorithm Character. Structure learning algorithm. Default `"hc"`
#'   (Hill-Climbing). See [bnlearn::boot.strength()] for options.
#' @param score Character. Scoring criterion. Default `"bic-g"` (BIC for
#'   Gaussian data).
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
#' The DAG is learned **with the presence variable included** in the network.
#' This ensures the learned structure reflects how environmental variables
#' influence species occurrence, not just inter-variable correlations.
#'
#' The bootstrap procedure (Friedman et al. 1999) provides confidence measures
#' for each edge:
#' - **Strength**: Proportion of bootstrap replicates containing the edge.
#' - **Direction**: Proportion of times the edge points from A -> B
#'   (vs B -> A).
#'
#' After learning, edges connecting to the `presence` node are excluded from
#' the returned DAG, so that only environmental causal structure is retained
#' for downstream feature engineering.
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
                     response = "presence",
                     R = 100L,
                     algorithm = "hc",
                     score = "bic-g",
                     strength_threshold = 0.7,
                     direction_threshold = 0.6,
                     max_rows = 8000L,
                     seed = NULL,
                     verbose = TRUE) {
  check_suggested("bnlearn", "for DAG structure learning")

  env_vars <- get_env_vars(data, response = response)
  if (length(env_vars) < 3) {
    cli::cli_abort("Need at least 3 environmental variables for DAG learning.")
  }

  # Build DAG data: env vars + presence, all numeric
  dag_cols <- c(env_vars, response)
  dag_df <- as.data.frame(data[, dag_cols, drop = FALSE])
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
      "Learning DAG: {length(env_vars)} vars, {nrow(dag_df)} obs, R={R}..."
    )
  }

  if (!is.null(seed)) set.seed(seed)
  boot_str <- bnlearn::boot.strength(
    dag_df, R = R, algorithm = algorithm,
    algorithm.args = list(score = score)
  )

  # Filter by thresholds
  strong <- boot_str[
    boot_str$strength >= strength_threshold &
      boot_str$direction >= direction_threshold, ,
    drop = FALSE
  ]

  # Exclude edges involving the response variable
  env_edges <- strong[
    strong$from != response & strong$to != response, ,
    drop = FALSE
  ]
  env_edges <- as.data.frame(env_edges)
  rownames(env_edges) <- NULL

  # Build nodes list (only env vars that appear in edges)
  nodes <- unique(c(env_edges$from, env_edges$to))
  nodes <- union(nodes, env_vars) # include all env vars as potential nodes

  if (verbose) {
    cli::cli_inform(
      "DAG: {nrow(env_edges)} edges retained (strength >= {strength_threshold})."
    )
  }

  new_cast_dag(
    edges = env_edges,
    nodes = env_vars,
    boot_R = R,
    strength_threshold = strength_threshold,
    direction_threshold = direction_threshold,
    score = score
  )
}
