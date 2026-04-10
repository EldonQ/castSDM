#' Causal Feature Engineering
#'
#' Constructs augmented feature matrices using ATE-weighted variables and
#' DAG-guided interaction terms. Three feature layers are produced:
#' 1. **ATE-weighted**: Base variables scaled by causal effect magnitude.
#' 2. **DAG interactions**: Pairwise products of the top-K strongest connected
#'    variable pairs, weighted by edge strength (CAST model only).
#' 3. **Raw**: Unmodified variables (for baseline models).
#'
#' @param data A `data.frame` of standardized environmental variables.
#' @param screen A [cast_screen] object.
#' @param dag A [cast_dag] object.
#' @param ate A [cast_ate] object.
#' @param model_type Character. Which feature set to build: `"cast"` (all
#'   layers), `"mlp_ate"` (ATE-weighted only), `"mlp"` (raw only). Default
#'   `"cast"`.
#' @param max_interactions Integer. Maximum number of DAG interaction features
#'   to add. Edges are ranked by bootstrap strength and only the top-K are
#'   used. This prevents dimension explosion when the DAG is dense. Default
#'   `15L`. Set to `Inf` to use all edges (original behaviour).
#'
#' @return A list with:
#' \describe{
#'   \item{`features`}{A `data.frame` of engineered features.}
#'   \item{`n_base`}{Number of base variables.}
#'   \item{`n_interactions`}{Number of DAG interaction features actually added.}
#'   \item{`n_total`}{Total number of features.}
#'   \item{`ate_weights`}{Named numeric vector of ATE weights applied.}
#'   \item{`interaction_names`}{Names of interaction columns (if any).}
#' }
#'
#' @seealso [cast_screen()], [cast_ate()], [cast_dag()], [cast_fit()]
#'
#' @export
cast_features <- function(data, screen, dag, ate,
                          model_type       = "cast",
                          max_interactions = 15L) {
  all_vars  <- dag$nodes
  cast_vars <- screen$selected
  edges     <- dag$edges
  est       <- ate$estimates

  X_base <- as.data.frame(data[, all_vars, drop = FALSE])
  p      <- length(all_vars)

  if (model_type == "mlp") {
    return(list(
      features          = X_base,
      n_base            = p,
      n_interactions    = 0L,
      n_total           = p,
      ate_weights       = stats::setNames(rep(1, p), all_vars),
      interaction_names = character(0)
    ))
  }

  # -- Layer 1: ATE weighting --
  # Cap the weight multiplier at 3.0 to prevent a single large ATE coefficient
  # from completely dominating the feature space.
  ate_weights <- stats::setNames(rep(1.0, p), all_vars)
  for (v in all_vars) {
    idx <- which(est$variable == v)
    if (length(idx) > 0 && isTRUE(est$significant[idx[1]])) {
      coef_val <- est$coef[idx[1]]
      if (is.finite(coef_val)) {
        ate_weights[v] <- min(1.0 + abs(coef_val), 3.0)
      }
    }
  }

  X_weighted <- X_base
  for (v in all_vars) {
    X_weighted[[v]] <- X_weighted[[v]] * ate_weights[v]
  }

  if (model_type == "mlp_ate") {
    return(list(
      features          = X_weighted,
      n_base            = p,
      n_interactions    = 0L,
      n_total           = p,
      ate_weights       = ate_weights,
      interaction_names = character(0)
    ))
  }

  # -- Layer 2: DAG interaction features (CAST only) --
  # Restrict to edges between screened variables, then rank by bootstrap
  # strength and take only the top max_interactions pairs.  This prevents
  # input-dimension explosion when the learned DAG is dense.
  int_cols  <- list()
  int_names <- character(0)

  if (nrow(edges) > 0 && length(cast_vars) > 0) {
    # Filter to edges whose both endpoints are in the screened set
    edge_mask <- edges$from %in% cast_vars &
                 edges$to   %in% cast_vars &
                 edges$from %in% all_vars  &
                 edges$to   %in% all_vars
    eligible_edges <- edges[edge_mask, , drop = FALSE]

    # Sort descending by bootstrap strength and cap
    if (nrow(eligible_edges) > 0) {
      eligible_edges <- eligible_edges[
        order(eligible_edges$strength, decreasing = TRUE), ,
        drop = FALSE
      ]
      n_use <- min(nrow(eligible_edges), as.integer(max_interactions))
      eligible_edges <- eligible_edges[seq_len(n_use), , drop = FALSE]

      for (k in seq_len(nrow(eligible_edges))) {
        from_v   <- eligible_edges$from[k]
        to_v     <- eligible_edges$to[k]
        col_name <- paste0("int_", from_v, "_", to_v)
        edge_w   <- eligible_edges$strength[k]
        int_cols[[col_name]] <- X_base[[from_v]] * X_base[[to_v]] * edge_w
        int_names <- c(int_names, col_name)
      }
    }
  }

  X_out <- if (length(int_cols) > 0) {
    cbind(X_weighted, as.data.frame(int_cols))
  } else {
    X_weighted
  }

  list(
    features          = X_out,
    n_base            = p,
    n_interactions    = length(int_cols),
    n_total           = ncol(X_out),
    ate_weights       = ate_weights,
    interaction_names = int_names
  )
}
