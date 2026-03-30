#' Causal Feature Engineering
#'
#' Constructs augmented feature matrices using ATE-weighted variables and
#' DAG-guided interaction terms. Three feature layers are produced:
#' 1. **ATE-weighted**: Base variables scaled by causal effect magnitude.
#' 2. **DAG interactions**: Pairwise products of connected variables weighted
#'    by edge strength (CAST model only).
#' 3. **Raw**: Unmodified variables (for baseline models).
#'
#' @param data A `data.frame` of standardized environmental variables.
#' @param screen A [cast_screen] object.
#' @param dag A [cast_dag] object.
#' @param ate A [cast_ate] object.
#' @param model_type Character. Which feature set to build: `"cast"` (all
#'   layers), `"mlp_ate"` (ATE-weighted only), `"mlp"` (raw only). Default
#'   `"cast"`.
#'
#' @return A list with:
#' \describe{
#'   \item{`features`}{A `data.frame` of engineered features.}
#'   \item{`n_base`}{Number of base variables.}
#'   \item{`n_interactions`}{Number of DAG interaction features.}
#'   \item{`n_total`}{Total number of features.}
#'   \item{`ate_weights`}{Named numeric vector of ATE weights applied.}
#'   \item{`interaction_names`}{Names of interaction columns (if any).}
#' }
#'
#' @seealso [cast_screen()], [cast_ate()], [cast_dag()], [cast_fit()]
#'
#' @export
cast_features <- function(data, screen, dag, ate,
                          model_type = "cast") {
  all_vars <- dag$nodes
  cast_vars <- screen$selected
  edges <- dag$edges
  est <- ate$estimates

  X_base <- as.data.frame(data[, all_vars, drop = FALSE])
  p <- length(all_vars)

  if (model_type == "mlp") {
    # Raw features, no transformation
    return(list(
      features = X_base,
      n_base = p, n_interactions = 0L,
      n_total = p,
      ate_weights = stats::setNames(rep(1, p), all_vars),
      interaction_names = character(0)
    ))
  }

  # -- Layer 1: ATE weighting --
  ate_weights <- stats::setNames(rep(1.0, p), all_vars)
  for (v in all_vars) {
    idx <- which(est$variable == v)
    if (length(idx) > 0 && isTRUE(est$significant[idx[1]])) {
      coef_val <- est$coef[idx[1]]
      if (is.finite(coef_val)) {
        ate_weights[v] <- 1.0 + abs(coef_val)
      }
    }
  }

  X_weighted <- X_base
  for (v in all_vars) {
    X_weighted[[v]] <- X_weighted[[v]] * ate_weights[v]
  }

  if (model_type == "mlp_ate") {
    return(list(
      features = X_weighted,
      n_base = p, n_interactions = 0L,
      n_total = p,
      ate_weights = ate_weights,
      interaction_names = character(0)
    ))
  }

  # -- Layer 2: DAG interaction features (CAST only) --
  int_cols <- list()
  int_names <- character(0)

  if (nrow(edges) > 0 && length(cast_vars) > 0) {
    for (k in seq_len(nrow(edges))) {
      from_v <- edges$from[k]
      to_v <- edges$to[k]
      if (from_v %in% cast_vars && to_v %in% cast_vars &&
          from_v %in% all_vars && to_v %in% all_vars) {
        col_name <- paste0("int_", from_v, "_", to_v)
        edge_w <- edges$strength[k]
        int_cols[[col_name]] <-
          X_base[[from_v]] * X_base[[to_v]] * edge_w
        int_names <- c(int_names, col_name)
      }
    }
  }

  if (length(int_cols) > 0) {
    X_out <- cbind(X_weighted, as.data.frame(int_cols))
  } else {
    X_out <- X_weighted
  }

  list(
    features = X_out,
    n_base = p,
    n_interactions = length(int_cols),
    n_total = ncol(X_out),
    ate_weights = ate_weights,
    interaction_names = int_names
  )
}
