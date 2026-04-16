#' Backdoor Criterion Check via dagitty
#'
#' Converts the learned DAG to a dagitty graph and identifies valid
#' adjustment sets (backdoor criterion) for each treatment variable
#' relative to the outcome (`presence`). This verifies whether the
#' causal effect of each variable can be identified given the DAG
#' structure and available covariates.
#'
#' @param dag A [cast_dag] object.
#' @param outcome Character. The outcome variable name. Default `"presence"`.
#'   If not among DAG nodes, added as an isolated node.
#' @param verbose Logical. Print summary. Default `TRUE`.
#'
#' @return A `data.frame` with columns:
#'   \describe{
#'     \item{variable}{Treatment variable name.}
#'     \item{identifiable}{Logical. Whether a valid adjustment set exists.}
#'     \item{adjustment_set}{Character. Comma-separated minimal adjustment set
#'       (or `NA` if not identifiable).}
#'     \item{n_paths}{Integer. Number of causal paths from this variable to
#'       the outcome.}
#'   }
#'
#' @details
#' The **backdoor criterion** (Pearl, 2009) states that the causal effect of
#' \eqn{X} on \eqn{Y} is identifiable if there exists a set \eqn{Z} that:
#' (1) blocks all non-causal (backdoor) paths from \eqn{X} to \eqn{Y}, and
#' (2) does not contain any descendant of \eqn{X}.
#'
#' This function wraps [dagitty::adjustmentSets()] to find minimal sufficient
#' adjustment sets. When a variable is **not identifiable**, it means the
#' DAG structure implies that unmeasured confounders may bias the estimate
#' (though the DML procedure may still provide useful estimates under
#' approximate conditions).
#'
#' @references
#' Pearl, J. (2009). *Causality*. Cambridge University Press.
#'
#' Textor, J., van der Zander, B., Gilthorpe, M.S., Liskiewicz, M., &
#' Ellison, G.T.H. (2016). Robust causal inference using directed acyclic
#' graphs: the R package 'dagitty'. *International Journal of Epidemiology*,
#' 45(6), 1887-1894.
#'
#' @seealso [cast_dag()], [cast_ate()]
#'
#' @export
cast_backdoor <- function(dag, outcome = "presence", verbose = TRUE) {
  if (!inherits(dag, "cast_dag")) {
    cli::cli_abort("{.arg dag} must be a {.cls cast_dag} object.")
  }
  check_suggested("dagitty", "for backdoor criterion analysis")

  edges <- dag$edges
  nodes <- dag$nodes

  # Build dagitty DAG string
  edge_strs <- if (nrow(edges) > 0) {
    paste(edges$from, "->", edges$to, collapse = "; ")
  } else {
    ""
  }

  # Ensure outcome is in the node set
  all_nodes <- union(nodes, outcome)
  node_str <- paste(all_nodes, collapse = "; ")
  dag_str <- paste0("dag { ", node_str, "; ", edge_strs, " }")

  dg <- dagitty::dagitty(dag_str)

  # For each env variable, check backdoor criterion to outcome
  env_nodes <- setdiff(all_nodes, outcome)
  results <- vector("list", length(env_nodes))

  for (i in seq_along(env_nodes)) {
    v <- env_nodes[i]
    adj <- tryCatch(
      dagitty::adjustmentSets(dg, exposure = v, outcome = outcome,
                              type = "minimal"),
      error = function(e) list()
    )

    paths <- tryCatch(
      dagitty::paths(dg, from = v, to = outcome, directed = TRUE),
      error = function(e) list(paths = character(0))
    )
    n_paths <- length(paths$paths)

    if (length(adj) > 0) {
      adj_set <- paste(sort(as.character(adj[[1]])), collapse = ", ")
      results[[i]] <- data.frame(
        variable = v,
        identifiable = TRUE,
        adjustment_set = adj_set,
        n_paths = n_paths,
        stringsAsFactors = FALSE
      )
    } else {
      results[[i]] <- data.frame(
        variable = v,
        identifiable = FALSE,
        adjustment_set = NA_character_,
        n_paths = n_paths,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL

  if (verbose) {
    n_id <- sum(out$identifiable)
    cli::cli_inform(c(
      "Backdoor criterion: {n_id}/{nrow(out)} variables identifiable",
      "i" = "Variables with valid adjustment sets can be causally interpreted."
    ))
  }

  out
}
