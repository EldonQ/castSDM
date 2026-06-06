#' DAG-Guided Variable Selection via Markov Blanket + RF Importance
#'
#' Selects predictor variables by fusing two signals:
#'
#' 1. **Markov Blanket (MB)** of the response node extracted from a causal
#'    DAG — theoretically the minimal sufficient feature set.
#' 2. **RF permutation importance** — data-driven predictive relevance.
#'
#' All MB variables are automatically selected.
#' Non-MB variables are added only if their RF importance exceeds the
#' median importance of the MB variables.
#'
#' Each selected variable is assigned a **screening role**:
#' - `"mb_direct"`: adjacent to the response in the response-focused MB graph.
#' - `"mb_associated"`: retained through non-adjacent Markov Blanket structure
#'   in an unconstrained local graph.
#' - `"importance_added"`: outside the MB but retained by RF importance.
#' - `"importance_screened"`: retained by RF importance when the MB is dense
#'   or otherwise non-selective.
#'
#' @param dag A [cast_dag] object, ideally learned with
#'   `include_response = TRUE`.
#' @param data A `data.frame` with species occurrence and environmental
#'   variables.
#' @param response Character. Name of the response column. Default
#'   `"presence"`.
#' @param num_trees Integer. Number of trees for the RF importance step.
#'   Default `300`.
#' @param min_vars Integer. Minimum variables to retain. Default `5`.
#' @param min_fraction Numeric in `[0, 1]`. Minimum fraction of candidate
#'   variables to retain. Default `0.3`.
#' @param seed Integer or `NULL`. Random seed. Default `NULL`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A `cast_select` object with components:
#' \describe{
#'   \item{selected}{Character vector of selected variable names.}
#'   \item{scores}{A `data.frame` with per-variable scores and metadata.}
#'   \item{roles}{A `data.frame` with columns `variable` and `role`.}
#' }
#'
#' @details
#' When the DAG was learned without the response node
#' (`include_response = FALSE`), the Markov Blanket cannot be computed.
#' In that case, variable selection falls back to pure RF importance
#' with a score-gap heuristic; all selected variables receive the
#' `"importance_screened"` role.
#'
#' @seealso [cast_dag()], [cast_fit()]
#'
#' @export
cast_select <- function(dag,
                        data,
                        response = "presence",
                        num_trees = 300L,
                        min_vars = 5L,
                        min_fraction = 0.3,
                        seed = NULL,
                        verbose = TRUE) {
  check_suggested("ranger", "for RF importance")

  # --- Determine candidate env vars (exclude response, lon, lat) -----------
  response_node <- dag$response_node
  env_vars <- setdiff(dag$nodes, c(response, "lon", "lat"))
  if (length(env_vars) < 3) {
    cli::cli_abort("Need at least 3 environmental variables for selection.")
  }

  # --- Extract Markov Blanket (if response was in DAG) ---------------------
  mb <- NULL
  if (!is.null(response_node)) {
    mb <- extract_markov_blanket(dag, response_node)
    stage1_mb <- dag$metadata$mb_vars_stage1 %||% character()
    if (length(stage1_mb) > 0) {
      mb$all <- unique(intersect(stage1_mb, env_vars))
      if (isTRUE(dag$metadata$response_as_sink)) {
        mb$parents <- mb$all
        mb$children <- character()
        mb$co_parents <- character()
      }
    }
    if (verbose && length(mb$all) > 0) {
      cli::cli_inform(
        "Markov Blanket of {.val {response_node}}: {length(mb$all)} variables ({length(mb$parents)} direct MB, {length(mb$children) + length(mb$co_parents)} associated MB)."
      )
    }
  }

  # --- RF permutation importance -------------------------------------------
  Y <- data[[response]]
  X <- as.data.frame(data[, env_vars, drop = FALSE])
  for (col in names(X)) X[[col]] <- as.numeric(X[[col]])
  X[is.na(X)] <- 0
  if (!is.null(seed)) set.seed(seed)
  rf_imp <- ranger::ranger(
    y ~ .,
    data = cbind(y = as.factor(Y), X),
    num.trees = num_trees,
    importance = "permutation",
    verbose = FALSE
  )$variable.importance

  # --- Build scoring data.frame --------------------------------------------
  deg_df <- compute_edge_degrees(dag$edges, env_vars)
  scores_df <- deg_df
  scores_df$importance <- rf_imp[scores_df$variable]
  scores_df$importance[is.na(scores_df$importance)] <- 0

  # Mark MB membership
  scores_df$in_mb <- scores_df$variable %in% (mb$all %||% character())
  scores_df$mb_role <- vapply(scores_df$variable, function(v) {
    if (v %in% (mb$parents %||% character())) return("mb_direct")
    if (v %in% (mb$children %||% character())) return("mb_associated")
    if (v %in% (mb$co_parents %||% character())) return("mb_associated")
    "none"
  }, character(1))

  # Normalize importance to [0, 1]
  scores_df$imp_norm <- normalize01(scores_df$importance)

  # --- Selection logic -----------------------------------------------------
  min_keep <- max(min_vars, ceiling(length(env_vars) * min_fraction))
  select_by_importance <- function(scores, min_keep) {
    scores <- scores[order(-scores$importance), ]
    imp_sorted <- scores$importance

    if (length(imp_sorted) <= min_keep) {
      return(scores$variable)
    }

    gaps <- imp_sorted[-length(imp_sorted)] - imp_sorted[-1L]
    valid_idx <- which(seq_along(gaps) >= min_keep &
                         seq_along(gaps) < length(imp_sorted))
    if (length(valid_idx) == 0) {
      scores$variable[seq_len(min_keep)]
    } else {
      cut_at <- valid_idx[which.max(gaps[valid_idx])]
      scores$variable[seq_len(cut_at)]
    }
  }
  mb_selective <- TRUE

  if (!is.null(mb) && length(mb$all) > 0) {
    # -- Fusion: MB vars auto-selected + importance-gated non-MB vars --
    mb_vars <- intersect(mb$all, env_vars)
    mb_dense <- length(mb_vars) >= ceiling(0.8 * length(env_vars))

    if (isTRUE(mb_dense)) {
      if (verbose) {
        cli::cli_warn(c(
          "Dense Markov Blanket: {length(mb_vars)}/{length(env_vars)} predictors are in MB.",
          "i" = "Treating MB membership as non-selective and falling back to RF importance screening."
        ))
      }
      mb_selective <- FALSE
      selected <- select_by_importance(scores_df, min_keep)
    } else {
      mb_importance <- scores_df$importance[scores_df$variable %in% mb_vars]
      imp_threshold <- if (length(mb_importance) > 0) {
        stats::median(mb_importance)
      } else {
        0
      }

      # Non-MB vars that pass the threshold
      non_mb <- scores_df[!scores_df$in_mb, , drop = FALSE]
      extra_vars <- non_mb$variable[non_mb$importance > imp_threshold]

      selected <- unique(c(mb_vars, extra_vars))

      # Enforce minimum
      if (length(selected) < min_keep) {
        remaining <- setdiff(env_vars, selected)
        rem_imp <- scores_df$importance[match(remaining, scores_df$variable)]
        top_rem <- remaining[order(-rem_imp)][seq_len(
          min(min_keep - length(selected), length(remaining))
        )]
        selected <- unique(c(selected, top_rem))
      }
    }
  } else {
    # -- Fallback: pure RF importance with score-gap heuristic --
    selected <- select_by_importance(scores_df, min_keep)
  }

  # --- Assign screening roles ----------------------------------------------
  roles_df <- data.frame(
    variable = selected,
    role = vapply(selected, function(v) {
      if (!isTRUE(mb_selective)) return("importance_screened")
      if (v %in% (mb$parents %||% character())) return("mb_direct")
      if (v %in% (mb$children %||% character())) return("mb_associated")
      if (v %in% (mb$co_parents %||% character())) return("mb_associated")
      "importance_added"
    }, character(1)),
    stringsAsFactors = FALSE
  )

  # Sort: MB roles first, then by importance
  role_order <- c(
    mb_direct = 1,
    mb_associated = 2,
    importance_added = 3,
    importance_screened = 4,
    parent = 1,
    child = 2,
    co_parent = 2,
    predictive = 3
  )
  roles_df$role_rank <- role_order[roles_df$role]
  roles_df$imp_val <- scores_df$importance[match(roles_df$variable,
                                                  scores_df$variable)]
  roles_df <- roles_df[order(roles_df$role_rank, -roles_df$imp_val), ]
  roles_df$role_rank <- NULL
  roles_df$imp_val <- NULL
  rownames(roles_df) <- NULL

  # Order scores_df: selected first, then by importance
  scores_df <- scores_df[order(-scores_df$in_mb, -scores_df$importance), ]
  rownames(scores_df) <- NULL

  if (verbose) {
    n_mb <- sum(roles_df$role %in% c("mb_direct", "mb_associated"))
    n_pred <- sum(roles_df$role %in% c("importance_added", "importance_screened"))
    cli::cli_inform(
      "Selection: {length(selected)}/{length(env_vars)} variables ({n_mb} from MB, {n_pred} importance-screened)."
    )
  }

  new_cast_select(
    selected = selected,
    scores = scores_df,
    roles = roles_df
  )
}
