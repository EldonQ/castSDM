#' Assign Causal Roles to Variables
#'
#' Builds a role table from a `cast_screen` object produced by
#' [cast_select()]. When the screen contains Markov Blanket roles
#' (`parent`, `child`, `co_parent`, `predictive`), those are used
#' directly. For backward compatibility with older screen objects
#' that lack roles, a degree-based heuristic is applied.
#'
#' @param screen A `cast_screen` object from [cast_select()].
#' @param dag A [cast_dag] object.
#'
#' @return A [cast_roles] object.
#'
#' @seealso [cast_select()], [cast_dag()]
#'
#' @export
cast_roles <- function(screen, dag) {
  selected_vars <- screen$selected
  edges <- dag$edges

  # Compute in/out-degree for selected vars (always useful metadata)
  deg_df <- compute_edge_degrees(edges, selected_vars)

  # Read roles from screen if available (cast_select produces these)
  if (!is.null(screen$roles) && is.data.frame(screen$roles) &&
      "role" %in% names(screen$roles)) {
    role_tbl <- screen$roles
    # Merge in/out-degree
    role_tbl <- merge(role_tbl, deg_df[, c("variable", "in_degree",
                                            "out_degree")],
                      by = "variable", all.x = TRUE)
    role_tbl$in_degree[is.na(role_tbl$in_degree)] <- 0L
    role_tbl$out_degree[is.na(role_tbl$out_degree)] <- 0L
    role_tbl$role_score <- ifelse(
      role_tbl$in_degree == 0,
      role_tbl$out_degree + 1,
      role_tbl$out_degree / (role_tbl$in_degree + 1)
    )
  } else {
    # Fallback: degree-based role assignment (backward compat)
    role_tbl <- deg_df
    role_tbl$role_score <- ifelse(
      role_tbl$in_degree == 0,
      role_tbl$out_degree + 1,
      role_tbl$out_degree / (role_tbl$in_degree + 1)
    )
    role_tbl <- role_tbl[order(-role_tbl$role_score), ]
    n_vars <- nrow(role_tbl)
    if (n_vars < 3) {
      labels <- c("parent", "child", "predictive")
      role_tbl$role <- labels[seq_len(n_vars)]
    } else {
      s1 <- ceiling(n_vars / 3)
      s2 <- ceiling((n_vars - s1) / 2)
      s3 <- max(0, n_vars - s1 - s2)
      role_tbl$role <- c(
        rep("parent", s1),
        rep("child", s2),
        rep("predictive", s3)
      )
    }
  }

  rownames(role_tbl) <- NULL
  new_cast_roles(role_tbl)
}
