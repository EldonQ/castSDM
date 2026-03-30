#' Assign Causal Roles to Variables
#'
#' Classifies screened variables into causal roles based on DAG topology:
#' **Root causes** (exogenous drivers), **Mediators** (intermediate pathways),
#' or **Terminal** (outcome proxies).
#'
#' @param screen A [cast_screen] object.
#' @param dag A [cast_dag] object.
#'
#' @return A [cast_roles] object.
#'
#' @details
#' Role assignment is based on in-degree and out-degree in the DAG:
#' - **Root**: Zero in-degree, high out-degree (exogenous causal drivers).
#' - **Mediator**: Both in-degree and out-degree > 0 (causal pathways).
#' - **Terminal**: Low out-degree or high in-degree (close to outcome).
#'
#' Role score: if `in_degree == 0` then `out_degree + 1`, else
#' `out_degree / (in_degree + 1)`.
#'
#' @seealso [cast_screen()], [cast_features()]
#'
#' @export
cast_roles <- function(screen, dag) {
  selected_vars <- screen$selected
  edges <- dag$edges

  # Compute out-degree and in-degree for selected vars
  if (nrow(edges) > 0) {
    out_agg <- stats::aggregate(
      to ~ from, data = edges, FUN = length
    )
    names(out_agg) <- c("variable", "out_degree")
    in_agg <- stats::aggregate(
      from ~ to, data = edges, FUN = length
    )
    names(in_agg) <- c("variable", "in_degree")
  } else {
    out_agg <- data.frame(
      variable = character(0), out_degree = integer(0)
    )
    in_agg <- data.frame(
      variable = character(0), in_degree = integer(0)
    )
  }

  role_df <- data.frame(
    variable = selected_vars, stringsAsFactors = FALSE
  )
  role_df <- merge(role_df, out_agg, by = "variable", all.x = TRUE)
  role_df <- merge(role_df, in_agg, by = "variable", all.x = TRUE)
  role_df$out_degree[is.na(role_df$out_degree)] <- 0
  role_df$in_degree[is.na(role_df$in_degree)] <- 0

  role_df$role_score <- ifelse(
    role_df$in_degree == 0,
    role_df$out_degree + 1,
    role_df$out_degree / (role_df$in_degree + 1)
  )
  role_df <- role_df[order(-role_df$role_score), ]

  # Assign roles by rank partitioning
  n_vars <- nrow(role_df)
  n_groups <- 3L
  if (n_vars < n_groups) {
    labels <- c("Root", "Mediator", "Terminal")
    role_df$role <- labels[seq_len(n_vars)]
  } else {
    s1 <- ceiling(n_vars / n_groups)
    s2 <- ceiling((n_vars - s1) / 2)
    s3 <- max(0, n_vars - s1 - s2)
    role_df$role <- c(
      rep("Root", s1),
      rep("Mediator", s2),
      rep("Terminal", s3)
    )
  }
  rownames(role_df) <- NULL

  new_cast_roles(role_df)
}
