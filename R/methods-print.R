# Print Methods ---------------------------------------------------------------

#' @export
print.cast_dag <- function(x, ...) {
  n_edges <- if (is.null(x$edges)) 0L else nrow(x$edges)
  n_nodes <- length(x$nodes)
  cli::cli_h1("CAST DAG")
  cli::cli_ul(c(
    "Nodes: {n_nodes}",
    "Edges: {n_edges} (strength >= {x$strength_threshold})",
    "Bootstrap replicates: {x$boot_R}",
    "Score: {x$score}"
  ))
  invisible(x)
}

#' @export
print.cast_ate <- function(x, ...) {
  n_total <- nrow(x$estimates)
  n_sig <- sum(x$estimates$significant, na.rm = TRUE)
  cli::cli_h1("CAST ATE (Double Machine Learning)")
  cli::cli_ul(c(
    "Variables tested: {n_total}",
    "Significant (p < {x$alpha}): {n_sig}",
    "Cross-fitting folds: {x$K}"
  ))
  invisible(x)
}

#' @export
print.cast_screen <- function(x, ...) {
  n_selected <- length(x$selected)
  cli::cli_h1("CAST Variable Screening")
  cli::cli_ul(c(
    "Selected variables: {n_selected}",
    "Weights: DAG={round(x$weights['w_dag'], 2)}, ATE={round(x$weights['w_ate'], 2)}, IMP={round(x$weights['w_imp'], 2)}"
  ))
  cli::cli_text("Variables: {.val {x$selected}}")
  invisible(x)
}

#' @export
print.cast_roles <- function(x, ...) {
  role_counts <- table(x$roles$role)
  cli::cli_h1("CAST Causal Roles")
  for (r in names(role_counts)) {
    cli::cli_li("{r}: {role_counts[r]}")
  }
  invisible(x)
}

#' @export
print.cast_fit <- function(x, ...) {
  model_names <- names(x$models)
  cli::cli_h1("CAST Model Fit")
  cli::cli_ul(c(
    "Models: {.val {model_names}}",
    "CAST variables: {length(x$cast_vars)}"
  ))
  invisible(x)
}

#' @export
print.cast_eval <- function(x, ...) {
  cli::cli_h1("CAST Model Evaluation")
  print(x$metrics)
  invisible(x)
}

#' @export
print.cast_predict <- function(x, ...) {
  n_sites <- nrow(x$predictions)
  cli::cli_h1("CAST Spatial Predictions")
  cli::cli_ul(c(
    "Sites: {n_sites}",
    "Models: {.val {x$models}}"
  ))
  invisible(x)
}

#' @export
print.cast_cate <- function(x, ...) {
  n_sites <- length(unique(paste(x$effects$lon, x$effects$lat)))
  cli::cli_h1("CAST Spatial CATE")
  cli::cli_ul(c(
    "Sites: {n_sites}",
    "Variables: {.val {x$variables}}",
    "Causal forest trees: {x$n_trees}"
  ))
  invisible(x)
}

#' @export
print.cast_result <- function(x, ...) {
  cli::cli_h1("CAST Pipeline Result")
  cli::cli_ul(c(
    "DAG: {nrow(x$dag$edges)} edges",
    "ATE: {sum(x$ate$estimates$significant)} significant variables",
    "Screened: {length(x$screen$selected)} variables",
    "Models: {.val {names(x$fit$models)}}",
    "Predictions: {if (!is.null(x$predict)) 'Yes' else 'No'}",
    "CATE: {if (!is.null(x$cate)) 'Yes' else 'No'}"
  ))
  invisible(x)
}
