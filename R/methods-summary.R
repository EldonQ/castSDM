# Summary Methods --------------------------------------------------------------

#' @export
summary.cast_dag <- function(object, ...) {
  cli::cli_h1("castSDM DAG Summary")
  cli::cli_ul(c(
    "Nodes: {length(object$nodes)}",
    "Edges: {nrow(object$edges)}",
    if (!is.null(object$response_node))
      "Response node: {object$response_node}" else NULL,
    "Structure method: {object$structure_method %||% 'pc'}",
    "Bootstrap R: {object$boot_R}",
    "Strength threshold: {object$strength_threshold}",
    "Direction threshold: {object$direction_threshold}"
  ))
  if (nrow(object$edges) > 0) {
    cli::cli_h2("Top edges by strength")
    top <- object$edges[order(-object$edges$strength), ]
    print(utils::head(top, 10))
  }
  invisible(object)
}

#' @export
summary.cast_eval <- function(object, ...) {
  cli::cli_h1("castSDM Model Evaluation Summary")
  print(object$metrics)
  if (nrow(object$metrics) > 1) {
    best_auc <- object$metrics[which.max(object$metrics$auc_mean), ]
    cli::cli_inform("Best AUC: {best_auc$model} = {round(best_auc$auc_mean, 4)}")
  }
  invisible(object)
}

#' @export
summary.cast_result <- function(object, ...) {
  cli::cli_h1("castSDM Pipeline Result Summary")
  cli::cli_h2("DAG")
  cli::cli_ul(c(
    "Nodes: {length(object$dag$nodes)} | Edges: {nrow(object$dag$edges)}",
    if (!is.null(object$dag$response_node))
      "Response node: {object$dag$response_node}" else NULL
  ))
  cli::cli_h2("Variable Selection")
  cli::cli_ul(c(
    "Selected: {length(object$screen$selected)} variables",
    "{paste(object$screen$selected, collapse = ', ')}"
  ))
  if (!is.null(object$screen$roles) && nrow(object$screen$roles) > 0) {
    role_tbl <- table(object$screen$roles$role)
    cli::cli_text("Roles: {paste(names(role_tbl), role_tbl, sep = '=', collapse = ', ')}")
  }
  cli::cli_h2("Models")
  if (!is.null(object$eval)) {
    print(object$eval$metrics)
  }
  if (!is.null(object$ensemble)) {
    cli::cli_h2("Ensemble")
    cli::cli_ul(c(
      "Method: {object$ensemble$method}",
      "Threshold: {round(object$ensemble$threshold, 3)}"
    ))
  }
  invisible(object)
}
