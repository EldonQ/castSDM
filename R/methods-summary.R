# Summary Methods --------------------------------------------------------------

#' @export
summary.cast_dag <- function(object, ...) {
  cat("castSDM DAG Summary\n")
  cat("===================\n")
  cat("Nodes:", length(object$nodes), "\n")
  cat("Edges:", nrow(object$edges), "\n")
  if (!is.null(object$response_node))
    cat("Response node:", object$response_node, "\n")
  cat("Structure method:", object$structure_method %||% "pc", "\n")
  cat("Bootstrap R:", object$boot_R, "\n")
  cat("Strength threshold:", object$strength_threshold, "\n")
  cat("Direction threshold:", object$direction_threshold, "\n")
  if (nrow(object$edges) > 0) {
    cat("\nTop edges by strength:\n")
    top <- object$edges[order(-object$edges$strength), ]
    print(utils::head(top, 10))
  }
  invisible(object)
}

#' @export
summary.cast_eval <- function(object, ...) {
  cat("castSDM Model Evaluation Summary\n")
  cat("================================\n")
  print(object$metrics)
  if (nrow(object$metrics) > 1) {
    best_auc <- object$metrics[which.max(object$metrics$auc_mean), ]
    cat("\nBest AUC:", best_auc$model, "=", round(best_auc$auc_mean, 4), "\n")
  }
  invisible(object)
}

#' @export
summary.cast_result <- function(object, ...) {
  cat("castSDM Pipeline Result Summary\n")
  cat("===============================\n\n")
  cat("-- DAG --\n")
  cat("  Nodes:", length(object$dag$nodes),
      "| Edges:", nrow(object$dag$edges), "\n")
  if (!is.null(object$dag$response_node))
    cat("  Response node:", object$dag$response_node, "\n")
  cat("\n")
  cat("-- Variable Selection --\n")
  cat("  Selected:", length(object$screen$selected), "variables\n")
  cat("  ", paste(object$screen$selected, collapse = ", "), "\n")
  if (!is.null(object$screen$roles) && nrow(object$screen$roles) > 0) {
    role_tbl <- table(object$screen$roles$role)
    cat("  Roles:", paste(names(role_tbl), role_tbl, sep = "=", collapse = ", "), "\n")
  }
  cat("\n")
  cat("-- Models --\n")
  if (!is.null(object$eval)) {
    print(object$eval$metrics)
  }
  if (!is.null(object$ensemble)) {
    cat("\n-- Ensemble --\n")
    cat("  Method:", object$ensemble$method, "\n")
    cat("  Threshold:", round(object$ensemble$threshold, 3), "\n")
  }
  invisible(object)
}
