# Summary Methods --------------------------------------------------------------

#' @export
summary.cast_dag <- function(object, ...) {
  cat("CAST DAG Summary\n")
  cat("================\n")
  cat("Nodes:", length(object$nodes), "\n")
  cat("Edges:", nrow(object$edges), "\n")
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
  cat("CAST Model Evaluation Summary\n")
  cat("=============================\n")
  print(object$metrics)
  if (nrow(object$metrics) > 1) {
    best_auc <- object$metrics[which.max(object$metrics$auc_mean), ]
    cat("\nBest AUC:", best_auc$model, "=", round(best_auc$auc_mean, 4), "\n")
  }
  invisible(object)
}

#' @export
summary.cast_result <- function(object, ...) {
  cat("CAST Pipeline Result Summary\n")
  cat("============================\n\n")
  cat("-- DAG --\n")
  cat("  Nodes:", length(object$dag$nodes),
      "| Edges:", nrow(object$dag$edges), "\n\n")
  cat("-- ATE --\n")
  cat("  Variables tested:", nrow(object$ate$estimates),
      "| Significant:", sum(object$ate$estimates$significant), "\n\n")
  cat("-- Screening --\n")
  cat("  Selected:", length(object$screen$selected), "variables\n")
  cat("  ", paste(object$screen$selected, collapse = ", "), "\n\n")
  cat("-- Models --\n")
  if (!is.null(object$eval)) {
    print(object$eval$metrics)
  }
  if (!is.null(object$shap)) {
    cat("\n-- SHAP --\n")
    shap_names <- names(Filter(Negate(is.null), object$shap))
    cat("  Explained models:", paste(shap_names, collapse = ", "), "\n")
    for (nm in shap_names) {
      sh <- object$shap[[nm]]
      cat(sprintf("  %s: %d features, %d rows, engine=%s\n",
                  nm, sh$n_features, nrow(sh$shap), sh$engine))
    }
  }
  invisible(object)
}
