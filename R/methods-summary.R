# Summary Methods -----------------------------------------------------------

#' @export
summary.cast_eval <- function(object, ...) {
  cli::cli_h1("castSDM Model Evaluation Summary")
  print(object$metrics)
  if (nrow(object$metrics) > 1) {
    best <- which.max(object$metrics$auc_mean)
    if (length(best)) {
      best_auc <- object$metrics[best, ]
      cli::cli_inform("Best AUC: {best_auc$model} = {round(best_auc$auc_mean, 4)}")
    } else {
      cli::cli_inform("No model has a finite AUC.")
    }
  }
  invisible(object)
}

#' @export
summary.cast_result <- function(object, ...) {
  cli::cli_h1("castSDM Pipeline Result Summary")
  cli::cli_h2("Variable Selection")
  cli::cli_ul(c(
    "Selected: {length(object$screen$selected)} variables",
    "{paste(object$screen$selected, collapse = ', ')}"
  ))
  if (!is.null(object$screen$method)) {
    cli::cli_text("Screening method: {object$screen$method}")
  }
  thr <- object$screen$diagnostics$null_threshold
  if (!is.null(thr) && is.finite(thr)) {
    cli::cli_text("Permutation-null threshold: {signif(thr, 3)}")
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
