# Summary Methods -----------------------------------------------------------

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
  cli::cli_h2("Variable Selection")
  cli::cli_ul(c(
    "Selected: {length(object$screen$selected)} variables",
    "{paste(object$screen$selected, collapse = ', ')}"
  ))
  if (!is.null(object$screen$diagnostics$engine)) {
    cli::cli_text(
      "Screening engine: {object$screen$diagnostics$engine}"
    )
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
