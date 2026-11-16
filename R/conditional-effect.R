# Descriptive reporting for the saved screen; no additional fitting.

#' @keywords internal
#' @noRd
.cast_extract_screen <- function(object) {
  screen <- if (inherits(object, "cast_select")) {
    object
  } else if (inherits(object, "cast_fit")) {
    object$screen
  } else if (inherits(object, "cast_result")) {
    object$screen
  } else {
    NULL
  }
  if (is.null(screen) || !inherits(screen, "cast_select")) {
    cli::cli_abort(
      "{.arg object} must be a {.cls cast_select}, {.cls cast_fit}, or {.cls cast_result} carrying a screen."
    )
  }
  screen
}

#' Describe the Fitted Variable Screen
#'
#' Returns selection status, qualification reasons, and a descriptive
#' training-scale group signal when available. Group signal is the standard
#' deviation of a variable's fitted additive component on the selector training
#' rows; it is not ecological importance, a causal effect or a p-value.
#' @param object A `cast_select`, `cast_fit`, or `cast_result` carrying a screen.
#' @return A `cast_importance` with the screen table and diagnostics.
#' @seealso [cast_select()], [cast_effect_table()]
#' @export
cast_importance <- function(object) {
  screen <- .cast_extract_screen(object)
  sc <- screen$scores
  cols <- intersect(c("variable", "role", "reason", "group_signal_sd", "freq",
    "step_added", "selected", "selected_reason"), names(sc))
  effects <- sc[, cols, drop = FALSE]
  metric <- if ("group_signal_sd" %in% names(sc)) "training_group_signal_sd" else
    if ("freq" %in% names(sc)) "fold_selection_frequency" else "selection_status"
  new_cast_importance(effects = effects, metric = metric, diagnostics = screen$diagnostics)
}

#' @keywords internal
#' @noRd
.cast_predict_matrix <- function(fit, X_raw, models, chunk_size = 200000L, clamp = FALSE) {
  n <- nrow(X_raw)
  predict_block <- function(Xi) {
    preds <- lapply(models, function(m) {
      value <- tryCatch(predict_single_model(fit$models[[m]], Xi, clamp = clamp,
        reference = fit$scaling$reference), error = function(e)
          cli::cli_abort("Prediction failed for {.val {m}}: {conditionMessage(e)}"))
      if (!is.numeric(value) || length(value) != nrow(Xi) || any(!is.finite(value)))
        cli::cli_abort("Model {.val {m}} returned missing or invalid scores on complete inputs.")
      value
    })
    do.call(cbind, preds)
  }
  # Predict in row chunks. A single full-grid predict() (especially mgcv GAM,
  # which allocates an n-by-ncoef basis matrix) can need many GB on national
  # 1km grids (~10M cells) and segfault the process. Chunking bounds peak
  # memory to one block at a time; results are identical (row-independent).
  if (n <= chunk_size) return(predict_block(X_raw))
  out <- matrix(NA_real_, nrow = n, ncol = length(models))
  for (s in seq.int(1L, n, by = chunk_size)) {
    idx <- seq.int(s, min(s + chunk_size - 1L, n))
    out[idx, ] <- predict_block(X_raw[idx, , drop = FALSE])
  }
  out
}
