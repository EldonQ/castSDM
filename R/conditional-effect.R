# Importance reporting ------------------------------------------------------
#
# cast_importance() tidies the screen's forward-selection path: each admitted
# predictor with the inner-CV loss improvement at admission. It reuses the
# screen so no new estimation engine or hand-tuned knob is introduced.

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

#' Forward-Selection Path from the Screen
#'
#' Turns a two-stage screen into a tidy per-predictor table carrying the
#' forward-selection path: each admitted predictor with its admission step
#' and the inner cross-validated loss improvement at admission (`loss_gain`).
#' Excluded survivors are listed with missing steps and gains.
#'
#' @section Interpretation (read before citing):
#' The table describes the \emph{fitted screening procedure}. It is not a
#' causal effect and carries no p-values: a causal reading additionally
#' requires a justified adjustment set, no uncontrolled confounding,
#' consistency, joint support and adequate response and observation models.
#' These assumptions are not verified by the ranking (Byrnes & Dee 2025).
#'
#' Gains carry no sign convention beyond "lower loss is better", so read the
#' path together with [cast_effect_table()] for response direction.
#'
#' @param object A `cast_select` from `method = "two_stage"`, or a
#'   `cast_fit` / `cast_result` that carries such a screen.
#'
#' @return A `cast_importance` object.
#' @references
#' Byrnes, J. E. K. & Dee, L. E. (2025). Causal inference with observational
#' data and unobserved confounding variables. *Ecology Letters*, 28(1), e70023.
#' @seealso [cast_select()], [cast_effect_table()]
#' @export
cast_importance <- function(object) {
  screen <- .cast_extract_screen(object)
  sc <- screen$scores
  needed <- c("step_added", "loss_gain")
  if (!all(needed %in% names(sc))) {
    cli::cli_abort(c(
      "{.fn cast_importance} needs a two-stage screen carrying a forward-selection path.",
      i = "Run {.code cast_select(..., method = \"two_stage\")} first."))
  }
  effects <- data.frame(
    variable              = sc$variable,
    step_added            = sc$step_added,
    loss_gain             = sc$loss_gain,
    selected              = sc$selected,
    kept_by_design        = sc$kept_by_design,
    selected_reason       = sc$selected_reason,
    stringsAsFactors = FALSE
  )
  effects <- effects[order(is.na(effects$step_added), effects$step_added), , drop = FALSE]
  rownames(effects) <- NULL

  diagnostics <- screen$diagnostics
  diagnostics$measure <- "loss_gain"
  new_cast_importance(
    effects = effects,
    metric = screen$diagnostics$metric %||% "brier",
    diagnostics = diagnostics
  )
}


#' @keywords internal
#' @noRd
.cast_predict_matrix <- function(fit, X_raw, models, chunk_size = 200000L) {
  n <- nrow(X_raw)
  predict_block <- function(Xi) {
    preds <- lapply(models, function(m) {
      tryCatch(predict_single_model(fit$models[[m]], Xi),
               error = function(e) rep(NA_real_, nrow(Xi)))
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

#' Sensitivity product (removed in 0.12.0)
#'
#' Removed: use [cast_effect_table()] / [cast_effect_map()].
#' @param ... Ignored.
#' @return Never returns; always aborts.
#' @export
cast_sensitivity <- function(...) {
  cli::cli_abort("`cast_sensitivity()` was removed in 0.12.0; use `cast_effect_table()` / `cast_effect_map()`.")
}
