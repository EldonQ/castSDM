# Importance reporting ------------------------------------------------------
#
# cast_importance() tidies the screen's single conditional-effect column: the
# shift effect that selected each predictor, with conditional-permutation
# p-values and null thresholds. It reuses the screen so no new estimation
# engine or hand-tuned knob is introduced.

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

#' Conditional Effect Table from the Screen
#'
#' Turns a two-stage screen into a tidy per-predictor table carrying the
#' single \strong{conditional effect} that selected the predictor (stage 2 of
#' [cast_select()]). Within-stratum permutation p-values refer only to this
#' statistic. Response-dependent stage-1 filtering is not repeated under the
#' null, so these are exploratory scores, not confirmatory tests with
#' guaranteed error control.
#'
#' @section Interpretation (read before citing):
#' The column describes the \emph{fitted model}. It is not a causal effect:
#' a causal reading additionally requires a justified adjustment set, no
#' uncontrolled confounding, consistency, joint support and adequate response
#' and observation models. These assumptions are not verified by the ranking
#' (Byrnes & Dee 2025).
#'
#' Importance carries no sign, so read it together with
#' [cast_effect_table()] for direction.
#'
#' @param object A `cast_select` from `method = "two_stage"`, or a
#'   `cast_fit` / `cast_result` that carries such a screen.
#'
#' @return A `cast_importance` object.
#' @references
#' Hooker, G., Mentch, L. & Zhou, S. (2021). Unrestricted permutation forces
#' extrapolation: variable importance requires at least one more model, or
#' there is no free variable importance. *Statistics and Computing*, 31, 82.
#' \doi{10.1007/s11222-021-10057-z}.
#' @seealso [cast_select()], [cast_effect_table()]
#' @export
cast_importance <- function(object) {
  screen <- .cast_extract_screen(object)
  sc <- screen$scores
  needed <- c("interventional_effect", "p_value")
  if (!all(needed %in% names(sc))) {
    cli::cli_abort(c(
      "{.fn cast_importance} needs a two-stage screen carrying a stage-2 statistic.",
      i = "Run {.code cast_select(..., method = \"two_stage\")} first."))
  }
  sc <- sc[is.finite(sc$interventional_effect), , drop = FALSE]
  if (!nrow(sc)) {
    cli::cli_abort(c(
      "The screen holds no finite conditional-effect estimates.",
      i = "{.code method = \"full\"} skips stage 2, so there is nothing to report."))
  }
  effects <- data.frame(
    variable              = sc$variable,
    interventional_effect = sc$interventional_effect,
    null_threshold        = sc$null_threshold,
    p_value               = sc$p_value,
    selected              = sc$selected,
    kept_by_design        = sc$kept_by_design,
    selected_reason       = sc$selected_reason,
    stringsAsFactors = FALSE
  )
  effects <- effects[order(-effects$interventional_effect), , drop = FALSE]
  rownames(effects) <- NULL

  diagnostics <- screen$diagnostics
  diagnostics$measure <- "interventional_effect"
  new_cast_importance(
    effects = effects,
    alpha = screen$diagnostics$alpha %||% 0.05,
    threshold = screen$diagnostics$null_threshold %||% NA_real_,
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
