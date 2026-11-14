# Consensus selection -------------------------------------------------------
#
# Fold inclusion frequencies describe selection stability, not causal or
# geographic transferability. Consensus is an optional descriptive subset.

#' Consensus Selection Across Spatial Folds
#'
#' Aggregates the fold-specific selections stored in a [cast_cv()] object
#' into a selection-frequency table and a consensus variable set: predictors
#' retained in at least `threshold` (a proportion) of the outer folds.
#'
#' Frequencies describe stability under the recorded fold design. Correlated
#' predictors may substitute for each other; frequency does not identify
#' ecological causes or prove geographic transferability. This optional
#' thresholded set is not used by the default automatic screen.
#'
#' @param cv A `cast_cv` object produced with a non-`NULL` `select_method`.
#' @param threshold Minimum fold-selection frequency (in `[0, 1]`) for
#'   membership in the consensus set. Default `0.5` (majority of folds).
#'
#' @return A `cast_select` object whose `selected` members form the
#'   consensus set and whose `scores` data.frame holds the per-predictor
#'   fold frequency (`variable`, `freq`, `selected`).
#' @seealso [cast_cv()], [cast_fit()]
#' @export
cast_consensus <- function(cv, threshold = 0.5) {
  if (!inherits(cv, "cast_cv")) {
    cli::cli_abort("{.arg cv} must be a {.cls cast_cv} object.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      threshold < 0 || threshold > 1) {
    cli::cli_abort("{.arg threshold} must be a single number in [0, 1].")
  }
  if (!is.null(cv$selection_freq) && nrow(cv$selection_freq)) {
    freq <- cv$selection_freq
  } else {
    # Fallback: rebuild frequency over ALL folds (denominator = total folds;
    # an empty/failed fold contributes zero, never shrinks the denominator).
    sets <- cv$selections %||% list()
    if (!length(sets)) {
      cli::cli_abort(c(
        "The CV object holds no fold selections.",
        i = "Run {.code cast_cv(..., select_method = \"vif\")} first."
      ))
    }
    all_vars <- unique(unlist(sets))
    f <- vapply(all_vars, function(v) {
      mean(vapply(sets, function(s) v %in% (s %||% character(0)), logical(1)))
    }, numeric(1))
    freq <- data.frame(variable = all_vars, freq = unname(f),
                       stringsAsFactors = FALSE)
    freq <- freq[order(-freq$freq), , drop = FALSE]
    rownames(freq) <- NULL
  }

  freq$selected <- freq$freq >= threshold
  scores <- freq[, c("variable", "freq", "selected"), drop = FALSE]
  if (!any(scores$selected)) {
    cli::cli_warn(c(
      "Empty consensus set: no predictor reached the fold-frequency threshold of {threshold}.",
      i = "Lower {.arg threshold} or inspect {.code cv$selection_freq}; {.fun cast_fit} aborts on an empty screen."
    ))
  }
  n_folds <- length(cv$selections)
  if (n_folds < 1L) n_folds <- cv$k

  new_cast_select(
    selected = scores$variable[scores$selected],
    scores = scores,
    method = "consensus",
    diagnostics = list(
      engine = sprintf("Consensus across %d spatial folds", n_folds),
      threshold = threshold,
      n_folds = n_folds
    )
  )
}
