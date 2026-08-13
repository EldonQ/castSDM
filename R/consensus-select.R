# Consensus selection -------------------------------------------------------
#
# Nested cross-validation re-runs the conditional screen inside every outer
# spatial fold (cast_cv()). A predictor that is retained in most folds is
# conditionally informative across regions -- the operational counterpart of
# "transferable driver". cast_consensus() aggregates the fold selections
# into a frequency table and a thresholded consensus set, and returns a
# cast_select object so the consensus set plugs straight into cast_fit().
#
# This answers the spatial-stability weakness honestly: instead of hiding
# fold disagreement, the package surfaces it (cast_cv()$selection_freq) and
# lets the user decide how stable the final set must be.

#' Consensus Selection Across Spatial Folds
#'
#' Aggregates the fold-specific selections stored in a [cast_cv()] object
#' into a selection-frequency table and a consensus variable set: predictors
#' retained in at least `threshold` (a proportion) of the outer folds.
#'
#' Nested selection is re-run in every spatial fold, so a predictor selected
#' in only one fold is likely a region-specific (or noise) association; a
#' predictor selected in most folds is a candidate driver that transfers
#' across regions. The frequency table makes that trade-off visible and the
#' consensus set operationalises it.
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
    sets <- Filter(function(s) length(s) > 0L, cv$selections %||% list())
    if (!length(sets)) {
      cli::cli_abort(c(
        "The CV object holds no fold selections.",
        i = "Run {.code cast_cv(..., select_method = \"cpi\")} first."
      ))
    }
    all_vars <- unique(unlist(sets))
    f <- vapply(all_vars, function(v) {
      mean(vapply(sets, function(s) v %in% s, logical(1)))
    }, numeric(1))
    freq <- data.frame(variable = all_vars, freq = unname(f),
                       stringsAsFactors = FALSE)
    freq <- freq[order(-freq$freq), , drop = FALSE]
    rownames(freq) <- NULL
  }

  freq$selected <- freq$freq >= threshold
  scores <- freq[, c("variable", "freq", "selected"), drop = FALSE]
  n_folds <- if (length(cv$selections)) {
    sum(vapply(cv$selections, function(s) length(s) > 0L, logical(1)))
  } else {
    cv$k
  }

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
