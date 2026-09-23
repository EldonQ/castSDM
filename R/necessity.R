# Necessity (knockout) diagnostic removed in 0.12.0.
# The 0.11.0 RF-knockout dAUC used a different estimator from the effect
# products, so paired reading could not distinguish substitutability from
# estimator switching. Use cast_effect_table()/cast_effect_map() only.

#' Necessity diagnostic (removed in 0.12.0)
#'
#' Removed: the knockout audit is discontinued. Read the single shift effect
#' from [cast_effect_table()] / [cast_effect_map()].
#' @param ... Ignored.
#' @return Never returns; always aborts.
#' @export
cast_necessity <- function(...) {
  cli::cli_abort("`cast_necessity()` was removed in 0.12.0; use `cast_effect_table()` / `cast_effect_map()`.")
}
