#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' \pkg{castSDM} explains how environmental changes alter predicted habitat
#' suitability and where those changes occur. [cast_attribution()] compares
#' current and future predictor values using frozen model weights, with maps
#' of overall change, driver-only effects and retained-driver comparisons.
#' Range coverage and model disagreement accompany the results.
#'
#' The package also provides data preparation, variable selection, RF/BRT/
#' MaxEnt/GAM fitting, nested spatial validation and raster projection. These
#' support the central environmental-attribution workflow; selection and
#' scenario comparisons are predictive diagnostics, not causal discovery.
#' Presence/background outputs remain relative habitat suitability.
#'
#' @section Quick Start:
#' ```
#' result <- cast(species_data, env_data)
#' summary(result)
#' plot(result)
#' ```
#'
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom stats predict cor sd median var na.omit lm
#' @importFrom utils head
NULL

# Suppress R CMD check NOTEs for non-standard evaluation variables
utils::globalVariables(c("self", ".data", "index"))
