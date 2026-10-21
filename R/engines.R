#' Learner Engine Registry
#'
#' Single source of truth for the per-engine facts that used to be scattered
#' across `if (nm == "rf")` chains: the backing Suggests package, the colour
#' used in figures, the tunable hyperparameter names consumed by
#' `.cast_tuned_args()`, and the display label.
#'
#' Adding a fifth learner previously required edits in four files and eleven
#' dispatch sites (`model-fit.R`, `model-evaluate.R`, `cast.R`,
#' `methods-plot.R`); a forgotten site produced a silent NA model or an `NA`
#' plot colour rather than an error. A new engine is now one entry here plus
#' its fit/predict functions.
#'
#' @format A named list; names are the internal engine keys (`"rf"`,
#'   `"maxent"`, `"gam"`). Each entry carries:
#'   \describe{
#'     \item{package}{Backing \pkg{Suggests} package.}
#'     \item{label}{Human-readable name for messages and figures.}
#'     \item{color}{Hex colour used in plots.}
#'     \item{tune_args}{Named character vector mapping a field of
#'       `fit$models[[x]]$tune$best` to the corresponding `cast_fit()`
#'       formal, so a full-data refit can reuse the tuned value.}
#'   }
#' @keywords internal
#' @noRd
.cast_engines <- function() {
list(
    rf = list(
      package = "ranger",
      label   = "Random forest",
      color   = "#4DBBD5",
      # names(tune$best) -> cast_fit() formal
      tune_args = c(mtry = "rf_mtry")
    ),
    maxent = list(
      package = "maxnet",
      label   = "MaxEnt",
      color   = "#B09C85",
      tune_args = c(classes = "maxent_classes", regmult = "maxent_regmult")
    ),
    gam = list(
      package = "mgcv",
      label   = "GAM",
      color   = "#00A087",
      tune_args = c()
    )
  )
}


#' Engine Keys
#'
#' @keywords internal
#' @noRd
.cast_engine_names <- function() names(.cast_engines())


#' Backing Package for an Engine
#'
#' @param name Engine key.
#' @return Package name, or `NULL` for an unknown engine.
#' @keywords internal
#' @noRd
.cast_engine_package <- function(name) {
  eng <- .cast_engines()[[name]]
  if (is.null(eng)) NULL else eng$package
}


#' Plot Colours for the Known Engines
#'
#' Returns one colour per requested engine. An engine with no registered
#' colour falls back to the next colour in the qualitative palette instead of
#' `NA`, so adding a learner cannot silently produce an invisible series.
#'
#' @param mdl_names Character vector of engine keys.
#' @return Named character vector of hex colours.
#' @keywords internal
#' @noRd
.cast_engine_colors <- function(mdl_names) {
  eng <- .cast_engines()
  fallback <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  out <- character(length(mdl_names))
  names(out) <- mdl_names
  for (i in seq_along(mdl_names)) {
    m <- mdl_names[i]
    col <- if (!is.null(eng[[m]])) eng[[m]]$color else NA_character_
    if (is.na(col)) col <- fallback[[((i - 1L) %% length(fallback)) + 1L]]
    out[[i]] <- col
  }
  out
}