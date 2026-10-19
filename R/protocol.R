# Analysis protocol ----------------------------------------------------------
#
# Single declaration point for the scientific settings that drive cast_run().
# Stage functions keep their standalone arguments (they remain usable on
# their own), but there is exactly ONE place where a pipeline default is
# declared; cast_run() validates the protocol, passes explicit values to
# every stage, echoes the effective values into the returned object and,
# when output_dir is given, into the run manifest.

#' Build a castSDM Analysis Protocol
#'
#' One validated object carrying every scientific setting of the pipeline:
#' engines and their hyperparameters, variable-selection rules, spatial
#' cross-validation design, threshold and ensemble policy, extrapolation
#' handling and reproducibility settings. Pass it to [cast_run()].
#'
#' Defaults are deliberate protocol choices, not recommendations for every
#' data set: in particular `prevalence_target = 0.5` calibrates every engine
#' (including MaxEnt, whose fit accepts no row weights, via a stored
#' logit-scale offset) to one common scale before ensembling, and
#' `clamp = TRUE` caps predictors to their training range at prediction
#' time, so all engines extrapolate under the same rule. MESS and AOA
#' diagnostics are reported alongside either way.
#'
#' @param models Engines to fit: any of `"rf"`, `"brt"`, `"maxent"`,
#'   `"gam"`. Default all four.
#' @param train_fraction Hold-out fraction rule: fraction of rows assigned
#'   to the training split (spatial block split). Default `0.7`.
#' @param select_method `"covsel"` (embedded selection, Adde et al. 2023) or
#'   `"full"` (keep every varying predictor). Default `"covsel"`.
#' @param select_keep Predictors forced through screening (see
#'   [cast_select()]).
#' @param select_corcut Absolute Pearson correlation threshold for the
#'   collinearity filter. Default `0.7`.
#' @param select_algorithms Embedding learners for `"covsel"`. Default all
#'   three (`"glm"`, `"gam"`, `"rf"`).
#' @param select_ncov Target covariate count for `"covsel"` (`NULL` =
#'   `ceiling(log2(n_presences))`).
#' @param select_maxncov Maximum covariate count for `"covsel"`. Default `12`.
#' @param rf_ntree RF trees. Default `500`.
#' @param rf_mtry RF `mtry` (`NULL` keeps the ranger default `sqrt(p)`).
#' @param brt_n_trees,brt_depth,brt_shrinkage BRT rounds, depth and learning
#'   rate. Defaults `2000`, `5`, `0.005`; depth and rate are replaced by the
#'   tuning grid when `tune = TRUE`.
#' @param maxent_classes,maxent_regmult MaxEnt feature classes and
#'   regularisation multiplier (`NULL` = maxnet defaults).
#' @param threshold_method Binary threshold rule (see [cast_threshold()]).
#'   Calculated on pooled out-of-fold predictions. Default `"max_tss"`.
#' @param tune Run the per-engine grid search inside every outer training
#'   fold. Default `FALSE`.
#' @param tune_folds Inner folds for the grid search. Default `3`.
#' @param prevalence_target Numeric in (0, 1), or `NULL` for uncalibrated
#'   scores. Default `0.5`.
#' @param cv_k Number of outer spatial folds. Default `5`.
#' @param n_repeat Number of independent spatial-blocking repeats. Default `1`.
#' @param block_method `"blockcv"` (default; blockCV spatial blocks with
#'   autocorrelation-range sizing), `"grid"` or `"env"` (see [cast_cv()]).
#' @param block_size Optional block edge length in metres for `"blockcv"`.
#' @param buffer,buffer_unit Fold-boundary exclusion band for training rows
#'   (see [cast_cv()]). Default `0`.
#' @param aoa Calibrate and report the area of applicability. Default `TRUE`.
#' @param ensemble_method,ensemble_min_score,ensemble_decay,ensemble_fallback,ensemble_min_metric_folds
#'   Ensemble policy, passed to [cast_ensemble()]. Defaults: weighted,
#'   `0.5`, `1`, `"best"`, `2`.
#' @param clamp Clamp predictors to the training range at prediction time.
#'   Default `TRUE` (applied uniformly to every engine; MESS/AOA flags are
#'   computed on the unclamped input).
#' @param num_threads Threads for the learners. Default `1`.
#' @param seed Base seed; every stage applies it locally without touching
#'   the caller's RNG state.
#'
#' @return A `cast_protocol` object (a validated list). Use `print()` to see
#'   the effective values.
#' @seealso [cast_run()]
#' @export
cast_protocol <- function(
  models = c("rf", "brt", "maxent", "gam"),
  train_fraction = 0.7,
  select_method = c("covsel", "full"),
  select_keep = character(0),
  select_corcut = 0.7,
  select_algorithms = c("glm", "gam", "rf"),
  select_ncov = NULL,
  select_maxncov = 12L,
  rf_ntree = 500L,
  rf_mtry = NULL,
  brt_n_trees = 2000L,
  brt_depth = 5L,
  brt_shrinkage = 0.005,
  maxent_classes = NULL,
  maxent_regmult = NULL,
  threshold_method = "max_tss",
  tune = FALSE,
  tune_folds = 3L,
  prevalence_target = 0.5,
  cv_k = 5L,
  n_repeat = 1L,
  block_method = c("blockcv", "grid", "env"),
  block_size = NULL,
  buffer = 0,
  buffer_unit = c("deg", "km"),
  aoa = TRUE,
  ensemble_method = c("weighted", "best", "equal", "median", "committee"),
  ensemble_min_score = 0.5,
  ensemble_decay = 1,
  ensemble_fallback = c("best", "equal", "error"),
  ensemble_min_metric_folds = 2L,
  clamp = TRUE,
  num_threads = 1L,
  seed = NULL
) {
  # Multi-value defaults document the allowed choices; the first element is
  # the effective default (match.arg semantics).
  first <- function(v) if (length(v) > 1L) v[[1L]] else v
  p <- list(
    models = models,
    train_fraction = train_fraction,
    select_method = first(select_method),
    select_keep = select_keep,
    select_corcut = select_corcut,
    select_algorithms = select_algorithms,
    select_ncov = select_ncov,
    select_maxncov = select_maxncov,
    rf_ntree = rf_ntree,
    rf_mtry = rf_mtry,
    brt_n_trees = brt_n_trees,
    brt_depth = brt_depth,
    brt_shrinkage = brt_shrinkage,
    maxent_classes = maxent_classes,
    maxent_regmult = maxent_regmult,
    threshold_method = threshold_method,
    tune = tune,
    tune_folds = tune_folds,
    prevalence_target = prevalence_target,
    cv_k = cv_k,
    n_repeat = n_repeat,
    block_method = first(block_method),
    block_size = block_size,
    buffer = buffer,
    buffer_unit = first(buffer_unit),
    aoa = aoa,
    ensemble_method = first(ensemble_method),
    ensemble_min_score = ensemble_min_score,
    ensemble_decay = ensemble_decay,
    ensemble_fallback = first(ensemble_fallback),
    ensemble_min_metric_folds = ensemble_min_metric_folds,
    clamp = clamp,
    num_threads = num_threads,
    seed = seed
  )
  class(p) <- "cast_protocol"
  .cast_validate_protocol(p)
}

#' Validate a Cast Protocol Object
#'
#' Normalises enum-valued fields and enforces every documented range, so a
#' malformed protocol fails once, at declaration time, instead of surfacing
#' deep inside a stage.
#'
#' @param p A `cast_protocol` list.
#' @return The validated protocol.
#' @keywords internal
#' @noRd
.cast_validate_protocol <- function(p) {
  if (!is.list(p)) {
    cli::cli_abort("{.arg protocol} must be a {.cls cast_protocol} object from {.fun cast_protocol}.")
  }
  required <- c(
    "models", "train_fraction", "select_method", "select_keep",
    "select_corcut", "select_algorithms", "select_ncov", "select_maxncov",
    "rf_ntree", "rf_mtry", "brt_n_trees", "brt_depth", "brt_shrinkage",
    "maxent_classes", "maxent_regmult", "threshold_method", "tune",
    "tune_folds", "prevalence_target", "cv_k", "n_repeat", "block_method",
    "block_size", "buffer", "buffer_unit", "aoa", "ensemble_method",
    "ensemble_min_score", "ensemble_decay", "ensemble_fallback",
    "ensemble_min_metric_folds", "clamp", "num_threads", "seed"
  )
  missing <- setdiff(required, names(p))
  if (length(missing)) {
    cli::cli_abort(c(
      "{.arg protocol} is missing field{?s}: {.val {missing}}.",
      i = "Build one with {.fun cast_protocol}."
    ))
  }
  p$models <- tolower(as.character(p$models))
  known_engines <- .cast_engine_names()
  bad <- setdiff(p$models, known_engines)
  if (length(bad) || !length(p$models) || anyNA(p$models)) {
    cli::cli_abort("{.arg protocol$models} must name engines from {.val {known_engines}}.")
  }
  p$train_fraction <- .cast_check_scalar(
    p$train_fraction, "protocol$train_fraction", lower = 0, upper = 1,
    min_inclusive = FALSE, max_inclusive = FALSE)
  enum_value <- function(value, field, allowed) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !value %in% allowed) {
      cli::cli_abort("{.arg protocol${field}} must be one of {.val {allowed}}.")
    }
    value
  }
  p$select_method <- enum_value(p$select_method, "select_method",
                                c("covsel", "full"))
  if (!is.character(p$select_keep) || anyNA(p$select_keep) ||
      anyDuplicated(p$select_keep)) {
    cli::cli_abort("{.arg protocol$select_keep} must be unique predictor names.")
  }
  p$select_corcut <- .cast_check_scalar(p$select_corcut,
    "protocol$select_corcut", lower = 0, upper = 1)
  p$select_algorithms <- tolower(as.character(p$select_algorithms))
  bad_alg <- setdiff(p$select_algorithms, c("glm", "gam", "rf"))
  if (length(bad_alg) || !length(p$select_algorithms)) {
    cli::cli_abort("{.arg protocol$select_algorithms} must contain one or more of 'glm', 'gam', 'rf'.")
  }
  if (!is.null(p$select_ncov)) {
    p$select_ncov <- .cast_check_scalar(p$select_ncov,
      "protocol$select_ncov", lower = 1, integer = TRUE)
  }
  p$select_maxncov <- .cast_check_scalar(p$select_maxncov,
    "protocol$select_maxncov", lower = 1, integer = TRUE)
  p$rf_ntree <- .cast_check_scalar(p$rf_ntree, "protocol$rf_ntree",
    lower = 1, integer = TRUE)
  if (!is.null(p$rf_mtry)) {
    p$rf_mtry <- .cast_check_scalar(p$rf_mtry, "protocol$rf_mtry",
      lower = 1, integer = TRUE)
  }
  p$brt_n_trees <- .cast_check_scalar(p$brt_n_trees, "protocol$brt_n_trees",
    lower = 1, integer = TRUE)
  p$brt_depth <- .cast_check_scalar(p$brt_depth, "protocol$brt_depth",
    lower = 1, integer = TRUE)
  p$brt_shrinkage <- .cast_check_scalar(p$brt_shrinkage,
    "protocol$brt_shrinkage", lower = 0, upper = 1, min_inclusive = FALSE)
  if (!is.null(p$maxent_classes) &&
      (!is.character(p$maxent_classes) || length(p$maxent_classes) != 1L ||
       !grepl("^[lqpht]+$", p$maxent_classes))) {
    cli::cli_abort("{.arg protocol$maxent_classes} must be a feature-class string like 'lqh'.")
  }
  if (!is.null(p$maxent_regmult)) {
    p$maxent_regmult <- .cast_check_scalar(p$maxent_regmult,
      "protocol$maxent_regmult", lower = 0, min_inclusive = FALSE)
  }
  p$threshold_method <- enum_value(
    as.character(p$threshold_method), "threshold_method",
    eval(formals(cast_threshold)$method))
  if (!is.logical(p$tune) || length(p$tune) != 1L || is.na(p$tune)) {
    cli::cli_abort("{.arg protocol$tune} must be TRUE or FALSE.")
  }
  p$tune_folds <- .cast_check_scalar(p$tune_folds, "protocol$tune_folds",
    lower = 2, integer = TRUE)
  if (!is.null(p$prevalence_target)) {
    p$prevalence_target <- .cast_check_scalar(p$prevalence_target,
      "protocol$prevalence_target", lower = 0, upper = 1,
      min_inclusive = FALSE, max_inclusive = FALSE)
  }
  p$cv_k <- .cast_check_scalar(p$cv_k, "protocol$cv_k", lower = 2, integer = TRUE)
  p$n_repeat <- .cast_check_scalar(p$n_repeat, "protocol$n_repeat",
    lower = 1, integer = TRUE)
  p$block_method <- enum_value(p$block_method, "block_method",
                               c("blockcv", "grid", "env"))
  if (!is.null(p$block_size)) {
    p$block_size <- .cast_check_scalar(p$block_size, "protocol$block_size",
      lower = 0, min_inclusive = FALSE)
  }
  p$buffer <- .cast_check_scalar(p$buffer, "protocol$buffer", lower = 0)
  p$buffer_unit <- enum_value(p$buffer_unit, "buffer_unit", c("deg", "km"))
  if (!is.logical(p$aoa) || length(p$aoa) != 1L || is.na(p$aoa)) {
    cli::cli_abort("{.arg protocol$aoa} must be TRUE or FALSE.")
  }
  p$ensemble_method <- enum_value(
    as.character(p$ensemble_method), "ensemble_method",
    c("weighted", "best", "equal", "median", "committee"))
  p$ensemble_min_score <- .cast_check_scalar(p$ensemble_min_score,
    "protocol$ensemble_min_score", lower = -Inf, upper = Inf)
  p$ensemble_decay <- .cast_check_decay(p$ensemble_decay)
  p$ensemble_fallback <- enum_value(p$ensemble_fallback, "ensemble_fallback",
                                    c("best", "equal", "error"))
  p$ensemble_min_metric_folds <- .cast_check_scalar(
    p$ensemble_min_metric_folds, "protocol$ensemble_min_metric_folds",
    lower = 1, integer = TRUE)
  if (!is.logical(p$clamp) || length(p$clamp) != 1L || is.na(p$clamp)) {
    cli::cli_abort("{.arg protocol$clamp} must be TRUE or FALSE.")
  }
  p$num_threads <- .cast_check_scalar(p$num_threads, "protocol$num_threads",
    lower = 1, integer = TRUE)
  if (!is.null(p$seed)) {
    p$seed <- .cast_check_scalar(p$seed, "protocol$seed", integer = TRUE)
  }
  p
}

#' @export
print.cast_protocol <- function(x, ...) {
  fmt <- function(v) {
    if (is.null(v)) return("(none)")
    paste(format(v), collapse = ", ")
  }
  cli::cli_h1("castSDM protocol (effective values)")
  cli::cli_ul(c(
    "Engines: {fmt(x$models)}",
    "Split: train_fraction = {x$train_fraction}",
    "Selection: {x$select_method} (corcut = {x$select_corcut}; algorithms = {fmt(x$select_algorithms)}; ncov = {fmt(x$select_ncov)}; maxncov = {x$select_maxncov})",
    "Engine settings: rf_ntree = {x$rf_ntree}; rf_mtry = {fmt(x$rf_mtry)}; brt_n_trees = {x$brt_n_trees}; brt_depth = {x$brt_depth}; brt_shrinkage = {x$brt_shrinkage}; maxent_classes = {fmt(x$maxent_classes)}; maxent_regmult = {fmt(x$maxent_regmult)}",
    "Tuning: tune = {x$tune} (tune_folds = {x$tune_folds})",
    "Threshold: {x$threshold_method}; prevalence_target = {fmt(x$prevalence_target)}",
    "Spatial CV: k = {x$cv_k}; repeats = {x$n_repeat}; blocks = {x$block_method}; block_size = {fmt(x$block_size)}; buffer = {x$buffer} {x$buffer_unit}; aoa = {x$aoa}",
    "Ensemble: {x$ensemble_method} (min_score = {x$ensemble_min_score}; decay = {x$ensemble_decay}; fallback = {x$ensemble_fallback}; min_metric_folds = {x$ensemble_min_metric_folds})",
    "Prediction: clamp = {x$clamp}",
    "Reproducibility: num_threads = {x$num_threads}; seed = {fmt(x$seed)}"
  ))
  invisible(x)
}
