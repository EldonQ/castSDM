#' Select Variables for Species Distribution Models
#'
#' Two-stage, literature-standard, zero-knob variable selection:
#' \enumerate{
#'   \item \strong{Stage 1 — collinearity thinning}: rank predictors by
#'     absolute marginal association with the response, then greedily keep
#'     predictors whose pairwise correlation with all kept predictors is
#'     \eqn{\le 0.7} (Dormann et al. 2013). This is the recipe used across
#'     conventional SDM pipelines (N-SDM covsel Stage 1, correlation
#'     filtering in biomod2/wallace), implemented natively.
#'   \item \strong{Stage 2 — importance filter}: fit a probability random
#'     forest on the stage-1 survivors and keep predictors with permutation
#'     importance > 0. Importance \eqn{\le 0} means the predictor does not
#'     improve the fitted model; the zero point is the natural, data-driven
#'     threshold (no tuning).
#' }
#' Variable selection serves parsimony and projection robustness, not causal
#' attribution. Causal questions belong to [cast_effect_table()] and
#' [cast_effect_map()].
#'
#' @param data Data frame with response and predictors (coordinates allowed;
#'   they are never selected).
#' @param response Binary response column. Default `"presence"`.
#' @param method `"two_stage"` (default) or `"full"` (keep every predictor).
#' @param num_trees Trees for the stage-2 forest. Default 300.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_select` object: `selected` (kept predictors), `scores`
#'   (per-predictor association, stage-1 status, permutation importance and
#'   `selected` flag).
#'
#' @references
#' Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal
#' with it in ecological studies. \emph{Ecography} 36: 27–46.
#' @export
cast_select <- function(data, response = "presence",
                        method = c("two_stage", "full"),
                        num_trees = 300L, seed = NULL, verbose = TRUE, ...) {
  if (is.null(method)) cli::cli_abort("{.arg method} must be specified explicitly.", i = "Pass one of {.val two_stage} or {.val full}.")
  ignored <- list(...)
  if (length(ignored)) cli::cli_warn(c(
    "Deprecated {.arg cast_select} arguments were ignored: {.val {names(ignored)}}.",
    "i" = "Selection is now the zero-knob two-stage procedure."))
  if (!is.null(method) && !method %in% c("two_stage", "full")) {
    cli::cli_warn("method {.val {method}} is retired; using {.val two_stage}.")
    method <- "two_stage"
  }
  method <- match.arg(method)
  env_vars <- get_env_vars(data, response)
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")
  num_trees <- as.integer(num_trees)

  if (identical(method, "full")) {
    scores <- data.frame(variable = env_vars, assoc = NA_real_,
                         collinear_thinned = FALSE, perm_importance = NA_real_,
                         selected = TRUE)
    return(new_cast_select(selected = env_vars, scores = scores,
                           method = "full", diagnostics = list()))
  }

  # ---- Stage 1: marginal ranking + greedy pairwise |r| <= 0.7 --------------
  assoc <- vapply(env_vars, function(v)
    abs(suppressWarnings(stats::cor(data[[v]], data[[response]],
        use = "pairwise.complete.obs"))), numeric(1))
  assoc[!is.finite(assoc)] <- 0
  ord <- names(sort(assoc, decreasing = TRUE))
  thinned <- stats::setNames(rep(FALSE, length(env_vars)), env_vars)
  kept <- character(0)
  for (v in ord) {
    if (!length(kept) || all(abs(vapply(kept, function(k)
        suppressWarnings(stats::cor(data[[v]], data[[k]],
            use = "pairwise.complete.obs")), numeric(1))) <= 0.7)) {
      kept <- c(kept, v)
    } else thinned[[v]] <- TRUE
  }
  if (verbose) cli::cli_inform("Stage 1: {length(env_vars)} -> {length(kept)} after collinearity thinning.")

  # ---- Stage 2: permutation importance > 0 on the survivors ----------------
  if (length(kept) >= 2L) {
    if (!is.null(seed)) set.seed(seed)
    fit <- ranger::ranger(
      x = data[, kept, drop = FALSE],
      y = factor(data[[response]]),
      probability = TRUE, num.trees = num_trees,
      importance = "permutation", seed = seed %||% 1L, num.threads = 1L)
    imp <- fit$variable.importance[kept]
    imp[!is.finite(imp)] <- 0
    selected <- kept[imp > 0]
  } else selected <- kept
  if (verbose) cli::cli_inform("Stage 2: {length(kept)} -> {length(selected)} after importance filter.")
  if (!length(selected)) {
    cli::cli_warn("No predictor carried positive permutation importance; keeping the stage-1 set.")
    selected <- kept
  }

  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    collinear_thinned = unname(thinned[env_vars]),
    perm_importance = as.numeric(NA),
    selected = env_vars %in% selected)
  scores$perm_importance[match(kept, scores$variable)] <-
    unname(if (exists("imp", inherits = FALSE)) imp[kept] else NA_real_)
  scores$selected[match(kept, scores$variable) & scores$collinear_thinned] <- FALSE

  new_cast_select(selected = selected, scores = scores, method = "two_stage",
                  diagnostics = list(stage1_kept = kept, num_trees = num_trees))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
