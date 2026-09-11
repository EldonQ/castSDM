#' Select Variables for Species Distribution Models
#'
#' Two-stage, literature-standard variable selection:
#' \enumerate{
#'   \item \strong{Stage 1 — collinearity thinning}: rank predictors by
#'     absolute marginal association with the response, then greedily keep
#'     predictors whose pairwise correlation with all kept predictors is
#'     \eqn{\le 0.7} (Dormann et al. 2013). This is the recipe used across
#'     conventional SDM pipelines (N-SDM covsel Stage 1, correlation
#'     filtering in biomod2/wallace), implemented natively.
#'   \item \strong{Stage 2 — importance filter against a permutation null}:
#'     fit a probability random forest on the stage-1 survivors, then refit
#'     it `n_perm` times on a permuted response to build the null
#'     distribution of permutation importance. Keep predictors whose
#'     importance exceeds the 95th percentile of that null (Altmann et al.
#'     2010). Permutation importance under the null is centred on zero, not
#'     bounded by it, so a bare `importance > 0` rule retains roughly half
#'     of all uninformative predictors; the null quantile is the calibrated
#'     replacement.
#' }
#' Variable selection serves parsimony and projection robustness, not causal
#' attribution. Causal questions belong to [cast_effect_table()],
#' [cast_effect_map()] and [cast_necessity()].
#'
#' @param data Data frame with response and predictors (coordinates allowed;
#'   they are never selected).
#' @param response Binary response column. Default `"presence"`.
#' @param method `"two_stage"` (default) or `"full"` (keep every predictor).
#' @param num_trees Trees per forest. Default 300.
#' @param n_perm Response permutations used to build the stage-2 null.
#' @param seed Random seed.
#' @param verbose Print progress.
#' @param ... Ignored, with a warning. Retired selection knobs land here.
#'
#' @return A `cast_select` object: `selected` (kept predictors), `scores`
#'   (per-predictor association, stage-1 status, permutation importance,
#'   permutation `p_value`, BH-adjusted `p_adjusted`, and the `selected`
#'   flag). Selection uses the raw permutation p-value, which is by
#'   construction the 95th-percentile null threshold; `p_adjusted` is
#'   reported for interpretation and is not the selection rule.
#'
#' @references
#' Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal
#' with it in ecological studies. \emph{Ecography} 36: 27–46.
#'
#' Altmann, A., Toloşi, L., Sander, O. & Lengauer, T. (2010). Permutation
#' importance: a corrected feature importance measure.
#' \emph{Bioinformatics} 26: 1340–1347.
#' @export
cast_select <- function(data, response = "presence",
                        method = c("two_stage", "full"),
                        num_trees = 300L, n_perm = 49L, seed = NULL,
                        verbose = TRUE, ...) {
  ignored <- list(...)
  if (length(ignored)) cli::cli_warn(c(
    "Deprecated {.arg cast_select} arguments were ignored: {.val {names(ignored)}}.",
    "i" = "Selection is now the two-stage procedure; see {.fn cast_select}."))
  if (is.null(method) || !length(method)) {
    cli::cli_abort(c("{.arg method} must be specified explicitly.",
                     i = "Pass one of {.val two_stage} or {.val full}."))
  }
  method <- as.character(method)[1L]
  if (!method %in% c("two_stage", "full")) {
    cli::cli_warn("method {.val {method}} is retired; using {.val two_stage}.")
    method <- "two_stage"
  }
  env_vars <- get_env_vars(data, response)
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")
  num_trees <- as.integer(num_trees)
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) {
    cli::cli_abort("{.arg n_perm} must be at least 1.")
  }

  if (identical(method, "full")) {
    scores <- data.frame(variable = env_vars, assoc = NA_real_,
                         collinear_thinned = FALSE, perm_importance = NA_real_,
                         p_value = NA_real_, p_adjusted = NA_real_,
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

  # ---- Stage 2: importance above the permutation null ----------------------
  imp <- stats::setNames(rep(NA_real_, length(kept)), kept)
  p_value <- imp
  threshold <- NA_real_
  if (length(kept) >= 2L) {
    check_suggested("ranger", "for stage-2 importance filtering")
    X <- data[, kept, drop = FALSE]
    y <- factor(data[[response]])
    imp_of <- function(y_use, s) {
      f <- ranger::ranger(x = X, y = y_use, probability = TRUE,
                          num.trees = num_trees, importance = "permutation",
                          seed = s, num.threads = 1L)
      out <- f$variable.importance[kept]
      out[!is.finite(out)] <- 0
      out
    }
    if (!is.null(seed)) set.seed(seed)
    imp <- imp_of(y, seed %||% 1L)
    if (verbose) cli::cli_inform("Stage 2: building the null from {n_perm} response permutation{?s}...")
    # Pool the null across predictors: importances share one scale, so the
    # pooled draws calibrate the threshold far better than n_perm draws per
    # predictor would at an affordable number of refits.
    null_draws <- unlist(lapply(seq_len(n_perm), function(i) {
      imp_of(sample(y), (seed %||% 1L) + i)
    }), use.names = FALSE)
    n_null <- length(null_draws)
    threshold <- unname(stats::quantile(null_draws, 0.95, names = FALSE))
    p_value <- vapply(kept, function(v)
      (1 + sum(null_draws >= imp[[v]])) / (1 + n_null), numeric(1))
    selected <- kept[p_value <= 0.05]
  } else {
    selected <- kept
  }
  p_adjusted <- if (length(kept) >= 2L) {
    stats::p.adjust(p_value, method = "BH")
  } else p_value
  if (verbose) cli::cli_inform("Stage 2: {length(kept)} -> {length(selected)} above the permutation null.")
  if (!length(selected)) {
    cli::cli_warn(c(
      "No predictor exceeded the permutation null; keeping the stage-1 set.",
      i = "Read this as weak evidence for any single predictor, not as a clean screen."))
    selected <- kept
  }

  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    collinear_thinned = unname(thinned[env_vars]),
    perm_importance = unname(imp[env_vars]),
    p_value = unname(p_value[env_vars]),
    p_adjusted = unname(p_adjusted[env_vars]),
    selected = env_vars %in% selected,
    stringsAsFactors = FALSE)

  new_cast_select(selected = selected, scores = scores, method = "two_stage",
                  diagnostics = list(stage1_kept = kept, num_trees = num_trees,
                                     n_perm = n_perm,
                                     null_quantile = 0.95,
                                     null_threshold = threshold,
                                     alpha = 0.05))
}
