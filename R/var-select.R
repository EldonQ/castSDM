# Variable selection ----------------------------------------------------------
#
# Two prespecified modes:
#   "screen" -- variance/unique-value filter, optional construct de-duplication,
#     correlation pruning, then an elastic-net logistic screen (glmnet).
#   "full"   -- keep every varying predictor (no data-driven screening).
#
# The screen is fitted on the supplied data; cast_cv() re-runs it inside every
# outer training fold, so held-out rows never influence predictor choice.

#' Select Variables for Species Distribution Models
#'
#' Screen predictors before model fitting. `method = "screen"` applies a
#' transparent sequence: low-information predictors are dropped, repeated
#' constructs can be de-duplicated with `var_id`, correlated survivors are
#' pruned at `corcut`, and the remaining set is screened with an L2-penalised
#' (elastic-net) logistic model. `method = "full"` keeps every varying
#' predictor and is the honest choice for a small prespecified set.
#'
#' Selection serves parsimony for interpretation and projection; it does not
#' identify a causal adjustment set.
#'
#' @param data Data frame with response and predictors (coordinates allowed;
#'   they are never selected).
#' @param response Binary response column. Default `"presence"`.
#' @param method `"screen"` (default) or `"full"` (keep every predictor).
#' @param corcut Absolute Pearson correlation threshold for pruning. Default
#'   `0.7` (Dormann et al. 2013).
#' @param var_id Optional character vector aligned with the predictors giving
#'   a construct identity (e.g. all years/focal scales of one covariate share
#'   an id). The strongest-ranked member of each construct is retained.
#' @param category Optional character vector aligned with predictors. Accepted
#'   for API compatibility; the screen does not use category quotas.
#' @param ncov Target number of covariates. Default `NULL` =
#'   `ceiling(log2(n_presences))`, capped at `maxncov`.
#' @param maxncov Maximum number of covariates. Default `12`.
#' @param weights Optional numeric row weights used by the screen.
#' @param ensure Optional named list of logical flag vectors aligned with the
#'   predictors, declaring coverage requirements. If the final set contains no
#'   flagged covariate, the strongest-ranked flagged survivor is admitted.
#' @param env_vars Optional explicit predictor set. `NULL` keeps auto-detection.
#' @param keep Character vector forced through the screen (expert knowledge).
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_select` object with `selected`, per-predictor `scores`, and
#'   diagnostics (filtered set, coefficients, lambda, coverage admissions).
#'
#' @references
#' Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal
#' with it in ecological studies. \emph{Ecography} 36: 27-46.
#'
#' Zou, H. & Hastie, T. (2005). Regularization and variable selection via the
#' elastic net. \emph{Journal of the Royal Statistical Society: Series B}
#' 67: 301-320.
#'
#' Friedman, J., Hastie, T. & Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent.
#' \emph{Journal of Statistical Software} 33: 1-22.
#'
#' @seealso [cast_effect_table()], [cast_importance()]
#' @export
cast_select <- function(data, response = "presence",
                        method = c("screen", "full"),
                        seed = NULL, verbose = TRUE, keep = character(0),
                        corcut = 0.7, var_id = NULL, category = NULL,
                        ncov = NULL, maxncov = 12L, weights = NULL,
                        ensure = NULL, env_vars = NULL) {
  if (!missing(method) && (!is.character(method) || length(method) != 1L ||
                          is.na(method) ||
                          !method %in% c("screen", "full"))) {
    cli::cli_abort("{.arg method} must be 'screen' or 'full'.")
  }
  method <- match.arg(method)
  if (is.null(env_vars)) {
    env_vars <- get_env_vars(data, response)
  } else {
    if (!is.character(env_vars) || anyNA(env_vars) || anyDuplicated(env_vars) ||
        !all(env_vars %in% setdiff(names(data), c(response, "lon", "lat")))) {
      cli::cli_abort("{.arg env_vars} must be unique predictor names in {.arg data}, not coordinates or the response.")
    }
    .cast_check_numeric_predictors(data[, env_vars, drop = FALSE], arg = "env_vars")
    supported <- vapply(env_vars, function(v) {
      z <- data[[v]][!is.na(data[[v]])]
      variance <- if (length(z) > 1L) stats::var(z) else NA_real_
      is.finite(variance) && variance > 1e-10
    }, logical(1))
    unsupported <- env_vars[!supported]
    if (length(unsupported)) {
      cli::cli_warn(c(
        "{length(unsupported)} variable{?s} in {.arg env_vars} lack supported variation in this sample and cannot be selected:",
        " " = "{.val {unsupported}}"))
      env_vars <- env_vars[supported]
    }
  }
  if (!is.character(keep) || anyNA(keep) || anyDuplicated(keep) ||
      !all(keep %in% env_vars)) {
    cli::cli_abort("{.arg keep} must contain unique numeric predictor names that vary in the training data, not coordinates or the response.")
  }
  .cast_check_response(data[[response]], response)
  .cast_check_numeric_predictors(data[, env_vars, drop = FALSE], arg = "data")
  if (!length(env_vars)) cli::cli_abort("Need at least one varying numeric predictor.")
  if (method == "screen" && length(env_vars) < 2L) {
    cli::cli_abort("Screening needs at least two predictors; use method = 'full' for a prespecified smaller set.")
  }
  if (!is.numeric(corcut) || length(corcut) != 1L || is.na(corcut) ||
      corcut < 0 || corcut > 1) {
    cli::cli_abort("{.arg corcut} must be a single number in [0, 1].")
  }
  align_metadata <- function(x, arg) {
    if (is.null(x)) return(NULL)
    if (!is.character(x) || anyNA(x)) {
      cli::cli_abort("{.arg {arg}} must be NULL or an NA-free character vector aligned with the predictors.")
    }
    if (!is.null(names(x))) {
      x <- x[match(env_vars, names(x))]
    } else if (length(x) != length(env_vars)) {
      cli::cli_abort("{.arg {arg}} must be NULL, named, or a character vector of length {length(env_vars)} aligned with the predictors.")
    }
    if (anyNA(x)) {
      cli::cli_abort("{.arg {arg}} does not cover all predictors: {.val {env_vars[is.na(x)]}}.")
    }
    x
  }
  var_id <- align_metadata(var_id, "var_id")
  category <- align_metadata(category, "category")
  if (!is.null(ncov) && (!is.numeric(ncov) || length(ncov) != 1L ||
                         is.na(ncov) || ncov < 1)) {
    cli::cli_abort("{.arg ncov} must be NULL or a single positive integer.")
  }
  maxncov <- as.integer(maxncov)
  if (is.na(maxncov) || maxncov < 1L) {
    cli::cli_abort("{.arg maxncov} must be a positive integer.")
  }
  if (!is.null(weights) && (!is.numeric(weights) || length(weights) != nrow(data))) {
    cli::cli_abort("{.arg weights} must be NULL or a numeric vector with one value per row of {.arg data}.")
  }
  if (!is.null(ensure)) {
    if (!is.list(ensure) || !length(ensure) || is.null(names(ensure)) ||
        any(!nzchar(names(ensure))) || anyDuplicated(names(ensure))) {
      cli::cli_abort("{.arg ensure} must be NULL or a uniquely named list of logical flag vectors.")
    }
    for (nm in names(ensure)) {
      f <- ensure[[nm]]
      if (!is.logical(f) || anyNA(f)) {
        cli::cli_abort("{.arg ensure}[[{.val {nm}}]] must be a NA-free logical vector aligned with the predictors.")
      }
      if (!is.null(names(f))) {
        f <- f[match(env_vars, names(f))]
      } else if (length(f) != length(env_vars)) {
        cli::cli_abort("{.arg ensure}[[{.val {nm}}]] must be a logical vector of length {length(env_vars)} aligned with the predictors.")
      }
      if (anyNA(f)) {
        cli::cli_abort("{.arg ensure}[[{.val {nm}}]] does not cover all predictors.")
      }
      ensure[[nm]] <- f
    }
  }
  if (!is.null(seed)) seed <- .cast_check_scalar(seed, "seed", integer = TRUE)

  if (identical(method, "full")) {
    scores <- data.frame(variable = env_vars, assoc = NA_real_,
                         stage1_p = NA_real_, stage1_rank = NA_integer_,
                         collinear_thinned = FALSE,
                         step_added = NA_integer_,
                         loss_gain = NA_real_,
                         kept_by_design = env_vars %in% keep,
                         selected_reason = "full",
                         selected = TRUE)
    return(new_cast_select(selected = env_vars, scores = scores,
                           method = "full", diagnostics = list(keep = keep)))
  }

  .cast_screen_select(data, response, env_vars, keep, corcut,
                      var_id, ncov, maxncov, weights, ensure, seed, verbose)
}

# ---- independent screen ----------------------------------------------------

#' Univariate quadratic-logistic signal for ranking
#'
#' Smaller p-value of the two `poly(x, 2)` terms in `y ~ poly(x, 2)`.
#' Returns `NA` when the fit fails.
#' @keywords internal
#' @noRd
.cast_stage1_p <- function(data, response, v) {
  tryCatch({
    x <- suppressWarnings(as.numeric(data[[v]]))
    y <- data[[response]]
    ok <- is.finite(x) & !is.na(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < 10L || stats::sd(x) <= 0 || length(unique(y)) < 2L) return(NA_real_)
    fit <- suppressWarnings(stats::glm(y ~ stats::poly(x, degree = 2),
                                       family = stats::binomial()))
    coefs <- suppressWarnings(summary(fit)$coefficients)
    if (is.null(coefs) || nrow(coefs) < 3L) return(NA_real_)
    p <- suppressWarnings(min(coefs[2:3, 4], na.rm = TRUE))
    if (!is.finite(p)) NA_real_ else p
  }, error = function(e) NA_real_)
}

#' Independent elastic-net screen
#'
#' The screen is deliberately independent of any external selection library:
#' it combines a variance/unique-value filter, optional construct de-duplication,
#' greedy correlation pruning, and an elastic-net logistic model fitted with
#' \pkg{glmnet}. The elastic-net fit is used only to rank and subset survivors;
#' final inference is performed by the fitted SDM engines.
#' @keywords internal
#' @noRd
.cast_screen_select <- function(data, response, env_vars, keep, corcut,
                                var_id, ncov, maxncov, weights, ensure,
                                seed, verbose) {
  check_suggested("glmnet", "for the elastic-net screen")
  y <- as.integer(data[[response]])
  if (!all(y %in% c(0L, 1L))) cli::cli_abort("The screen response must be binary 0/1.")
  X <- as.data.frame(data[, env_vars, drop = FALSE], check.names = FALSE)
  for (col in names(X)) X[[col]] <- as.numeric(X[[col]])

  assoc <- vapply(env_vars, function(v) {
    z <- suppressWarnings(stats::cor(X[[v]], y, use = "pairwise.complete.obs"))
    if (is.finite(z)) abs(z) else 0
  }, numeric(1))
  stage1_p <- vapply(env_vars, function(v) .cast_stage1_p(data, response, v),
                     numeric(1))
  names(assoc) <- names(stage1_p) <- env_vars
  ord <- order(is.na(stage1_p), stage1_p, -assoc, env_vars)

  unique_n <- vapply(X, function(v) length(unique(v[is.finite(v)])), integer(1))
  too_few <- names(unique_n)[unique_n < 5L]
  forced_too_few <- intersect(too_few, keep)
  if (length(forced_too_few)) {
    cli::cli_abort(c(
      "Forced covariate(s) have fewer than five unique finite values: {.val {forced_too_few}}.",
      i = "The elastic-net screen cannot rank a near-constant forced layer."
    ))
  }
  dropped <- setdiff(too_few, keep)
  survivors <- setdiff(env_vars, dropped)

  # One member per construct, strongest univariate signal first.
  kept <- character(0)
  seen_id <- character(0)
  thinned <- stats::setNames(rep(FALSE, length(env_vars)), env_vars)
  thinned[dropped] <- TRUE
  if (!is.null(var_id)) names(var_id) <- env_vars
  for (v in survivors[order(match(survivors, env_vars[ord]))]) {
    id <- if (is.null(var_id)) v else unname(var_id[[v]])
    if (id %in% seen_id) {
      thinned[[v]] <- TRUE
      next
    }
    if (length(kept) && any(abs(stats::cor(X[[v]], X[, kept, drop = FALSE],
        use = "pairwise.complete.obs")) > corcut, na.rm = TRUE)) {
      thinned[[v]] <- TRUE
      next
    }
    kept <- c(kept, v)
    seen_id <- c(seen_id, id)
  }
  kept <- union(keep, kept)
  # Forced predictors are admitted by design even when they duplicate or
  # correlate with a screened survivor; the exclusion flag must describe the
  # final candidate set, not the order of admission.
  thinned[keep] <- FALSE
  if (!length(kept)) {
    cli::cli_abort("No predictor survived the variance/correlation screen.")
  }
  n_pres <- sum(y == 1L)
  target0 <- if (is.null(ncov)) max(1L, ceiling(log2(n_pres))) else as.integer(ncov)
  target <- as.integer(min(target0, maxncov))
  cap_binding <- target0 > maxncov
  if (isTRUE(cap_binding)) {
    cli::cli_warn(c(
      "Covariate target exceeds maxncov ({maxncov}); the screen is truncated at {maxncov}.",
      i = "This is a user-set complexity guardrail, not a data-derived optimum."
    ))
  }

  selected <- kept
  coefs <- stats::setNames(rep(NA_real_, length(kept)), kept)
  lambda <- NA_real_
  if (length(kept) > 1L) {
    x <- as.matrix(X[, kept, drop = FALSE])
    storage.mode(x) <- "double"
    w <- if (is.null(weights)) NULL else as.numeric(weights) / mean(weights)
    nfolds <- max(2L, min(5L, n_pres, sum(y == 0L), 10L))
    foldid <- .cast_tune_folds(y, nfolds, seed)
    cv <- tryCatch(
      glmnet::cv.glmnet(x, y, family = "binomial", alpha = 0.5,
                        weights = w, foldid = foldid,
                        type.measure = "deviance", standardize = TRUE),
      error = function(e) NULL)
    if (!is.null(cv)) {
      b <- as.matrix(stats::coef(cv, s = "lambda.1se"))
      nz <- setdiff(rownames(b)[b[, 1L] != 0], "(Intercept)")
      if (!length(nz)) {
        b <- as.matrix(stats::coef(cv, s = "lambda.min"))
        nz <- setdiff(rownames(b)[b[, 1L] != 0], "(Intercept)")
      }
      if (length(nz)) {
        selected <- union(keep, nz)
        coefs <- stats::setNames(as.numeric(b[, 1L]), rownames(b))
        lambda <- cv$lambda.1se
      }
    }
    if (!length(setdiff(selected, keep))) {
      best <- kept[order(is.na(stage1_p[kept]), stage1_p[kept],
                         -assoc[kept])][1]
      selected <- union(keep, best)
    }
    if (length(selected) > target) {
      rank <- abs(coefs[selected])
      rank[is.na(rank)] <- -Inf
      selected <- union(keep, selected[order(-rank)][seq_len(max(0L, target - length(keep)))])
    }
  }

  # Coverage requirements. Only survivors of the variance/collinearity filter
  # can satisfy a declared class; resurrecting a thinned layer would contradict
  # the filter that is reported alongside the selection.
  ensured <- character(0)
  if (!is.null(ensure)) {
    for (nm in names(ensure)) {
      flag <- stats::setNames(ensure[[nm]], env_vars)
      if (any(flag[selected])) next
      cand <- intersect(kept, env_vars[flag])
      if (!length(cand)) {
        cli::cli_warn("ensure[[{.val {nm}}]] has no survivor after filtering; coverage cannot be satisfied.")
        next
      }
      pick <- cand[order(is.na(stage1_p[cand]), stage1_p[cand], cand)][1]
      selected <- union(selected, pick)
      ensured <- c(ensured, pick)
    }
  }
  selected <- union(keep, selected)
  kept <- union(kept, ensured)
  thinned[ensured] <- FALSE
  if (verbose) {
    cli::cli_inform(
      "Screen: {length(env_vars)} -> {length(kept)} candidates after variance/collinearity; {length(selected)} selected (target {target}, maxncov {maxncov})."
    )
  }

  step_rank <- stats::setNames(rep(NA_integer_, length(env_vars)), env_vars)
  if (any(is.finite(coefs))) {
    ord_coef <- names(sort(abs(coefs[is.finite(coefs)]), decreasing = TRUE))
    step_rank[ord_coef] <- seq_along(ord_coef)
  }
  step_rank[keep] <- 0L
  final_ranks <- data.frame(
    variable = env_vars[ord],
    stage1_rank = seq_along(env_vars),
    coefficient = unname(coefs[env_vars[ord]]),
    selected = env_vars[ord] %in% selected,
    stringsAsFactors = FALSE
  )
  reason <- stats::setNames(rep("excluded", length(env_vars)), env_vars)
  reason[selected] <- "screen"
  reason[keep] <- "prespecified"
  reason[ensured] <- "ensured"
  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    stage1_p = unname(stage1_p[env_vars]),
    stage1_rank = match(env_vars, env_vars[ord]),
    collinear_thinned = unname(thinned[env_vars]),
    step_added = unname(step_rank[env_vars]),
    loss_gain = NA_real_,
    kept_by_design = env_vars %in% keep,
    selected_reason = unname(reason[env_vars]),
    selected = env_vars %in% selected,
    stringsAsFactors = FALSE)
  new_cast_select(
    selected = selected,
    scores = scores,
    method = "screen",
    diagnostics = list(
      keep = keep, corcut = corcut, ncov_target = target0, ncov = target,
      ncov_use = length(selected), maxncov = maxncov,
      cap_binding = cap_binding, ensured = ensured, filtered = kept,
      lambda = lambda, coefficients = coefs, n_pres = n_pres,
      algorithms = "elastic_net", final_ranks = final_ranks))
}
