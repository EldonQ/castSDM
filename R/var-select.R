# Variable selection ----------------------------------------------------------
#
# Two prespecified modes:
#   "covsel" -- the two-step embedded covariate selection of Adde et al.
#     (2023): hierarchical collinearity filtering (Step A) followed by
#     embedded regularization with GLM/GAM/guided-RF (Step B). Runs the
#     vendored reference implementation (R/covsel-core.R) inside a local
#     seed; no GitHub-only dependency.
#   "full"  -- keep every varying predictor (no data-driven screening).
#
# The screen is always fitted on the supplied data; cast_cv() re-runs the
# same call inside every outer training fold, so held-out rows never
# influence predictor choice.

#' Select Variables for Species Distribution Models
#'
#' Screen predictor variables before model fitting. `method = "covsel"`
#' implements the two-step embedded covariate selection procedure of Adde
#' et al. (2023) through the vendored reference implementation of the
#' \pkg{covsel} package (Step A: hierarchical collinearity filtering; Step
#' B: elastic-net GLM, null-space-penalized GAM and guided regularized
#' random forest, aggregated into one ranking and truncated to `ncov`).
#' `method = "full"` keeps every varying predictor and is the honest choice
#' for a small prespecified set.
#'
#' Selection serves parsimony for interpretation and projection; it does not
#' identify a causal adjustment set.
#'
#' @param data Data frame with response and predictors (coordinates allowed;
#'   they are never selected).
#' @param response Binary response column. Default `"presence"`.
#' @param method `"covsel"` (default) or `"full"` (keep every predictor).
#' @param corcut Absolute Pearson correlation threshold for the collinearity
#'   filter levels. Default `0.7` (Dormann et al. 2013).
#' @param var_id Character vector aligned with the predictors giving the
#'   variable-level identity (e.g. all years and focal scales of one
#'   covariate share an id). At the variable level the filter keeps exactly
#'   one covariate per id. `NULL` skips that level.
#' @param category Character vector aligned with the predictors giving the
#'   category-level identity (e.g. "climate", "terrain"). `NULL` skips it.
#' @param algorithms Embedding learners for `"covsel"`: any of `"glm"`
#'   (elastic net), `"gam"` (null-space penalisation), `"rf"` (guided
#'   regularised random forest). Default all three. Requires \pkg{glmnet},
#'   \pkg{mgcv} and \pkg{RRF} respectively.
#' @param ncov Target number of covariates for `"covsel"`. Default `NULL` =
#'   `ceiling(log2(n_presences))` (the covsel default), capped at `maxncov`.
#' @param maxncov Maximum number of covariates for `"covsel"`. Default `12`.
#' @param weights Optional numeric row weights (e.g. prevalence targets)
#'   used by the `"covsel"` filter and embedding steps.
#' @param ensure Optional named list of logical flag vectors aligned with
#'   the predictors, declaring coverage requirements the data-driven search
#'   cannot guarantee. For each flag (e.g. `list(climate = c(TRUE, ...))`),
#'   if the final set contains no flagged covariate, the strongest-ranked
#'   flagged survivor is admitted (reason `"ensured"`). This is the generic
#'   mechanism behind "a scenario-projection model must contain at least one
#'   time-varying driver": the flag comes from covariate provenance, never a
#'   hardcoded variable name.
#' @param env_vars Optional explicit predictor set (a contract, as in
#'   [cast_prepare()]). When supplied, selection runs on `env_vars`
#'   intersected with what the data actually supports; `var_id`, `category`
#'   and `keep` are re-indexed to that intersection by name. `NULL` keeps
#'   auto-detection.
#' @param keep Character vector of predictors specified before screening,
#'   such as an exposure and a scientifically justified adjustment set.
#'   These are forced through the collinearity filter and the embedding
#'   (`force`): expert knowledge outranks the data-driven search. Candidate
#'   covariates must still be scientifically admissible (not colliders or
#'   mediators of a total effect); retention does not verify a sufficient
#'   adjustment set.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_select` object with `selected` (kept predictors, in
#'   admission order), `scores` (per-predictor marginal `assoc`, the
#'   univariate `stage1_p` ranking signal, the `collinear_thinned` flag, the
#'   admission `step_added` (`0` for prespecified predictors, `NA` for
#'   never-admitted ones; the covsel aggregate rank when embedded), the
#'   `loss_gain` (always `NA` for covsel: an aggregate ranking has no loss
#'   path), the `selected_reason` (`"prespecified"`, `"embed"`, `"filter"`,
#'   `"ensured"` or `"excluded"`), the `kept_by_design` indicator for
#'   `keep`, and the `selected` flag), and `diagnostics` (filtered set,
#'   per-algorithm ranks, final aggregate ranking, coverage admissions and
#'   effective settings).
#'
#' @references
#' Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal
#' with it in ecological studies. \emph{Ecography} 36: 27-46.
#'
#' Adde, A. et al. (2023). Too many candidates: embedded covariate selection
#' procedure for species distribution modelling with the covsel R package.
#' \emph{Ecological Informatics} 75: 102080.
#' \doi{10.1016/j.ecoinf.2023.102080}.
#'
#' @seealso [cast_effect_table()], [cast_importance()]
#' @export
cast_select <- function(data, response = "presence",
                        method = c("covsel", "full"),
                        seed = NULL, verbose = TRUE, keep = character(0),
                        corcut = 0.7, var_id = NULL, category = NULL,
                        algorithms = c("glm", "gam", "rf"),
                        ncov = NULL, maxncov = 12L, weights = NULL,
                        ensure = NULL, env_vars = NULL) {
  if (!missing(method) && (!is.character(method) || length(method) != 1L ||
                          is.na(method) ||
                          !method %in% c("covsel", "full"))) {
    cli::cli_abort("{.arg method} must be 'covsel' or 'full'.")
  }
  method <- match.arg(method)
  if (is.null(env_vars)) {
    env_vars <- get_env_vars(data, response)
  } else {
    # Validate the caller's set before dropping fold-constant columns. Explicit
    # names may override metadata auto-detection, but never response/coordinates.
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
  if (method == "covsel" && length(env_vars) < 3L) {
    cli::cli_abort("Screening needs at least three predictors; use method = 'full' for a prespecified smaller set.")
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
    # Named vectors are re-indexed to the supported predictors; unnamed vectors
    # must already match the (possibly intersected) predictor set in order.
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
      # Named flags are re-indexed to the supported predictors; unnamed flags
      # must already match the (possibly intersected) predictor set in order.
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
  if (!is.null(seed)) {
    seed <- .cast_check_scalar(seed, "seed", integer = TRUE)
  }

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

  .cast_covsel_select(data, response, env_vars, keep, corcut,
                      var_id, category, algorithms, ncov, maxncov,
                      weights, ensure, seed, verbose)
}

# ---- internal helpers -----------------------------------------------------

#' Univariate quadratic-logistic signal for ranking
#'
#' Smaller p-value of the two `poly(x, 2)` terms in `y ~ poly(x, 2)`.
#' Returns `NA` when the fit fails. Used for the per-predictor ranking
#' column and for the `ensure` coverage rule.
#' @keywords internal
#' @noRd
.cast_stage1_p <- function(data, response, v) {
  tryCatch({
    x <- suppressWarnings(as.numeric(data[[v]]))
    y <- data[[response]]
    ok <- is.finite(x) & !is.na(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) < 10L || stats::sd(x) <= 0 || length(unique(y)) < 2L) return(NA_real_)
    fit <- suppressWarnings(stats::glm(y ~ stats::poly(x, degree = 2),
                                       family = stats::binomial()))
    coefs <- suppressWarnings(summary(fit)$coefficients)
    if (is.null(coefs) || nrow(coefs) < 3L) return(NA_real_)
    p <- suppressWarnings(min(coefs[2:3, 4], na.rm = TRUE))
    if (!is.finite(p)) NA_real_ else p
  }, error = function(e) NA_real_)
}

#' covsel screen: vendored Step A + Step B with castSDM extensions
#'
#' Runs the vendored reference implementation in a token space (`v1`...`vp`
#' column names), so covariates with dots, spaces or other non-syntactic
#' characters survive the formula construction and prefix stripping intact;
#' tokens map back to the original names afterwards. Extensions on top of
#' the published procedure: row `weights`, forced `keep` covariates, the
#' `ensure` coverage rule, and the diagnostic `scores` schema.
#'
#' @keywords internal
#' @noRd
.cast_covsel_select <- function(data, response, env_vars, keep, corcut,
                                var_id, category, algorithms, ncov, maxncov,
                                weights, ensure, seed, verbose) {
  for (pkg in c("glmnet", "RRF", "mgcv")) {
    check_suggested(pkg, "for the covsel embedded selection")
  }
  algorithms <- tolower(algorithms)
  bad <- setdiff(algorithms, c("glm", "gam", "rf"))
  if (length(bad) || !length(algorithms)) {
    cli::cli_abort("{.arg algorithms} must contain one or more of 'glm', 'gam', 'rf'.")
  }
  pa <- as.integer(data[[response]])
  w <- if (is.null(weights)) NULL else as.numeric(weights)
  X <- as.data.frame(data[, env_vars, drop = FALSE], check.names = FALSE)
  for (col in names(X)) X[[col]] <- as.numeric(X[[col]])

  assoc <- vapply(env_vars, function(v)
    abs(suppressWarnings(stats::cor(data[[v]], data[[response]],
        use = "pairwise.complete.obs"))), numeric(1))
  assoc[!is.finite(assoc)] <- 0
  stage1_p <- vapply(env_vars, function(v) .cast_stage1_p(data, response, v),
                     numeric(1))

  # Token space: the vendored algorithm manipulates column names with
  # formula construction and prefix stripping, so non-syntactic or
  # dot-containing predictor names must not leak into it.
  tok <- stats::setNames(paste0("v", seq_along(env_vars)), env_vars)
  Xt <- X
  names(Xt) <- unname(tok[env_vars])
  keep_tok <- unname(tok[keep])

  # ---- Step A: hierarchical collinearity filter ---------------------------
  varg <- if (is.null(var_id)) NULL else unname(var_id)
  catg <- if (is.null(category)) NULL else unname(category)
  filtered_tok <- tryCatch(
    colnames(.covsel_filter(Xt, pa, weights = w, force = keep_tok,
                            corcut = corcut, variables = varg,
                            categories = catg, verbose = verbose)),
    error = function(e) {
      cli::cli_abort(c(
        "covsel collinearity filtering (Step A) failed.",
        "x" = conditionMessage(e)
      ))
    }
  )
  if (!length(filtered_tok) || anyNA(match(filtered_tok, unname(tok)))) {
    cli::cli_abort("covsel Step A returned no usable covariates; inspect the candidate set.")
  }
  filtered <- env_vars[match(filtered_tok, unname(tok[env_vars]))]
  if (verbose) {
    cli::cli_inform("covsel Step A: {length(env_vars)} -> {length(filtered)} non-collinear covariate{?s}.")
  }

  # ---- Step B: embedded selection -----------------------------------------
  n_pres <- sum(pa == 1L)
  ncov_use <- if (is.null(ncov)) max(1L, ceiling(log2(n_pres))) else as.integer(ncov)
  ncov_use <- as.integer(min(ncov_use, maxncov))
  Xf <- Xt[, filtered_tok, drop = FALSE]
  embed <- .cast_with_seed(
    seed,
    .covsel_embed(Xf, pa, weights = w,
                  force = intersect(keep_tok, filtered_tok),
                  algorithms = algorithms, ncov = ncov_use,
                  maxncov = maxncov, nthreads = 1L, verbose = verbose)
  )
  if (is.null(embed)) {
    cli::cli_warn(c(
      "No covariate survived the embedded selection; keeping the filtered set.",
      i = "The seeds are unchanged; this is a property of the data and algorithms."
    ))
    selected <- filtered
    reason <- stats::setNames(rep("filter", length(selected)), selected)
    ranks_2 <- data.frame(covariate = unname(tok[filtered]),
                          rank.f = seq_along(filtered),
                          stringsAsFactors = FALSE)
  } else {
    sel_tok <- embed$selected
    selected <- env_vars[match(sel_tok, unname(tok[env_vars]))]
    reason <- stats::setNames(rep("embed", length(selected)), selected)
    ranks_2 <- embed$ranks_2
  }
  reason[keep] <- "prespecified"

  # ---- Coverage requirements (ensure) --------------------------------------
  # A declared flag (e.g. climate, from covariate provenance) that the final
  # set does not cover is topped up with its strongest-ranked survivor.
  # The rule names no variable: the admitted covariate is whatever the data
  # ranks strongest, so the mechanism is species-agnostic while the outcome
  # adapts. Coverage guarantees the projection is answerable; it is not a
  # claim that the ensured covariate matters (importance reporting decides).
  ensured <- character(0)
  if (!is.null(ensure)) {
    for (nm in names(ensure)) {
      flag <- stats::setNames(ensure[[nm]], env_vars)
      covered <- any(flag[union(keep, selected)])
      if (covered) next
      cand <- intersect(filtered, env_vars[flag])
      if (!length(cand)) cand <- env_vars[flag]
      if (!length(cand)) {
        cli::cli_warn("ensure[[{.val {nm}}]] matches no predictor; coverage cannot be satisfied.")
        next
      }
      p <- stage1_p[match(cand, env_vars)]
      pick <- cand[order(is.na(p), p, cand)][1]
      selected <- c(selected, pick)
      ensured <- c(ensured, pick)
      reason[pick] <- "ensured"
      if (verbose) {
        cli::cli_inform("covsel ensure: no {.val {nm}} covariate admitted; adding the strongest-ranked survivor {.val {pick}}.")
      }
    }
  }
  selected <- union(keep, selected)

  # ---- Diagnostics and scores ---------------------------------------------
  orig_of_rank <- env_vars[match(ranks_2$covariate, unname(tok[env_vars]))]
  final_ranks <- data.frame(
    covariate = ifelse(is.na(orig_of_rank), ranks_2$covariate, orig_of_rank),
    rank.f = ranks_2$rank.f,
    stringsAsFactors = FALSE
  )
  step_added <- stats::setNames(rep(NA_integer_, length(env_vars)), env_vars)
  rk <- suppressWarnings(as.integer(ranks_2$rank.f))
  names(rk) <- orig_of_rank
  rk <- rk[!is.na(names(rk))]
  common <- intersect(names(rk), env_vars)
  step_added[common] <- rk[common]
  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    stage1_p = unname(stage1_p[env_vars]),
    stage1_rank = NA_integer_,
    collinear_thinned = !(env_vars %in% filtered) & !(env_vars %in% keep),
    step_added = unname(step_added[env_vars]),
    loss_gain = NA_real_,
    kept_by_design = env_vars %in% keep,
    selected_reason = unname(ifelse(
      env_vars %in% keep, "prespecified",
      ifelse(env_vars %in% selected, unname(reason[env_vars]), "excluded"))),
    selected = env_vars %in% selected,
    stringsAsFactors = FALSE
  )
  scores$selected_reason[is.na(scores$selected_reason)] <- "excluded"
  new_cast_select(selected = selected, scores = scores, method = "covsel",
                  diagnostics = list(
                    keep = keep, corcut = corcut, algorithms = algorithms,
                    ncov = ncov_use, maxncov = maxncov, ensured = ensured,
                    filtered = filtered,
                    embed_ranks = if (is.null(embed)) NULL else {
                      er <- embed$ranks_1
                      er$covariate <- env_vars[match(er$covariate,
                                                     unname(tok[env_vars]))]
                      er
                    },
                    final_ranks = final_ranks,
                    n_pres = n_pres))
}
