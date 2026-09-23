#' Select Variables for Species Distribution Models
#'
#' Two-stage variable selection whose second stage is calibrated against a
#' \emph{conditional} permutation null:
#' \enumerate{
#'   \item \strong{Stage 1 — collinearity thinning}: rank predictors by a
#'     univariate quadratic-logistic signal (the smaller p-value of the two
#'     `poly(x, 2)` terms, which sees U-shaped responses a linear correlation
#'     misses), then greedily keep predictors whose pairwise correlation with
#'     all kept predictors is \eqn{\le 0.7} (Dormann et al. 2013). Predictors
#'     whose GLM fails to fit fall back to the marginal-correlation rank.
#'   \item \strong{Stage 2 — conditional effect above a conditional null}:
#'     fit a probability random forest on the stage-1 survivors, then measure
#'     how far each predictor moves the fitted probability when it is
#'     \strong{shifted while every other predictor is held at its observed
#'     value} (a g-computation contrast; see [cast_effect_table()]). The same
#'     statistic is recomputed on forests refitted to within-stratum
#'     permuted predictors (strata from k-means on the survivors), and
#'     predictors whose Monte Carlo tail probability is at most 0.05 are
#'     kept. The kept set is capped at `ncov` predictors (default
#'     `ceiling(log2(n_presence))`, at most `maxncov`), ordered by the
#'     conditional effect, so model complexity stays tied to the number of
#'     presences (Adde et al. 2023).
#' }
#'
#' @section Why the second stage is conditional:
#' Marginal permutation breaks a predictor's correlation
#' with every other predictor. For collinear predictors the permuted rows leave
#' the observed data support, so the score is governed by the model's
#' extrapolation behaviour rather than by the predictor's influence
#' (Hooker, Mentch & Zhou 2021). A within-stratum permutation instead shuffles
#' each predictor only among rows with similar values on the other survivors
#' (k-means strata), so the null respects the observed joint distribution.
#' Selection does not find an adjustment set;
#' use a scientifically justified set directly when estimating causal effects.
#'
#' @param data Data frame with response and predictors (coordinates allowed;
#'   they are never selected).
#' @param response Binary response column. Default `"presence"`.
#' @param method `"two_stage"` (default) or `"full"` (keep every predictor).
#' @param num_trees Trees per forest. Default 300.
#' @param n_perm Conditional permutations used to build the stage-2 null.
#' @param shift_size Shift applied to the intervened predictor, in training
#'   standard deviations. Default `1`. The statistic averages the absolute
#'   change over `+shift_size` and `-shift_size`, so it does not depend on a
#'   sign convention.
#' @param max_rows Rows sampled (evenly, preserving order) for the
#'   g-computation contrast. Bounds stage-2 cost on large data sets.
#'   Default 2000.
#' @param ncov Maximum predictors retained, ordered by the interventional
#'   effect. Default `NULL` selects `ceiling(log2(n_presence))` where
#'   `n_presence` is the number of presences. Use `Inf` for no cap.
#' @param maxncov Upper bound applied to the automatic `ncov`. Default `12`.
#' @param seed Random seed.
#' @param verbose Print progress.
#' @param keep Character vector of predictors specified before screening, such
#'   as an exposure and a scientifically justified adjustment set. These bypass
#'   collinearity thinning and the null threshold, count toward `ncov`, and
#'   cannot be removed by the cap. The cap must accommodate them. This does not
#'   identify a sufficient adjustment set; candidate covariates must also be
#'   scientifically admissible (not colliders or mediators of a total effect).
#'
#' @return A `cast_select` object with `selected` (kept predictors), `scores`
#'   (per-predictor marginal `assoc`, the `stage1_p` ranking signal with its
#'   `stage1_rank`, the stage-1 `collinear_thinned` flag,
#'   `interventional_effect` -- the stage-2 selection statistic -- with its
#'   `effect_rank`, its conditional-permutation `p_value`,
#'   the feature-wise `null_threshold`,
#'   the `selected_reason` (`"null"`, `"null+top-ncov"`,
#'   `"fallback-top-ncov"`, `"prespecified"`, `"full"` or `"excluded"`), the
#'   `kept_by_design` indicator for `keep`, and the `selected` flag).
#'   Prespecified retention is not evidence of an effect. Selection uses
#'   `(1 + sum(null >= observed)) / (1 + n_perm)` for each predictor
#'   separately. These are screening diagnostics: the
#'   response-dependent stage-1 screen is held fixed during permutation, so
#'   they do not establish confirmatory p-values or FDR control.
#'   `passed_null` distinguishes evidence from a capped or fallback set when
#'   the null is not the binding constraint. Fewer than 19 permutations
#'   cannot resolve p <= 0.05.
#'
#' @references
#' Adde, A. et al. (2023). Too many candidates: embedded covariate selection
#' procedure for species distribution modelling with the covsel R package.
#' \emph{Ecological Informatics} 75: 102080.
#'
#' Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal
#' with it in ecological studies. \emph{Ecography} 36: 27-46.
#'
#' Hooker, G., Mentch, L. & Zhou, S. (2021). Unrestricted permutation forces
#' extrapolation: variable importance requires at least one more model, or
#' there is no free variable importance. \emph{Statistics and Computing} 31: 82.
#' \doi{10.1007/s11222-021-10057-z}.
#' @seealso [cast_effect_table()], [cast_importance()]
#' @export
cast_select <- function(data, response = "presence",
                         method = c("two_stage", "full"),
                         num_trees = 300L, n_perm = 49L, shift_size = 1,
                         max_rows = 2000L, ncov = NULL, maxncov = 12L,
                         seed = NULL, verbose = TRUE, keep = character(0)) {
  if (!missing(method) && (!is.character(method) || length(method) != 1L ||
                          is.na(method) || !method %in% c("two_stage", "full"))) {
    cli::cli_abort("{.arg method} must be one of 'two_stage' or 'full'.")
  }
  method <- match.arg(method)
  env_vars <- get_env_vars(data, response)
  if (!is.character(keep) || anyNA(keep) || anyDuplicated(keep) ||
      !all(keep %in% env_vars)) {
    cli::cli_abort("{.arg keep} must contain unique numeric predictor names that vary in the training data, not coordinates or the response.")
  }
  .cast_check_response(data[[response]], response)
  .cast_check_numeric_predictors(data[, env_vars, drop = FALSE], arg = "data")
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")
  num_trees <- as.integer(num_trees)
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) {
    cli::cli_abort("{.arg n_perm} must be at least 1.")
  }
  if (!is.numeric(shift_size) || length(shift_size) != 1L ||
      !is.finite(shift_size) || shift_size <= 0) {
    cli::cli_abort("{.arg shift_size} must be one positive finite number.")
  }
  max_rows <- as.integer(max_rows)
  if (is.na(max_rows) || max_rows < 2L) {
    cli::cli_abort("{.arg max_rows} must be at least 2.")
  }
  maxncov <- suppressWarnings(as.integer(maxncov))
  if (length(maxncov) != 1L || is.na(maxncov) || maxncov < 1L) {
    cli::cli_abort("{.arg maxncov} must be one positive integer.")
  }
  n_presence <- sum(data[[response]] == 1, na.rm = TRUE)
  ncov_auto <- max(1L, ceiling(log2(max(n_presence, 2L))))
  ncov_auto <- min(ncov_auto, maxncov)
  if (is.null(ncov)) {
    ncov <- ncov_auto
  } else if (length(ncov) == 1L && is.infinite(ncov) && ncov > 0) {
    ncov <- .Machine$integer.max
  } else {
    ncov <- suppressWarnings(as.integer(ncov))
    if (length(ncov) != 1L || is.na(ncov) || ncov < 1L) {
      cli::cli_abort("{.arg ncov} must be one positive integer, {.code NULL} or {.code Inf}.")
    }
  }

  if (identical(method, "full")) {
    scores <- data.frame(variable = env_vars, assoc = NA_real_,
                         stage1_p = NA_real_, stage1_rank = NA_integer_,
                         collinear_thinned = FALSE,
                         interventional_effect = NA_real_,
                         effect_rank = NA_integer_,
                         p_value = NA_real_,
                         null_threshold = NA_real_,
                         passed_null = FALSE, fallback = FALSE,
                         kept_by_design = env_vars %in% keep,
                         selected_reason = "full",
                         selected = TRUE)
    return(new_cast_select(selected = env_vars, scores = scores,
                           method = "full", diagnostics = list(keep = keep)))
  }
  if (length(keep) > ncov) {
    cli::cli_abort("{.arg ncov} must be at least length(keep); increase the cap rather than dropping prespecified predictors.")
  }

  # ---- Stage 1: quadratic-logistic ranking + greedy pairwise |r| <= 0.7 ----
  assoc <- vapply(env_vars, function(v)
    abs(suppressWarnings(stats::cor(data[[v]], data[[response]],
        use = "pairwise.complete.obs"))), numeric(1))
  assoc[!is.finite(assoc)] <- 0
  stage1_p <- vapply(env_vars, function(v) .cast_stage1_p(data, response, v),
                     numeric(1))
  varies <- vapply(env_vars, function(v) {
    z <- data[[v]][is.finite(data[[v]])]
    length(z) > 1L && stats::sd(z) > 0
  }, logical(1))
  rankable <- env_vars[varies]
  # Smaller p ranks first; predictors whose GLM failed (NA) fall back to the
  # marginal-correlation order. Ties break on the marginal association.
  ord <- rankable[order(is.na(stage1_p[rankable]),
                        stage1_p[rankable],
                        -assoc[rankable], rankable)]
  stage1_rank <- stats::setNames(rep(NA_integer_, length(env_vars)), env_vars)
  stage1_rank[ord] <- seq_along(ord)
  thinned <- stats::setNames(rep(FALSE, length(env_vars)), env_vars)
  thinned[!varies] <- TRUE
  kept <- keep
  for (v in setdiff(ord, keep)) {
    if (!length(kept) || all(abs(vapply(kept, function(k)
        suppressWarnings(stats::cor(data[[v]], data[[k]],
            use = "pairwise.complete.obs")), numeric(1))) <= 0.7)) {
      kept <- c(kept, v)
    } else thinned[[v]] <- TRUE
  }
  if (verbose) cli::cli_inform("Stage 1: {length(env_vars)} -> {length(kept)} after collinearity thinning.")

  # ---- Stage 2: conditional effect above a conditional null ---------------
  effect <- stats::setNames(rep(NA_real_, length(kept)), kept)
  p_value <- effect
  threshold <- NA_real_
  n_strata <- NA_integer_
  if (length(kept) >= 1L) {
    check_suggested("ranger", "for stage-2 conditional screening")
    X <- data[, kept, drop = FALSE]
    for (col in names(X)) X[[col]] <- as.numeric(X[[col]])
    ok <- rowSums(!is.finite(as.matrix(X))) == 0L
    X <- X[ok, , drop = FALSE]
    y <- factor(data[[response]][ok])
    if (!nrow(X) || length(unique(y)) < 2L) {
      cli::cli_abort(c(
        "Stage 2 needs complete predictor rows and two response classes.",
        "i" = "Check for missing predictor values or a single-class response."))
    }
    sds <- vapply(X, stats::sd, numeric(1))
    sds[!is.finite(sds) | sds <= 0] <- 1
    # Bound the g-computation cost; even steps keep the row distribution.
    if (nrow(X) > max_rows) {
      idx <- unique(round(seq(1, nrow(X), length.out = max_rows)))
      X <- X[idx, , drop = FALSE]
      y <- y[idx]
    }

    if (!is.null(seed)) set.seed(seed)
    base_seed <- seed %||% 1L
    fit_obs <- .cast_importance_fit(X, y, num_trees, base_seed)
    effect <- .cast_shift_effect(fit_obs$model, X, sds, shift_size)

    if (verbose) cli::cli_inform("Stage 2: building the conditional null from {n_perm} within-stratum permutation{?s}...")
    # Strata group rows with similar values on the OTHER survivors, so each
    # null replicate shuffles a predictor only among rows that match on the
    # rest. The null forests see the same joint X support as the observed
    # forest; only the conditional links to the response are broken. With one
    # stratum this degrades gracefully to a full permutation.
    strata_list <- .cast_perm_strata_list(X, base_seed)
    n_strata <- max(vapply(strata_list, function(s) length(unique(s)),
                           integer(1)))
    # Each feature has its own null scale; unrelated predictors are not
    # exchangeable null replicates for this predictor.
    null_draws <- vapply(seq_len(n_perm), function(i) {
      X_p <- .cast_stratum_permute_list(X, strata_list)
      m <- ranger::ranger(x = X_p, y = y, probability = TRUE,
                          num.trees = num_trees, seed = base_seed + i,
                          num.threads = 1L, write.forest = TRUE)
      .cast_shift_effect(m, X_p, sds, shift_size)
    }, numeric(length(kept)))
    if (is.null(dim(null_draws))) {
      null_draws <- matrix(null_draws, nrow = length(kept))
    }
    rownames(null_draws) <- kept
    threshold <- apply(null_draws, 1, stats::quantile, probs = 0.95, names = FALSE)
    p_value <- vapply(kept, function(v)
      (1 + sum(null_draws[v, ] >= effect[[v]])) / (1 + n_perm), numeric(1))
  }
  effect_rank <- stats::setNames(rep(NA_integer_, length(kept)), kept)
  effect_rank[names(sort(effect, decreasing = TRUE))] <- seq_along(kept)
  cap <- min(ncov, length(kept))
  pass <- kept[is.finite(p_value) & p_value <= 0.05]
  pass <- pass[order(-effect[pass])]
  optional_pass <- setdiff(pass, keep)
  optional_cap <- cap - length(keep)
  capped <- length(optional_pass) > optional_cap
  fallback <- !length(pass)
  reason <- stats::setNames(rep("excluded", length(kept)), kept)
  if (!length(kept)) {
    selected <- kept
    reason <- stats::setNames(character(0), character(0))
    capped <- FALSE
    fallback <- TRUE
    cli::cli_warn("No varying predictor survived stage 1; returning an empty set.")
  } else if (!length(pass)) {
    optional <- setdiff(names(sort(effect, decreasing = TRUE)), keep)
    selected <- c(keep, utils::head(optional, optional_cap))
    reason[selected] <- "fallback-top-ncov"
    cli::cli_warn(c(
      "No predictor exceeded the conditional null; retaining {length(keep)} prespecified and keeping the top {optional_cap} optional predictors.",
      "i" = "Read this as weak evidence for any single predictor, not as a clean screen."))
  } else {
    selected <- c(keep, utils::head(optional_pass, optional_cap))
    reason[selected] <- if (capped) "null+top-ncov" else "null"
    if (verbose) cli::cli_inform("Stage 2: retaining {length(keep)} prespecified and {length(selected) - length(keep)} null-screened predictors (cap = {cap}).")
  }
  reason[keep] <- "prespecified"
  passed_null <- pass

  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    stage1_p = unname(stage1_p[env_vars]),
    stage1_rank = unname(stage1_rank[env_vars]),
    collinear_thinned = unname(thinned[env_vars]),
    interventional_effect = unname(effect[env_vars]),
    effect_rank = unname(effect_rank[env_vars]),
    p_value = unname(p_value[env_vars]),
    null_threshold = unname(threshold[match(env_vars, names(threshold))]),
    passed_null = env_vars %in% passed_null,
    fallback = fallback & env_vars %in% selected & !env_vars %in% keep,
    kept_by_design = env_vars %in% keep,
    selected_reason = ifelse(env_vars %in% selected,
                             unname(reason[env_vars]), "excluded"),
    selected = env_vars %in% selected,
    stringsAsFactors = FALSE)
  scores$selected_reason[scores$selected_reason %in% c(NA, "NA")] <- "excluded"

  new_cast_select(selected = selected, scores = scores, method = "two_stage",
                  diagnostics = list(
                    stage1_kept = kept, stage1_metric = "poly2-glm-min-p",
                    num_trees = num_trees,
                    n_perm = n_perm, shift_size = shift_size,
                    max_rows = max_rows, null_quantile = 0.95,
                    null_threshold = threshold,
                    null_method = paste("feature-wise within-stratum",
                                        "permutation of the shift effect"),
                    null_strata = n_strata,
                    statistic = "interventional_effect",
                    n_presence = n_presence, ncov = cap, maxncov = maxncov,
                    keep = keep, capped = capped,
                    fallback = fallback, alpha = 0.05))
}

# ---- internal helpers -----------------------------------------------------

#' Univariate quadratic-logistic signal for stage-1 ranking
#'
#' Smaller p-value of the two `poly(x, 2)` terms in `y ~ poly(x, 2)`.
#' Returns `NA` when the fit fails so the caller can fall back to the
#' marginal-correlation order.
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

#' Fit a probability forest (no permutation importance)
#'
#' Stage 2 reports a single conditional-effect statistic, so the forest is
#' fit without the permutation-importance bookkeeping.
#' @keywords internal
#' @noRd
.cast_importance_fit <- function(X, y, num_trees, seed) {
  m <- ranger::ranger(x = X, y = y, probability = TRUE,
                      num.trees = num_trees, importance = "none",
                      seed = seed, num.threads = 1L)
  list(model = m, importance = stats::setNames(rep(NA_real_, ncol(X)), colnames(X)))
}

#' Strata for the conditional permutation null
#'
#' One stratification per survivor, each built from the OTHER survivors, so a
#' predictor's own signal cannot leak into its null through the clustering.
#' Few rows collapse to a single stratum (a full permutation); k-means
#' failure also falls back to one stratum rather than aborting. A lone
#' survivor has nothing to condition on and permutes fully.
#' @keywords internal
#' @noRd
.cast_perm_strata_list <- function(X, seed, k = 5L) {
  n <- nrow(X)
  out <- lapply(names(X), function(v) {
    others <- setdiff(names(X), v)
    if (!length(others)) return(rep(1L, n))
    .cast_perm_strata(X[, others, drop = FALSE], seed)
  })
  stats::setNames(out, names(X))
}

#' @keywords internal
#' @noRd
.cast_perm_strata <- function(X, seed, k = 5L) {
  n <- nrow(X)
  k_use <- min(as.integer(k), n)
  if (k_use < 2L || !ncol(X)) return(rep(1L, n))
  Xs <- scale(as.matrix(X))
  Xs[!is.finite(Xs)] <- 0
  km <- tryCatch(suppressWarnings(stats::kmeans(Xs, centers = k_use, nstart = 10L)),
                 error = function(e) NULL)
  if (is.null(km)) return(rep(1L, n))
  as.integer(km$cluster)
}

#' Permute each predictor independently within its own strata
#'
#' Single-row strata are left untouched (`sample()` on one value would sample
#' from `1:x` instead of permuting).
#' @keywords internal
#' @noRd
.cast_stratum_permute_list <- function(X, strata_list) {
  X_p <- X
  for (v in names(X_p)) {
    strata <- strata_list[[v]]
    for (s in unique(strata)) {
      idx <- which(strata == s)
      if (length(idx) < 2L) next
      X_p[idx, v] <- sample(X_p[idx, v])
    }
  }
  X_p
}

#' Interventional effect of shifting each predictor by +/- shift_size SD
#'
#' One g-computation contrast per predictor: shift it, hold every other
#' predictor at its observed value, and average the absolute change in the
#' fitted probability over both shift directions.
#' @keywords internal
#' @noRd
.cast_shift_effect <- function(model, X, sds, shift_size) {
  base <- stats::predict(model, data = X)$predictions[, "1"]
  out <- stats::setNames(numeric(ncol(X)), colnames(X))
  for (v in colnames(X)) {
    step <- shift_size * sds[[v]]
    total <- 0
    for (sgn in c(1, -1)) {
      Xs <- X
      Xs[[v]] <- Xs[[v]] + sgn * step
      cf <- stats::predict(model, data = Xs)$predictions[, "1"]
      d <- abs(cf - base)
      d[!is.finite(d)] <- NA_real_
      total <- total + mean(d, na.rm = TRUE)
    }
    out[[v]] <- total / 2
  }
  out
}
