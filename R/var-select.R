#' Select Variables for Species Distribution Models
#'
#' Two-stage variable selection whose second stage is calibrated against an
#' \emph{interventional} null rather than a permutation-importance null:
#' \enumerate{
#'   \item \strong{Stage 1 — collinearity thinning}: rank predictors by a
#'     univariate quadratic-logistic signal (the smaller p-value of the two
#'     `poly(x, 2)` terms, which sees U-shaped responses a linear correlation
#'     misses), then greedily keep predictors whose pairwise correlation with
#'     all kept predictors is \eqn{\le 0.7} (Dormann et al. 2013). Predictors
#'     whose GLM fails to fit fall back to the marginal-correlation rank.
#'   \item \strong{Stage 2 — interventional effect above a permutation null}:
#'     fit a probability random forest on the stage-1 survivors, then measure
#'     how far each predictor moves the fitted probability when it is
#'     \strong{shifted while every other predictor is held at its observed
#'     value} (a g-computation contrast; see [cast_effect_table()]). The same
#'     statistic is recomputed on forests refitted to a permuted response to
#'     build a feature-wise null, and predictors whose Monte Carlo tail
#'     probability is at most 0.05 are kept (Altmann et al. 2010). The kept
#'     set is capped at `ncov` predictors (default
#'     `ceiling(log2(n_presence))`, at most `maxncov`), ordered by the
#'     interventional effect, so model complexity stays tied to the number of
#'     presences (Adde et al. 2023).
#' }
#'
#' @section Why the second stage is interventional:
#' Permutation importance permutes a predictor, which breaks its correlation
#' with every other predictor. For collinear predictors the permuted rows leave
#' the observed data support, so the score is governed by the model's
#' extrapolation behaviour rather than by the predictor's influence
#' (Hooker, Mentch & Zhou 2021). An additive shift also can leave the joint
#' data support; holding other predictors fixed does not prevent extrapolation.
#' Its benefit is an explicitly defined contrast, not guaranteed causal
#' identification or superiority. Selection does not find an adjustment set;
#' use a scientifically justified set directly when estimating causal effects.
#'
#' @section Optional invariant selection:
#' `method = "tramicp"` delegates binary-logistic invariant selection to
#' `tramicp::glmICP()`. It requires scientifically justified environment labels,
#' a correctly specified response model, adequate measured causes/confounders,
#' and the method's sampling assumptions. Spatial CV does not establish these
#' assumptions. The resulting intersection is neither a verified list of
#' ecological causes nor a sufficient adjustment set. Presence/background
#' labels do not become true occurrence observations through this method.
#'
#' @param data Data frame with response and predictors (coordinates allowed;
#'   they are never selected).
#' @param response Binary response column. Default `"presence"`.
#' @param method `"two_stage"` (default), `"full"` (keep every predictor), or
#'   `"tramicp"` (opt-in full-subset binomial invariant causal prediction).
#' @param environment For `"tramicp"`, the name of a supplied environment column
#'   in `data`, never a vector or automatically generated CV grouping. Excluded
#'   from predictors; must have at least two levels after shared complete-row
#'   filtering. Not supported for other methods.
#' @param icp_alpha Invariance-test level for `"tramicp"`. Sets are accepted only
#'   when their p-value is strictly greater than this value. Default `0.05`.
#' @param icp_max_predictors Computational guard for the full `2^p` subset search.
#'   Default `12L`; exceeding it stops, without prescreening or top-K truncation.
#'   `"tramicp"` rejects nonempty `keep` and non-NULL `ncov`, and never augments or
#'   caps the accepted-set intersection. Scores report selection, not effects;
#'   diagnostics retain all subset tests and distinguish empty-result statuses.
#' @param num_trees Trees per forest. Default 300.
#' @param n_perm Response permutations used to build the stage-2 null.
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
#' @return A `cast_select` object. For `"tramicp"`, `selected` is the unaltered
#'   accepted-set intersection; `scores` contains `variable`, `selected` and
#'   `selected_reason`. Diagnostics retain `tested_sets`, `accepted_sets`,
#'   `intersection`, retained/excluded row indices and `status`: `"selected"`,
#'   `"empty_set_accepted"`, `"empty_intersection"` or `"no_accepted_sets"`.
#'   Invalid subset tests raise an error rather than imply nonselection.
#'   For the two-stage method: `selected` (kept predictors), `scores`
#'   (per-predictor marginal `assoc`, the `stage1_p` ranking signal with its
#'   `stage1_rank`, the stage-1 `collinear_thinned` flag,
#'   `interventional_effect` -- the stage-2 selection statistic -- with its
#'   `effect_rank`, its permutation `p_value`, the BH-adjusted `p_adjusted`,
#'   the feature-wise `null_threshold`, the retained `perm_importance`
#'   diagnostic, the `selected_reason` (`"null"`, `"null+top-ncov"`,
#'   `"fallback-top-ncov"`, `"prespecified"`, `"full"` or `"excluded"`), the
#'   `kept_by_design` indicator for `keep`, and the `selected` flag).
#'   Prespecified retention is not evidence of an effect. Optional selection uses
#'   `(1 + sum(null >= observed)) / (1 + n_perm)` for each predictor
#'   separately. `p_adjusted` adjusts over stage-1 survivors only and is not
#'   the selection rule. These are screening diagnostics: the
#'   response-dependent stage-1 screen is held fixed during permutation, so
#'   they do not establish confirmatory p-values or FDR control.
#'   `passed_null` distinguishes evidence from a capped or fallback set when
#'   the null is not the binding constraint. Fewer than 19 permutations
#'   cannot resolve p <= 0.05.
#'
#' @references
#' Kook, L. et al. (2024). Model-based causal feature selection for general
#' response types. \emph{Journal of the American Statistical Association}.
#' \doi{10.1080/01621459.2024.2395588}.
#'
#' Adde, A. et al. (2023). Too many candidates: embedded covariate selection
#' procedure for species distribution modelling with the covsel R package.
#' \emph{Ecological Informatics} 75: 102080.
#'
#' Altmann, A., Toloşi, L., Sander, O. & Lengauer, T. (2010). Permutation
#' importance: a corrected feature importance measure.
#' \emph{Bioinformatics} 26: 1340-1347.
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
                         method = c("two_stage", "full", "tramicp"),
                         num_trees = 300L, n_perm = 49L, shift_size = 1,
                         max_rows = 2000L, ncov = NULL, maxncov = 12L,
                         seed = NULL, verbose = TRUE, keep = character(0),
                         environment = NULL, icp_alpha = 0.05,
                         icp_max_predictors = 12L) {
  if (!missing(method) && (!is.character(method) || length(method) != 1L ||
                          is.na(method) || !method %in% c("two_stage", "full", "tramicp"))) {
    cli::cli_abort("{.arg method} must be one of 'two_stage', 'full' or 'tramicp'.")
  }
  method <- match.arg(method)
  if (identical(method, "tramicp")) {
    return(.cast_select_tramicp(data, response, environment, icp_alpha,
                               icp_max_predictors, keep, ncov, seed, verbose))
  }
  if (!is.null(environment)) {
    cli::cli_abort("{.arg environment} is only supported for method = 'tramicp'; remove environment metadata from data for other methods.")
  }
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
                         perm_importance = NA_real_,
                         p_value = NA_real_, p_adjusted = NA_real_,
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

  # ---- Stage 2: interventional effect above the permutation null -----------
  imp <- stats::setNames(rep(NA_real_, length(kept)), kept)
  effect <- imp
  p_value <- imp
  threshold <- NA_real_
  if (length(kept) >= 1L) {
    check_suggested("ranger", "for stage-2 interventional screening")
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
    imp <- fit_obs$importance
    effect <- .cast_shift_effect(fit_obs$model, X, sds, shift_size)

    if (verbose) cli::cli_inform("Stage 2: building the interventional null from {n_perm} response permutation{?s}...")
    # Each feature has its own null scale; unrelated predictors are not
    # exchangeable null replicates for this predictor. The null forests carry
    # no importance, which the statistic does not use.
    null_draws <- vapply(seq_len(n_perm), function(i) {
      y_p <- sample(y)
      m <- ranger::ranger(x = X, y = y_p, probability = TRUE,
                          num.trees = num_trees, seed = base_seed + i,
                          num.threads = 1L, write.forest = TRUE)
      .cast_shift_effect(m, X, sds, shift_size)
    }, numeric(length(kept)))
    if (is.null(dim(null_draws))) {
      null_draws <- matrix(null_draws, nrow = length(kept))
    }
    rownames(null_draws) <- kept
    threshold <- apply(null_draws, 1, stats::quantile, probs = 0.95, names = FALSE)
    p_value <- vapply(kept, function(v)
      (1 + sum(null_draws[v, ] >= effect[[v]])) / (1 + n_perm), numeric(1))
  }
  p_adjusted <- if (length(kept) >= 2L) {
    stats::p.adjust(p_value, method = "BH")
  } else p_value
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
      "No predictor exceeded the interventional null; retaining {length(keep)} prespecified and keeping the top {optional_cap} optional predictors.",
      "i" = "Read this as weak evidence for any single predictor, not as a clean screen."))
  } else {
    selected <- c(keep, utils::head(optional_pass, optional_cap))
    reason[selected] <- if (capped) "null+top-ncov" else "null"
    if (verbose) cli::cli_inform("Stage 2: retaining {length(keep)} prespecified and {length(selected) - length(keep)} null-screened predictors (cap = {cap}).")
  }
  reason[keep] <- "prespecified"
  passed_null <- pass

  agreement <- if (length(kept) >= 3L &&
                   all(is.finite(imp)) && all(is.finite(effect))) {
    suppressWarnings(stats::cor(imp, effect, method = "spearman"))
  } else NA_real_

  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    stage1_p = unname(stage1_p[env_vars]),
    stage1_rank = unname(stage1_rank[env_vars]),
    collinear_thinned = unname(thinned[env_vars]),
    interventional_effect = unname(effect[env_vars]),
    effect_rank = unname(effect_rank[env_vars]),
    perm_importance = unname(imp[env_vars]),
    p_value = unname(p_value[env_vars]),
    p_adjusted = unname(p_adjusted[env_vars]),
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
                    null_method = paste("feature-wise response permutation of",
                                        "the interventional effect"),
                    statistic = "interventional_effect",
                    importance_agreement = agreement,
                    n_presence = n_presence, ncov = cap, maxncov = maxncov,
                    keep = keep, capped = capped,
                    fallback = fallback, alpha = 0.05))
}

# ---- internal helpers -----------------------------------------------------

# Full-subset ICP: no response-dependent screening or post-selection cap.
.cast_select_tramicp <- function(data, response, environment, icp_alpha,
                                 icp_max_predictors, keep, ncov, seed, verbose) {
  if (!is.data.frame(data) || is.null(names(data)) || anyNA(names(data)) ||
      any(!nzchar(trimws(names(data)))) || anyDuplicated(names(data))) {
    cli::cli_abort("{.arg data} must be a data frame with unique, nonempty column names.")
  }
  if (!is.character(response) || length(response) != 1L || is.na(response) ||
      !response %in% names(data)) {
    cli::cli_abort("{.arg response} must name one column in data.")
  }
  if (!is.character(environment) || length(environment) != 1L ||
      is.na(environment) || !environment %in% names(data) ||
      environment == response) {
    cli::cli_abort("{.arg environment} must name one supplied column in data, distinct from the response; vectors and automatic CV environments are not supported.")
  }
  if (length(keep)) {
    cli::cli_abort("Nonempty {.arg keep} is not supported for tramicp: the invariant intersection must not be altered.")
  }
  if (!is.null(ncov)) {
    cli::cli_abort("{.arg ncov} must be NULL for tramicp: the invariant intersection must not be capped.")
  }
  if (!is.numeric(icp_alpha) || length(icp_alpha) != 1L ||
      !is.finite(icp_alpha) || icp_alpha <= 0 || icp_alpha >= 1) {
    cli::cli_abort("{.arg icp_alpha} must be one finite number strictly between 0 and 1.")
  }
  if (!is.numeric(icp_max_predictors) || length(icp_max_predictors) != 1L ||
      !is.finite(icp_max_predictors) || icp_max_predictors < 1 ||
      icp_max_predictors != floor(icp_max_predictors)) {
    cli::cli_abort("{.arg icp_max_predictors} must be one positive finite integer.")
  }
  # Match the package's coordinate/metadata exclusions, but do not use
  # get_env_vars(): it silently screens nonnumeric and low-variance columns.
  metadata <- c("lon", "lat", "HID", "species", "sid", "family", "category",
                "fraction", "id", "ID", "site", "cell_id", "grid_id", "group",
                "spid", "siteid", "occ", "fold")
  predictors <- setdiff(names(data), c(response, environment, metadata))
  if (!length(predictors)) cli::cli_abort("tramicp needs at least one numeric predictor.")
  if (length(predictors) > icp_max_predictors) {
    cli::cli_abort("Full subset search needs 2^{length(predictors)} tests, exceeding {.arg icp_max_predictors} = {icp_max_predictors}. Supply a scientifically prespecified smaller candidate set or explicitly increase the computational budget; no top-K truncation is performed.")
  }
  X <- data[, predictors, drop = FALSE]
  if (!all(vapply(X, function(x) is.numeric(x) && is.null(dim(x)), logical(1)))) {
    cli::cli_abort("tramicp requires numeric predictor columns (not factors, characters or matrices).")
  }
  y <- data[[response]]
  if ((!is.numeric(y) && !is.logical(y)) || !is.null(dim(y)) ||
      !all(y[!is.na(y)] %in% c(0, 1))) {
    cli::cli_abort("tramicp requires a binary 0/1 response.")
  }
  e <- data[[environment]]
  if ((!is.factor(e) && !is.character(e) && !is.numeric(e) && !is.logical(e)) ||
      !is.null(dim(e)) || (is.numeric(e) && any(!is.finite(e) & !is.na(e)))) {
    cli::cli_abort("The environment column must contain finite categorical labels.")
  }
  ok <- stats::complete.cases(X, y, e)
  X <- X[ok, , drop = FALSE]
  y <- as.numeric(y[ok])
  e <- droplevels(factor(e[ok]))
  if (!nrow(X) || length(unique(y)) != 2L) {
    cli::cli_abort("tramicp needs complete rows with two response classes after shared complete-row filtering.")
  }
  if (any(!is.finite(as.matrix(X)))) {
    cli::cli_abort("tramicp requires finite predictors after shared complete-row filtering.")
  }
  if (nlevels(e) < 2L) {
    cli::cli_abort("tramicp needs at least two environment levels after shared complete-row filtering.")
  }
  aliases <- paste0("X", seq_along(predictors))
  names(X) <- aliases
  X$Y <- y
  X$E <- e
  if (!is.null(seed)) set.seed(seed)
  backend <- tryCatch(
    .cast_run_tramicp(X, stats::reformulate(aliases, response = "Y"),
                      icp_alpha, verbose),
    error = function(e) cli::cli_abort("tramicp backend failed: {conditionMessage(e)}", parent = e))
  parsed <- .cast_parse_tramicp(backend, aliases, predictors, icp_alpha)
  selected <- parsed$intersection
  scores <- data.frame(variable = predictors, selected = predictors %in% selected,
                       selected_reason = ifelse(predictors %in% selected,
                         "invariant_intersection", "not_identified"),
                       stringsAsFactors = FALSE)
  new_cast_select(selected, scores, method = "tramicp", diagnostics = c(parsed,
    list(environment = environment, row_ids = which(ok),
         excluded_row_ids = which(!ok), environment_counts = table(e, dnn = environment),
         alpha = icp_alpha, icp_max_predictors = icp_max_predictors,
         alias_map = stats::setNames(predictors, aliases),
         backend_version = tryCatch(as.character(utils::packageVersion("tramicp")),
                                    error = function(e) NA_character_),
         backend_result = backend,
         interpretation = paste(
           "Invariant-intersection identification assumes a correctly specified logit model,",
           "valid environments, no hidden confounding and iid observations.",
           "It establishes neither adjustment sufficiency nor ecological causality",
           "from presence/background labels; subset p-values are not effect importance."))))
}

# Keep dependency calls isolated so the wrapper contract can be tested without
# installing tramicp. Resolve its unexported binary-GLM residual method exactly.
.cast_run_tramicp <- function(data, formula, icp_alpha, verbose) {
  check_suggested("tramicp", "for invariant causal prediction")
  controls <- tramicp::dicp_controls(
    type = "residual", test = "gcm.test", alpha = icp_alpha,
    residuals = utils::getFromNamespace("residuals.binglm", "tramicp"),
    crossfit = FALSE, stop_if_empty_set_invariant = FALSE)
  tramicp::glmICP(formula = formula, data = data, env = ~ E,
                  family = stats::binomial(link = "logit"), verbose = verbose,
                  type = "residual", test = "gcm.test", controls = controls,
                  alpha = icp_alpha, greedy = FALSE, max_size = NULL,
                  mandatory = NULL)
}

.cast_parse_tramicp <- function(backend, aliases, predictors, alpha) {
  fail <- function() cli::cli_abort(paste(
    "tramicp returned failed, nonfinite or incomplete subset tests;",
    "every one of the 2^p subsets must have a valid p-value."))
  if (!is.list(backend) || !is.list(backend$tests) ||
      length(backend$tests) != 2^length(aliases)) fail()
  subsets <- lapply(backend$tests, function(x) {
    if (!is.list(x) || inherits(x, "try-error") || !is.character(x$set)) fail()
    s <- x$set
    if (identical(s, "Empty")) s <- character(0)
    if (anyNA(s) || anyDuplicated(s) || !all(s %in% aliases)) fail()
    aliases[aliases %in% s]
  })
  keys <- vapply(subsets, function(s) if (!length(s)) "Empty" else
    paste(s, collapse = "+"), character(1))
  # Unique valid subsets, with cardinality 2^p, prove full enumeration.
  if (anyDuplicated(keys)) fail()
  pvalues <- vapply(backend$tests, function(x) {
    if (!is.list(x$test) || inherits(x$test, "try-error")) fail()
    p <- x$test$p.value
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 0 || p > 1) fail()
    unname(p)
  }, numeric(1))
  sp <- backend$set_pvals
  if (!is.numeric(sp) || length(sp) != length(keys) || is.null(names(sp)) ||
      anyNA(names(sp)) || anyDuplicated(names(sp)) || any(!is.finite(sp))) fail()
  # Backend labels use their original term order, not necessarily X1,...,Xp.
  labels <- vapply(backend$tests, function(x) paste(x$set, collapse = "+"), character(1))
  if (!setequal(names(sp), labels) || !isTRUE(all.equal(
      unname(sp[labels]), unname(pvalues), tolerance = 0))) fail()
  subsets <- lapply(subsets, function(s) predictors[match(s, aliases)])
  accepted <- pvalues > alpha
  accepted_sets <- subsets[accepted]
  intersection <- if (length(accepted_sets)) Reduce(intersect, accepted_sets) else character(0)
  status <- if (!length(accepted_sets)) "no_accepted_sets" else
    if (any(lengths(accepted_sets) == 0L)) "empty_set_accepted" else
      if (!length(intersection)) "empty_intersection" else "selected"
  tested_sets <- data.frame(p_value = unname(pvalues), accepted = unname(accepted))
  tested_sets$set <- I(unname(subsets))
  list(status = status, accepted_sets = unname(accepted_sets),
       tested_sets = tested_sets, intersection = intersection)
}

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

#' Fit a probability forest and its permutation importance
#' @keywords internal
#' @noRd
.cast_importance_fit <- function(X, y, num_trees, seed) {
  m <- ranger::ranger(x = X, y = y, probability = TRUE,
                      num.trees = num_trees, importance = "permutation",
                      seed = seed, num.threads = 1L)
  out <- m$variable.importance[colnames(X)]
  out[!is.finite(out)] <- 0
  list(model = m, importance = out)
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
