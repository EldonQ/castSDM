#' Select Variables for Species Distribution Models
#'
#' `method = "cpi"` (default) is the castSDM causal selector. It tests each
#' predictor's conditional independence from occurrence given every other
#' predictor via the conditional predictive impact (Watson & Wright 2021):
#' a random forest is cross-fitted once, each predictor is replaced by a
#' Gaussian knockoff, and the increase in log-loss measures its conditional
#' contribution. Because the learner is a random forest, symmetric unimodal
#' (bell-shaped) niche responses are detected. Predictors surviving
#' Benjamini-Hochberg FDR control are kept.
#'
#' `method = "dml"` is an alternative causal selector estimating each
#' predictor's Neyman-orthogonal partially linear effect with
#' [DoubleML][DoubleML::DoubleMLPLR]; it assumes an approximately linear
#' conditional effect and is best paired with [cast_effect()] for a signed,
#' interpretable coefficient. `method = "rf"` is retained as a conventional
#' (associational) permutation-importance benchmark for comparison.
#'
#' @param data Data frame with response, coordinates, and predictors.
#' @param response Binary response column.
#' @param method `"cpi"` (default), `"dml"`, or the conventional `"rf"`
#'   benchmark.
#' @param alpha FDR level for the CPI/DML selectors. Default `0.05`.
#' @param max_candidates Predictors tested with CPI/DML; larger sets are
#'   pre-screened by RF importance for feasibility. Also the output ceiling for
#'   the RF benchmark. Default `30`.
#' @param dml_folds Cross-fitting folds for the CPI/DML selectors. Default `5`.
#' @param num_trees Trees for the RF nuisance/benchmark forests. Default `300`.
#'   The CPI and DML selectors cap this at `200` with an informational message
#'   (the nuisance forest keeps fitting tractable on large candidate sets).
#' @param min_vars Minimum retained variables. Default `3`.
#' @param force_include Character vector of predictor names to always retain,
#'   regardless of the selector's decision. Use for known niche-defining axes
#'   (e.g. temperature, elevation) that a within-calibration conditional screen
#'   may drop as redundant but that are essential for spatial transfer; the
#'   selector still runs normally on every other predictor. Names absent from
#'   the predictor pool are skipped. Default none.
#' @param cor_threshold Absolute-correlation threshold for the RF benchmark.
#' @param num_threads Threads for the ranger learners/benchmark. Default `1`.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_select` object. Its `scores` data.frame marks each
#'   predictor with `selected` (retained), `fallback` (kept only through the
#'   `min_vars` floor, i.e. it did not pass FDR), and `forced` (added via
#'   `force_include`).
#' @seealso [cast_effect()], [cast_counterfactual()]
#' @export
cast_select <- function(data,
                        response = "presence",
                        method = c("cpi", "dml", "rf"),
                        alpha = 0.05,
                        max_candidates = 30L,
                        dml_folds = 5L,
                        num_trees = 300L,
                        min_vars = 3L,
                        force_include = character(),
                        cor_threshold = 0.8,
                        num_threads = 1L,
                        seed = NULL,
                        verbose = TRUE) {
  method <- match.arg(method)
  env_vars <- get_env_vars(data, response)
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")

  num_trees <- as.integer(num_trees)
  if (method %in% c("cpi", "dml") && num_trees > 200L) {
    if (verbose) {
      cli::cli_inform(
        "CPI/DML: capping {.arg num_trees} at 200 (nuisance-forest tractability)."
      )
    }
    num_trees <- 200L
  }

  if (identical(method, "cpi")) {
    out <- .cast_select_cpi(
      data = data, env_vars = env_vars, response = response,
      alpha = alpha, max_candidates = max_candidates, n_folds = dml_folds,
      num_trees = num_trees, min_vars = min_vars,
      seed = seed, verbose = verbose, num_threads = num_threads
    )
    res <- .cast_apply_force_include(out$selected, out$scores, env_vars,
                                     force_include, verbose)
    return(new_cast_select(
      selected = res$selected, scores = res$scores, method = method,
      diagnostics = c(out$diagnostics,
                      list(forced = intersect(force_include, env_vars)))
    ))
  }

  if (identical(method, "dml")) {
    out <- .cast_select_dml(
      data = data, env_vars = env_vars, response = response,
      alpha = alpha, max_candidates = max_candidates, n_folds = dml_folds,
      num_trees = num_trees, min_vars = min_vars,
      seed = seed, verbose = verbose, num_threads = num_threads
    )
    res <- .cast_apply_force_include(out$selected, out$scores, env_vars,
                                     force_include, verbose)
    return(new_cast_select(
      selected = res$selected, scores = res$scores, method = method,
      diagnostics = c(out$diagnostics,
                      list(forced = intersect(force_include, env_vars)))
    ))
  }

  # -- Conventional RF permutation-importance benchmark ---------------------
  check_suggested("ranger", "for RF permutation selection")
  x <- as.data.frame(.cast_numeric_matrix(data, env_vars), check.names = FALSE)
  names(x) <- make.names(env_vars, unique = TRUE)
  x$.response <- factor(data[[response]])
  rf <- ranger::ranger(
    .response ~ ., data = x, importance = "permutation",
    num.trees = as.integer(num_trees), seed = seed %||% 42L,
    num.threads = as.integer(num_threads), verbose = FALSE
  )
  imp <- rf$variable.importance
  original <- env_vars[match(names(imp), make.names(env_vars, unique = TRUE))]
  imp <- stats::setNames(as.numeric(imp), original)
  imp <- imp[env_vars]
  imp[!is.finite(imp) | is.na(imp)] <- 0
  max_vars <- min(max(as.integer(min_vars), as.integer(max_candidates)),
                  length(env_vars))
  ranked <- env_vars[order(imp, decreasing = TRUE)]

  selected <- character()
  for (v in ranked) {
    if (length(selected) >= max_vars) break
    if (!length(selected)) {
      selected <- v
    } else {
      cors <- abs(stats::cor(
        .cast_numeric_matrix(data, v),
        .cast_numeric_matrix(data, selected), use = "pairwise.complete.obs"
      ))
      if (!any(cors >= cor_threshold, na.rm = TRUE)) selected <- c(selected, v)
    }
  }
  loop_selected <- selected
  if (length(selected) < min_vars) {
    selected <- unique(c(selected, ranked))[seq_len(min(min_vars, length(ranked)))]
  }
  fallback_added <- setdiff(selected, loop_selected)
  scores <- data.frame(
    variable = env_vars,
    rf_importance = unname(imp),
    selected = env_vars %in% selected,
    fallback = env_vars %in% fallback_added,
    forced = FALSE,
    exclusion_reason = ifelse(env_vars %in% selected, "selected", "lower_or_redundant"),
    stringsAsFactors = FALSE
  )
  scores <- scores[order(!scores$selected, -scores$rf_importance), , drop = FALSE]
  res <- .cast_apply_force_include(selected, scores, env_vars,
                                   force_include, verbose)
  new_cast_select(
    selected = res$selected, scores = res$scores, method = "rf",
    diagnostics = list(engine = "RF permutation importance",
                       cor_threshold = cor_threshold,
                       forced = intersect(force_include, env_vars))
  )
}

# Post-hoc guarantee that ecologically mandatory predictors survive selection.
# Applied by cast_select() after any selector (cpi/dml/rf) so the behaviour is
# identical across methods: forced predictors are unioned into `selected` and
# flagged in the score table. Names absent from the candidate pool are skipped.
.cast_apply_force_include <- function(selected, scores, env_vars,
                                      force_include, verbose = TRUE) {
  force_include <- intersect(unique(force_include), env_vars)
  if (!length(force_include)) return(list(selected = selected, scores = scores))
  newly <- setdiff(force_include, selected)
  selected <- unique(c(selected, force_include))
  if (!is.null(scores) && "variable" %in% names(scores)) {
    hit <- scores$variable %in% force_include
    if ("selected" %in% names(scores)) scores$selected[hit] <- TRUE
    if ("forced" %in% names(scores)) scores$forced[hit] <- TRUE
    if ("exclusion_reason" %in% names(scores)) {
      scores$exclusion_reason[hit] <- "force_included"
    }
  }
  if (verbose && length(newly)) {
    cli::cli_inform(
      "Force-included {length(newly)} predictor{?s} the selector dropped: {.val {newly}}."
    )
  }
  list(selected = selected, scores = scores)
}

#' Compare castSDM Conditional Screening Against Conventional Baselines
#'
#' Runs the predictor-selection strategies most researchers actually use --
#' a correlation filter (Dormann et al. 2013), stepwise VIF (Naimi & Araujo
#' 2016), a univariate marginal screen (Guisan & Zimmermann 2000), and
#' random-forest permutation importance (Cutler et al. 2007) -- alongside the
#' castSDM conditional predictive impact screen, on the *same* data. The
#' associational baselines cannot separate a genuine driver from a collinear
#' bystander, whereas CPI tests each predictor given all others; the returned
#' object makes that contrast explicit for the Results narrative.
#'
#' Because stepwise VIF and pairwise correlation scale poorly, the predictor
# pool is first reduced to the `pool` most marginally associated predictors
# (this mirrors practice -- nobody runs VIF on hundreds of raw layers) and
# every method, including CPI, is compared on that shared pool. The RF
# importance baseline keeps the same number of predictors as CPI so the
# comparison is budget-matched.
#'
#' @param data Data frame with response, coordinates, and predictors.
#' @param response Binary response column. Default `"presence"`.
#' @param cpi Optional precomputed `cast_select` object (`method = "cpi"`) to
#'   reuse instead of recomputing the conditional screen.
#' @param pool Integer. Size of the shared candidate pool. Default `40`.
#' @param cor_threshold Correlation-filter cutoff. Default `0.7`.
#' @param vif_threshold Stepwise-VIF cutoff. Default `10`.
#' @param alpha FDR level for the marginal/associational screens and CPI.
#'   Default `0.05`.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_screen_comparison` object.
#' @seealso [cast_select()], [plot.cast_screen_comparison()]
#' @export
cast_screen_comparison <- function(data, response = "presence", cpi = NULL,
                                   pool = 40L, cor_threshold = 0.7,
                                   vif_threshold = 10, alpha = 0.05,
                                   seed = NULL, verbose = TRUE) {
  env_vars <- get_env_vars(data, response)
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")

  y <- as.integer(data[[response]])
  xm <- .cast_numeric_matrix(data, env_vars)

  # castSDM conditional predictive impact (reuse if supplied).
  if (is.null(cpi)) {
    cpi <- cast_select(data, response = response, method = "cpi",
                       alpha = alpha, seed = seed, verbose = verbose)
  }

  # Shared candidate pool: the most marginally associated predictors, always
  # augmented with the CPI selection so conditionally-important but marginally
  # weak drivers stay visible in the comparison.
  assoc <- abs(suppressWarnings(stats::cor(xm, y, use = "pairwise.complete.obs")))
  assoc <- stats::setNames(as.numeric(assoc), env_vars)
  assoc[!is.finite(assoc)] <- 0
  pool <- min(as.integer(pool), length(env_vars))
  cand <- union(env_vars[order(assoc, decreasing = TRUE)][seq_len(pool)],
                intersect(env_vars, cpi$selected))

  df <- as.data.frame(xm[, cand, drop = FALSE], check.names = FALSE)
  df[[response]] <- y

  assoc_p <- function(vars) vapply(vars, function(v)
    suppressWarnings(stats::cor.test(df[[v]], df[[response]])$p.value), numeric(1))

  # 1. Correlation filter
  a <- abs(suppressWarnings(stats::cor(df[, cand, drop = FALSE], df[[response]])))
  a <- stats::setNames(as.numeric(a), cand); a[!is.finite(a)] <- 0
  kept <- character(0)
  for (v in cand[order(-a)]) {
    if (!length(kept) ||
        max(abs(stats::cor(df[[v]], df[, kept, drop = FALSE]))) < cor_threshold)
      kept <- c(kept, v)
  }
  cor_keep <- kept[stats::p.adjust(assoc_p(kept), "BH") < alpha]

  # 2. Stepwise VIF
  keep <- cand
  repeat {
    if (length(keep) < 2L) break
    vifs <- vapply(keep, function(v) {
      r2 <- summary(stats::lm(
        stats::reformulate(setdiff(keep, v), response = v), data = df))$r.squared
      if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
    }, numeric(1))
    if (max(vifs) <= vif_threshold) break
    keep <- setdiff(keep, names(which.max(vifs)))
  }
  vif_keep <- keep[stats::p.adjust(assoc_p(keep), "BH") < alpha]

  # 3. Univariate marginal screen
  uni_keep <- cand[stats::p.adjust(assoc_p(cand), "BH") < alpha]

  # CPI retentions restricted to the shared pool.
  cpi_keep <- intersect(cand, cpi$selected)

  # 4. RF permutation importance: budget-matched to the CPI retention count
  #    so the associational benchmark and the conditional screen compete on
  #    the same number of predictors.
  rf_keep <- character(0)
  if (requireNamespace("ranger", quietly = TRUE)) {
    rdf <- df[, c(response, cand)]
    rdf[[response]] <- factor(rdf[[response]])
    rf <- ranger::ranger(
      stats::reformulate(make.names(cand), response = response),
      data = stats::setNames(rdf, c(response, make.names(cand))),
      num.trees = 300L, importance = "permutation",
      probability = TRUE, seed = seed %||% 42L, num.threads = 1L
    )
    imp <- rf$variable.importance
    imp_orig <- stats::setNames(as.numeric(imp),
                                cand[match(names(imp), make.names(cand))])
    ranked_imp <- sort(imp_orig[imp_orig > 0], decreasing = TRUE)
    n_keep <- min(length(ranked_imp), max(1L, length(cpi_keep)))
    rf_keep <- names(ranked_imp)[seq_len(n_keep)]
  }

  membership <- data.frame(
    variable = cand,
    correlation = cand %in% cor_keep,
    vif = cand %in% vif_keep,
    univariate = cand %in% uni_keep,
    rf = cand %in% rf_keep,
    cpi = cand %in% cpi_keep,
    stringsAsFactors = FALSE
  )
  methods <- c("correlation", "vif", "univariate", "rf", "cpi")
  membership$n_methods <- rowSums(membership[methods])
  membership <- membership[membership$n_methods > 0, , drop = FALSE]

  if (verbose) {
    cli::cli_inform(c(
      "Screen comparison on {pool} shared candidates:",
      i = "correlation {sum(membership$correlation)} | vif {sum(membership$vif)} | univariate {sum(membership$univariate)} | rf {sum(membership$rf)} | CPI {sum(membership$cpi)}"
    ))
  }

  new_cast_screen_comparison(
    membership = membership, methods = methods, cpi_method = "cpi",
    diagnostics = list(pool = pool, cor_threshold = cor_threshold,
                       vif_threshold = vif_threshold, alpha = alpha,
                       n_predictors = length(env_vars))
  )
}
