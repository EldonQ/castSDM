# Conditional Predictive Impact (CPI) variable selection --------------------
#
# For each candidate predictor we test its conditional independence from
# occurrence given every other predictor: a random forest is cross-fitted once,
# each predictor is replaced by a second-order Gaussian knockoff, and the
# increase in held-out log-loss measures its conditional contribution
# (Watson & Wright 2021). Predictors whose one-sided test survives
# Benjamini-Hochberg FDR control are kept.

#' @keywords internal
#' @noRd
.cast_quiet_mlr3 <- function() {
  if (requireNamespace("lgr", quietly = TRUE)) {
    suppressWarnings({
      try(lgr::get_logger("mlr3")$set_threshold("off"), silent = TRUE)
      try(lgr::get_logger("bbotk")$set_threshold("off"), silent = TRUE)
    })
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.cast_prescreen_importance <- function(x_std, y, keep, num_trees, seed,
                                       num_threads = 1L) {
  rf_dat <- as.data.frame(x_std, check.names = FALSE)
  rf_dat$.response <- factor(y)
  rf <- ranger::ranger(
    .response ~ ., data = rf_dat, probability = TRUE,
    num.trees = as.integer(num_trees), importance = "permutation",
    seed = seed %||% 42L, num.threads = as.integer(num_threads),
    verbose = FALSE
  )
  imp <- rf$variable.importance
  imp[!is.finite(imp)] <- 0
  names(sort(imp, decreasing = TRUE))[seq_len(min(keep, length(imp)))]
}

#' Response-Stratified Fold Assignment
#'
#' Assigns each class of the binary response evenly across folds so every
#' fold keeps both classes (rare-species guard for the CPI cross-fit).
#'
#' @param y Integer 0/1 response.
#' @param k Number of folds.
#' @param seed Random seed.
#' @return Integer fold ids in `1:k`.
#' @keywords internal
#' @noRd
.cast_stratified_folds <- function(y, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fold_id <- integer(length(y))
  for (cl in sort(unique(y))) {
    idx <- which(y == cl)
    fold_id[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  as.integer(fold_id)
}

# Assign observations to a fine spatial grid of (up to) `n_blocks` cells.
# Each grid cell is one inference block: per-observation log-loss differences
# are averaged within a block and the block means feed the one-sided t-test,
# so spatially autocorrelated observations inside a block contribute one
# effective unit instead of being treated as independent (cluster-robust
# variance on spatial blocks). Bins are quantile-based (balanced counts) but
# produce contiguous lon x lat rectangles, exactly like the grid-block CV.
.cast_spatial_blocks <- function(lon, lat, n_blocks, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(lon)
  if (n < 2L) return(rep(1L, n))
  side <- max(2L, ceiling(sqrt(max(2L, as.integer(n_blocks)))))
  xb <- .cast_bin(lon, side)
  yb <- .cast_bin(lat, side)
  as.integer(factor(interaction(xb, yb, drop = TRUE)))
}


#' In-Package CPI Core Estimator
#'
#' Follows the \pkg{cpi} (0.1.5) reference flow: one cross-fitted learner
#' (`resample(..., store_models = TRUE)`), second-order Gaussian knockoffs from
#' [knockoff::create.second_order()], and the per-observation log-loss
#' difference between full and knockoff-reduced predictions (with the reference
#' `1e-15` probability clamp). Two departures from the reference:
#'
#' - predictors are first **Gaussian-quantile-transformed** (rank -> `qnorm`,
#'   a monotone map) so `create.second_order`'s multivariate-Gaussian
#'   assumption holds for skewed bioclimatic predictors; a tree learner is
#'   invariant to this transform, so the CPI is unchanged, and
#' - a **single** knockoff draw is used (the reference uses one draw; repeated
#'   draws only vary the knockoff Monte-Carlo and their median p-value is not a
#'   valid pooled p-value).
#'
#' Folds are response-stratified (rare-species safety). `inference = "block"`
#' (default) averages the per-observation differences to one mean per spatial
#' grid block and runs the one-sided t-test on those block means
#' (`df = n_blocks - 1`) -- a cluster-robust variance estimate on spatial
#' blocks. `inference = "fold"` aggregates to one mean per cross-fitting fold
#' (`df = folds - 1`); `inference = "observation"` (the \pkg{cpi} default)
#' treats n observations as i.i.d. and is anti-conservative under spatial
#' autocorrelation; both are kept for comparability.
#'
#' @param x Standardized numeric predictor matrix (column names used).
#' @param y Integer 0/1 response.
#' @param n_folds Cross-fitting folds.
#' @param num_trees Trees for the random-forest learner.
#' @param num_threads Threads for the ranger learner.
#' @param inference `"block"` (default), `"fold"`, or `"observation"`.
#' @param lon,lat Coordinate vectors for `inference = "block"`.
#' @param n_blocks Target number of spatial blocks (auto if `NULL`).
#' @param seed Random seed.
#'
#' @return A `data.frame` with `Variable`, `CPI`, `SE`, `statistic`,
#'   `p.value` (one-sided).
#' @keywords internal
#' @noRd
.cast_cpi_core <- function(x, y, n_folds, num_trees, num_threads = 1L,
                           inference = c("block", "fold", "observation"),
                           seed = NULL, lon = NULL, lat = NULL,
                           n_blocks = NULL) {
  inference <- match.arg(inference)
  check_suggested("mlr3", "for the CPI learner backend")
  check_suggested("mlr3learners", "for the ranger learner")
  check_suggested("ranger", "for the ranger learner")
  check_suggested("knockoff", "for second-order Gaussian knockoffs")
  .cast_quiet_mlr3()

  vars <- colnames(x)
  n <- nrow(x)
  k <- as.integer(n_folds)
  fold_id <- .cast_stratified_folds(y, k, seed)
  if (identical(inference, "block")) {
    if (is.null(lon) || is.null(lat) || length(lon) != n || length(lat) != n) {
      cli::cli_abort(c(
        "Block inference needs {.arg lon}/{.arg lat} for every row.",
        i = "Pass coordinate columns to {.fn cast_select}, or use {.code inference = \"fold\"}."
      ))
    }
    if (is.null(n_blocks)) n_blocks <- min(50L, max(20L, floor(n / 20L)))
    blk <- .cast_spatial_blocks(lon, lat, n_blocks, seed)
  }
  test_sets <- lapply(seq_len(k), function(f) which(fold_id == f))
  train_sets <- lapply(seq_len(k), function(f) which(fold_id != f))

  # Gaussian quantile transform (monotone per column): create.second_order
  # assumes ~multivariate Gaussian, which skewed environmental predictors are
  # not; the rank->qnorm map restores the assumption without changing a tree
  # learner's fit, so the knockoff exchangeability (and hence FDR) holds.
  xg <- apply(x, 2L, function(z) {
    stats::qnorm(rank(z, ties.method = "average") / (length(z) + 1))
  })
  colnames(xg) <- vars

  task_df <- data.frame(.presence = factor(y), xg, check.names = FALSE)
  task <- mlr3::TaskClassif$new("cast_cpi", task_df, target = ".presence")
  learner <- mlr3::lrn("classif.ranger", predict_type = "prob",
                       num.trees = as.integer(num_trees),
                       seed = seed %||% 42L,
                       num.threads = as.integer(num_threads), verbose = FALSE)
  # Stratified custom resampling: same cross-fit object the cpi reference
  # builds from rsmp("cv"), but with both classes guaranteed per fold.
  resampling <- mlr3::rsmp("custom")
  resampling$instantiate(task, train_sets = train_sets, test_sets = test_sets)
  fit <- mlr3::resample(task, learner, resampling, store_models = TRUE)

  # Per-observation log-loss, identical to cpi:::compute_loss (eps clamp).
  logloss <- function(truth, prob) {
    eps <- 1e-15
    ii <- match(as.character(truth), colnames(prob))
    p <- prob[cbind(seq_len(nrow(prob)), ii)]
    -log(pmax(eps, pmin(1 - eps, p)))
  }

  err_full <- numeric(n)
  for (f in seq_len(k)) {
    pred <- fit$learners[[f]]$predict(task, row_ids = test_sets[[f]])
    err_full[test_sets[[f]]] <- logloss(pred$truth, pred$prob)
  }

  # Single knockoff draw (the reference flow).
  if (!is.null(seed)) set.seed(seed + 1000L)
  x_tilde <- knockoff::create.second_order(xg)

  stat <- lapply(seq_along(vars), function(j) {
    reduced_df <- task_df
    reduced_df[[vars[j]]] <- as.numeric(x_tilde[, j])
    reduced_task <- mlr3::TaskClassif$new(
      sprintf("cast_cpi_reduced_%d", j), reduced_df, target = ".presence"
    )
    err_reduced <- numeric(n)
    for (f in seq_len(k)) {
      pred <- fit$learners[[f]]$predict(reduced_task,
                                        row_ids = test_sets[[f]])
      err_reduced[test_sets[[f]]] <- logloss(pred$truth, pred$prob)
    }
    dif <- err_reduced - err_full
    dif_test <- if (identical(inference, "fold")) {
      vapply(test_sets, function(idx) mean(dif[idx]), numeric(1))
    } else if (identical(inference, "block")) {
      as.numeric(tapply(dif, blk, mean))
    } else {
      dif
    }
    cpi_j <- mean(dif)
    se_j <- stats::sd(dif_test) / sqrt(length(dif_test))
    stat_j <- NA_real_
    p_j <- NA_real_
    if (length(unique(dif_test)) < 2L) {
      # Degenerate: the aggregated differences carry no variance.
      stat_j <- 0
      p_j <- 1
    } else {
      tt <- tryCatch(
        stats::t.test(dif_test, alternative = "greater"),
        error = function(e) NULL
      )
      if (!is.null(tt)) {
        stat_j <- unname(tt$statistic)
        p_j <- tt$p.value
      }
    }
    c(CPI = cpi_j, SE = se_j, statistic = stat_j, p.value = p_j)
  })
  res <- do.call(rbind, stat)
  data.frame(
    Variable = vars,
    CPI = res[, "CPI"],
    SE = res[, "SE"],
    statistic = res[, "statistic"],
    p.value = res[, "p.value"],
    stringsAsFactors = FALSE
  )
}


#' Conditional Predictive Impact (CPI) Variable Selection
#'
#' Tests each predictor's conditional independence from occurrence given every
#' other predictor using the conditional predictive impact (Watson & Wright
#' 2021). A random-forest learner is cross-fitted once; each predictor is then
#' replaced by a Gaussian knockoff and the increase in log-loss measures its
#' conditional contribution. Because the learner is a random forest, symmetric
#' unimodal (bell-shaped) niche responses are detected - unlike a scalar
#' partially linear coefficient. Predictors whose one-sided CPI test survives
#' Benjamini-Hochberg FDR control are retained.
#'
#' The estimator is an in-package reimplementation of the \pkg{cpi} (0.1.5)
#' reference flow (cross-fit learner, `create.second_order` knockoffs,
#' per-observation log-loss differences) with response-stratified folds and a
#' spatial-inference layer: with the default `inference = "block"` each spatial
#' grid block contributes one mean log-loss difference and the one-sided t-test
#' runs on the block means (cluster-robust, `df = n_blocks - 1`).
#' `inference = "fold"` aggregates to one mean per cross-fitting fold
#' (`df = folds - 1`); `inference = "observation"` reproduces the reference
#' per-observation t-test; that test assumes i.i.d. observations and is
#' **anti-conservative on spatially autocorrelated data** (effective sample
#' size << n), so p-values and the FDR screen are then over-optimistic. Even
#' block-level inference only mitigates the test layer - the cross-fit folds
#' remain random, so train/test are spatially adjacent within the cross-fit.
#'
#' @param data,env_vars,response Prepared inputs from [cast_select()].
#' @param alpha FDR level for the adjusted p-values. Default `0.05`.
#' @param max_candidates Optional candidate ceiling. `NULL` (default) tests
#'   every non-constant predictor conditionally; a positive integer triggers a
#'   single random-forest importance pre-screen to this many candidates (a
#'   marginal approximation used only for computational feasibility).
#' @param n_folds Cross-fitting folds. Default `10`; fold-level inference
#'   needs enough folds for its t-test (`df = folds - 1`).
#' @param num_trees Trees for the random-forest learner. Default `100`.
#' @param min_vars Fallback floor if nothing passes FDR. Default `0` (allow an
#'   empty selection); set a positive integer only to enforce an ecological
#'   prior, in which case the topped-up predictors are flagged `fallback`.
#' @param inference `"block"` (default), `"fold"`, or `"observation"`; see
#'   [.cast_cpi_core()].
#' @param lon,lat Coordinate vectors forwarded to the block-inference layer.
#' @param n_blocks Target number of spatial blocks (auto if `NULL`).
#' @param seed,verbose Standard control arguments.
#'
#' @return A list with `selected`, `scores`, and `diagnostics`.
#' @keywords internal
#' @noRd
.cast_select_cpi <- function(data, env_vars, response, alpha = 0.05,
                             max_candidates = NULL, n_folds = 10L,
                             num_trees = 100L, min_vars = 0L,
                             inference = c("block", "fold", "observation"),
                             seed = NULL, verbose = TRUE,
                             num_threads = 1L,
                             lon = NULL, lat = NULL, n_blocks = NULL) {
  check_suggested("mlr3", "for the CPI learner backend")
  check_suggested("mlr3learners", "for the ranger learner")
  check_suggested("ranger", "for the ranger learner")
  check_suggested("knockoff", "for second-order Gaussian knockoffs")
  .cast_quiet_mlr3()
  inference <- match.arg(inference)

  y <- as.integer(data[[response]])
  x <- .cast_numeric_matrix(data, env_vars)
  keep_var <- apply(x, 2L, stats::sd) > 1e-10
  active <- env_vars[keep_var]
  if (length(active) < 2L) {
    cli::cli_abort("CPI selection needs at least two non-constant predictors.")
  }
  x <- x[, keep_var, drop = FALSE]
  x_std <- scale(x)
  x_std[!is.finite(x_std)] <- 0
  safe_names <- make.names(active, unique = TRUE)
  colnames(x_std) <- safe_names

  candidates <- active
  prescreen <- "none"
  if (!is.null(max_candidates) && length(active) > max_candidates) {
    top_safe <- .cast_prescreen_importance(x_std, y, max_candidates,
                                           num_trees, seed, num_threads)
    candidates <- active[match(top_safe, safe_names)]
    prescreen <- "RF importance (top candidates)"
    if (verbose) {
      cli::cli_inform(
        "CPI: pre-screened {length(active)} -> {length(candidates)} candidates by RF importance."
      )
    }
  }
  cand_safe <- safe_names[match(candidates, active)]

  # Rare-species guard: stratified folds can only keep both classes per fold
  # when the rarer class has at least n_folds members.
  min_class <- min(table(y))
  if (min_class < n_folds) {
    cli::cli_abort(c(
      "CPI needs at least {n_folds} presence{?s} and background{?s} each for {n_folds} folds; the rarer class has {min_class}.",
      i = "Reduce {.arg n_folds} (e.g. {max(2L, min_class)}) or use {.code method = \"rf\"} for rare species."
    ))
  }
  if (identical(inference, "block")) {
    if (is.null(lon) || is.null(lat) || length(lon) != nrow(x_std)) {
      cli::cli_abort(c(
        "Block inference needs {.arg lon}/{.arg lat} for every row.",
        i = "Pass coordinate columns to {.fn cast_select}, or use {.code inference = \"fold\"}."
      ))
    }
    if (is.null(n_blocks)) n_blocks <- min(50L, max(20L, floor(nrow(x_std) / 20L)))
    blk <- .cast_spatial_blocks(lon, lat, n_blocks, seed)
  }
  res <- tryCatch(
    .cast_cpi_core(
      x_std[, cand_safe, drop = FALSE], y,
      n_folds = n_folds, num_trees = num_trees, num_threads = num_threads,
      inference = inference, seed = seed,
      lon = lon, lat = lat, n_blocks = n_blocks
    ),
    error = function(e) {
      cli::cli_abort(c(
        "CPI selection failed: {conditionMessage(e)}",
        i = "On rare or highly clustered species try fewer folds ({.arg n_folds}, currently {n_folds}) or {.code method = \"rf\"}."
      ), call = NULL)
    }
  )

  # CPI rows are already in candidate order (safe names used internally).
  orig <- candidates
  cpi_val <- res$CPI
  se_val  <- res$SE
  stat    <- res$statistic
  p_raw   <- res$p.value

  p_adj <- rep(NA_real_, length(p_raw))
  ok <- is.finite(p_raw)
  if (any(ok)) p_adj[ok] <- stats::p.adjust(p_raw[ok], method = "BH")

  selected <- orig[is.finite(p_adj) & p_adj < alpha & is.finite(cpi_val) &
                     cpi_val > 0]
  fdr_passed <- selected
  if (length(selected) < min_vars) {
    # Union semantics (matching the RF benchmark): keep every FDR passer and
    # top up with the strongest remaining candidates to reach the floor.
    # Only candidates with a finite, positive CPI may enter the fallback --
    # a candidate whose CPI estimation failed (NA) must never be force-included.
    finite_cpi <- is.finite(cpi_val) & cpi_val > 0
    ord <- order(cpi_val[finite_cpi], decreasing = TRUE)
    selected <- unique(c(selected, orig[finite_cpi][ord]))
    selected <- selected[seq_len(min(min_vars, length(selected)))]
    selected <- selected[!is.na(selected)]
  }
  fallback_added <- setdiff(selected, fdr_passed)

  scores <- data.frame(
    variable = orig,
    cpi = cpi_val,
    std_error = se_val,
    statistic = stat,
    abs_statistic = abs(stat),
    p_value = p_raw,
    p_adjusted = p_adj,
    selected = orig %in% selected,
    fallback = orig %in% fallback_added,
    forced = FALSE,
    stringsAsFactors = FALSE
  )
  not_tested <- setdiff(env_vars, orig)
  if (length(not_tested)) {
    scores <- rbind(scores, data.frame(
      variable = not_tested, cpi = NA_real_, std_error = NA_real_,
      statistic = NA_real_, abs_statistic = NA_real_, p_value = NA_real_,
      p_adjusted = NA_real_, selected = FALSE, fallback = FALSE,
      forced = FALSE, stringsAsFactors = FALSE
    ))
  }
  scores <- scores[order(!scores$selected, -scores$cpi), , drop = FALSE]
  rownames(scores) <- NULL

  diagnostics <- list(
    engine = "Conditional Predictive Impact (random-forest, log-loss knockoff)",
    alpha = alpha,
    fdr_method = "BH",
    test = if (identical(inference, "block")) {
      "one-sided t on spatial block means (cluster-robust)"
    } else if (identical(inference, "fold")) {
      "one-sided t on fold means (df = folds - 1)"
    } else {
      "one-sided t on per-observation differences (assumes i.i.d.)"
    },
    inference = inference,
    n_blocks = if (identical(inference, "block")) {
      length(unique(blk))
    } else NULL,
    n_folds = as.integer(n_folds),
    n_candidates = length(candidates),
    n_predictors = length(env_vars),
    prescreen = prescreen
  )
  if (verbose) {
    cli::cli_inform(
      "CPI selection: {length(selected)}/{length(candidates)} predictors significant (FDR < {alpha})."
    )
  }
  list(selected = selected, scores = scores, diagnostics = diagnostics)
}

