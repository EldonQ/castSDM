# Double Machine Learning (DML) variable selection ---------------------------
#
# For each candidate predictor we estimate its Neyman-orthogonal partially
# linear effect on occurrence while flexibly controlling for every other
# predictor (Chernozhukov et al. 2018; DoubleML, Bach et al. 2022). Predictors
# whose effect is significant after Benjamini-Hochberg FDR control are kept.
# The only user choice is the FDR level; the cross-fitting folds and nuisance
# learner are method defaults rather than hand-tuned ecological thresholds.

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

#' @keywords internal
#' @noRd
.cast_dml_effect <- function(df, treat, controls, n_folds, num_trees, seed,
                             num_threads = 1L) {
  dml_data <- DoubleML::DoubleMLData$new(
    df, y_col = ".presence", d_cols = treat, x_cols = controls
  )
  ml_l <- mlr3::lrn("regr.ranger", num.trees = as.integer(num_trees),
                    num.threads = as.integer(num_threads), verbose = FALSE)
  ml_m <- mlr3::lrn("regr.ranger", num.trees = as.integer(num_trees),
                    num.threads = as.integer(num_threads), verbose = FALSE)
  obj <- DoubleML::DoubleMLPLR$new(
    dml_data, ml_l = ml_l, ml_m = ml_m, n_folds = as.integer(n_folds)
  )
  if (!is.null(seed)) set.seed(seed)
  obj$fit()
  c(estimate = unname(obj$coef[treat]),
    std_error = unname(obj$se[treat]),
    statistic = unname(obj$t_stat[treat]),
    p_value = unname(obj$pval[treat]))
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
#' @param data,env_vars,response Prepared inputs from [cast_select()].
#' @param alpha FDR level for the adjusted p-values. Default `0.05`.
#' @param max_candidates When more predictors are supplied, a single
#'   random-forest importance pass pre-screens this many for feasibility; the
#'   CPI conditioning set is the pre-screened set.
#' @param n_folds Cross-fitting folds. Default `5`.
#' @param num_trees Trees for the random-forest learner. Default `100`.
#' @param min_vars Fallback floor if nothing passes FDR. Default `1`.
#' @param seed,verbose Standard control arguments.
#'
#' @return A list with `selected`, `scores`, and `diagnostics`.
#' @keywords internal
#' @noRd
.cast_select_cpi <- function(data, env_vars, response, alpha = 0.05,
                             max_candidates = 30L, n_folds = 5L,
                             num_trees = 100L, min_vars = 1L,
                             seed = NULL, verbose = TRUE,
                             num_threads = 1L) {
  check_suggested("cpi", "for CPI variable selection")
  check_suggested("mlr3", "for the CPI learner backend")
  check_suggested("mlr3learners", "for the ranger learner")
  check_suggested("ranger", "for the ranger learner")
  .cast_quiet_mlr3()

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
  if (length(active) > max_candidates) {
    top_safe <- .cast_prescreen_importance(x_std, y, max_candidates,
                                           num_trees, seed, num_threads)
    candidates <- active[match(top_safe, safe_names)]
    if (verbose) {
      cli::cli_inform(
        "CPI: pre-screened {length(active)} -> {length(candidates)} candidates by RF importance."
      )
    }
  }
  cand_safe <- safe_names[match(candidates, active)]

  task_df <- data.frame(.presence = factor(y), x_std[, cand_safe, drop = FALSE],
                        check.names = FALSE)
  task <- mlr3::TaskClassif$new("cast_cpi", task_df, target = ".presence")
  learner <- mlr3::lrn("classif.ranger", predict_type = "prob",
                       num.trees = as.integer(num_trees),
                       num.threads = as.integer(num_threads), verbose = FALSE)
  resampling <- mlr3::rsmp("cv", folds = as.integer(n_folds))

  if (!is.null(seed)) set.seed(seed)
  res <- cpi::cpi(
    task = task, learner = learner, resampling = resampling,
    measure = "classif.logloss", test = "t", verbose = FALSE
  )

  # Map CPI rows (safe names) back to original predictor names.
  orig <- candidates[match(res$Variable, cand_safe)]
  cpi_val <- res$CPI
  se_val  <- res$SE
  stat    <- res$statistic
  p_raw   <- res$p.value

  p_adj <- rep(NA_real_, length(p_raw))
  ok <- is.finite(p_raw)
  if (any(ok)) p_adj[ok] <- stats::p.adjust(p_raw[ok], method = "BH")

  selected <- orig[is.finite(p_adj) & p_adj < alpha & is.finite(cpi_val) &
                     cpi_val > 0]
  if (length(selected) < min_vars) {
    ord <- order(cpi_val, decreasing = TRUE)
    selected <- orig[ord][seq_len(min(min_vars, length(orig)))]
    selected <- selected[!is.na(selected)]
  }

  scores <- data.frame(
    variable = orig,
    cpi = cpi_val,
    std_error = se_val,
    statistic = stat,
    abs_statistic = abs(stat),
    p_value = p_raw,
    p_adjusted = p_adj,
    selected = orig %in% selected,
    stringsAsFactors = FALSE
  )
  not_tested <- setdiff(env_vars, orig)
  if (length(not_tested)) {
    scores <- rbind(scores, data.frame(
      variable = not_tested, cpi = NA_real_, std_error = NA_real_,
      statistic = NA_real_, abs_statistic = NA_real_, p_value = NA_real_,
      p_adjusted = NA_real_, selected = FALSE, stringsAsFactors = FALSE
    ))
  }
  scores <- scores[order(!scores$selected, -scores$cpi), , drop = FALSE]
  rownames(scores) <- NULL

  diagnostics <- list(
    engine = "Conditional Predictive Impact (random-forest, log-loss knockoff)",
    alpha = alpha,
    fdr_method = "BH",
    test = "one-sided t",
    n_folds = as.integer(n_folds),
    n_candidates = length(candidates),
    n_predictors = length(env_vars)
  )
  if (verbose) {
    cli::cli_inform(
      "CPI selection: {length(selected)}/{length(candidates)} predictors significant (FDR < {alpha})."
    )
  }
  list(selected = selected, scores = scores, diagnostics = diagnostics)
}

#' Double Machine Learning Variable Selection
#'
#' Estimates each predictor's orthogonalized partially linear effect on the
#' binary response while controlling for all other predictors, then retains the
#' predictors whose effect survives Benjamini-Hochberg FDR control.
#'
#' @param data,env_vars,response Prepared inputs from [cast_select()].
#' @param alpha FDR level for the adjusted p-values. Default `0.05`.
#' @param max_candidates Predictors tested with DML. When more predictors are
#'   supplied, a single random-forest importance pass selects this many for
#'   feasibility; every selected predictor is still adjusted for the full set.
#' @param n_folds Cross-fitting folds. Default `5`.
#' @param num_trees Trees per random-forest nuisance learner. Default `100`.
#' @param min_vars Fallback floor if nothing passes FDR. Default `1`.
#' @param seed,verbose Standard control arguments.
#'
#' @return A list with `selected`, `scores`, and `diagnostics`.
#' @keywords internal
#' @noRd
.cast_select_dml <- function(data, env_vars, response, alpha = 0.05,
                             max_candidates = 30L, n_folds = 5L,
                             num_trees = 100L, min_vars = 1L,
                             seed = NULL, verbose = TRUE,
                             num_threads = 1L) {
  check_suggested("DoubleML", "for DML variable selection")
  check_suggested("mlr3", "for DML nuisance learners")
  check_suggested("mlr3learners", "for the ranger nuisance learner")
  check_suggested("ranger", "for the ranger nuisance learner")
  .cast_quiet_mlr3()

  y <- as.integer(data[[response]])
  x <- .cast_numeric_matrix(data, env_vars)
  keep_var <- apply(x, 2L, stats::sd) > 1e-10
  active <- env_vars[keep_var]
  if (length(active) < 2L) {
    cli::cli_abort("DML selection needs at least two non-constant predictors.")
  }
  x <- x[, keep_var, drop = FALSE]
  x_std <- scale(x)
  x_std[!is.finite(x_std)] <- 0
  safe_names <- make.names(active, unique = TRUE)
  colnames(x_std) <- safe_names

  candidates <- active
  if (length(active) > max_candidates) {
    top_safe <- .cast_prescreen_importance(x_std, y, max_candidates,
                                           num_trees, seed, num_threads)
    candidates <- active[match(top_safe, safe_names)]
    if (verbose) {
      cli::cli_inform(
        "DML: pre-screened {length(active)} -> {length(candidates)} candidates by RF importance."
      )
    }
  }

  df <- data.frame(.presence = y, x_std, check.names = FALSE)
  cand_safe <- safe_names[match(candidates, active)]

  stats_list <- lapply(seq_along(candidates), function(i) {
    treat <- cand_safe[i]
    controls <- setdiff(safe_names, treat)
    out <- tryCatch(
      .cast_dml_effect(df, treat, controls, n_folds, num_trees,
                       if (is.null(seed)) NULL else seed + i, num_threads),
      error = function(e) {
        if (verbose) cli::cli_warn("DML failed for {candidates[i]}: {e$message}")
        c(estimate = NA_real_, std_error = NA_real_,
          statistic = NA_real_, p_value = NA_real_)
      }
    )
    out
  })
  est <- do.call(rbind, stats_list)

  p_raw <- est[, "p_value"]
  p_adj <- rep(NA_real_, length(p_raw))
  ok <- is.finite(p_raw)
  if (any(ok)) p_adj[ok] <- stats::p.adjust(p_raw[ok], method = "BH")

  selected <- candidates[is.finite(p_adj) & p_adj < alpha]
  if (length(selected) < min_vars) {
    ord <- order(abs(est[, "statistic"]), decreasing = TRUE)
    selected <- candidates[ord][seq_len(min(min_vars, length(candidates)))]
    selected <- selected[!is.na(selected)]
  }

  scores <- data.frame(
    variable = candidates,
    estimate = est[, "estimate"],
    std_error = est[, "std_error"],
    statistic = est[, "statistic"],
    abs_statistic = abs(est[, "statistic"]),
    p_value = p_raw,
    p_adjusted = p_adj,
    selected = candidates %in% selected,
    stringsAsFactors = FALSE
  )
  not_tested <- setdiff(env_vars, candidates)
  if (length(not_tested)) {
    scores <- rbind(scores, data.frame(
      variable = not_tested, estimate = NA_real_, std_error = NA_real_,
      statistic = NA_real_, abs_statistic = NA_real_, p_value = NA_real_,
      p_adjusted = NA_real_, selected = FALSE, stringsAsFactors = FALSE
    ))
  }
  scores <- scores[order(!scores$selected, -scores$abs_statistic), , drop = FALSE]
  rownames(scores) <- NULL

  diagnostics <- list(
    engine = "DoubleML PLR (random-forest nuisance)",
    alpha = alpha,
    fdr_method = "BH",
    n_folds = as.integer(n_folds),
    n_candidates = length(candidates),
    n_predictors = length(env_vars)
  )
  if (verbose) {
    cli::cli_inform(
      "DML selection: {length(selected)}/{length(candidates)} predictors significant (FDR < {alpha})."
    )
  }
  list(selected = selected, scores = scores, diagnostics = diagnostics)
}
