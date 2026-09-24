#' Select Variables for Species Distribution Models
#'
#' Two-stage variable selection:
#' \enumerate{
#'   \item \strong{Stage 1 — collinearity thinning}: rank predictors by a
#'     univariate quadratic-logistic signal (the smaller p-value of the two
#'     `poly(x, 2)` terms, which sees U-shaped responses a linear correlation
#'     misses), then greedily keep predictors whose pairwise correlation with
#'     all kept predictors is \eqn{\le 0.7} (Dormann et al. 2013). Predictors
#'     whose GLM fails to fit fall back to the marginal-correlation rank.
#'   \item \strong{Stage 2 — spatial forward selection}: starting from the
#'     prespecified predictors (if any), repeatedly add the stage-1 survivor
#'     that most improves the inner cross-validated loss of a probability
#'     random forest, and stop when nothing improves it further. Folds are
#'     spatial when coordinates are available, random otherwise. There is no
#'     predictor-count cap and no fallback set: if no candidate improves on
#'     the intercept-only model, the selection is empty. Stopping is driven
#'     by predictive performance (Meyer et al. 2018, 2019), not by an
#'     arbitrary count.
#' }
#'
#' @section Why forward selection on spatial-CV loss:
#' Importance scores rank predictors by how much the fitted model leans on
#' them, which rewards collinear stand-ins for the correlation they borrow
#' from the true drivers (Hooker, Mentch & Zhou 2021). A forward search asks
#' the opposite question — does adding this predictor improve predictions on
#' held-out spatial folds? — so a redundant proxy adds nothing once its
#' parents are in the model, while a count cap is never needed: the search
#' stops by itself. Selection serves parsimony for interpretation and
#' projection; it does not identify a causal adjustment set.
#'
#' @param data Data frame with response and predictors (coordinates allowed;
#'   they are never selected, but `lon`/`lat` define the inner spatial folds
#'   when both are present and finite).
#' @param response Binary response column. Default `"presence"`.
#' @param method `"two_stage"` (default) or `"full"` (keep every predictor).
#' @param num_trees Trees per forest in the stage-2 forward search.
#'   Default 300.
#' @param max_rows Rows sampled (evenly, preserving order) before the
#'   stage-2 search. Bounds stage-2 cost on large data sets.
#'   Default 2000.
#' @param metric Inner-CV loss minimised by the forward search: `"brier"`
#'   (default, mean squared error of the predicted probability, sensitive to
#'   calibration) or `"auc"` (rank discrimination, insensitive to calibration
#'   shifts).
#' @param tolerance Non-negative number. A candidate is admitted only if its
#'   paired inner-CV loss improvement exceeds both this absolute floor and
#'   two standard errors of the fold differences (a 2-SE guard in the spirit
#'   of the `glmnet`/`rpart` one-standard-error rule, set a notch stricter
#'   because greedy search always takes the best of many candidates).
#'   Default `0`, i.e. the 2-SE guard alone stops the search; outer nested
#'   spatial CV in [cast_cv()] remains the honest performance estimate.
#' @param n_folds Number of inner folds for the stage-2 search. Default `3`.
#' @param seed Random seed.
#' @param verbose Print progress.
#' @param keep Character vector of predictors specified before screening, such
#'   as an exposure and a scientifically justified adjustment set. These bypass
#'   collinearity thinning and are always retained: expert knowledge outranks
#'   the data-driven search. Candidate covariates must still be scientifically
#'   admissible (not colliders or mediators of a total effect); retention does
#'   not verify a sufficient adjustment set.
#'
#' @return A `cast_select` object with `selected` (kept predictors, in
#'   admission order), `scores` (per-predictor marginal `assoc`, the
#'   `stage1_p` ranking signal with its `stage1_rank`, the stage-1
#'   `collinear_thinned` flag, the admission `step_added` (`0` for
#'   prespecified predictors, `NA` for never-admitted ones), the
#'   `loss_gain` (inner-CV loss improvement at admission, `NA` otherwise),
#'   the `selected_reason` (`"prespecified"`, `"forward"`, `"full"` or
#'   `"excluded"`), the `kept_by_design` indicator for `keep`, and the
#'   `selected` flag). An empty `selected` set is a valid answer: it means
#'   no candidate improved on the intercept-only model.
#'
#' @references
#' Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal
#' with it in ecological studies. \emph{Ecography} 36: 27-46.
#'
#' Meyer, H. et al. (2018). Improving performance of spatio-temporal machine
#' learning models using forward feature selection and target-oriented
#' validation. \emph{Environmental Modelling & Software} 101: 1-9.
#' \doi{10.1016/j.envsoft.2017.12.001}.
#'
#' Meyer, H. et al. (2019). Importance of spatial predictor variable selection
#' in machine learning applications — moving from data reproduction to spatial
#' prediction. \emph{Ecological Modelling} 411: 108815.
#' \doi{10.1016/j.ecolmodel.2019.108815}.
#' @seealso [cast_effect_table()], [cast_importance()]
#' @export
cast_select <- function(data, response = "presence",
                         method = c("two_stage", "full"),
                         num_trees = 300L, max_rows = 2000L,
                         metric = c("brier", "auc"),
                         tolerance = 0, n_folds = 3L,
                         seed = NULL, verbose = TRUE, keep = character(0)) {
  if (!missing(method) && (!is.character(method) || length(method) != 1L ||
                          is.na(method) || !method %in% c("two_stage", "full"))) {
    cli::cli_abort("{.arg method} must be one of 'two_stage' or 'full'.")
  }
  method <- match.arg(method)
  if (missing(metric)) metric <- "brier"
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
      !metric %in% c("brier", "auc")) {
    cli::cli_abort("{.arg metric} must be {.val brier} or {.val auc}.")
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
  if (is.na(num_trees) || num_trees < 1L) {
    cli::cli_abort("{.arg num_trees} must be at least 1.")
  }
  max_rows <- as.integer(max_rows)
  if (is.na(max_rows) || max_rows < 2L) {
    cli::cli_abort("{.arg max_rows} must be at least 2.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    cli::cli_abort("{.arg tolerance} must be one finite non-negative number.")
  }
  n_folds <- as.integer(n_folds)
  if (is.na(n_folds) || n_folds < 2L) {
    cli::cli_abort("{.arg n_folds} must be at least 2.")
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

  # ---- Stage 2: spatial forward selection on inner-CV loss -----------------
  step_added <- stats::setNames(rep(NA_integer_, length(env_vars)), env_vars)
  loss_gain <- stats::setNames(rep(NA_real_, length(env_vars)), env_vars)
  reason <- stats::setNames(rep("excluded", length(env_vars)), env_vars)
  path <- data.frame(step = integer(0), added = character(0),
                     loss = numeric(0), gain = numeric(0),
                     stringsAsFactors = FALSE)
  diagnostics <- list(stage1_kept = kept, stage1_metric = "poly2-glm-min-p",
                      num_trees = num_trees, max_rows = max_rows,
                      metric = metric, tolerance = tolerance,
                      n_folds = n_folds, keep = keep)
  if (!length(kept)) {
    cli::cli_warn("No varying predictor survived stage 1; returning an empty set.")
    selected <- kept
  } else {
    fwd <- .cast_forward_search(data, response, kept, keep, num_trees,
                                max_rows, metric, tolerance, n_folds,
                                seed, verbose)
    selected <- fwd$selected
    step_added[names(fwd$step_added)] <- fwd$step_added
    loss_gain[names(fwd$loss_gain)] <- fwd$loss_gain
    reason[selected] <- fwd$reason[selected]
    reason[keep] <- "prespecified"
    path <- fwd$path
    diagnostics <- c(diagnostics, fwd$diagnostics)
    diagnostics$path <- fwd$path
  }
  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    stage1_p = unname(stage1_p[env_vars]),
    stage1_rank = unname(stage1_rank[env_vars]),
    collinear_thinned = unname(thinned[env_vars]),
    step_added = unname(step_added[env_vars]),
    loss_gain = unname(loss_gain[env_vars]),
    kept_by_design = env_vars %in% keep,
    selected_reason = unname(reason[env_vars]),
    selected = env_vars %in% selected,
    stringsAsFactors = FALSE)
  scores$selected_reason[scores$selected_reason %in% c(NA, "NA")] <- "excluded"

  new_cast_select(selected = selected, scores = scores, method = "two_stage",
                  diagnostics = diagnostics)
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

#' Inner folds for the stage-2 forward search
#'
#' Spatial grid folds when finite `lon`/`lat` columns exist, otherwise plain
#' random folds. Never aborts on degenerate coordinates: [make_spatial_folds()]
#' already collapses those.
#' @keywords internal
#' @noRd
.cast_inner_folds <- function(data, n_folds) {
  if (all(c("lon", "lat") %in% names(data))) {
    lon <- suppressWarnings(as.numeric(data$lon))
    lat <- suppressWarnings(as.numeric(data$lat))
    if (all(is.finite(lon)) && all(is.finite(lat))) {
      return(list(folds = make_spatial_folds(lon, lat, k = n_folds,
                                             method = "grid", seed = NULL),
                  method = "spatial"))
    }
  }
  list(folds = sample(rep(seq_len(n_folds), length.out = nrow(data))),
       method = "random")
}

#' Inner-CV loss of one candidate predictor set
#'
#' Fits a probability forest per inner fold and returns the mean held-out
#' loss together with the per-fold losses (for paired comparison). Folds
#' without two response classes (or failed predictions) contribute `NA`; a
#' set with no evaluable fold scores `NA`.
#' @keywords internal
#' @noRd
.cast_inner_loss <- function(X, y_num, folds, vars, metric, num_trees,
                             seed_base, counter) {
  lv <- sort(unique(folds))
  fl <- vapply(lv, function(f) {
    tr <- which(folds != f)
    te <- which(folds == f)
    if (!length(tr) || !length(te) ||
        length(unique(y_num[tr])) < 2L) return(NA_real_)
    counter$seed <- counter$seed + 1L
    m <- tryCatch(
      ranger::ranger(x = X[tr, vars, drop = FALSE], y = factor(y_num[tr]),
                     probability = TRUE, num.trees = num_trees,
                     seed = seed_base + counter$seed, num.threads = 1L),
      error = function(e) NULL)
    if (is.null(m)) return(NA_real_)
    p <- tryCatch(
      suppressWarnings(stats::predict(m, data = X[te, vars, drop = FALSE])$predictions[, "1"]),
      error = function(e) NULL)
    if (is.null(p)) return(NA_real_)
    if (identical(metric, "brier")) {
      ok <- is.finite(p)
      if (!any(ok)) return(NA_real_)
      mean((p[ok] - y_num[te][ok])^2)
    } else {
      compute_auc(as.integer(y_num[te]), as.numeric(p))
    }
  }, numeric(1))
  list(loss = if (all(!is.finite(fl))) NA_real_ else mean(fl[is.finite(fl)]),
       folds = stats::setNames(fl, lv), fits = length(lv))
}

#' Intercept-only inner-CV loss (the honest null model)
#'
#' For `"brier"`, each fold predicts its own training prevalence; for `"auc"`,
#' a constant scores 0.5 by definition (reported only where the test fold
#' carries two classes, keeping the pairing honest).
#' @keywords internal
#' @noRd
.cast_null_loss <- function(y_num, folds, metric) {
  lv <- sort(unique(folds))
  if (identical(metric, "auc")) {
    fl <- vapply(lv, function(f) {
      if (length(unique(y_num[folds == f])) < 2L) return(NA_real_)
      0.5
    }, numeric(1))
    return(list(loss = 0.5, folds = stats::setNames(fl, lv)))
  }
  fl <- vapply(lv, function(f) {
    tr <- which(folds != f)
    te <- which(folds == f)
    if (!length(tr) || !length(te)) return(NA_real_)
    p0 <- mean(y_num[tr])
    if (!is.finite(p0)) return(NA_real_)
    mean((p0 - y_num[te])^2)
  }, numeric(1))
  list(loss = if (all(!is.finite(fl))) NA_real_ else mean(fl[is.finite(fl)]),
       folds = stats::setNames(fl, lv))
}

#' Paired admission test: 2-SE rule with an absolute floor
#'
#' Gains are signed so positive means better (`ref - cand` for Brier,
#' `cand - ref` for AUC), paired by inner fold. Admits iff the mean gain
#' over finite paired folds exceeds both `tolerance` and two standard errors
#' of the gains. The 2-SE bar (rather than 1-SE) guards against the
#' best-of-many luck intrinsic to greedy search over correlated folds, whose
#' overlap makes the naive SE optimistic. Returns the mean gain (or `NA`
#' when no paired fold exists).
#' @keywords internal
#' @noRd
.cast_admits <- function(ref_folds, cand_folds, tolerance, metric) {
  common <- intersect(names(ref_folds), names(cand_folds))
  g <- if (identical(metric, "brier")) {
    ref_folds[common] - cand_folds[common]
  } else {
    cand_folds[common] - ref_folds[common]
  }
  g <- unname(g[is.finite(g)])
  if (!length(g)) return(list(admit = FALSE, gain = NA_real_))
  se <- if (length(g) >= 2L) 2 * stats::sd(g) / sqrt(length(g)) else 0
  gain <- mean(g)
  list(admit = is.finite(gain) && gain > max(tolerance, se), gain = gain)
}

#' Spatial forward selection on inner-CV loss
#'
#' Starts from the prespecified set (or the best pair when nothing is
#' prespecified), then greedily admits the candidate with the largest loss
#' improvement while it exceeds `tolerance`. Returns an empty selection —
#' with a warning, not a fallback set — when nothing improves on the
#' intercept-only model.
#' @keywords internal
#' @noRd
.cast_forward_search <- function(data, response, kept, keep, num_trees,
                                 max_rows, metric, tolerance, n_folds,
                                 seed, verbose) {
  check_suggested("ranger", "for stage-2 forward selection")
  X <- data[, kept, drop = FALSE]
  for (col in names(X)) X[[col]] <- as.numeric(X[[col]])
  y_raw <- data[[response]]
  y_num <- suppressWarnings(as.numeric(y_raw))
  ok <- rowSums(!is.finite(as.matrix(X))) == 0L &
    is.finite(y_num) & y_num %in% c(0, 1)
  rows <- which(ok)
  if (!length(rows)) {
    cli::cli_abort(c(
      "Stage 2 needs complete predictor rows and two response classes.",
      "i" = "Check for missing predictor values or a single-class response."))
  }
  # Bound the search cost; even steps keep the row distribution.
  if (length(rows) > max_rows) {
    rows <- rows[unique(round(seq(1, length(rows), length.out = max_rows)))]
  }
  X <- X[rows, , drop = FALSE]
  y_num <- y_num[rows]
  geo <- data[rows, , drop = FALSE]
  if (length(unique(y_num)) < 2L) {
    cli::cli_abort(c(
      "Stage 2 needs complete predictor rows and two response classes.",
      "i" = "Check for missing predictor values or a single-class response."))
  }
  if (!is.null(seed)) set.seed(seed)
  seed_base <- seed %||% 1L
  inner <- .cast_inner_folds(geo, n_folds)
  folds <- inner$folds
  if (length(unique(folds)) < 2L) {
    cli::cli_abort(c(
      "Stage 2 needs at least two non-empty inner folds.",
      "i" = "Provide more rows or fewer {.arg n_folds}."))
  }
  counter <- new.env(parent = emptyenv())
  counter$seed <- 0L
  eval_loss <- function(vars) {
    .cast_inner_loss(X, y_num, folds, vars, metric, num_trees,
                     seed_base, counter)
  }
  null_res <- .cast_null_loss(y_num, folds, metric)
  n_fits <- 0L
  path <- data.frame(step = integer(0), added = character(0),
                     loss = numeric(0), gain = numeric(0),
                     stringsAsFactors = FALSE)
  step_added <- stats::setNames(rep(NA_integer_, length(kept)), kept)
  loss_gain <- stats::setNames(rep(NA_real_, length(kept)), kept)
  status <- "selected"

  current <- keep
  pool <- setdiff(kept, keep)
  cur <- if (length(current)) {
    r <- eval_loss(current)
    n_fits <- n_fits + r$fits
    r
  } else {
    list(loss = NA_real_, folds = stats::setNames(numeric(0), character(0)), fits = 0L)
  }
  if (length(keep)) {
    step_added[keep] <- 0L
    g0 <- .cast_admits(null_res$folds, cur$folds, tolerance, metric)
    loss_gain[keep] <- g0$gain
    path <- rbind(path, data.frame(step = 0L, added = paste(keep, collapse = "+"),
                                   loss = cur$loss,
                                   gain = g0$gain,
                                   stringsAsFactors = FALSE))
    if (verbose) cli::cli_inform("Stage 2: starting from {length(keep)} prespecified predictor{?s}.")
  } else if (length(pool) >= 2L) {
    # CAST-style start: the best pair must first beat the intercept model,
    # otherwise the honest answer is an empty set, not a fallback set.
    pairs <- utils::combn(pool, 2, simplify = FALSE)
    if (verbose) cli::cli_inform("Stage 2: testing {length(pairs)} candidate pairs in {n_folds} inner {inner$method} folds...")
    cand <- lapply(pairs, function(pr) {
      r <- eval_loss(pr)
      n_fits <<- n_fits + r$fits
      r
    })
    verdict <- lapply(cand, function(r) .cast_admits(null_res$folds, r$folds, tolerance, metric))
    gains <- vapply(verdict, function(v) ifelse(is.finite(v$gain), v$gain, -Inf), numeric(1))
    best <- which.max(gains)
    if (!length(best) || !isTRUE(verdict[[best]]$admit)) {
      cli::cli_warn(c(
        "No candidate pair improves on the intercept-only model; returning an empty set.",
        "i" = "This is a finding (no detectable signal), not a failure. Prespecify {.arg keep} to force predictors in."))
      return(list(selected = character(0), step_added = step_added,
                  loss_gain = loss_gain, reason = stats::setNames(character(0), character(0)),
                  path = path,
                  diagnostics = list(inner_method = inner$method, null_loss = null_res$loss,
                                     final_loss = NA_real_, n_fits = n_fits,
                                     status = "empty_selection")))
    }
    current <- pairs[[best]]
    cur <- cand[[best]]
    g <- verdict[[best]]$gain
    step_added[current] <- 1L
    loss_gain[current] <- g
    path <- rbind(path, data.frame(step = 1L, added = paste(current, collapse = "+"),
                                   loss = cur$loss, gain = g,
                                   stringsAsFactors = FALSE))
    pool <- setdiff(pool, current)
    if (verbose) cli::cli_inform("Stage 2: starting pair {.val {current}}.")
  } else {
    # A single candidate: admit it only if it beats the intercept model.
    r <- eval_loss(pool)
    n_fits <- n_fits + r$fits
    v <- .cast_admits(null_res$folds, r$folds, tolerance, metric)
    if (!isTRUE(v$admit)) {
      cli::cli_warn(c(
        "The single candidate does not improve on the intercept-only model; returning an empty set.",
        "i" = "This is a finding (no detectable signal), not a failure. Prespecify {.arg keep} to force predictors in."))
      return(list(selected = character(0), step_added = step_added,
                  loss_gain = loss_gain, reason = stats::setNames(character(0), character(0)),
                  path = path,
                  diagnostics = list(inner_method = inner$method, null_loss = null_res$loss,
                                     final_loss = NA_real_, n_fits = n_fits,
                                     status = "empty_selection")))
    }
    current <- pool
    cur <- r
    g <- v$gain
    step_added[current] <- 1L
    loss_gain[current] <- g
    path <- rbind(path, data.frame(step = 1L, added = current,
                                   loss = cur$loss, gain = g,
                                   stringsAsFactors = FALSE))
    pool <- character(0)
  }

  # Greedy additions while any candidate passes the paired 1-SE admission test.
  step <- max(path$step)
  repeat {
    if (!length(pool)) break
    cand <- lapply(pool, function(v) {
      r <- eval_loss(c(current, v))
      n_fits <<- n_fits + r$fits
      r
    })
    verdict <- lapply(cand, function(r) .cast_admits(cur$folds, r$folds, tolerance, metric))
    gains <- vapply(verdict, function(x) ifelse(is.finite(x$gain), x$gain, -Inf), numeric(1))
    best <- which.max(gains)
    if (!length(best) || !isTRUE(verdict[[best]]$admit)) break
    step <- step + 1L
    v <- pool[best]
    current <- c(current, v)
    cur <- cand[[best]]
    step_added[v] <- step
    loss_gain[v] <- verdict[[best]]$gain
    path <- rbind(path, data.frame(step = step, added = v,
                                   loss = cur$loss, gain = verdict[[best]]$gain,
                                   stringsAsFactors = FALSE))
    pool <- setdiff(pool, v)
    if (verbose) cli::cli_inform("Stage 2: step {step} admits {.val {v}}.")
  }
  if (verbose) cli::cli_inform("Stage 2: retaining {length(current)} predictor{?s} after {step} forward step{?s}.")

  reason <- stats::setNames(rep("forward", length(current)), current)
  list(selected = current, step_added = step_added, loss_gain = loss_gain,
       reason = reason, path = path,
       diagnostics = list(inner_method = inner$method, null_loss = null_res$loss,
                          final_loss = cur$loss, n_fits = n_fits,
                          status = "selected"))
}
