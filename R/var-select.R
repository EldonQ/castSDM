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
#' @param method `"two_stage"` (default), `"full"` (keep every predictor) or
#'   `"covsel"` (the two-step embedded selection of Adde et al. 2023: a
#'   hierarchical collinearity filter followed by elastic-net / GAM / guided-RF
#'   embedding with a `ceiling(log2(n_presences)) - 1` count rule).
#' @param corcut Correlation threshold for the `"covsel"` filter levels.
#'   Default `0.7` (the community standard; Dormann et al. 2013).
#' @param var_id Character vector aligned with the predictors giving the
#'   variable-level identity (e.g. all years and focal scales of one covariate
#'   share an id). At the variable level the filter keeps exactly one
#'   covariate per id. `NULL` skips that level.
#' @param category Character vector aligned with the predictors giving the
#'   category-level identity (e.g. "climate", "terrain"). `NULL` skips it.
#' @param algorithms Embedding learners for `"covsel"`: any of `"glm"`
#'   (elastic net), `"gam"` (null-space penalisation), `"rf"` (guided
#'   regularised random forest). Default all three.
#' @param ncov Target number of covariates for `"covsel"`. Default `NULL` =
#'   `ceiling(log2(n_presences)) - 1` (the Adde et al. rule), capped at
#'   `maxncov`.
#' @param maxncov Maximum number of covariates for `"covsel"`. Default `12`.
#' @param weights Optional numeric row weights (e.g. prevalence targets) used
#'   by the `"covsel"` filter and embedding steps.
#' @param ensure Optional named list of logical flag vectors aligned with the
#'   predictors, declaring coverage requirements the data-driven search cannot
#'   guarantee. For each flag (e.g. `list(climate = c(TRUE, FALSE, ...))`), if
#'   the final set contains no flagged covariate, the strongest Step-A-ranked
#'   flagged survivor is admitted (reason `"ensured"`). This is the generic
#'   mechanism behind "a scenario-projection model must contain at least one
#'   time-varying driver": the flag comes from covariate provenance (e.g. the
#'   pool manifest's climate declaration), never a hardcoded variable name, so
#'   the rule is identical for every species while the admitted covariate is
#'   whatever the data ranks strongest. Coverage is a purpose-adequacy
#'   requirement, not an ecological claim; importance reporting stays honest
#'   about how much the ensured covariate contributes.
#' @param env_vars Optional explicit predictor set (a contract, as in
#'   [cast_prepare()]). When supplied, selection runs on
#'   `env_vars` intersected with what the data actually supports (a layer can
#'   be constant inside a CV fold's training sample even when it varies in the
#'   full data); `var_id`, `category` and `keep` are re-indexed to that
#'   intersection by name. `NULL` keeps auto-detection.
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
#'   two standard errors of the fold differences, with at least two finite
#'   paired folds. This heuristic does not control a false-positive rate.
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
#'   `selected` flag). An empty `selected` set is a valid answer: no
#'   candidate passed admission, including when too few paired folds exist.
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
                        method = c("two_stage", "full", "covsel"),
                        num_trees = 300L, max_rows = 2000L,
                        metric = c("brier", "auc"),
                        tolerance = 0, n_folds = 3L,
                        seed = NULL, verbose = TRUE, keep = character(0),
                        corcut = 0.7, var_id = NULL, category = NULL,
                        algorithms = c("glm", "gam", "rf"),
                        ncov = NULL, maxncov = 12L, weights = NULL,
                        ensure = NULL, env_vars = NULL) {
  if (!missing(method) && (!is.character(method) || length(method) != 1L ||
                          is.na(method) ||
                          !method %in% c("two_stage", "full", "covsel"))) {
    cli::cli_abort("{.arg method} must be one of 'two_stage', 'full' or 'covsel'.")
  }
  method <- match.arg(method)
  if (missing(metric)) metric <- "brier"
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
      !metric %in% c("brier", "auc")) {
    cli::cli_abort("{.arg metric} must be {.val brier} or {.val auc}.")
  }
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
  if (method != "full" && length(env_vars) < 3L) {
    cli::cli_abort("Screening needs at least three predictors; use method = 'full' for a prespecified smaller set.")
  }
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

  if (identical(method, "covsel")) {
    return(.cast_covsel_select(data, response, env_vars, keep, corcut,
                               var_id, category, algorithms, ncov, maxncov,
                               weights, ensure, seed, verbose))
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
#' over at least two finite paired folds exceeds both `tolerance` and two
#' standard errors of the gains. This is a heuristic, not a significance test.
#' Returns the mean gain (or `NA` when no paired fold exists).
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
  gain <- mean(g)
  if (length(g) < 2L) return(list(admit = FALSE, gain = gain))
  se <- 2 * stats::sd(g) / sqrt(length(g))
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
  # Bound the search cost with a response-stratified subsample: every
  # presence row is kept and only the background rows are evenly thinned
  # (even spacing keeps their geographic spread). A plain systematic
  # subsample would starve the minority class on presence-background data
  # with a large background (e.g. ~20 presences within 10,000 rows) and
  # starve the inner folds of evaluable paired folds.
  if (length(rows) > max_rows) {
    pres <- rows[y_num[rows] == 1]
    bg <- rows[y_num[rows] == 0]
    if (length(pres) >= max_rows) {
      pres <- pres[unique(round(seq(1, length(pres),
                                    length.out = max_rows)))]
      rows <- pres
    } else {
      n_bg <- max(0L, max_rows - length(pres))
      bg <- bg[unique(round(seq(1, length(bg), length.out = n_bg)))]
      rows <- sort(c(pres, bg))
    }
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
    gains <- vapply(verdict, function(v) if (isTRUE(v$admit)) v$gain else -Inf, numeric(1))
    best <- which.max(gains)
    if (!length(best) || !isTRUE(verdict[[best]]$admit)) {
      cli::cli_warn(c(
        "No candidate pair passes the paired admission guard; returning an empty set.",
        "i" = "Admission requires at least two finite paired folds and sufficient gain. An empty set does not establish absence of signal."))
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
        "The single candidate does not pass the paired admission guard; returning an empty set.",
        "i" = "Admission requires at least two finite paired folds and sufficient gain. An empty set does not establish absence of signal."))
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

  step <- max(path$step)
  repeat {
    if (!length(pool)) break
    cand <- lapply(pool, function(v) {
      r <- eval_loss(c(current, v))
      n_fits <<- n_fits + r$fits
      r
    })
    verdict <- lapply(cand, function(r) .cast_admits(cur$folds, r$folds, tolerance, metric))
    gains <- vapply(verdict, function(v) if (isTRUE(v$admit)) v$gain else -Inf, numeric(1))
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

# ========================================================================
# method = "covsel": the two-step embedded selection of Adde et al. (2023)
# Step A (filter): hierarchical collinearity thinning -- one covariate per
#   variable (corcut = 0), then per category (corcut), then globally (corcut),
#   ranking by the smaller p of the two poly(x, 2) logistic terms. Forced
#   covariates (`keep`) bypass thinning and drop their collinear partners.
# Step B (embed): model-specific embedded selection -- GLM elastic net,
#   GAM null-space penalisation, guided regularised RF -- each ranking the
#   survivors; the average-rank table keeps covariates selected by every
#   algorithm first, then the union, subset to ncov = ceiling(log2(n_pres)) - 1
#   (capped at maxncov). If no algorithm retains anything, the filtered set is
#   returned (the sabinaNSDM "nocorr" default).
# ========================================================================

.cast_covsel_select <- function(data, response, env_vars, keep, corcut,
                                var_id, category, algorithms, ncov, maxncov,
                                weights, ensure, seed, verbose) {
  pa <- as.integer(data[[response]])
  w <- if (is.null(weights)) rep(1, length(pa)) else as.numeric(weights)
  assoc <- vapply(env_vars, function(v)
    abs(suppressWarnings(stats::cor(data[[v]], data[[response]],
        use = "pairwise.complete.obs"))), numeric(1))
  assoc[!is.finite(assoc)] <- 0
  stage1_p <- vapply(env_vars, function(v) .cast_stage1_p(data, response, v),
                     numeric(1))

  # ---- Step A: hierarchical collinearity filter --------------------------
  X <- as.data.frame(data[, env_vars, drop = FALSE], check.names = FALSE)
  for (col in names(X)) X[[col]] <- as.numeric(X[[col]])
  thinned <- stats::setNames(rep(FALSE, length(env_vars)), env_vars)
  filter_level <- stats::setNames(rep(NA_character_, length(env_vars)), env_vars)
  pool <- env_vars
  # Level 1 -- one covariate per variable (corcut = 0 is structural: a
  # "representative" is exactly one measurement of the same construct).
  if (!is.null(var_id)) {
    grp <- split(pool, var_id[match(pool, env_vars)])
    out <- .cast_filter_group(grp, X, pa, w, keep, corcut = 0)
    pool <- unlist(out$kept, use.names = FALSE)
    dropped <- setdiff(unlist(grp, use.names = FALSE), pool)
    dropped <- setdiff(dropped, keep)
    thinned[dropped] <- TRUE
    filter_level[dropped] <- "variable"
    if (verbose) {
      cli::cli_inform("covsel Step A (variable level): {length(unlist(grp, use.names = FALSE))} -> {length(pool)} covariates (one per variable).")
    }
  }
  # Level 2 -- per category
  if (!is.null(category)) {
    grp <- split(pool, category[match(pool, env_vars)])
    out <- .cast_filter_group(grp, X, pa, w, keep, corcut = corcut)
    pool <- unlist(out$kept, use.names = FALSE)
    dropped <- setdiff(unlist(grp, use.names = FALSE), pool)
    dropped <- setdiff(dropped, keep)
    thinned[dropped] <- TRUE
    filter_level[dropped] <- "category"
    if (verbose) {
      cli::cli_inform("covsel Step A (category level): -> {length(pool)} covariates (|r| <= {corcut}).")
    }
  }
  # Level 3 -- global remainder
  if (length(pool) > 1L) {
    kept <- .cast_filteralgo(X[, pool, drop = FALSE], pa, w, keep, corcut)
    dropped <- setdiff(pool, kept)
    dropped <- setdiff(dropped, keep)
    thinned[dropped] <- TRUE
    filter_level[dropped] <- "global"
    pool <- c(setdiff(kept, keep), keep)
  }
  filtered <- pool
  if (verbose) {
    cli::cli_inform("covsel Step A done: {length(env_vars)} -> {length(filtered)} non-collinear covariates.")
  }

  # ---- Step B: embedded selection ----------------------------------------
  n_pres <- sum(pa == 1L)
  ncov_use <- if (is.null(ncov)) max(1L, ceiling(log2(n_pres)) - 1L) else as.integer(ncov)
  ncov_use <- min(ncov_use, maxncov)
  embed <- .cast_embed_select(X[, filtered, drop = FALSE], pa, w,
                              algorithms = algorithms, ncov = ncov_use,
                              seed = seed, verbose = verbose)
  if (is.null(embed)) {
    cli::cli_warn(c(
      "No covariate survived the embedded selection; keeping the filtered set.",
      "i" = "This mirrors the sabinaNSDM default (maxncov = 'nocorr')."))
    selected <- filtered
    reason <- stats::setNames(rep("filter", length(selected)), selected)
    ranks_2 <- data.frame(covariate = filtered, rank.f = seq_along(filtered))
  } else {
    selected <- embed$selected
    reason <- stats::setNames(rep("embed", length(selected)), selected)
    ranks_2 <- embed$ranks_2
  }
  reason[keep] <- "prespecified"
  # ---- Coverage requirements (ensure) --------------------------------------
  # A declared flag (e.g. climate, from covariate provenance) that the final
  # set does not cover is topped up with its strongest Step-A-ranked survivor.
  # The rule names no variable: the admitted covariate is whatever the data
  # ranks strongest, so the mechanism is species-agnostic while the outcome
  # adapts. Coverage guarantees the projection is answerable; it is not a
  # claim that the ensured covariate matters (importance reporting decides).
  ensured <- character(0)
  if (!is.null(ensure)) {
    for (nm in names(ensure)) {
      flag <- ensure[[nm]]
      covered <- any(flag[match(union(keep, selected), env_vars)])
      if (covered) next
      cand <- intersect(filtered, env_vars[flag])
      if (!length(cand)) cand <- env_vars[flag]
      if (!length(cand)) {
        cli::cli_warn("ensure[[{.val {nm}]] matches no predictor; coverage cannot be satisfied.")
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

  step_added <- stats::setNames(rep(NA_integer_, length(env_vars)), env_vars)
  rk <- ranks_2$rank.f
  names(rk) <- ranks_2$covariate
  common <- intersect(names(rk), env_vars)
  step_added[common] <- as.integer(rk[common])
  scores <- data.frame(
    variable = env_vars,
    assoc = unname(assoc[env_vars]),
    stage1_p = unname(stage1_p[env_vars]),
    stage1_rank = NA_integer_,
    collinear_thinned = unname(thinned[env_vars]),
    step_added = unname(step_added[env_vars]),
    loss_gain = NA_real_,
    kept_by_design = env_vars %in% keep,
    selected_reason = unname(ifelse(
      env_vars %in% keep, "prespecified",
      ifelse(env_vars %in% selected, unname(reason[env_vars]), "excluded"))),
    selected = env_vars %in% selected,
    stringsAsFactors = FALSE)
  scores$selected_reason[is.na(scores$selected_reason)] <- "excluded"
  new_cast_select(selected = selected, scores = scores, method = "covsel",
                  diagnostics = list(
                    keep = keep, corcut = corcut, algorithms = algorithms,
                    ncov = ncov_use, maxncov = maxncov, ensured = ensured,
                    filter_level = filter_level, filtered = filtered,
                    embed_ranks = if (is.null(embed)) NULL else embed$ranks_1,
                    final_ranks = ranks_2, n_pres = n_pres))
}

# ---- covsel Step A: greedy collinearity thinning (covsel.filteralgo) ------
# Ranks candidates by the smaller p of the two poly(x, 2) logistic terms,
# then greedily keeps the top-ranked covariate and drops everything
# correlated with it beyond `corcut`, repeating on the remainder. Forced
# covariates are kept and their collinear partners removed up front.
.cast_filteralgo <- function(covdata, pa, weights = NULL, force = NULL,
                             corcut = 0.7) {
  candidates <- names(covdata)
  if (ncol(covdata) <= 1L) return(candidates)
  # Drop covariates with fewer than 10 unique values (unless forced): the
  # embedding step needs variation to fit.
  nuniq <- vapply(covdata, function(x) length(unique(stats::na.omit(x))), integer(1))
  drop <- setdiff(candidates[nuniq < 10L], force)
  if (length(drop)) covdata <- covdata[, setdiff(candidates, drop), drop = FALSE]
  if (ncol(covdata) <= 1L) return(names(covdata))
  # Forced covariates: keep, and remove their collinear partners.
  force_in <- intersect(force, names(covdata))
  if (length(force_in)) {
    cm <- abs(stats::cor(covdata, use = "pairwise.complete.obs"))
    partners <- unique(unlist(lapply(force_in, function(f)
      setdiff(names(which(cm[, f] > corcut)), force_in))))
    if (length(partners)) covdata <- covdata[, setdiff(names(covdata), partners), drop = FALSE]
  }
  if (ncol(covdata) <= 1L) return(names(covdata))
  # Rank by the quadratic-logistic signal (same statistic as stage 1).
  res <- sort(vapply(names(covdata), function(v) {
    x <- suppressWarnings(as.numeric(covdata[[v]]))
    ok <- is.finite(x) & !is.na(pa)
    if (sum(ok) < 10L || stats::sd(x[ok]) <= 0) return(NA_real_)
    fit <- suppressWarnings(stats::glm(pa[ok] ~ stats::poly(x[ok], degree = 2),
                                       family = stats::binomial(),
                                       weights = if (is.null(weights)) NULL else weights[ok]))
    co <- suppressWarnings(summary(fit)$coefficients)
    if (is.null(co) || nrow(co) < 3L) return(NA_real_)
    p <- suppressWarnings(min(co[2:3, 4], na.rm = TRUE))
    if (!is.finite(p)) NA_real_ else p
  }, numeric(1)))
  ord <- names(res)   # sort() already orders by p (NA last), ties by name
  covdata <- covdata[, ord, drop = FALSE]
  # Greedy thinning: keep the head, drop its correlated partners, repeat.
  sanctuaire <- character(0)
  cm <- abs(stats::cor(covdata, use = "pairwise.complete.obs"))
  if (all(cm[cm != 1] < corcut)) {
    sanctuaire <- names(covdata)
  } else {
    while (ncol(cm) > 1L) {
      var.name <- colnames(cm)[1]
      sanctuaire <- c(sanctuaire, var.name)
      ix <- which(cm[, 1] <= corcut)
      cm <- cm[ix, ix, drop = FALSE]
      if (all(cm[cm != 1] < corcut)) {
        sanctuaire <- c(sanctuaire, colnames(cm))
        break
      }
    }
    if (!length(sanctuaire)) sanctuaire <- colnames(cm)
  }
  # Forced covariates are always retained.
  c(setdiff(sanctuaire, force), force_in)
}

# Apply .cast_filteralgo per group and re-label the results.
.cast_filter_group <- function(grp, X, pa, weights, force, corcut) {
  kept <- lapply(grp, function(vs) {
    vs <- vs[vs %in% names(X)]
    .cast_filteralgo(X[, vs, drop = FALSE], pa, weights, force, corcut)
  })
  list(kept = kept)
}

# ---- covsel Step B: model-specific embedded selection --------------------
# Three regularised learners rank the filtered covariates; the aggregate
# table prefers covariates selected by every learner, then the union.
.cast_embed_select <- function(covdata, pa, weights = NULL,
                               algorithms = c("glm", "gam", "rf"), ncov = 10L,
                               seed = NULL, verbose = FALSE) {
  algorithms <- tolower(algorithms)
  # glmnet and mgcv are deterministic; the guided RRF is not. Seed once here
  # so the whole embedding (and therefore the selection) is reproducible.
  if (!is.null(seed)) set.seed(seed)
  ranks_1 <- data.frame(covariate = character(0), rank = integer(0),
                        model = character(0), stringsAsFactors = FALSE)
  if ("glm" %in% algorithms) {
    r <- try({
      check_suggested("glmnet", "for the covsel elastic-net embedding")
      form <- stats::as.formula(paste0(
        "as.factor(pa) ~ ",
        paste(paste0("poly(`", names(covdata), "`, 2)"), collapse = " + "), "-1"))
      x <- stats::model.matrix(form, cbind(covdata, pa = pa))
      mdl <- suppressWarnings(glmnet::cv.glmnet(
        x, as.factor(pa), alpha = 0.5, weights = weights, family = "binomial",
        type.measure = "deviance"))
      b <- as.matrix(stats::coef(mdl, s = mdl$lambda.1se))
      b <- data.frame(covariate = rownames(b), coef = abs(as.numeric(b[, 1])))
      b <- b[b$coef > 0 & b$covariate != "(Intercept)", , drop = FALSE]
      if (nrow(b) < 1L) {
        b <- as.matrix(stats::coef(mdl, s = mdl$lambda.min))
        b <- data.frame(covariate = rownames(b), coef = abs(as.numeric(b[, 1])))
        b <- b[b$coef > 0 & b$covariate != "(Intercept)", , drop = FALSE]
      }
      if (nrow(b) >= 1L) {
        b$covariate <- sub("^poly\\(`?(.*)`?,\\s*2\\)[0-9]+$", "\\1", b$covariate)
        b <- b[order(-b$coef), , drop = FALSE]
        b <- b[!duplicated(b$covariate), , drop = FALSE]
        data.frame(covariate = b$covariate, rank = seq_len(nrow(b)),
                   model = "glm", stringsAsFactors = FALSE)
      } else NULL
    }, silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r)) {
      ranks_1 <- rbind(ranks_1, r)
    } else if (verbose) {
      cli::cli_inform("covsel embed: glm branch produced no selection.")
    }
  }
  if ("gam" %in% algorithms) {
    r <- try({
      check_suggested("mgcv", "for the covsel GAM embedding")
      form <- stats::as.formula(paste0(
        "pa ~ ", paste(paste0("s(`", names(covdata), "`, bs = 'cr')"),
                       collapse = " + ")))
      mdl <- suppressWarnings(mgcv::bam(
        form, data = cbind(covdata, pa = pa), weights = weights,
        family = stats::binomial(), method = "fREML", select = TRUE,
        discrete = TRUE))
      tb <- summary(mdl)$s.table
      keep <- rownames(tb)[tb[, "p-value"] < 0.9]
      if (length(keep)) {
        ord <- keep[order(-abs(tb[keep, "Chi.sq"]))]
        ord <- sub("^s\\(`?(.*)`?\\)$", "\\1", ord)
        data.frame(covariate = ord, rank = seq_along(ord), model = "gam",
                   stringsAsFactors = FALSE)
      } else NULL
    }, silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r)) {
      ranks_1 <- rbind(ranks_1, r)
    } else if (verbose) {
      cli::cli_inform("covsel embed: gam branch produced no selection.")
    }
  }
  if ("rf" %in% algorithms) {
    r <- try({
      check_suggested("RRF", "for the covsel guided-RF embedding")
      w <- if (is.null(weights)) rep(1, length(pa)) else weights
      rf0 <- RRF::RRF(covdata, as.factor(pa), flagReg = 0)
      imp <- rf0$importance[, "MeanDecreaseGini"]
      gamma <- 0.5
      coefReg <- (1 - gamma) + gamma * (imp / max(imp))
      mdl <- RRF::RRF(covdata, as.factor(pa),
                      classwt = c("0" = min(w), "1" = max(w)),
                      coefReg = coefReg, flagReg = 1)
      imp2 <- mdl$importance[, "MeanDecreaseGini"]
      keep <- names(imp2)[imp2 > 0]
      if (length(keep)) {
        ord <- keep[order(-imp2[keep])]
        data.frame(covariate = ord, rank = seq_along(ord), model = "rf",
                   stringsAsFactors = FALSE)
      } else NULL
    }, silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r)) {
      ranks_1 <- rbind(ranks_1, r)
    } else if (verbose) {
      cli::cli_inform("covsel embed: rf branch produced no selection.")
    }
  }
  if (!nrow(ranks_1)) return(NULL)
  # Aggregate: covariates selected by every learner first (by summed rank),
  # then the union-ranked remainder.
  n_models <- length(unique(ranks_1$model))
  in_all <- names(which(table(ranks_1$covariate) == n_models))
  parts <- list()
  if (length(in_all)) {
    inter <- ranks_1[ranks_1$covariate %in% in_all, , drop = FALSE]
    inter <- stats::aggregate(rank ~ covariate, data = inter, FUN = sum)
    inter <- inter[order(inter$rank), , drop = FALSE]
    inter$rank.f <- seq_len(nrow(inter))
    parts[["inter"]] <- inter[, c("covariate", "rank.f")]
  }
  rest <- ranks_1[!ranks_1$covariate %in% in_all, , drop = FALSE]
  if (nrow(rest)) {
    rest <- stats::aggregate(rank ~ covariate, data = rest, FUN = sum)
    rest <- rest[order(rest$rank), , drop = FALSE]
    start <- if (length(parts)) max(parts[["inter"]]$rank.f) + 1L else 1L
    rest$rank.f <- seq.int(start, length.out = nrow(rest))
    parts[["rest"]] <- rest[, c("covariate", "rank.f")]
  }
  if (!length(parts)) return(NULL)
  ranks_2 <- do.call(rbind, parts)
  ncov <- min(ncov, nrow(ranks_2))
  sel <- ranks_2$covariate[seq_len(ncov)]
  # Keep only covariates that exist in the input set (label hygiene).
  sel <- intersect(sel, names(covdata))
  list(selected = sel, ranks_1 = ranks_1, ranks_2 = ranks_2)
}
