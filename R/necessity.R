# ==========================================================================
# Necessity (knockout) diagnostic.
# cast_effect_table() answers "does the model respond when I intervene on
# driver X?" (sensitivity). This file answers the complementary question
# "does the model still perform without driver X?" (necessity). A driver
# with a large interventional effect but no knockout cost is replaceable by
# a collinear partner, so attribution to it is not identified. Report both.
# ==========================================================================

.necessity_vars <- function(data, screen, variables, response) {
  if (!is.null(variables)) return(as.character(variables))
  if (is.null(screen)) return(get_env_vars(data, response))
  if (is.character(screen)) return(screen)
  sel <- if (inherits(screen, "cast_select")) {
    screen$selected
  } else if (inherits(screen, "cast_fit")) {
    screen$env_vars
  } else if (inherits(screen, "cast_result")) {
    screen$fit$env_vars %||% screen$screen$selected
  } else {
    NULL
  }
  if (is.null(sel)) {
    cli::cli_abort(c(
      "{.arg screen} must be a {.cls cast_select}, {.cls cast_fit}, {.cls cast_result}, or a character vector.",
      i = "Pass {.code screen = NULL} to use every predictor in {.arg data}."))
  }
  as.character(sel)
}

.necessity_auc <- function(train, test, vars, response, num_trees, seed,
                           num_threads = 1L) {
  if (!length(vars)) return(NA_real_)
  if (length(unique(train[[response]])) < 2L) return(NA_real_)
  rf <- tryCatch(
    ranger::ranger(x = train[, vars, drop = FALSE],
                   y = factor(train[[response]]),
                   probability = TRUE, num.trees = num_trees,
                   seed = seed, num.threads = num_threads),
    error = function(e) NULL)
  if (is.null(rf)) return(NA_real_)
  pp <- stats::predict(rf, test[, vars, drop = FALSE])$predictions
  p <- if (is.matrix(pp)) {
    if ("1" %in% colnames(pp)) pp[, "1"] else pp[, ncol(pp)]
  } else as.numeric(pp)
  compute_auc(as.integer(test[[response]]), as.numeric(p))
}

#' Necessity of Each Driver by Spatial Knockout
#'
#' Refits the model without each driver in turn and measures how much
#' held-out discrimination is lost. Folds are spatial, and both the full and
#' the knocked-out model are refitted inside every training fold, so the
#' loss is an out-of-sample cost, not an in-sample importance score.
#'
#' @section Read with the sensitivity diagnostic:
#' Necessity and sensitivity answer different questions and can disagree.
#' This diagnostic always refits a random forest, even when the supplied
#' screen came from a multi-engine fit. It measures RF predictive necessity,
#' not necessity for that ensemble. A small cost can reflect redundancy,
#' limited power, estimator choice or the discrimination metric. A positive
#' cost does not identify a causal driver. Report effects and costs together,
#' with their estimator and uncertainty. A fixed predictor set selected using
#' these same rows makes this a conditional, post-selection diagnostic; use
#' an independently specified set for confirmatory evaluation.
#'
#' @param data Data frame with `lon`, `lat`, the binary response and the
#'   predictors.
#' @param screen Which predictors to test: a `cast_select` object (its
#'   `selected` set), a `cast_fit` / `cast_result` (its fitted predictors),
#'   a character vector, or `NULL` to use every predictor in `data`.
#' @param variables Optional character vector overriding `screen`.
#' @param response Binary response column. Default `"presence"`.
#' @param k Number of spatial folds. Default 5.
#' @param block_method Spatial blocking passed to the fold builder:
#'   `"grid"` (default), `"grid_random"`, or `"cluster"`. Ignored when
#'   `folds` is supplied.
#' @param folds Optional integer vector of length `nrow(data)` assigning each
#'   row to a fold. When supplied, it overrides internal fold construction and
#'   `block_method`, so a caller with pre-computed or frozen spatial folds
#'   (e.g. a pre-registered protocol) can knock out on exactly those folds.
#' @param num_trees Trees per random forest. Default 300.
#' @param num_threads Threads passed to `ranger` for each fit. Default 1,
#'   which keeps results bit-for-bit reproducible. Raise it to parallelise
#'   the many refits (one full model plus one per driver, per fold); with a
#'   fixed `seed`, `ranger` stays deterministic across thread counts.
#' @param seed Random seed.
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return A `cast_necessity` object. `necessity` is a data.frame with one
#'   row per driver: `mean_dAUC` (mean held-out AUC lost by dropping it),
#'   `sd_dAUC`, `min_dAUC`, `max_dAUC`, `pct_folds_positive`, `n_folds`,
#'   and `necessary` (`mean_dAUC > 0`). `pct_folds_positive` is there so a
#'   stricter rule can be applied without refitting.
#'   `necessary` is a legacy descriptive positive-mean flag, not a test of
#'   statistical significance or causal identification.
#' @seealso [cast_effect_table()], [cast_effect_map()]
#' @export
cast_necessity <- function(data, screen = NULL, variables = NULL,
                           response = "presence", k = 5L,
                           block_method = c("grid", "grid_random", "cluster"),
                           folds = NULL,
                           num_trees = 300L, num_threads = 1L,
                           seed = NULL, verbose = TRUE) {
  block_method <- match.arg(block_method)
  check_suggested("ranger", "for the necessity knockout diagnostic")
  validate_species_data(data, required_cols = c("lon", "lat", response),
                        response = response)
  vars <- .necessity_vars(data, screen, variables, response)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    cli::cli_abort("{.arg data} is missing predictor{?s}: {.val {missing_vars}}.")
  }
  if (length(vars) < 2L) {
    cli::cli_abort("Knockout needs at least two predictors; got {length(vars)}.")
  }
  num_trees <- as.integer(num_trees)
  num_threads <- as.integer(num_threads)

  if (!is.null(folds)) {
    if (length(folds) != nrow(data)) {
      cli::cli_abort("{.arg folds} must have one entry per row of {.arg data} ({nrow(data)}).")
    }
    folds <- as.integer(folds)
    if (anyNA(folds) || length(unique(folds)) < 2L) {
      cli::cli_abort("{.arg folds} must assign rows to at least two folds without missing values.")
    }
    block_method <- "custom"
  } else {
    k <- as.integer(k)
    if (k < 2L) cli::cli_abort("{.arg k} must be at least 2.")
    folds <- make_spatial_folds(data$lon, data$lat, k, block_method, seed)
  }

  fold_ids <- sort(unique(folds))
  auc_full <- stats::setNames(rep(NA_real_, length(fold_ids)), fold_ids)
  dauc <- matrix(NA_real_, nrow = length(vars), ncol = length(fold_ids),
                 dimnames = list(vars, as.character(fold_ids)))
  skipped <- integer(0)

  for (j in seq_along(fold_ids)) {
    f <- fold_ids[[j]]
    train <- data[folds != f, , drop = FALSE]
    test <- data[folds == f, , drop = FALSE]
    if (length(unique(train[[response]])) < 2L ||
        length(unique(test[[response]])) < 2L) {
      skipped <- c(skipped, f)
      next
    }
    fold_seed <- if (is.null(seed)) NULL else seed + f
    if (verbose) cli::cli_inform("fold {j}/{length(fold_ids)}: full model + {length(vars)} knockouts...")
    auc_full[[j]] <- .necessity_auc(train, test, vars, response, num_trees,
                                   fold_seed %||% 1L, num_threads = num_threads)
    if (!is.finite(auc_full[[j]])) next
    for (v in vars) {
      a <- .necessity_auc(train, test, setdiff(vars, v), response, num_trees,
                          fold_seed %||% 1L, num_threads = num_threads)
      dauc[v, j] <- auc_full[[j]] - a
    }
  }
  if (length(skipped)) {
    n_skipped <- length(skipped)
    cli::cli_warn("{n_skipped} fold{?s} had a single response class and {?was/were} skipped: {.val {skipped}}.")
  }
  if (all(!is.finite(dauc))) {
    cli::cli_abort("No fold produced a usable knockout comparison; the diagnostic is undefined here.")
  }

  necessity <- do.call(rbind, lapply(vars, function(v) {
    d <- dauc[v, ]
    d <- d[is.finite(d)]
    if (!length(d)) {
      return(data.frame(variable = v, mean_dAUC = NA_real_, sd_dAUC = NA_real_,
                        min_dAUC = NA_real_, max_dAUC = NA_real_,
                        pct_folds_positive = NA_real_, n_folds = 0L,
                        necessary = NA, stringsAsFactors = FALSE))
    }
    data.frame(variable = v, mean_dAUC = mean(d),
               sd_dAUC = if (length(d) > 1L) stats::sd(d) else NA_real_,
               min_dAUC = min(d), max_dAUC = max(d),
               pct_folds_positive = mean(d > 0), n_folds = length(d),
               necessary = mean(d) > 0, stringsAsFactors = FALSE)
  }))
  necessity <- necessity[order(-necessity$mean_dAUC), , drop = FALSE]
  rownames(necessity) <- NULL

  new_cast_necessity(
    necessity = necessity,
    fold_dauc = dauc,
    auc_full = auc_full,
    folds = folds,
    k = length(fold_ids),
    block_method = block_method,
    diagnostics = list(variables = vars, num_trees = num_trees,
                       skipped_folds = skipped)
  )
}

#' @export
print.cast_necessity <- function(x, ...) {
  cli::cli_text("{.strong Necessity (knockout) diagnostic}: {x$k} spatial fold{?s}, {.val {x$block_method}} blocking")
  cli::cli_text("dAUC = held-out AUC of the full model minus the model without that driver.")
  print(as.data.frame(x$necessity))
  cli::cli_text("RF predictive diagnostic; positive dAUC is not a causal identification test.")
  invisible(x)
}
