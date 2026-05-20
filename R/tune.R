#' Unified Per-Algorithm Hyperparameter Grid Search
#'
#' Runs split-sample evaluation over a user-supplied (or default) grid for
#' each requested algorithm and returns a `cast_tune` object that
#' [cast_fit()] consumes via its `tune_result` argument. Selection is
#' scored on `AUC + TSS - 1` (clamped to `[0, 2]`) averaged across splits;
#' ties broken by AUC.
#'
#' Supported algorithms (`models`):
#' \itemize{
#'   \item `"cast"` (CI-MLP): tunes `hidden_size`, `dropout`, `lr`.
#'   \item `"rf"`: tunes `rf_ntree`, plus internally `mtry` defaults via
#'     `ranger`.
#'   \item `"brt"`: tunes `brt_n_trees`, `brt_depth`, learning rate.
#'   \item `"gam"`: no hyper-grid (single fit, included for completeness).
#'   \item `"maxent"`: no hyper-grid via `maxnet` defaults.
#' }
#'
#' Each grid point is evaluated on `n_splits` random train/validation
#' splits of `val_fraction` size. Evaluation uses the same internal AUC
#' + TSS scorer as [cast_evaluate()], skipping CBI to keep tuning fast.
#'
#' @param data Training `data.frame`. Must contain `response`.
#' @param models Character vector. Algorithms to tune. Default
#'   `c("cast", "rf", "brt")`.
#' @param screen,dag,ate As in [cast_fit()] (only required when tuning
#'   `"cast"`).
#' @param grid Optional named list of per-model parameter grids. Each
#'   element is a `data.frame` (or list-of-vectors that gets passed to
#'   `expand.grid`). If `NULL`, uses [cast_default_grid()].
#' @param n_splits Integer. Random splits per grid point. Default `30L`
#'   (raise toward `100` for paper-grade tuning).
#' @param val_fraction Numeric. Validation hold-out per split. Default
#'   `0.25`.
#' @param response Character. Default `"presence"`.
#' @param seed Integer or `NULL`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A `cast_tune` object: list of best parameters per model plus
#'   the full grid scores.
#'
#' @seealso [cast_fit()], [cast_default_grid()]
#' @export
cast_tune <- function(data,
                      models       = c("cast", "rf", "brt"),
                      screen       = NULL,
                      dag          = NULL,
                      ate          = NULL,
                      grid         = NULL,
                      n_splits     = 30L,
                      val_fraction = 0.25,
                      response     = "presence",
                      seed         = NULL,
                      verbose      = TRUE) {
  models <- tolower(models)
  if (!is.null(seed)) set.seed(seed)
  grid <- grid %||% cast_default_grid()
  Y <- data[[response]]
  n  <- nrow(data)

  # Helper: produce n_splits stratified train/val index pairs.
  pos <- which(Y == 1L); neg <- which(Y == 0L)
  splits <- lapply(seq_len(n_splits), function(i) {
    vp <- sample(pos, max(1L, round(val_fraction * length(pos))))
    vn <- sample(neg, max(1L, round(val_fraction * length(neg))))
    list(val = c(vp, vn))
  })

  best <- list(); all_scores <- list()

  for (mdl in models) {
    if (verbose) cli::cli_h2("Tuning {.val {mdl}}")
    g <- grid[[mdl]]
    if (is.null(g)) {
      if (verbose) cli::cli_inform("No grid for {.val {mdl}}; using cast_fit defaults.")
      best[[mdl]] <- list()
      next
    }
    if (is.list(g) && !is.data.frame(g)) g <- expand.grid(g, stringsAsFactors = FALSE)

    scores <- numeric(nrow(g))
    aucs   <- numeric(nrow(g))
    for (i in seq_len(nrow(g))) {
      params <- as.list(g[i, , drop = FALSE])
      ev <- vapply(splits, function(sp) {
        tr_idx <- setdiff(seq_len(n), sp$val)
        tr <- data[tr_idx, , drop = FALSE]
        va <- data[sp$val, , drop = FALSE]
        y_va <- Y[sp$val]
        pv <- tryCatch(
          .cast_tune_fit_predict(mdl, tr, va, params, screen, dag, ate,
                                 response = response),
          error = function(e) rep(NA_real_, length(y_va))
        )
        if (all(is.na(pv))) return(c(NA_real_, NA_real_))
        m <- evaluate_model_full(pv, y_va)
        c(m["auc"], m["auc"] + m["tss"] - 1)
      }, numeric(2))
      aucs[i]   <- mean(ev[1, ], na.rm = TRUE)
      scores[i] <- mean(ev[2, ], na.rm = TRUE)
      if (verbose)
        cli::cli_inform(
          "  [{i}/{nrow(g)}] {.val {paste(names(params), unlist(params), sep='=', collapse=', ')}}: score={signif(scores[i],3)} (AUC={signif(aucs[i],3)})"
        )
    }
    g$score <- scores; g$auc <- aucs
    all_scores[[mdl]] <- g

    ord <- order(scores, aucs, decreasing = TRUE, na.last = TRUE)
    best_row <- g[ord[1L], , drop = FALSE]
    best[[mdl]] <- as.list(best_row[setdiff(names(best_row),
                                            c("score", "auc"))])
    if (verbose) cli::cli_inform(
      "  -> best {.val {mdl}}: {.val {paste(names(best[[mdl]]), unlist(best[[mdl]]), sep='=', collapse=', ')}}"
    )
  }

  out <- list(best = best, scores = all_scores,
              n_splits = n_splits, models = models)
  class(out) <- "cast_tune"
  out
}


#' Default Hyperparameter Grids for [cast_tune()]
#'
#' Reasonable defaults that finish in minutes for typical SDM datasets.
#' Override by passing your own grid to `cast_tune(grid = ...)`.
#'
#' @return Named list of `data.frame`s (or list-of-vectors).
#' @export
cast_default_grid <- function() {
  list(
    cast = expand.grid(
      hidden_size = c(32L, 64L, 128L),
      dropout     = c(0.1, 0.2, 0.3),
      lr          = c(5e-4, 1e-3, 2e-3),
      stringsAsFactors = FALSE
    ),
    rf = expand.grid(
      rf_ntree = c(200L, 500L, 1000L),
      stringsAsFactors = FALSE
    ),
    brt = expand.grid(
      brt_n_trees = c(300L, 700L, 1500L),
      brt_depth   = c(3L, 5L, 7L),
      stringsAsFactors = FALSE
    )
  )
}


#' @export
print.cast_tune <- function(x, ...) {
  cat("<cast_tune>\n")
  cat("  models   :", paste(x$models, collapse = ", "), "\n")
  cat("  splits   :", x$n_splits, "\n")
  cat("  best     :\n")
  for (m in names(x$best)) {
    p <- x$best[[m]]
    if (length(p) == 0) {
      cat("    -", m, ": (defaults)\n")
    } else {
      cat("    -", m, ":",
          paste(names(p), unlist(p), sep = "=", collapse = ", "), "\n")
    }
  }
  invisible(x)
}


# Internal: fit a single algo with given params, return val predictions.
.cast_tune_fit_predict <- function(mdl, tr, va, params,
                                   screen, dag, ate, response) {
  fit_args <- list(
    data     = tr,
    models   = mdl,
    response = response,
    n_runs   = 1L,
    n_epochs = 100L,
    verbose  = FALSE
  )
  if (mdl == "cast") {
    fit_args$screen <- screen
    fit_args$dag    <- dag
    fit_args$ate    <- ate
  }
  fit_args <- utils::modifyList(fit_args, params)
  f <- do.call(cast_fit, fit_args)
  ev <- cast_evaluate(f, va, response = response)
  # cast_evaluate returns a data.frame; we re-predict directly to get raw probs:
  mdl_info <- f$models[[mdl]]
  X_raw <- as.data.frame(va[, f$env_vars, drop = FALSE])
  X_raw[is.na(X_raw)] <- 0
  X_sc <- as.data.frame(
    scale(X_raw, center = f$scaling$means, scale = f$scaling$sds)
  )
  X_sc[is.na(X_sc)] <- 0
  predict_single_model(mdl_info, X_raw, X_sc, f$screen, f$dag, f$ate)
}
