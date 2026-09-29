#' Response Curves and Partial Dependence for Fitted SDMs
#'
#' Single-variable response curves for a fitted [cast_fit()] object, either
#' with the other predictors held fixed at their training median, or as
#' partial dependence averaged over the training background (Friedman 2001).
#' An optional two-variable grid exposes pairwise interaction surfaces.
#' The result is a tidy long-format data frame with a [plot()] method.
#'
#' @param fit A `cast_fit` object from [cast_fit()].
#' @param variables Character. Predictor names to profile. Default `NULL`:
#'   every predictor stored in the fit.
#' @param data Optional `data.frame` supplying the background distribution
#'   used for the fixed values (medians) and for `type = "pdp"`
#'   marginalization. Default `NULL`: the training data stored in the fit.
#'   Only the fitted predictors are used; missing values are imputed with
#'   the training medians.
#' @param type Character. `"fixed"` (default): all other predictors held at
#'   their training median. `"pdp"`: predictions averaged over the training
#'   background of the other predictors (partial dependence; more honest
#'   about correlations, and more expensive).
#' @param grid_size Integer >= 2. Points per variable curve. Default `50`.
#' @param pair Optional length-2 character vector of distinct predictors.
#'   When supplied, a two-variable interaction grid is produced instead of
#'   individual curves.
#' @param pair_size Integer >= 2. Grid points per axis of the bivariate
#'   grid. Default `25`.
#' @param models Character. Subset of fitted engines to profile. Default
#'   `NULL`: all fitted engines.
#' @param clamp Logical. Clamp predictions to the training environmental
#'   range (see [cast_predict()]). Default `TRUE`.
#'
#' @return A `cast_response` object (long-format `data.frame`).
#'   Univariate: columns `variable`, `value`, `model`, `prediction`.
#'   Bivariate: `var1`, `value1`, `var2`, `value2`, `model`, `prediction`.
#'
#' @details
#' Grids span the observed training range of each predictor, so with the
#' default `clamp = TRUE` no curve section relies on extrapolation.
#' Partial-dependence curves average over the joint training distribution,
#' which for correlated predictors yields more realistic statements than
#' the fixed-profile variant; neither establishes causality (see the
#' package-level note on interventional effects).
#'
#' @references
#' Friedman, J. H. (2001). Greedy function approximation: a gradient
#' boosting machine. *Annals of Statistics*, 29(5), 1189-1232.
#'
#' @seealso [cast_fit()], [cast_predict()]
#'
#' @export
cast_response_curves <- function(fit,
                                 variables = NULL,
                                 data = NULL,
                                 type = c("fixed", "pdp"),
                                 grid_size = 50L,
                                 pair = NULL,
                                 pair_size = 25L,
                                 models = NULL,
                                 clamp = TRUE) {
  if (!inherits(fit, "cast_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls cast_fit} object from {.fn cast_fit}.")
  }
  type <- match.arg(type)
  grid_size <- as.integer(grid_size)
  pair_size <- as.integer(pair_size)
  if (is.na(grid_size) || grid_size < 2L) {
    cli::cli_abort("{.arg grid_size} must be an integer >= 2.")
  }
  if (is.na(pair_size) || pair_size < 2L) {
    cli::cli_abort("{.arg pair_size} must be an integer >= 2.")
  }

  env_vars <- fit$env_vars
  if (is.null(variables)) variables <- env_vars
  if (!is.character(variables) || !all(variables %in% env_vars)) {
    cli::cli_abort(c(
      "All {.arg variables} must be predictors of the fitted model.",
      x = "Unknown: {.val {setdiff(variables, env_vars)}}."
    ))
  }
  if (!is.null(pair)) {
    if (!is.character(pair) || length(pair) != 2L || pair[1] == pair[2] ||
        !all(pair %in% env_vars)) {
      cli::cli_abort(
        "{.arg pair} must be two distinct predictors of the fitted model."
      )
    }
  }

  # Background rows: the imputed training data stored in the fit, or a
  # caller-supplied frame imputed on the same training statistics.
  bg <- if (is.null(data)) {
    fit$scaling$reference
  } else {
    X <- as.data.frame(data[, env_vars, drop = FALSE], check.names = FALSE)
    .cast_check_numeric_predictors(X, arg = "data")
    for (col in names(X)) X[[col]] <- as.numeric(X[[col]])
    .cast_impute(X, fit$scaling$impute)
  }
  # Deterministic cap on the pdp marginalization sample: equidistant rows
  # keep the average reproducible without a seed.
  pdp_cap <- 2000L
  if (type == "pdp" && nrow(bg) > pdp_cap) {
    idx <- unique(floor(seq(1, nrow(bg), length.out = pdp_cap)))
    bg <- bg[idx, , drop = FALSE]
  }

  available <- names(fit$models)
  if (is.null(models)) models <- available
  bad_m <- setdiff(models, available)
  if (length(bad_m)) {
    cli::cli_abort(c(
      "{.arg models} must be a subset of the fitted engines.",
      x = "Unknown: {.val {bad_m}}; fitted: {.val {available}}."
    ))
  }
  models <- models[vapply(models, function(m) {
    !is.null(fit$models[[m]]$model)
  }, logical(1))]

  .predict_mean_or_vec <- function(m, X) {
    tryCatch(
      predict_single_model(fit$models[[m]], X, clamp = isTRUE(clamp)),
      error = function(e) rep(NA_real_, nrow(X))
    )
  }

  if (!is.null(pair)) {
    g1 <- seq(min(bg[[pair[1]]], na.rm = TRUE),
              max(bg[[pair[1]]], na.rm = TRUE), length.out = pair_size)
    g2 <- seq(min(bg[[pair[2]]], na.rm = TRUE),
              max(bg[[pair[2]]], na.rm = TRUE), length.out = pair_size)
    grid <- expand.grid(value1 = g1, value2 = g2)
    X <- bg[rep(1L, nrow(grid)), , drop = FALSE]
    # Fixed base: every non-pair predictor at its training median.
    for (o in setdiff(env_vars, pair)) {
      X[[o]] <- stats::median(bg[[o]], na.rm = TRUE)
    }
    X[[pair[1]]] <- grid$value1
    X[[pair[2]]] <- grid$value2
    out_list <- lapply(models, function(m) {
      p <- .predict_mean_or_vec(m, X)
      data.frame(var1 = pair[1], value1 = grid$value1,
                 var2 = pair[2], value2 = grid$value2,
                 model = m, prediction = p, stringsAsFactors = FALSE)
    })
  } else {
    out_list <- lapply(variables, function(v) {
      vals <- bg[[v]]
      g <- seq(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE),
               length.out = grid_size)
      other <- setdiff(env_vars, v)
      base <- if (type == "fixed") {
        # Fixed profile: others at their training medians.
        X <- bg[rep(1L, length(g)), , drop = FALSE]
        for (o in other) X[[o]] <- stats::median(bg[[o]], na.rm = TRUE)
        X
      } else {
        NULL
      }
      preds <- lapply(g, function(val) {
        if (type == "fixed") {
          X <- base
          X[[v]] <- val
          p <- vapply(models, function(m) {
            mean(.predict_mean_or_vec(m, X))
          }, numeric(1))
          stats::setNames(p, models)
        } else {
          # Partial dependence: swap the profiled column into the whole
          # background and average predictions per engine.
          X <- bg
          X[[v]] <- val
          vapply(models, function(m) {
            mean(.predict_mean_or_vec(m, X))
          }, numeric(1))
        }
      })
      pmat <- do.call(rbind, preds)  # grid_size x models
      data.frame(
        variable = v, value = g,
        model = rep(models, each = length(g)),
        prediction = as.numeric(t(pmat)),
        stringsAsFactors = FALSE
      )
    })
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  class(out) <- c("cast_response", "data.frame")
  out
}

#' Plot Response Curves
#'
#' Univariate curves are drawn per variable (colour = engine, free y
#' scales); bivariate grids are drawn as raster heatmaps per engine.
#'
#' @param x A `cast_response` object from [cast_response_curves()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.cast_response <- function(x, ...) {
  check_suggested("ggplot2", "for plotting")
  if ("variable" %in% names(x)) {
    ggplot2::ggplot(
      x,
      ggplot2::aes(x = .data$value, y = .data$prediction,
                   color = .data$model)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~ .data$variable, scales = "free") +
      ggplot2::labs(
        x = NULL, y = "Predicted suitability",
        color = "Model", title = "Response curves"
      ) +
      ggplot2::theme_bw()
  } else {
    ggplot2::ggplot(
      x,
      ggplot2::aes(x = .data$value1, y = .data$value2, fill = .data$prediction)
    ) +
      ggplot2::geom_raster() +
      ggplot2::facet_wrap(~ .data$model) +
      ggplot2::scale_fill_viridis_c(name = "Suitability") +
      ggplot2::labs(
        x = x$var1[1], y = x$var2[1], title = "Interaction surface"
      ) +
      ggplot2::theme_bw()
  }
}
