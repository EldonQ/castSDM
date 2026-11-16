# Exact additive GAM evaluation reuses one-variable term values. Other engines
# and unsupported GAM formulas use their ordinary prediction methods.
.cast_response_predictor <- function(info, reference, clamp, training_reference = reference) {
  terms <- .cast_gam_terms(info)
  if (is.null(terms)) return(function(newdata)
    predict_single_model(info, newdata, clamp = isTRUE(clamp), reference = training_reference))
  if (isTRUE(clamp)) reference <- .cast_clamp(reference, training_reference)
  a <- .cast_predict_gam(info$model, reference, type = "terms")
  constant <- attr(a, "constant") %||% 0
  lookup <- lapply(names(terms), function(v) list(value = reference[[v]],
    contribution = a[, terms[[v]]]))
  names(lookup) <- names(terms)
  function(newdata) {
    if (isTRUE(clamp)) newdata <- .cast_clamp(newdata, training_reference)
    eta <- rep(constant, nrow(newdata))
    for (v in names(terms)) {
      item <- lookup[[v]]
      at <- match(newdata[[v]], item$value)
      missing <- which(is.na(at))
      if (length(missing)) {
        rows <- missing[!duplicated(newdata[[v]][missing])]
        added <- .cast_predict_gam(info$model, newdata[rows, , drop = FALSE],
          type = "terms", terms = terms[[v]])
        item$value <- c(item$value, newdata[[v]][rows])
        item$contribution <- c(item$contribution, added[, terms[[v]]])
        at[missing] <- match(newdata[[v]][missing], item$value)
        lookup[[v]] <<- item
      }
      eta <- eta+item$contribution[at]
    }
    .cast_apply_offset(info$model$family$linkinv(eta), info$offset)
  }
}

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
#' @param type Character. `"ale"` uses the published [iml::FeatureEffect]
#'   implementation of centred accumulated local effects. `"fixed"` (default): all other predictors held at
#'   their training median. `"pdp"`: predictions averaged over the training
#'   background of the other predictors (partial dependence; more expensive).
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
#' Grids span each predictor's marginal range in the supplied background.
#' All supplied reference rows are used; no hidden subsampling is applied.
#' Clamping to training ranges does not guarantee joint environmental support.
#' Partial dependence replaces the profiled predictor while retaining the
#' background values of the others, which can create unsupported combinations
#' when predictors are correlated. Neither profile establishes causality (see
#' the package-level note on interventional effects).
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
                                 type = c("fixed", "ale", "pdp"),
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
  if (!is.character(variables) || !length(variables) || anyNA(variables) ||
      anyDuplicated(variables) || !all(variables %in% env_vars)) {
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
  if (!length(models)) cli::cli_abort("No fitted model is available for response curves.")
  # A fixed profile needs only one median row to initialise its term cache.
  reference <- if (type == "fixed") as.data.frame(lapply(bg, stats::median, na.rm = TRUE)) else bg
  predictors <- stats::setNames(lapply(models, function(m)
    .cast_response_predictor(fit$models[[m]], reference, clamp,
      training_reference = fit$scaling$reference)), models)

  .predict_mean_or_vec <- function(m, X) {
    value <- tryCatch(
      predictors[[m]](X),
      error = function(e) cli::cli_abort("Response prediction failed for {.val {m}}: {conditionMessage(e)}")
    )
    if (length(value) != nrow(X) || any(!is.finite(value)))
      cli::cli_abort("Response prediction for {.val {m}} returned missing or invalid scores.")
    value
  }

  if (type == "ale" && !is.null(pair)) cli::cli_abort("ALE currently supports one-variable curves; pair grids use fixed profiles.")
  if (type == "ale") {
    check_suggested("iml", "for the published accumulated local effects implementation")
    ale_predictors <- lapply(models, function(m) iml::Predictor$new(
      model = fit$models[[m]], data = bg,
      predict.function = function(model, newdata) .predict_mean_or_vec(m, newdata)))
    names(ale_predictors) <- models
    out_list <- lapply(variables, function(v) {
      do.call(rbind, lapply(models, function(m) {
        if (!.cast_is_varying(bg[[v]])) return(data.frame(variable = v,
          value = bg[[v]][1L], model = m, prediction = NA_real_))
        effect <- iml::FeatureEffect$new(ale_predictors[[m]], feature = v,
          method = "ale", grid.size = grid_size)$results
        data.frame(variable = v, value = effect[[v]], model = m,
          prediction = effect$.value, stringsAsFactors = FALSE)
      }))
    })
  } else if (!is.null(pair)) {
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
      if (type == "fixed") {
        # Fixed profile: others at their training medians.
        X <- bg[rep(1L, length(g)), , drop = FALSE]
        for (o in other) X[[o]] <- stats::median(bg[[o]], na.rm = TRUE)
        X[[v]] <- g
        # Predict the whole grid once per engine instead of repeating each
        # value across a grid-sized frame and then averaging identical rows.
        pmat <- do.call(cbind, lapply(models, function(m) {
          .predict_mean_or_vec(m, X)
        }))
      } else {
        # Partial dependence: replace the profiled column in the background
        # and average predictions per engine at each grid value.
        preds <- lapply(g, function(val) {
          X <- bg
          X[[v]] <- val
          vapply(models, function(m) {
            mean(.predict_mean_or_vec(m, X))
          }, numeric(1))
        })
        pmat <- do.call(rbind, preds)
      }
      data.frame(
        variable = v, value = g,
        model = rep(models, each = length(g)),
        prediction = as.numeric(pmat),  # model columns, each in grid order
        stringsAsFactors = FALSE
      )
    })
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  class(out) <- c("cast_response", "data.frame")
  attr(out, "response_type") <- type
  attr(out, "reference_n") <- nrow(bg)
  out
}

#' Plot Response Curves
#'
#' Univariate curves are drawn per variable (colour = engine, free y
#' scales); bivariate grids are drawn as raster heatmaps per engine.
#'
#' @param x A `cast_response` object from [cast_response_curves()].
#' @param ... Unused.
#' @param ensemble Whether to show an equal-weight mean and cross-model range.
#'   Default `FALSE` draws individual model curves. The range is not a confidence interval.
#' @return A `ggplot` object.
#' @export
plot.cast_response <- function(x, ..., ensemble = FALSE) {
  check_suggested("ggplot2", "for plotting")
  ale <- identical(attr(x, "response_type"), "ale")
  y_label <- if (ale) "Accumulated local effect (centred score)" else "Predicted suitability"
  if ("variable" %in% names(x)) {
    if (ensemble) {
      d <- do.call(rbind, lapply(split(x, x$variable), function(z) {
        values <- sort(unique(z$value))
        do.call(rbind, lapply(values, function(v) {
          p <- z$prediction[z$value == v]; p <- p[is.finite(p)]
          data.frame(variable = z$variable[1], value = v,
            mean = if (length(p)) mean(p) else NA_real_,
            lo = if (length(p)) min(p) else NA_real_, hi = if (length(p)) max(p) else NA_real_)
        }))
      }))
      return(ggplot2::ggplot(d, ggplot2::aes(.data$value, .data$mean)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), fill = "#9FD3DB", alpha = .4) +
        ggplot2::geom_line(colour = "#256E92", linewidth = .5) +
        ggplot2::facet_wrap(~ .data$variable, scales = "free_x", ncol = 2,
          labeller = ggplot2::as_labeller(.cast_var_label)) +
        ggplot2::labs(x = "Predictor value", y = y_label, title = if (ale) "Accumulated local effects" else "Model response curves",
          caption = "Equal-weight model mean; band = model range, not a confidence interval. No smoothing.") +
        theme_cast(base_size = 8))
    }
    ggplot2::ggplot(
      x,
      ggplot2::aes(x = .data$value, y = .data$prediction,
                   color = .data$model)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~ .data$variable, scales = "free") +
      ggplot2::labs(
        x = NULL, y = y_label,
        color = "Model", title = if (ale) "Accumulated local effects" else "Response curves"
      ) +
      ggplot2::scale_colour_manual(values = .cast_engine_colors(unique(x$model))) +
      theme_cast(base_size = 8)
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
