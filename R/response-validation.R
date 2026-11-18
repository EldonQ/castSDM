# Reuse each spatial training fit while it is live. No additional learner or
# fitted model is introduced; evaluation locations and scenario values are fixed.
.cast_response_validation_data <- function(x, candidates) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || !all(c("current", "future", "drivers", "scenario") %in% names(x)))
    cli::cli_abort("response_data needs current, future, drivers and scenario.")
  if (!is.data.frame(x$current) || !is.data.frame(x$future) ||
      !nrow(x$current) || nrow(x$current) != nrow(x$future) ||
      !all(candidates %in% names(x$current)) || !all(candidates %in% names(x$future)))
    cli::cli_abort("response_data must pair the same locations and include every candidate predictor, including predictors selected only inside a fold.")
  .cast_check_numeric_predictors(x$current[, candidates, drop = FALSE], "response_data$current")
  .cast_check_numeric_predictors(x$future[, candidates, drop = FALSE], "response_data$future")
  if (!is.character(x$drivers) || anyNA(x$drivers) || anyDuplicated(x$drivers) ||
      !all(x$drivers %in% candidates) || !is.character(x$scenario) ||
      length(x$scenario) != 1L || is.na(x$scenario) || !nzchar(x$scenario))
    cli::cli_abort("response_data needs unique candidate driver names and one scenario identifier.")
  x$support_probs <- x$support_probs %||% c(.01, .99)
  if (!is.numeric(x$support_probs) || length(x$support_probs) != 2L ||
      any(!is.finite(x$support_probs)) || any(x$support_probs < 0 | x$support_probs > 1) ||
      x$support_probs[1L] >= x$support_probs[2L])
    cli::cli_abort("response_data$support_probs must be two increasing probabilities.")
  x$current <- x$current[, candidates, drop = FALSE]
  x$future <- x$future[, candidates, drop = FALSE]
  x
}

.cast_response_validation_fold <- function(fit, x, models, clamp) {
  a <- x$current[, fit$env_vars, drop = FALSE]
  b <- x$future[, fit$env_vars, drop = FALSE]
  ok <- .cast_finite_rows(a) & .cast_finite_rows(b)
  contrasts <- c("total", if (length(x$drivers)) paste0("only:", x$drivers),
    if (length(x$drivers)) paste0("retained:", x$drivers))
  values <- array(NA_real_, c(nrow(a), length(models), length(contrasts)),
    dimnames = list(NULL, models, contrasts))
  support <- rep(FALSE, nrow(a))
  observed_range <- support
  if (any(ok)) {
    a <- a[ok, , drop = FALSE]; b <- b[ok, , drop = FALSE]
    bounds <- .cast_support_bounds(fit$scaling$reference, x$support_probs)
    support[ok] <- .cast_attribution_inside(a, bounds) & .cast_attribution_inside(b, bounds)
    range_bounds <- .cast_support_bounds(fit$scaling$reference, c(0, 1))
    observed_range[ok] <- .cast_attribution_inside(a, range_bounds) & .cast_attribution_inside(b, range_bounds)
    # Missing engines are errors here, never silently reweighted.
    context <- .cast_prediction_context(fit, stats::setNames(rep(1, length(models)), models), a, b, clamp)
    base <- context$predict(); future <- context$predict(names(b))
    values[ok, , "total"] <- future-base
    for (v in x$drivers) {
      selected <- v %in% fit$env_vars
      values[ok, , paste0("only:", v)] <- if (selected) context$predict(v)-base else 0
      values[ok, , paste0("retained:", v)] <- if (selected)
        future-context$predict(setdiff(names(b), v)) else 0
    }
  }
  list(values = values, supported = support, observed_range = observed_range,
    selected = x$drivers %in% fit$env_vars)
}

.cast_response_validation_summary <- function(cv, weights) {
  validation <- cv$response_validation
  if (is.null(validation)) return(NULL)
  runs <- validation$runs
  expected <- cv$k * cv$n_repeat
  if (length(runs) != expected || any(vapply(runs, is.null, logical(1))))
    return(data.frame(scenario = validation$scenario, status = "incomplete_spatial_fits",
      n_scheduled = expected, n_available = sum(!vapply(runs, is.null, logical(1)))))
  models <- names(weights)
  if (!length(models) || !all(models %in% dimnames(runs[[1L]]$values)[[2L]])) return(NULL)
  common <- Reduce(`&`, lapply(runs, `[[`, "supported"))
  scopes <- list(marginal_training_quantiles = common)
  if (all(vapply(runs, function(r) !is.null(r$observed_range), logical(1))))
    scopes$observed_training_range <- Reduce(`&`, lapply(runs, `[[`, "observed_range"))
  contrasts <- dimnames(runs[[1L]]$values)[[3L]]
  rows <- list()
  for (criterion in names(scopes)) for (contrast in contrasts) {
    driver <- if (contrast == "total") NA_character_ else sub("^(only|retained):", "", contrast)
    selected <- if (is.na(driver)) rep(TRUE, length(runs)) else
      vapply(runs, function(r) isTRUE(r$selected[match(driver, validation$drivers)]), logical(1))
    matrices <- lapply(runs, function(run) {
      p <- matrix(run$values[, models, contrast, drop = FALSE], ncol = length(models))
      cbind(p, ensemble = as.numeric(p %*% weights))
    })
    for (j in seq_len(length(models)+1L)) {
      z <- do.call(cbind, lapply(matrices, function(p) p[, j]))
      keep <- scopes[[criterion]] & .cast_finite_rows(z)
      y <- z[keep, , drop = FALSE]
      means <- if (nrow(y)) colMeans(y) else rep(NA_real_, expected)
      signs <- .cast_response_sign(y)
      unanimous <- if (nrow(y)) apply(signs, 1L, function(s) length(unique(s)) == 1L) else logical()
      selected_agree <- if (!is.na(driver) && sum(selected) >= 2L && nrow(y)) {
        selected_signs <- signs[, selected, drop = FALSE]
        sum(rowSums(selected_signs > 0) == sum(selected) |
          rowSums(selected_signs < 0) == sum(selected))
      } else NA_integer_
      rows[[length(rows)+1L]] <- data.frame(scenario = validation$scenario,
        status = if (any(keep)) "descriptive_sensitivity" else "no_common_supported_locations",
        criterion = criterion,
        lower_probability = if (criterion == "observed_training_range") 0 else validation$support_probs[1L],
        upper_probability = if (criterion == "observed_training_range") 1 else validation$support_probs[2L],
        contrast = sub(":.*", "", contrast), driver = if (contrast == "total") NA_character_ else driver,
        model = c(models, "ensemble")[j], n_scheduled = expected, n_available = length(runs),
        n_locations = length(keep), n_common = sum(keep), common_fraction = mean(keep),
        n_common_presence = if (length(validation$observations) == length(keep))
          sum(keep & validation$observations == 1L) else NA_integer_,
        n_common_background = if (length(validation$observations) == length(keep))
          sum(keep & validation$observations == 0L) else NA_integer_,
        n_driver_selected = if (!is.na(driver)) sum(selected) else NA_integer_,
        n_direction_agree_when_selected = selected_agree,
        direction_agree_when_selected = if (!is.na(selected_agree)) selected_agree/nrow(y) else NA_real_,
        mean_response_min = if (nrow(y)) min(means) else NA_real_,
        mean_response_max = if (nrow(y)) max(means) else NA_real_,
        mean_response_sd = if (nrow(y)) stats::sd(means) else NA_real_,
        direction_unanimous_fraction = if (nrow(y)) mean(unanimous) else NA_real_,
        all_zero_fraction = if (nrow(y)) mean(rowSums(signs != 0) == 0) else NA_real_,
        interpretation = "Common locations across all spatial fits. Selected-fit agreement counts a shared positive or negative sign in at least two fits retaining the driver. Range/quantile boxes describe marginal coverage. Small or background-only intersections limit inference.")
    }
  }
  do.call(rbind, rows)
}
