#' Use Externally Fitted SDMs for Scenario Comparisons
#'
#' Wraps mature modelling tools without fitting or changing their models.
#' The callback owns any preprocessing and returns relative suitability on
#' one declared, comparable scale. This adapter is intended for prediction,
#' response diagnostics and [cast_attribution()], not internal model tuning.
#'
#' @param models Named list of externally fitted model objects. Names must be
#'   unique and non-empty. Use `list(gam = fitted_gam)` for one model.
#' @param reference Non-empty numeric data frame of training predictors, in
#'   the same units and column order expected by the callback. Each predictor
#'   must contain finite training values. Used for range diagnostics only;
#'   castSDM cannot infer how the external tool fitted or transformed data.
#' @param predict_function Function with arguments `(model, newdata)` returning
#'   one numeric suitability in `[0, 1]` per row, or a named list of such
#'   functions with exactly the same names as `models`. Missing predictions
#'   are allowed; they remain missing in paired attribution.
#' @return A [cast_fit] wrapper. External objects and callbacks are retained
#'   together; load their required packages in a new R session. No response
#'   labels or validation scores are inferred. Supply externally validated
#'   named `weights` or use equal weights in [cast_attribution()].
#' @examples
#' set.seed(1)
#' reference <- data.frame(temperature = rnorm(100), rain = runif(100))
#' response <- rbinom(100, 1, plogis(reference$temperature - reference$rain))
#' model <- stats::glm(response ~ ., reference, family = stats::binomial())
#' fit <- cast_as_fit(list(glm = model), reference,
#'   function(model, newdata) stats::predict(model, newdata, type = "response"))
#' future <- reference
#' future$temperature <- future$temperature + 0.2
#' attribution <- cast_attribution(fit, reference, future, verbose = FALSE)
#' @export
cast_as_fit <- function(models, reference, predict_function) {
  if (!is.list(models) || !length(models) || is.null(names(models)) ||
      anyNA(names(models)) || any(!nzchar(names(models))) ||
      anyDuplicated(names(models))) {
    cli::cli_abort("{.arg models} must be a non-empty list with unique non-empty names.")
  }
  if (!is.data.frame(reference) || !nrow(reference) || !ncol(reference) ||
      anyNA(names(reference)) || any(!nzchar(names(reference))) ||
      anyDuplicated(names(reference))) {
    cli::cli_abort("{.arg reference} must be a non-empty predictor data frame with unique names.")
  }
  .cast_check_numeric_predictors(reference, arg = "reference")
  if (any(!vapply(reference, function(x) any(is.finite(x)), logical(1)))) {
    cli::cli_abort("Every {.arg reference} predictor needs finite training values.")
  }
  callbacks <- if (is.function(predict_function)) {
    stats::setNames(rep(list(predict_function), length(models)), names(models))
  } else predict_function
  if (!is.list(callbacks) || length(callbacks) != length(models) ||
      is.null(names(callbacks)) || anyDuplicated(names(callbacks)) ||
      !setequal(names(callbacks), names(models)) ||
      !all(vapply(callbacks, is.function, logical(1)))) {
    cli::cli_abort("{.arg predict_function} must be a function or matching named list of functions.")
  }
  wrapped <- lapply(names(models), function(n) list(type = "external",
    name = n, model = models[[n]], predict_function = callbacks[[n]]))
  names(wrapped) <- names(models)
  impute <- vapply(reference, function(x) stats::median(x[is.finite(x)]), numeric(1))
  means <- vapply(reference, function(x) mean(x[is.finite(x)]), numeric(1))
  sds <- vapply(reference, function(x) stats::sd(x[is.finite(x)]), numeric(1))
  sds[!is.finite(sds) | sds < 1e-10] <- 1
  new_cast_fit(wrapped, names(reference), names(reference),
    scaling = list(reference = reference, impute = impute, means = means, sds = sds))
}
