#' Build Data-Derived Environmental Domains
#'
#' Creates domain labels used by the CAST stable selector. Domains are learned
#' from coordinates or environmental covariates only; the response is used
#' solely to reject domains that cannot support binary SDM fitting.
#'
#' @param data A data frame containing `lon`, `lat`, `presence`, and predictors.
#' @param method Domain construction method: `"spatial"` (default) or
#'   `"environment"`.
#' @param k Requested number of domains.
#' @param response Response column used for minimum class-count checks.
#' @param env_vars Environmental columns for `method = "environment"`.
#' @param min_n Minimum observations per domain.
#' @param min_class Minimum presences and backgrounds per domain.
#' @param seed Random seed.
#'
#' @return A factor with one domain label per row.
#' @export
cast_domains <- function(data,
                         method = c("spatial", "environment"),
                         k = 4L,
                         response = "presence",
                         env_vars = NULL,
                         min_n = 40L,
                         min_class = 5L,
                         seed = NULL) {
  method <- match.arg(method)
  if (!is.data.frame(data) || !response %in% names(data)) {
    cli::cli_abort("{.arg data} must contain response column {.val {response}}.")
  }
  y <- as.integer(data[[response]])
  if (!all(y %in% c(0L, 1L))) {
    cli::cli_abort("{.arg response} must contain only 0/1 values.")
  }

  if (method == "spatial") {
    if (!all(c("lon", "lat") %in% names(data))) {
      cli::cli_abort("Spatial domains require {.field lon} and {.field lat}.")
    }
    z <- scale(cbind(lon = data$lon, lat = data$lat))
  } else {
    env_vars <- env_vars %||% get_env_vars(data, response)
    if (length(env_vars) < 2L) {
      cli::cli_abort("Environmental domains require at least two predictors.")
    }
    x <- .cast_numeric_matrix(data, env_vars)
    keep <- apply(x, 2L, stats::sd) > 1e-10
    x <- x[, keep, drop = FALSE]
    if (ncol(x) < 2L) cli::cli_abort("Too few non-constant predictors for domains.")
    pc <- stats::prcomp(x, center = TRUE, scale. = TRUE, rank. = min(5L, ncol(x)))
    z <- pc$x[, seq_len(min(3L, ncol(pc$x))), drop = FALSE]
  }

  z[!is.finite(z)] <- 0
  k <- max(2L, min(as.integer(k), nrow(data) %/% max(1L, as.integer(min_n))))
  if (k < 2L) cli::cli_abort("Insufficient rows to construct two valid domains.")
  base_seed <- seed %||% 1L

  while (k >= 2L) {
    for (attempt in seq_len(30L)) {
      set.seed(base_seed + 100L * k + attempt)
      km <- stats::kmeans(z, centers = k, nstart = 1L, iter.max = 100L)
      d <- factor(km$cluster)
      ok <- vapply(levels(d), function(level) {
        idx <- d == level
        sum(idx) >= min_n && sum(y[idx] == 1L) >= min_class &&
          sum(y[idx] == 0L) >= min_class
      }, logical(1))
      if (all(ok)) return(d)
    }
    k <- k - 1L
  }
  cli::cli_abort(
    "Could not construct at least two domains with {min_class} observations per class."
  )
}

#' @keywords internal
#' @noRd
.cast_numeric_matrix <- function(data, vars) {
  x <- as.data.frame(data[, vars, drop = FALSE])
  for (nm in names(x)) {
    x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
    med <- stats::median(x[[nm]], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    x[[nm]][!is.finite(x[[nm]]) | is.na(x[[nm]])] <- med
  }
  as.matrix(x)
}
