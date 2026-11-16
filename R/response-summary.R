# Unweighted member responses on the same paired domain as the ensemble.
.cast_response_summary <- function(sums, drivers) {
  divide <- function(x, n) ifelse(n > 0, x/n, NA_real_)
  at <- seq_len(1L+2L*length(drivers))
  data.frame(driver = c(NA_character_, drivers, drivers, NA_character_),
    contrast = c("total", rep("only", length(drivers)), rep("retained", length(drivers)), "remainder"),
    n = sums[, "n"], n_supported = sums[, "ns"], support = divide(sums[, "ns"], sums[, "n"]),
    mean_signed = divide(sums[, "sum"], sums[, "n"]), mean_abs = divide(sums[, "abs"], sums[, "n"]),
    mean_signed_supported = divide(sums[, "sums"], sums[, "ns"]), mean_abs_supported = divide(sums[, "abss"], sums[, "ns"]),
    mean_model_sd = c(divide(sums[at, "sd"], sums[at, "nsd"]), NA_real_), row.names = NULL)
}

.cast_model_response_sums <- function(value, support, weights, driver, contrast) {
  valid <- .cast_finite_rows(value)
  supported <- valid & !is.na(support) & support
  assessed <- !length(support) || any(!is.na(support))
  do.call(rbind, lapply(seq_along(weights), function(j) {
    x <- value[, j]; sign <- .cast_response_sign(x[valid])
    data.frame(model = names(weights)[j], driver = driver, contrast = contrast,
      n = sum(valid), n_supported = if (assessed) sum(supported) else NA_real_, sum = sum(x[valid]),
      sum_abs = sum(abs(x[valid])), sum_supported = if (assessed) sum(x[supported]) else NA_real_,
      sum_abs_supported = if (assessed) sum(abs(x[supported])) else NA_real_, negative = sum(sign < 0),
      zero = sum(sign == 0), positive = sum(sign > 0), row.names = NULL)
  }))
}

.cast_model_response_add <- function(total, block) {
  if (is.null(total)) return(block)
  if (is.null(block)) return(total)
  if (!identical(total[, 1:3], block[, 1:3]))
    cli::cli_abort("Per-model response summary keys differ between raster batches.")
  total[, -(1:3)] <- total[, -(1:3)] + block[, -(1:3)]
  total
}

.cast_model_response_summary <- function(sums, weights) {
  if (is.null(sums)) return(NULL)
  divide <- function(x, n) ifelse(n > 0, x/n, NA_real_)
  data.frame(sums[, c("model", "driver", "contrast", "n", "n_supported")],
    ensemble_weight = unname(weights[sums$model]),
    support = divide(sums$n_supported, sums$n),
    mean_signed = divide(sums$sum, sums$n), mean_abs = divide(sums$sum_abs, sums$n),
    mean_signed_supported = divide(sums$sum_supported, sums$n_supported),
    mean_abs_supported = divide(sums$sum_abs_supported, sums$n_supported),
    negative_fraction = divide(sums$negative, sums$n), zero_fraction = divide(sums$zero, sums$n),
    positive_fraction = divide(sums$positive, sums$n), row.names = NULL)
}
