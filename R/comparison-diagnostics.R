# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# Reductions performed on the prediction block, never by re-predicting.
.cast_check_thresholds <- function(thresholds) {
  if (is.null(thresholds)) return(numeric())
  if (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds)) || anyDuplicated(thresholds))
    cli::cli_abort("thresholds must be unique finite numbers.")
  thresholds
}

.cast_threshold_accumulator <- function(thresholds, n_faces = 4L) {
  e <- new.env(parent = emptyenv())
  e$thresholds <- thresholds
  e$counts <- array(0, c(length(thresholds), n_faces, 3L))
  e$n <- numeric(3L)
  e
}

.cast_threshold_accumulate <- function(e, faces, support, common_support = NULL) {
  valid <- .cast_finite_rows(faces)
  domains <- list(valid, valid & !is.na(support) & support,
    if (is.null(common_support)) rep(FALSE, nrow(faces)) else
      valid & !is.na(common_support) & common_support)
  for (d in seq_along(domains)) {
    keep <- domains[[d]]
    e$n[d] <- e$n[d] + sum(keep)
    if (!any(keep)) next
    x <- faces[keep, , drop = FALSE]
    for (i in seq_along(e$thresholds))
      e$counts[i, , d] <- e$counts[i, , d] + colSums(x >= e$thresholds[i])
  }
  invisible(e)
}

.cast_threshold_summary <- function(e, has_support, has_common) {
  rows <- list()
  domains <- c(1L, if (has_support) 2L, if (has_common) 3L)
  for (i in seq_along(e$thresholds)) for (d in domains) {
    b <- .cast_factorial_binary(e$counts[i, , d], e$n[d])
    b$threshold <- e$thresholds[i]
    b$domain <- c("full", "supported", "common_supported")[d]
    rows[[length(rows)+1L]] <- b
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

.cast_sign_values <- function(per_model, tolerance = 1e-12) {
  n <- nrow(per_model[[1L]]); m <- length(per_model)
  k_count <- ncol(per_model[[1L]])
  agreement <- conflict <- positive <- negative <- zero <- matrix(NA_real_, n, k_count)
  if (m > 1L) for (k in seq_len(k_count)) {
    x <- matrix(vapply(per_model, function(p) p[, k], numeric(n)), n, m)
    pos <- rowSums(x > tolerance); neg <- rowSums(x < -tolerance)
    neutral <- m-pos-neg
    agreement[, k] <- pmax(pos, neg, neutral)/m
    conflict[, k] <- (pos > 0 & neg > 0)*1
    positive[, k] <- (pos == m)*1; negative[, k] <- (neg == m)*1
    zero[, k] <- (neutral == m)*1
  }
  list(agreement = agreement, conflict = conflict, positive = positive,
       negative = negative, zero = zero)
}

.cast_sign_accumulator <- function(n_models, components = .cast_factorial_names()) {
  e <- new.env(parent = emptyenv())
  e$n_models <- n_models
  e$components <- components
  e$n <- numeric(3L)
  e$sums <- array(0, c(length(components), 5L, 3L))
  e
}

.cast_sign_accumulate <- function(e, sign, support, common_support) {
  domains <- list(rep(TRUE, nrow(sign$agreement)), !is.na(support) & support,
    if (is.null(common_support)) rep(FALSE, nrow(sign$agreement)) else
      !is.na(common_support) & common_support)
  for (d in seq_along(domains)) {
    keep <- domains[[d]]
    e$n[d] <- e$n[d] + sum(keep)
    if (e$n_models < 2L || !any(keep)) next
    for (j in seq_along(sign)) e$sums[, j, d] <- e$sums[, j, d] +
      colSums(sign[[j]][keep, , drop = FALSE])
  }
  invisible(e)
}

.cast_sign_summary <- function(e, has_support, has_common) {
  do.call(rbind, lapply(c(1L, if (has_support) 2L, if (has_common) 3L), function(d) {
    s <- if (e$n_models > 1L && e$n[d] > 0) matrix(e$sums[, , d], nrow = length(e$components))/e$n[d] else matrix(NA_real_, length(e$components), 5L)
    data.frame(component = e$components, domain = c("full", "supported", "common_supported")[d],
      n_valid = e$n[d], n_models = e$n_models, mean_sign_agreement = s[, 1L],
      mixed_sign_share = s[, 2L], unanimous_positive_share = s[, 3L],
      unanimous_negative_share = s[, 4L], unanimous_zero_share = s[, 5L],
      sign_tolerance = 1e-12)
  }))
}

# A persistent reader avoids opening/decompressing the same GeoTIFF tiles for
# every short row strip. Align blocks to GeoTIFF tile height (writer: 256 rows).
.cast_threshold_scan <- function(result, thresholds, has_support, has_common,
                                  memory_mb = 256) {
  raster <- !is.null(result$raster_file)
  face_names <- if (identical(result$design, "paired")) c("current", "joint") else c("current", "joint", "climate", "lulc")
  n_faces <- length(face_names)
  cols <- c(face_names, "support")
  if (raster) {
    input <- terra::rast(result$raster_file)
    has_common <- has_common && "common_support" %in% names(input)
    if (has_common) cols <- c(cols, "common_support")
    input <- input[[cols]]
    terra::readStart(input)
    on.exit(terra::readStop(input), add = TRUE)
    target <- max(1, floor(memory_mb*1e6/(8*terra::ncol(input)*length(cols)*4)))
    rows <- if (target >= 256) max(256L, floor(target/256)*256L) else target
    rows <- min(rows, terra::nrow(input))
    starts <- seq.int(1L, terra::nrow(input), by = rows)
  } else {
    has_common <- has_common && "common_support" %in% names(result$data)
    if (has_common) cols <- c(cols, "common_support")
    input <- as.matrix(result$data[, cols, drop = FALSE])
    starts <- seq.int(1L, nrow(input), by = 100000L)
  }
  acc <- .cast_threshold_accumulator(thresholds, n_faces)
  for (first in starts) {
    x <- if (raster) terra::readValues(input, row = first,
      nrows = min(rows, terra::nrow(input)-first+1L), mat = TRUE) else
        input[seq.int(first, min(nrow(input), first+99999L)), , drop = FALSE]
    .cast_threshold_accumulate(acc, x[, seq_len(n_faces), drop = FALSE],
      is.finite(x[, n_faces+1L]) & x[, n_faces+1L] != 0,
      if (has_common) is.finite(x[, n_faces+2L]) & x[, n_faces+2L] != 0)
  }
  .cast_threshold_summary(acc, has_support, has_common)
}
