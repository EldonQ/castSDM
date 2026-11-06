# Shared four-corner arithmetic and online reductions.
.cast_factorial_names <- function() c("climate_main", "lulc_main", "joint_change", "interaction",
  "climate_conditional", "lulc_conditional", "climate_average", "lulc_average")

.cast_factorial_values <- function(x) {
  a <- x[, 3L] - x[, 1L]; b <- x[, 4L] - x[, 1L]; j <- x[, 2L] - x[, 1L]
  i <- j - a - b
  z <- cbind(a, b, j, i, x[, 2L] - x[, 4L], x[, 2L] - x[, 3L], a + i/2, b + i/2)
  colnames(z) <- .cast_factorial_names()
  z
}

.cast_factorial_accumulator <- function() {
  e <- new.env(parent = emptyenv())
  e$moments <- matrix(0, 8L, 10L, dimnames = list(.cast_factorial_names(),
    c("n", "ns", "sum", "abs", "sums", "abss", "sd", "nsd", "surface", "nsurface")))
  e$binary <- e$binary_supported <- numeric(4L)
  e$n_binary <- e$n_binary_supported <- 0
  e
}

.cast_factorial_accumulate <- function(e, faces, effects, supported, binary = NULL,
                                        model_sd = NULL, surface_sd = NULL) {
  for (k in seq_len(ncol(effects))) {
    d <- effects[, k]; ok <- is.finite(d); s <- ok & !is.na(supported) & supported
    sd <- if (is.null(model_sd)) rep(NA_real_, length(d)) else model_sd[, k]
    su <- if (is.null(surface_sd)) rep(NA_real_, length(d)) else surface_sd[, k]
    ks <- ok & is.finite(sd); ku <- ok & is.finite(su)
    e$moments[k, ] <- e$moments[k, ] + c(sum(ok), sum(s), sum(d[ok]), sum(abs(d[ok])),
      sum(d[s]), sum(abs(d[s])), sum(sd[ks]), sum(ks), sum(su[ku]), sum(ku))
  }
  if (!is.null(binary)) {
    if (any(is.finite(binary) & binary != 0 & binary != 1) || any(is.infinite(binary)))
      cli::cli_abort("Binary inputs must contain only 0, 1 or NA.")
    valid <- rowSums(!is.finite(faces)) == 0L & rowSums(!is.finite(binary)) == 0L
    s <- valid & !is.na(supported) & supported
    e$binary <- e$binary + colSums(binary[valid, , drop = FALSE])
    e$binary_supported <- e$binary_supported + colSums(binary[s, , drop = FALSE])
    e$n_binary <- e$n_binary + sum(valid); e$n_binary_supported <- e$n_binary_supported + sum(s)
  }
  invisible(e)
}

.cast_factorial_summary <- function(e, has_support) {
  m <- e$moments; div <- function(x, n) ifelse(n > 0, x/n, NA_real_)
  data.frame(component = rownames(m), scale = "continuous", n_valid = m[, "n"],
    n_supported = if (has_support) m[, "ns"] else NA_real_,
    support_share = if (has_support) div(m[, "ns"], m[, "n"]) else NA_real_,
    mean = div(m[, "sum"], m[, "n"]), mean_abs = div(m[, "abs"], m[, "n"]),
    mean_supported = if (has_support) div(m[, "sums"], m[, "ns"]) else NA_real_,
    mean_abs_supported = if (has_support) div(m[, "abss"], m[, "ns"]) else NA_real_,
    mean_model_sd = div(m[, "sd"], m[, "nsd"]),
    mean_surface_sd = div(m[, "surface"], m[, "nsurface"]), row.names = NULL)
}

.cast_factorial_binary <- function(counts, n) {
  estimate <- as.numeric(.cast_factorial_values(matrix(counts, nrow = 1L)))
  data.frame(component = .cast_factorial_names(), n_valid = n, n_current_range = counts[1L],
    n_future_range = c(counts[3L], counts[4L], counts[2L], NA_real_, counts[2L], counts[2L], NA_real_, NA_real_),
    estimate = estimate, pct_of_current = if (counts[1L] > 0) 100*estimate/counts[1L] else NA_real_, row.names = NULL)
}

.cast_factorial_contract <- function() {
  data.frame(key = c("claim_level", "design", .cast_factorial_names(), "support", "model_sd", "surface_sd", "not_identified"),
    value = c("model_conditional_scenario_contrasts", "2x2: current, joint, first-group-only, second-group-only",
      "climate - current (baseline simple effect; legacy column name)",
      "lulc - current (baseline simple effect; legacy column name)", "joint - current",
      "joint - climate - lulc + current", "joint - lulc", "joint - climate",
      "[(climate - current) + (joint - lulc)] / 2", "[(lulc - current) + (joint - climate)] / 2",
      "Common finite domain; supported subset separately, not a confidence interval",
      "Across models after paired contrasts; not sampling uncertainty; NA for one model",
      "Optional independent-surface sqrt-sum diagnostic, not contrast SD or a confidence interval",
      "Not identified: real-world causal effects, occupancy, absolute occurrence probability or realised distribution"),
    stringsAsFactors = FALSE)
}

.cast_single_raster <- function(x, arg, reference = NULL) {
  if (is.character(x) && length(x) == 1L && !is.na(x) && file.exists(x)) x <- terra::rast(x)
  if (!inherits(x, "SpatRaster") || terra::nlyr(x) != 1L)
    cli::cli_abort("{.arg {arg}} must be a single-layer SpatRaster or readable raster path.")
  if (!is.null(reference) && !isTRUE(terra::compareGeom(reference, x,
      lyrs = FALSE, stopOnError = FALSE, messages = FALSE)))
    cli::cli_abort("{.arg {arg}} is not aligned with current.")
  x
}

.cast_factorial_surface_sd <- function(x) {
  if (is.null(x)) return(NULL)
  if (any(x < 0, na.rm = TRUE) || any(is.infinite(x)))
    cli::cli_abort("Surface SD inputs must be non-negative finite values or NA.")
  a <- x[, 1L]^2; b <- x[, 2L]^2; c <- x[, 3L]^2; d <- x[, 4L]^2
  z <- cbind(sqrt(a+c), sqrt(a+d), sqrt(a+b), sqrt(a+b+c+d), sqrt(b+d), sqrt(b+c),
    sqrt(a+b+c+d)/2, sqrt(a+b+c+d)/2)
  colnames(z) <- .cast_factorial_names()
  z
}

#' Compare Four Aligned Factorial Suitability Surfaces
#'
#' Computes baseline simple effects, conditional and average marginal effects,
#' and non-additivity on the declared suitability scale. Legacy `climate_main`
#' and `lulc_main` columns are baseline simple effects, not averaged effects.
#' These are fitted-model scenario contrasts, not identified causal effects.
#' Use [cast_compare()] to construct all four corners from an existing model.
#' @param current,joint,climate,lulc Equally sized numeric vectors or one-column
#'   frames, or four aligned single-layer `terra::SpatRaster`s/paths.
#' @param current_binary,joint_binary,climate_binary,lulc_binary Optional 0/1/NA
#'   inputs of matching size/geometry. Supply all four or none.
#' @param threshold Optional finite threshold. Binary values use `>=`.
#' @param mask Matching input; only finite nonzero values participate. All
#'   contrasts and counts share the four-surface finite domain.
#' @param support Matching input; finite nonzero values mark training support.
#'   Unsupported values remain in raw contrasts and are reported separately.
#' @param current_sd,joint_sd,climate_sd,lulc_sd Optional non-negative surface
#'   SDs, all four or none. The independent-surface approximation is named
#'   `surface_sd_*`; it is not paired contrast SD.
#' @param scenario,ssp,period Optional metadata.
#' @param filename Optional GeoTIFF path. Raster results otherwise use a
#'   temporary file to keep computation bounded in memory.
#' @param overwrite,compression Raster publication options.
#' @param verbose Print progress.
#' @param block_size Maximum cells per processing block (approximately, by row).
#' @return A `cast_factorial` object with `summary`, optional `binary` and
#'   `binary_supported`, contrast `raster` or table `data`, and `contract`.
#'   `mean_surface_sd` is unpaired; `mean_model_sd` needs model-level inputs.
#'   Range sizes are cell counts, never automatically square kilometres.
#' @examples
#' fa <- cast_factorial_attribution(c(.2, .3), c(.8, .7),
#'   c(.4, .5), c(.5, .4), threshold = .5, verbose = FALSE)
#' fa$summary
#' @export
cast_factorial_attribution <- function(current, joint, climate, lulc,
    current_binary = NULL, joint_binary = NULL, climate_binary = NULL, lulc_binary = NULL,
    threshold = NULL, mask = NULL, support = NULL,
    current_sd = NULL, joint_sd = NULL, climate_sd = NULL, lulc_sd = NULL,
    scenario = NULL, ssp = NULL, period = NULL, filename = NULL, overwrite = FALSE,
    compression = "LZW", verbose = TRUE, block_size = 100000L) {
  block_size <- .cast_check_scalar(block_size, "block_size", lower = 1, integer = TRUE)
  if (!is.null(threshold)) threshold <- .cast_check_scalar(threshold, "threshold")
  inputs <- list(current = current, joint = joint, climate = climate, lulc = lulc)
  binary <- list(current_binary, joint_binary, climate_binary, lulc_binary)
  sds <- list(current_sd, joint_sd, climate_sd, lulc_sd)
  complete <- function(x, label) {
    n <- sum(!vapply(x, is.null, logical(1)))
    if (!n %in% c(0L, 4L)) cli::cli_abort("Supply all four {label} inputs, or none.")
    n == 4L
  }
  has_binary <- complete(binary, "binary"); has_sd <- complete(sds, "SD")
  has_support <- !is.null(support)
  raster <- inherits(current, "SpatRaster") || is.character(current)
  if (raster) {
    check_suggested("terra", "for raster factorial comparisons")
    ref <- .cast_single_raster(current, "current")
    inputs <- Map(function(x, n) .cast_single_raster(x, n, ref), inputs, names(inputs))
    if (has_binary) binary <- lapply(binary, .cast_single_raster, arg = "binary", reference = ref)
    if (has_sd) sds <- lapply(sds, .cast_single_raster, arg = "SD", reference = ref)
    if (!is.null(mask)) mask <- .cast_single_raster(mask, "mask", ref)
    if (has_support) support <- .cast_single_raster(support, "support", ref)
    stack <- do.call(c, unname(inputs))
    if (has_binary) stack <- c(stack, do.call(c, unname(binary)))
    if (has_sd) stack <- c(stack, do.call(c, unname(sds)))
    if (!is.null(mask)) stack <- c(stack, mask)
    if (has_support) stack <- c(stack, support)
    rows <- max(1L, as.integer(block_size/terra::ncol(ref)))
    starts <- seq.int(1L, terra::nrow(ref), by = rows)
    filename <- filename %||% tempfile("cast_factorial_", fileext = ".tif")
    layers <- c(names(inputs), .cast_factorial_names(), if (has_sd) paste0("surface_sd_", .cast_factorial_names()))
    writer <- .cast_raster_writer(ref, filename, layers, datatype = "FLT8S", compression = compression, overwrite = overwrite)
    on.exit(.cast_writer_abort(writer), add = TRUE)
  } else {
    one <- function(x) {
      if (is.data.frame(x) && ncol(x) == 1L) x <- x[[1L]]
      if (!is.numeric(x) && !is.logical(x)) cli::cli_abort("Table inputs must be numeric vectors or one-column frames.")
      as.numeric(x)
    }
    inputs <- lapply(inputs, one); n <- length(inputs[[1L]])
    if (!n || any(lengths(inputs) != n)) cli::cli_abort("The four table inputs must have equal nonzero length.")
    extras <- c(if (has_binary) binary, if (has_sd) sds,
      if (!is.null(mask)) list(mask), if (has_support) list(support))
    extras <- lapply(extras, one)
    if (any(lengths(extras) != n)) cli::cli_abort("All optional table inputs must match current's length.")
    stack <- do.call(cbind, c(inputs, extras)); starts <- seq.int(1L, n, by = block_size)
    if (!is.null(filename)) cli::cli_abort("filename requires raster inputs.")
    data <- matrix(NA_real_, n, 12L + if (has_sd) 8L else 0L)
    colnames(data) <- c(names(inputs), .cast_factorial_names(), if (has_sd) paste0("surface_sd_", .cast_factorial_names()))
  }
  acc <- .cast_factorial_accumulator()
  for (i in seq_along(starts)) {
    first <- starts[i]
    if (raster) {
      count <- min(rows, terra::nrow(ref)-first+1L)
      values <- terra::values(stack, row = first, nrows = count, mat = TRUE)
    } else {
      index <- seq.int(first, min(n, first+block_size-1L)); values <- stack[index, , drop = FALSE]
    }
    faces <- values[, 1:4, drop = FALSE]; pos <- 4L
    bins <- if (has_binary) values[, pos+1:4, drop = FALSE] else NULL
    if (has_binary) pos <- pos+4L
    sd <- if (has_sd) values[, pos+1:4, drop = FALSE] else NULL
    if (has_sd) pos <- pos+4L
    domain <- rowSums(!is.finite(faces)) == 0L
    if (!is.null(mask)) {
      pos <- pos+1L; domain <- domain & is.finite(values[, pos]) & values[, pos] != 0
    }
    supported <- if (has_support) is.finite(values[, pos+1L]) & values[, pos+1L] != 0 else rep(NA, nrow(faces))
    faces[!domain, ] <- NA_real_; effects <- .cast_factorial_values(faces)
    if (!has_binary && !is.null(threshold)) bins <- (faces >= threshold)*1
    if (!is.null(bins)) bins[!domain, ] <- NA_real_
    surface <- .cast_factorial_surface_sd(sd)
    if (!is.null(surface)) surface[!domain, ] <- NA_real_
    .cast_factorial_accumulate(acc, faces, effects, supported, bins, surface_sd = surface)
    z <- cbind(faces, effects, surface)
    if (raster) .cast_writer_block(writer, z, first, count) else data[index, ] <- z
  }
  if (raster) .cast_writer_finish(writer)
  result <- structure(list(summary = .cast_factorial_summary(acc, has_support),
    binary = if (has_binary || !is.null(threshold)) .cast_factorial_binary(acc$binary, acc$n_binary) else NULL,
    binary_supported = if (has_support && (has_binary || !is.null(threshold)))
      .cast_factorial_binary(acc$binary_supported, acc$n_binary_supported) else NULL,
    raster = if (raster) terra::rast(filename) else NULL,
    raster_file = if (raster) filename else NA_character_, data = if (raster) NULL else as.data.frame(data),
    contract = .cast_factorial_contract(), scenario = scenario, ssp = ssp, period = period,
    threshold = threshold), class = "cast_factorial")
  if (verbose) cli::cli_inform("Computed eight model-conditional contrasts on a common finite domain.")
  result
}

#' @export
print.cast_factorial <- function(x, ...) {
  cli::cli_h1("castSDM factorial scenario contrasts")
  cli::cli_inform("Legacy *_main columns are baseline simple effects; *_average averages both levels.")
  print(x$summary, row.names = FALSE)
  if (!is.null(x$binary)) print(x$binary, row.names = FALSE)
  invisible(x)
}
#' @export
summary.cast_factorial <- function(object, ...) { print(object); invisible(object) }
#' @export
plot.cast_factorial <- function(x, layers = NULL, ...) {
  layers <- layers %||% c("climate_main", "lulc_main", "joint_change", "interaction")
  if (is.null(x$raster)) {
    s <- x$summary[x$summary$component %in% layers, ]
    if (!nrow(s)) cli::cli_abort("No requested contrast is present.")
    graphics::barplot(stats::setNames(s$mean, s$component), ylab = "Mean suitability contrast", ...)
  } else {
    layers <- intersect(layers, names(x$raster))
    if (!length(layers)) cli::cli_abort("No requested contrast layer is present.")
    terra::plot(x$raster[[layers]], ...)
  }
  invisible(x)
}
