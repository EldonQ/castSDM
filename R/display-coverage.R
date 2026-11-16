# Report representation of the original footprint without changing native data.
.cast_display_coverage <- function(native, display, resolution,
    references = stats::setNames("full_domain", names(native)[1L])) {
  native <- native[[names(references)]]; display <- display[[names(references)]]
  rows <- max(1L, min(256L, as.integer(1000000L/terra::ncol(native))))
  terra::readStart(native)
  on.exit(terra::readStop(native), add = TRUE)
  valid <- missing <- integer(terra::nlyr(native))
  for (first in seq.int(1L, terra::nrow(native), by = rows)) {
    count <- min(rows, terra::nrow(native)-first+1L)
    x <- terra::readValues(native, row = first, nrows = count, mat = TRUE)
    for (i in seq_len(ncol(x))) {
      at <- which(is.finite(x[, i]))
      if (!length(at)) next
      cells <- (first-1L)*terra::ncol(native)+at
      xy <- terra::project(terra::xyFromCell(native, cells), terra::crs(native), terra::crs(display))
      y <- terra::extract(display[[i]], xy)[[1L]]
      valid[i] <- valid[i]+length(at); missing[i] <- missing[i]+sum(!is.finite(y))
    }
  }
  data.frame(reference_layer = names(native), scope = unname(references),
    mask_order = ifelse(references == "full_domain", "not applicable", "native mask before aggregation"),
    native_grid_cells = terra::ncell(native),
    native_valid_cells = valid, display_valid_cells = terra::global(is.finite(display), "sum", na.rm = TRUE)[[1L]],
    native_centres_missing_in_display = missing,
    native_centres_represented_fraction = ifelse(valid > 0L, 1-missing/valid, NA_real_),
    display_resolution_degrees = resolution,
    resampling = "area-weighted average; categorical mode; support 0/1/2 (outside/inside/mixed); novelty maximum",
    interpretation = "Display cells summarise overlapping native predictions; they do not extend the analysis domain.")
}
