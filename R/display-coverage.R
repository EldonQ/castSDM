# Report representation of the original footprint without changing native data.
.cast_display_coverage <- function(native, display, resolution) {
  native <- native[[1L]]; display <- display[[1L]]
  rows <- max(1L, as.integer(100000L/terra::ncol(native)))
  valid <- missing <- 0L
  for (first in seq.int(1L, terra::nrow(native), by = rows)) {
    count <- min(rows, terra::nrow(native)-first+1L)
    x <- .cast_raster_block(native, first, count)[[1L]]
    at <- which(is.finite(x))
    if (!length(at)) next
    cells <- (first-1L)*terra::ncol(native)+at
    xy <- terra::project(terra::xyFromCell(native, cells), terra::crs(native), terra::crs(display))
    y <- terra::extract(display, xy)[[1L]]
    valid <- valid+length(at); missing <- missing+sum(!is.finite(y))
  }
  data.frame(reference_layer = names(native), native_grid_cells = terra::ncell(native),
    native_valid_cells = valid, display_valid_cells = terra::global(is.finite(display), "sum", na.rm = TRUE)[[1L]],
    native_centres_missing_in_display = missing,
    native_centres_represented_fraction = if (valid) 1-missing/valid else NA_real_,
    display_resolution_degrees = resolution,
    resampling = "area-weighted average; categorical mode; binary support minimum; novelty maximum",
    interpretation = "Display cells summarise overlapping native predictions; they do not extend the analysis domain.")
}
