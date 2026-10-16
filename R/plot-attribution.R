#' Plot Environmental Attribution of Habitat Change
#'
#' @param x A [cast_attribution] result.
#' @param type `"summary"` or `"map"`.
#' @param contrast `"only"` (driver-only future change) or `"retained"`
#'   (full future minus future with that driver held current).
#' @param driver Fitted driver name for a map. `NULL` maps total change.
#' @param basemap `"none"`, `"china"` or `"world"`. Geographic basemaps
#'   require longitude/latitude coordinates. Raster display is reprojected
#'   to WGS84 when needed; source rasters remain unchanged.
#' @param limits Optional symmetric suitability-change colour limits.
#' @param supported_only Whether to hide cells outside the training range.
#' @param max_cells Maximum display cells (regular aggregation for rasters;
#'   table maps retain all observations).
#' @param ... Ignored.
#' @return A `ggplot` object using the package's cartographic assets.
#' @export
plot.cast_attribution <- function(x, type = c("summary", "map"),
                                  contrast = c("only", "retained"), driver = NULL,
                                  basemap = "none", limits = NULL,
                                  supported_only = FALSE, max_cells = 100000L, ...) {
  check_suggested("ggplot2", "for attribution plots")
  .cast_check_scalar(max_cells, "max_cells", lower = 1, integer = TRUE)
  basemap <- match.arg(basemap, c("none", "china", "world"))
  type <- match.arg(type); contrast <- match.arg(contrast)
  if (type == "summary") {
    d <- x$summary[x$summary$contrast == contrast, , drop = FALSE]
    d$label <- .cast_var_label(d$driver)
    d$label <- factor(d$label, levels = d$label[order(d$mean_abs)])
    return(ggplot2::ggplot(d, ggplot2::aes(.data$mean_signed, .data$label)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey65") +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$mean_signed,
        yend = .data$label), colour = "grey65", linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$mean_signed), size = 2.5) +
      ggplot2::scale_colour_gradient2(low = "#B35806", mid = "#F7F7F7", high = "#2166AC",
        midpoint = 0, guide = "none") +
      ggplot2::labs(x = "Mean suitability change", y = NULL,
        title = if (contrast == "only") "Driver-only environmental change" else "Contribution of future driver changes",
        caption = "Frozen model weights; scenario contrasts are non-additive. Range coverage is reported in the summary table.") +
      theme_cast(base_size = 8))
  }
  if (!is.null(driver) && (!is.character(driver) || length(driver) != 1L ||
      !driver %in% x$drivers)) cli::cli_abort("{.arg driver} must name one mapped driver.")
  layer <- if (is.null(driver)) "total_change" else
    paste0(if (contrast == "only") "effect_" else "retained_", driver)
  support <- if (is.null(driver)) "support_total" else paste0("support_", layer)
  raster <- !is.null(x$raster_file)
  if (raster) {
    check_suggested("terra", "for attribution map display")
    r <- terra::rast(x$raster_file)[[c(layer, support)]]
    if (supported_only) r[[1]] <- terra::ifel(r[[2]] == 1, r[[1]], NA)
    factor <- max(1L, ceiling(sqrt(terra::ncell(r) / max_cells)))
    if (factor > 1L) r <- terra::aggregate(r[[1]], factor, fun = "mean", na.rm = TRUE)
    else r <- r[[1]]
    if (basemap != "none" && !terra::is.lonlat(r)) r <- terra::project(r, "EPSG:4326")
    d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(d) <- c("lon", "lat", "value")
    dx <- terra::res(r)[1]; dy <- terra::res(r)[2]
  } else {
    if (!all(c("lon", "lat") %in% names(x$data))) cli::cli_abort("Map tables need lon/lat coordinates.")
    d <- data.frame(lon = x$data$lon, lat = x$data$lat, value = x$data[[layer]])
    if (supported_only) d$value[x$data[[support]] != 1 | is.na(x$data[[support]])] <- NA
  }
  if (!nrow(d) || !any(is.finite(d$value))) cli::cli_abort("No finite values for this attribution map.")
  if (basemap != "none" && (any(abs(d$lon) > 180) || any(abs(d$lat) > 90))) {
    cli::cli_abort("Geographic basemaps require longitude/latitude coordinates.")
  }
  if (is.null(limits)) {
    hi <- max(abs(d$value), na.rm = TRUE)
    if (hi == 0) hi <- 0.01
    limits <- c(-hi, hi)
  }
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) ||
      limits[1] >= 0 || limits[2] <= 0 || abs(sum(limits)) > 1e-10) {
    cli::cli_abort("{.arg limits} must be two symmetric finite limits around zero.")
  }
  scale <- function(guide = "colourbar") ggplot2::scale_fill_gradient2(
    low = "#B35806", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
    limits = limits, oob = function(v, range) pmax(range[1], pmin(range[2], v)),
    na.value = NA, name = "Suitability\nchange", guide = guide)
  title <- if (is.null(driver)) "Overall habitat suitability change" else
    paste(if (contrast == "only") "Driver-only effect:" else "Future contribution:", .cast_var_label(driver))
  subtitle <- if (supported_only) "Within training ranges" else "Raw scenario comparison; inspect range coverage"
  if (basemap == "china" && raster) {
    return(.cast_map_surface(d, "value", function(guide) scale(if (guide == "none") "none" else "colourbar"),
      title = title, subtitle = subtitle, bg = "#F7F7F7", cell_size = c(dx, dy)))
  }
  p <- ggplot2::ggplot()
  if (basemap != "none") {
    check_suggested("sf", "for geographic mapping")
    bm <- load_basemap(basemap)
    if (!is.null(bm)) p <- p + ggplot2::geom_sf(data = bm, fill = "#F7F7F7",
      colour = "grey70", linewidth = 0.2)
  }
  if (raster) {
    p <- p + ggplot2::geom_tile(data = d, ggplot2::aes(.data$lon, .data$lat,
      fill = .data$value), width = dx, height = dy, linewidth = 0)
  } else {
    p <- p + ggplot2::geom_point(data = d, ggplot2::aes(.data$lon, .data$lat,
      fill = .data$value), shape = 22, size = 1, stroke = 0)
  }
  title <- if (is.null(driver)) "Overall habitat suitability change" else
    paste(if (contrast == "only") "Driver-only effect:" else "Future contribution:", .cast_var_label(driver))
  p <- p + scale() + ggplot2::labs(title = title,
    subtitle = if (supported_only) "Within training ranges" else "Raw scenario comparison; inspect range coverage",
    x = NULL, y = NULL) + theme_cast(base_size = 8) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")
  if (basemap == "none") return(p + ggplot2::coord_equal())
  p <- .add_china_dashline(p, basemap)
  p <- .add_china_outline(p, basemap) + .coord_for_map(basemap)
  .add_china_south_sea_inset(p, basemap, data = d, value_var = "value", raster = raster,
    scale = scale("none"), bg_fill = "#F7F7F7",
    cell_size = if (raster) c(dx, dy) else NULL)
}
