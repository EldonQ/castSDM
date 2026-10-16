.cast_display_raster <- function(path, layer = 1L, res_deg = 0.05,
                          method = "bilinear", mask = NULL,
                          keep_outside = FALSE, feather_deg = 0) {
  rd <- terra::project(terra::rast(path)[[layer]], "EPSG:4326",
                       res = res_deg, method = method)
  df <- terra::as.data.frame(rd, xy = TRUE, na.rm = TRUE)
  names(df) <- c("lon", "lat", "val")
  if (!is.null(mask)) {
    inside <- try({
      m <- sf::st_within(
        sf::st_as_sf(data.frame(lon = df$lon, lat = df$lat),
                     coords = c("lon", "lat"), crs = 4326), mask, sparse = FALSE)
      as.logical(m[, 1L])
    }, silent = TRUE)
    if (!inherits(inside, "try-error") && length(inside) == nrow(df)) {
      df$.inside <- inside
      if (!keep_outside) df <- df[inside, , drop = FALSE]
      if (!is.null(feather_deg) && is.numeric(feather_deg) && feather_deg > 0) {
        pts <- sf::st_as_sf(data.frame(lon = df$lon, lat = df$lat),
                            coords = c("lon", "lat"), crs = 4326)
        bnd <- sf::st_cast(mask, "MULTILINESTRING")
        df$.dist_deg <- as.numeric(sf::st_distance(pts, bnd)) / 111320
      }
    } else {
      warning("study-area clip could not be applied to ", basename(path),
              " (inside length ", length(inside), " vs ", nrow(df),
              " cells); map shows the full prediction grid", call. = FALSE)
    }
  }
  df
}

#' Plot a Saved Single-Shift Effect Map
#' @param path GeoTIFF path from [cast_effect_map()].
#' @param layer Raster layer name or index.
#' @param title,legend Display title and legend title.
#' @param signed Whether the contrast is signed.
#' @param limits,breaks Optional display scale limits and breaks.
#' @param mask Optional WGS84 sf reporting polygon for display clipping.
#' @param boundary,bounds,feather Optional display boundary, bounds and feathering.
#' @param method Raster display resampling: `"bilinear"` for continuous values
#'   or `"near"` for categorical layers and binary range flags. Unsigned binary
#'   layers with `limits = c(0, 1)` use a discrete outside/inside legend.
#' @return A ggplot drawn using the package cartographic assets.
#' @export
cast_plot_effect_map <- function(path, layer, title, legend, signed = FALSE,
                           limits = NULL, breaks = NULL, mask = NULL,
                           boundary = NULL, bounds = NULL, feather = 0,
                           method = c("bilinear", "near")) {
  method <- match.arg(method)
  if (is.null(path) || !file.exists(path)) stop("Missing raster for ", layer)
  df <- .cast_display_raster(path, layer = layer, method = method, mask = mask, keep_outside = TRUE,
                      feather_deg = feather)
  finite <- df$val[is.finite(df$val)]
  if (!length(finite)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, lineheight = 0.95, size = 2.4,
        colour = "grey30",
        label = paste0("Masked: no supported cells\n", layer,
                       "\nsingle-shift support = 0 after hard masking")) +
      ggplot2::xlim(-0.5, 0.5) + ggplot2::ylim(-0.5, 0.5) +
      ggplot2::labs(title = title,
        subtitle = "Effect undefined outside the training box (all cells NA)") +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 8),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 6))
    return(p)
  }
  if (is.null(limits)) {
    bound <- unname(stats::quantile(abs(finite), 0.90, na.rm = TRUE))
    if (!is.finite(bound) || bound <= 0) bound <- max(abs(finite))
    if (!is.finite(bound) || bound == 0) bound <- 1e-6
    limits <- if (signed) c(-bound, bound) else c(0, bound)
    stretch_note <- sprintf("display clipped at p90 of |response| (%s)",
                            format(bound, digits = 3))
  } else {
    stretch_note <- NULL
  }
  binary <- identical(method, "near") && !signed &&
    identical(as.numeric(limits), c(0, 1)) && all(finite %in% c(0, 1))
  if (binary) df$val <- factor(df$val, levels = c(0, 1))
  scale_fn <- if (binary) {
    function(guide) ggplot2::scale_fill_manual(
      values = c(`0` = "#EDF2F4", `1` = "#0B3D5C"),
      breaks = c("0", "1"), labels = c("Outside", "Inside"), drop = FALSE,
      na.value = NA, name = legend,
      guide = if (identical(guide, "none")) "none" else "legend")
  } else if (signed) {
    function(guide) ggplot2::scale_fill_gradient2(low = "#B35806", mid = "#F2F2F2",
      high = "#2166AC", midpoint = 0, limits = limits, oob = scales::squish,
      na.value = NA, name = legend,
      guide = if (identical(guide, "none")) "none" else "colourbar")
  } else {
    function(guide) ggplot2::scale_fill_gradientn(colours = c("#EDF2F4", "#9FD3DB", "#4FA3BE", "#256E92", "#0B3D5C"),
      limits = limits, breaks = breaks, oob = scales::squish, na.value = NA,
      name = legend, guide = if (identical(guide, "none")) "none" else "colourbar")
  }
  .cast_map_surface(df, "val", scale_fn, title,
    subtitle = paste(c(if (length(unique(finite)) == 1L)
      sprintf("Constant finite value: %g", finite[1]) else NULL, stretch_note),
      collapse = " | "),
    legend_name = legend, boundary = boundary, bounds = bounds, feather = feather)
}
