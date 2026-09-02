# ==========================================================================
# Standards-compliant China map rendering (geo-viz engine, ported).
# Reference: F:/LLQ/Very/skills/geo-viz (china.shp + nine-dash line + SCS
# inset + Nature-grade palettes + 600 dpi transparent PNG + cairo PDF).
# User-facing entry point: cast_plot_china().
# ==========================================================================

.cast_china_asset <- function(file) {
  system.file("china_basemap", file, package = "castSDM")
}

.cast_china_load_vec <- function(path) {
  v <- sf::st_read(path, quiet = TRUE)
  if (is.na(sf::st_crs(v)) || (sf::st_crs(v)$epsg %||% 0L) != 4326L) {
    v <- suppressWarnings(sf::st_transform(v, 4326))
  }
  sf::st_make_valid(v)
}

.cast_china_seq <- function(name, n = 256, rev = FALSE) {
  cols <- grDevices::hcl.colors(n, name); if (rev) rev(cols) else cols
}

.cast_china_themes <- function() list(
  blues   = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"),
  greens  = c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),
  purples = c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D"),
  oranges = c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"),
  reds    = c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
  teal    = c("#F7FCFD", "#E5F5F9", "#CCECE6", "#99D8C9", "#66C2A4", "#41AE76", "#238B45", "#006D2C", "#00441B"),
  ylgnbu  = c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"),
  ylorrd  = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
  mako    = c("#D5F0EA", "#9FDCCB", "#5FB3A3", "#2E8289", "#1E4F6B", "#16263F"),
  rocket  = c("#FBE6C8", "#F6A97A", "#E85D5D", "#B5305F", "#6E1A55", "#2C0B2E"),
  rdbu     = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
  rdylbu   = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
  spectral = c("#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"),
  brbg     = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30"),
  puor     = c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B"),
  prgn     = c("#40004B", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B"))

.cast_china_diverging <- function() c("rdbu", "rdylbu", "spectral", "brbg", "puor", "prgn")

.cast_china_palette <- function(name, reverse = FALSE) {
  th <- .cast_china_themes()
  perc <- c(viridis = "Viridis", cividis = "Cividis", inferno = "Inferno", plasma = "Plasma")
  cols <- if (!is.null(th[[name]])) th[[name]]
    else if (name %in% names(perc)) .cast_china_seq(perc[[name]])
    else if (name == "hillshade") grDevices::grey.colors(256, 0.05, 0.95)
    else cli::cli_abort("Unknown palette {.val {name}}. Sequential: {.val {setdiff(names(th), .cast_china_diverging())}}; diverging: {.val {.cast_china_diverging()}}.")
  if (isTRUE(reverse)) rev(cols) else cols
}

.cast_china_build_map <- function(df, china_sf, dashline_sf, title, legend_title,
                                  palette, limits, midpoint = NA_real_) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  BORDER_COL <- "#4E5963"
  INSET_XLIM <- c(105, 126); INSET_YLIM <- c(2, 26)
  bounds <- as.numeric(sf::st_bbox(china_sf))
  rescaler <- if (!is.na(midpoint)) {
    function(x, to = c(0, 1), from = range(x, na.rm = TRUE))
      scales::rescale_mid(x, to, from, mid = midpoint)
  } else scales::rescale
  fill_scale <- ggplot2::scale_fill_gradientn(
    colours = palette, name = legend_title, limits = limits,
    rescaler = rescaler, oob = scales::oob_squish, na.value = "transparent")
  main_theme <- ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12,
                                         margin = ggplot2::margin(b = 3)),
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.position = "right",
      legend.key.width = ggplot2::unit(0.42, "cm"),
      legend.key.height = ggplot2::unit(1.4, "cm"),
      plot.margin = ggplot2::margin(4, 4, 4, 4))
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)) +
    fill_scale +
    ggplot2::geom_sf(data = china_sf, fill = NA, colour = BORDER_COL, linewidth = 0.22)
  if (!is.null(dashline_sf))
    p <- p + ggplot2::geom_sf(data = dashline_sf, fill = NA, colour = BORDER_COL,
                              linewidth = 0.30)
  p <- p + ggplot2::coord_sf(xlim = bounds[c(1, 3)], ylim = bounds[c(2, 4)],
                             expand = FALSE, datum = NA) +
    ggplot2::labs(title = title) + main_theme
  sub <- df[df$x >= INSET_XLIM[1] & df$x <= INSET_XLIM[2] &
              df$y >= INSET_YLIM[1] & df$y <= INSET_YLIM[2], , drop = FALSE]
  inset <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = china_sf, fill = "#F5F6F7", colour = NA)
  if (nrow(sub) > 0)
    inset <- inset + ggplot2::geom_raster(
      data = sub, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)) + fill_scale
  inset <- inset +
    ggplot2::geom_sf(data = china_sf, fill = NA, colour = BORDER_COL, linewidth = 0.18)
  if (!is.null(dashline_sf))
    inset <- inset + ggplot2::geom_sf(data = dashline_sf, fill = NA, colour = BORDER_COL,
                                      linewidth = 0.25)
  inset <- inset +
    ggplot2::coord_sf(xlim = INSET_XLIM, ylim = INSET_YLIM, expand = FALSE, datum = NA) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none",
                   plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   panel.border = ggplot2::element_rect(fill = NA, colour = BORDER_COL, linewidth = 0.4),
                   plot.margin = ggplot2::margin(0, 0, 0, 0))
  xr <- bounds[3] - bounds[1]; yr <- bounds[4] - bounds[2]
  ins_w <- xr * 0.16; ins_h <- yr * 0.30
  xmax <- bounds[3] - xr * 0.012; xmin <- xmax - ins_w
  ymin <- bounds[2] + yr * 0.012; ymax <- ymin + ins_h
  p + ggplot2::annotation_custom(grob = ggplot2::ggplotGrob(inset),
                                 xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

#' Render a Raster as a Standards-Compliant China Map
#'
#' Wraps the geo-viz engine: China boundary + South China Sea nine-dash line
#' + inset, Nature-grade colour-blind-friendly palettes, transparent
#' background, 600 dpi PNG and cairo PDF export. The input raster is
#' reprojected to a unified WGS84 display grid with GDAL (not terra, which
#' segfaults on some machines) and masked to China.
#'
#' @param raster A single-layer `SpatRaster`.
#' @param output Output path without extension (writes `.png` and `.pdf`).
#' @param title Map title. Default empty.
#' @param legend Legend title. Default empty.
#' @param palette Palette name: sequential themes (`blues`, `greens`,
#'   `teal`, `mako`, `rocket`, ...), diverging themes (`rdbu`, `rdylbu`,
#'   `spectral`, `brbg`, `puor`, `prgn`; pair with `midpoint`), perceptual
#'   (`viridis`, `cividis`, `inferno`, `plasma`), or `hillshade`.
#' @param midpoint Diverging midpoint. Default `NA` (no mid-rescaling).
#' @param clamp Logical. Clamp colour limits to the 2-98th percentiles.
#' @param reverse Reverse the palette.
#' @param width,height,dpi Figure size and resolution.
#'
#' @return Invisibly, the output path prefix (files are written).
#' @export
cast_plot_china <- function(raster, output, title = "", legend = "",
                            palette = "viridis", midpoint = NA_real_,
                            clamp = FALSE, reverse = FALSE,
                            width = 7, height = 6.2, dpi = 600) {
  if (!requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    cli::cli_abort("Packages {.pkg terra} and {.pkg sf} are required.")
  }
  sf::sf_use_s2(FALSE)
  if (terra::nlyr(raster) != 1L) {
    cli::cli_abort("{.arg raster} must have exactly one layer; subset it first.")
  }
  china_shp <- .cast_china_asset("china.shp")
  if (!nzchar(china_shp) || !file.exists(china_shp)) {
    cli::cli_abort("Basemap assets missing from the installation.")
  }
  china_sf <- .cast_china_load_vec(china_shp)
  dash_shp <- .cast_china_asset("dashline.shp")
  dashline_sf <- if (nzchar(dash_shp) && file.exists(dash_shp))
    .cast_china_load_vec(dash_shp) else NULL
  china_vect <- terra::vect(china_shp)
  if (!is.na(terra::crs(china_vect)) &&
      (terra::crs(china_vect, describe = TRUE)$code %||% "") != "4326")
    china_vect <- terra::project(china_vect, "EPSG:4326")

  cache <- file.path(tempdir(), "cast_china_cache")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  dst <- file.path(cache, "display.tif")
  wopt <- c("-t_srs", "EPSG:4326",
            "-te", "73.40", "18.10", "135.15", "53.60",
            "-tr", "0.05", "0.05", "-r", "average",
            "-overwrite", "-of", "GTiff", "-co", "COMPRESS=DEFLATE",
            "-multi", "-ovr", "AUTO")
  sf::gdal_utils("warp", source = terra::sources(raster)[[1]], destination = dst,
                 options = wopt, quiet = TRUE)
  r <- terra::mask(terra::rast(dst)[[1]], china_vect)
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "value")
  if (!nrow(df)) cli::cli_abort("No valid data inside China after masking.")
  lim <- if (isTRUE(clamp)) {
    l <- as.numeric(stats::quantile(df$value, c(0.02, 0.98), na.rm = TRUE))
    if (!all(is.finite(l)) || diff(l) == 0) range(df$value, na.rm = TRUE) else l
  } else NULL
  pal <- .cast_china_palette(palette, reverse = reverse)
  p <- .cast_china_build_map(df, china_sf, dashline_sf, title, legend, pal, lim, midpoint)
  outdir <- dirname(normalizePath(output, mustWork = FALSE))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(output, ".png"), p, width = width, height = height,
                  dpi = dpi, bg = "transparent", limitsize = FALSE)
  ggplot2::ggsave(paste0(output, ".pdf"), p, width = width, height = height,
                  device = grDevices::cairo_pdf, bg = "transparent", limitsize = FALSE)
  invisible(output)
}
