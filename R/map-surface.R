# Shared cartographic renderer; preserves explicit cell size and package assets.
.cast_regular_raster_grid <- function(df, cell_size) {
  if (is.null(cell_size) || nrow(df) < 2L ||
      any(!is.finite(df$lon) | !is.finite(df$lat))) return(NULL)
  dx <- cell_size[1]; dy <- utils::tail(cell_size, 1)
  regular <- function(x, step) {
    d <- diff(sort(unique(x)))
    length(d) && all(abs(d/step-round(d/step)) < 1e-5)
  }
  if (!regular(df$lon, dx) || !regular(df$lat, dy)) return(NULL)
  # Transparent anchors retain fully absent grid rows/columns without
  # inventing values or allocating a second complete Cartesian table.
  x0 <- min(df$lon); y0 <- min(df$lat)
  ix <- round((df$lon-x0)/dx); iy <- round((df$lat-y0)/dy)
  missing_x <- setdiff(seq.int(0L, max(ix)), unique(ix))
  missing_y <- setdiff(seq.int(0L, max(iy)), unique(iy))
  n <- length(missing_x)+length(missing_y)
  if (n) {
    anchors <- df[rep(NA_integer_, n), , drop = FALSE]
    anchors$lon <- c(x0+missing_x*dx, rep(x0, length(missing_y)))
    anchors$lat <- c(rep(y0, length(missing_x)), y0+missing_y*dy)
    df <- rbind(df, anchors)
  }
  df
}

.cast_map_surface <- function(df, value_col, scale_fn, title, subtitle = NULL,
                             legend_name = NULL,
                             bg = .cast_map_bg(), boundary = NULL,
                             bounds = NULL, context = 0, feather = 0,
                             alpha_col = NULL, cell_size = 0.02,
                             basemap = "china", raster = FALSE) {
  outside <- NULL
  if (".inside" %in% names(df)) {
    outside <- df[!df$.inside, , drop = FALSE]
    df <- df[df$.inside, , drop = FALSE]
  }
  bm <- if (basemap == "none") NULL else load_basemap(basemap)
  if (!is.null(bm) && !is.null(boundary)) {
    boundary <- sf::st_transform(boundary, 4326)
    bm <- suppressWarnings(sf::st_intersection(bm, boundary))
  }
  p <- ggplot2::ggplot()
  if (!is.null(bm)) {
    p <- p + ggplot2::geom_sf(data = bm, fill = bg, color = "#bdc3c7", linewidth = 0.2)
  }
  if (!is.null(outside) && nrow(outside) && context > 0) {
    p <- p + ggplot2::geom_tile(data = outside,
      ggplot2::aes(x = .data$lon, y = .data$lat), fill = "#FFFFFF",
      alpha = context, width = cell_size[1], height = utils::tail(cell_size, 1))
  }
  map_aes <- if (!is.null(alpha_col) && alpha_col %in% names(df)) {
    ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data[[value_col]],
                 alpha = .data[[alpha_col]])
  } else if (!is.null(feather) && feather > 0 && ".dist_deg" %in% names(df)) {
    ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data[[value_col]],
                 alpha = .data$.dist_deg)
  } else {
    ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data[[value_col]])
  }
  if (!is.null(alpha_col) && alpha_col %in% names(df)) {
    p <- p + ggplot2::scale_alpha_identity(guide = "none")
  }
  # Dense regular display grids are embedded as a raster in PDF; boundaries,
  # labels and legends remain vector. Sparse/nonregular inputs retain tiles.
  raster_df <- if (isTRUE(raster)) .cast_regular_raster_grid(df, cell_size)
  p <- p + if (!is.null(raster_df)) ggplot2::geom_raster(data = raster_df, mapping = map_aes,
    na.rm = TRUE) else
    ggplot2::geom_tile(data = df, mapping = map_aes,
      width = cell_size[1], height = utils::tail(cell_size, 1))
  p <- p + scale_fn("main")
  if (!is.null(boundary)) {
    p <- p + ggplot2::geom_sf(data = boundary, fill = NA,
      color = "#4D4D4D", linewidth = 0.3, linetype = "dashed")
  }
  if ((is.null(alpha_col) || !alpha_col %in% names(df)) &&
      !is.null(feather) && feather > 0 && ".dist_deg" %in% names(df)) {
    p <- p + ggplot2::scale_alpha(range = c(0, 1), limits = c(0, feather),
      guide = "none", na.value = 1, oob = scales::squish)
  }
  p <- p +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL,
                  fill = legend_name) +
    ggplot2::theme_void(base_size = 10,
      base_family = getOption("castSDM.font_family", "sans")) +
    ggplot2::theme(
      text = ggplot2::element_text(family = getOption("castSDM.font_family", "sans")),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 8),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 6.5),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7),
      legend.position = "right",
      legend.key.width = ggplot2::unit(0.5, "cm"),
      legend.key.height = ggplot2::unit(1.5, "cm"),
      plot.margin = ggplot2::margin(4, 12, 4, 6),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank())
  p <- .add_china_dashline(p, basemap)
  p <- .add_china_outline(p, basemap)
  p <- p + if (!is.null(bounds)) {
    ggplot2::coord_sf(xlim = bounds[c(1L, 3L)], ylim = bounds[c(2L, 4L)],
                      expand = FALSE, label_graticule = "")
  } else {
    bb <- if (basemap == "china") .china_main_bounds() else NULL
    if (is.null(bb)) {
      ggplot2::coord_sf(expand = FALSE, label_graticule = "")
    } else {
      ggplot2::coord_sf(xlim = bb[c(1L, 3L)], ylim = bb[c(2L, 4L)],
                        expand = FALSE, label_graticule = "")
    }
  }
  wrap_label <- function(x, width) {
    if (is.null(x)) return(NULL)
    lines <- strsplit(gsub("_", " ", x, fixed = TRUE), "\n", fixed = TRUE)[[1L]]
    paste(vapply(lines, function(line) paste(strwrap(line, width = width),
      collapse = "\n"), character(1)), collapse = "\n")
  }
  p <- p + ggplot2::labs(title = wrap_label(title, 48),
                         subtitle = wrap_label(subtitle, 72))
  .add_china_south_sea_inset(p, basemap, data = df, value_var = value_col,
    bg_fill = bg, raster = TRUE, scale = scale_fn("none"),
    cell_size = cell_size, alpha_col = alpha_col)
}
