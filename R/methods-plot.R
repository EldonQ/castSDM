# Plot Methods ----------------------------------------------------------------

#' Plot Variable Selection
#'
#' Bar chart showing predictor score and selection status.
#'
#' @param x A `cast_select` object (from [cast_select()]).
#' @param var_labels Optional named character vector for display labels.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_select <- function(x, var_labels = NULL, ...) {
  check_suggested("ggplot2", "for plotting")

  scr <- x$scores
  scr$is_selected <- scr$variable %in% x$selected

  scr$status <- ifelse(scr$is_selected, "selected", "not selected")

  imp_candidates <- c("combined_score", "rf_importance")
  imp_col <- imp_candidates[imp_candidates %in% names(scr)][1]
  if (is.na(imp_col) || !length(imp_col)) {
    scr$importance_plot <- as.numeric(scr$is_selected)
    imp_col <- "importance_plot"
  }
  scr[[imp_col]] <- suppressWarnings(as.numeric(scr[[imp_col]]))
  scr[[imp_col]][!is.finite(scr[[imp_col]])] <- 0
  scr <- scr[order(-scr[[imp_col]]), ]

  if (!is.null(var_labels)) {
    scr$display <- ifelse(
      scr$variable %in% names(var_labels),
      var_labels[scr$variable], scr$variable
    )
  } else {
    scr$display <- scr$variable
  }
  scr$display <- factor(scr$display, levels = rev(scr$display))

  status_colors <- c(selected = "#2166AC", `not selected` = "grey75")

  n_sel <- sum(scr$is_selected)
  sub_txt <- sprintf("%d / %d variables selected", n_sel, nrow(scr))
  sub_txt <- paste0(sub_txt, " | method = ", x$method %||% "unknown")

  p <- ggplot2::ggplot(scr, ggplot2::aes(
    x = .data[[imp_col]], y = .data$display,
    fill = .data$status, alpha = .data$is_selected
  )) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = status_colors, name = "Status") +
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 0.95, "FALSE" = 0.35), guide = "none"
    ) +
    ggplot2::labs(
      title = "Variable screening",
      subtitle = sub_txt,
      x = if (identical(imp_col, "importance_plot")) "Selection indicator" else "Predictor score",
      y = ""
    ) +
    theme_cast(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_line(
        color = "grey93", linewidth = 0.3
      )
    )
  p
}


#' Plot Spatial Habitat Suitability Map
#'
#' @param x A `cast_predict` object.
#' @param model Character. Which model's HSS to plot. Default is the first
#'   model.
#' @param basemap Character. `"world"`, `"china"`, or `"none"`.
#' @param title Optional character string for plot title.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_predict <- function(x, model = NULL, basemap = "world",
                              title = NULL, ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("sf", "for geographic mapping")

  pred <- x$predictions
  model <- model %||% x$models[1]
  hss_col <- paste0("HSS_", model)

  if (!hss_col %in% names(pred)) {
    cli::cli_abort(c(
      "Model {.val {model}} not found in predictions.",
      i = "Available models: {.val {x$models}}.",
      i = "Use {.code plot(pred, model = \"{x$models[1]}\")}."
    ))
  }

  title <- title %||% sprintf("Habitat Suitability - %s", toupper(model))

  p <- ggplot2::ggplot()
  if (basemap != "none") {
    basemap_sf <- load_basemap(basemap)
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf,
        fill = "#f4f6f7", color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  p <- p +
    ggplot2::geom_point(
      data = pred,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data[[hss_col]]),
      size = 0.4, alpha = 0.85
    ) +
    ggplot2::scale_color_viridis_c(option = "turbo", name = "HSS", limits = c(0, 1)) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(
      base_size = 10,
      base_family = getOption("castSDM.font_family", "Arial")
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = getOption("castSDM.font_family", "Arial")),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.position = "right",
      legend.key.width = ggplot2::unit(0.5, "cm"),
      legend.key.height = ggplot2::unit(1.5, "cm")
    )

  p <- .add_china_dashline(p, basemap)
  p <- .add_china_outline(p, basemap) + .coord_for_basemap(basemap)
  .add_china_south_sea_inset(p, basemap)
}


#' Plot Evaluation Metrics Comparison
#'
#' Multi-panel bar chart comparing AUC, TSS, and CBI across fitted models.
#'
#' @param x A `cast_eval` object.
#' @param metrics Character vector. Which metrics to show. Default
#'   `c("auc","tss","cbi")`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (faceted).
#' @export
plot.cast_eval <- function(x, metrics = c("auc", "tss", "cbi"), ...) {
  check_suggested("ggplot2", "for plotting")

  m <- x$metrics
  model_colors <- c(
    rf = "#4DBBD5", brt = "#3C5488", maxent = "#B09C85",
    gam = "#00A087", esm = "#E64B35"
  )

  src <- if (isTRUE(x$cv_source)) "Spatial CV" else "Hold-out test"

  mean_cols <- paste0(metrics, "_mean")
  present   <- intersect(mean_cols, names(m))
  if (length(present) == 0) {
    present <- intersect(c("auc_mean", "tss_mean"), names(m))
  }

  rows <- list()
  for (col in present) {
    metric_nm <- sub("_mean$", "", col)
    rows[[col]] <- data.frame(
      model  = m$model,
      metric = toupper(metric_nm),
      value  = m[[col]],
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(rbind, rows)
  long$metric <- factor(long$metric, levels = toupper(sub("_mean$", "", present)))
  long$model  <- factor(long$model, levels = names(model_colors))

  ggplot2::ggplot(long, ggplot2::aes(
    x = .data$model, y = .data$value, fill = .data$model
  )) +
    ggplot2::geom_col(width = 0.65, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$value)),
      vjust = -0.4, size = 2.8,
      family = getOption("castSDM.font_family", "Arial"),
      fontface = "bold"
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 2) +
    ggplot2::scale_fill_manual(values = model_colors, guide = "none") +
    ggplot2::labs(title = "Model Performance Comparison", subtitle = src,
                  x = "", y = "Score") +
    theme_cast() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}


#' Plot Spatial CV Fold Map and Metrics
#'
#' Two-panel figure: (left) geographic fold assignment map; (right) per-fold
#' metric box/dot plot.
#'
#' @param x A `cast_cv` object.
#' @param lon Numeric vector. Longitudes of the data used in [cast_cv()].
#' @param lat Numeric vector. Latitudes of the data used in [cast_cv()].
#' @param metric Character. Metric to show in right panel. Default `"auc"`.
#' @param basemap Character. `"world"`, `"china"`, or `"none"`.
#' @param ... Ignored.
#'
#' @return A `patchwork` combined plot, or single ggplot if patchwork absent.
#' @export
plot.cast_cv <- function(x, lon = NULL, lat = NULL,
                         metric = "auc", basemap = "world", ...) {
  check_suggested("ggplot2", "for plotting")

  model_colors <- c(
    rf = "#4DBBD5", brt = "#3C5488", maxent = "#B09C85",
    gam = "#00A087", esm = "#E64B35"
  )
  fold_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62"
  )

  fd   <- x$fold_metrics
  mcol <- metric
  if (!mcol %in% names(fd)) {
    cli::cli_warn("Metric '{metric}' not in fold_metrics. Using 'auc'.")
    mcol <- "auc"
  }

  fd$model <- factor(fd$model, levels = names(model_colors))
  fd$fold  <- factor(fd$fold)

  p_metric <- ggplot2::ggplot(
    fd,
    ggplot2::aes(x = .data$fold, y = .data[[mcol]],
                 color = .data$model, group = .data$model)
  ) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = model_colors, name = "Model") +
    ggplot2::labs(
      title = sprintf("Per-fold %s", toupper(mcol)),
      subtitle = sprintf("%d-fold spatial %s CV", x$k, x$block_method),
      x = "Fold", y = toupper(mcol)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_cast()

  if (is.null(lon) || is.null(lat)) return(p_metric)

  map_df <- data.frame(lon = lon, lat = lat, fold = factor(x$folds))
  p_map <- ggplot2::ggplot()
  if (basemap != "none") {
    bm <- load_basemap(basemap)
    if (!is.null(bm)) {
      p_map <- p_map + ggplot2::geom_sf(
        data = bm, fill = "#f4f6f7",
        color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  p_map <- p_map +
    ggplot2::geom_point(
      data = map_df,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$fold),
      size = 0.8, alpha = 0.7
    ) +
    ggplot2::scale_color_manual(values = fold_colors[seq_len(x$k)], name = "Fold") +
    .coord_for_basemap(basemap) +
    ggplot2::labs(title = "Spatial fold assignment") +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  p_map <- .add_china_dashline(p_map, basemap)
  p_map <- .add_china_outline(p_map, basemap)
  p_map <- .add_china_south_sea_inset(p_map, basemap)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_map + p_metric + patchwork::plot_layout(widths = c(1.4, 1))
  } else {
    p_metric
  }
}


#' Plot castSDM Pipeline Result
#'
#' Currently renders the variable-selection panel.
#'
#' @param x A `cast_result` object.
#' @param var_labels Optional named character vector for display labels.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_result <- function(x, var_labels = NULL, ...) {
  check_suggested("ggplot2", "for plotting")
  plot.cast_select(x$screen, var_labels = var_labels)
}


#' Plot Ensemble Prediction Map
#'
#' @param x A `cast_ensemble` object.
#' @param basemap Character. `"world"`, `"china"`, or `"none"`.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
plot.cast_ensemble <- function(x, basemap = "world", ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("sf", "for geographic mapping")

  pred <- x$predictions
  if (!all(c("lon", "lat", "hss_ensemble") %in% names(pred))) {
    cli::cli_abort("Ensemble predictions must contain lon, lat, hss_ensemble.")
  }

  p <- ggplot2::ggplot()
  if (basemap != "none") {
    basemap_sf <- load_basemap(basemap)
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf,
        fill = "#f4f6f7", color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  if (nrow(pred) > 2e5) {
    p <- p +
      ggplot2::geom_raster(
        data = pred,
        ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$hss_ensemble)
      ) +
      ggplot2::scale_fill_viridis_c(option = "turbo", name = "Ensemble\nHSS", limits = c(0, 1))
  } else {
    p <- p +
      ggplot2::geom_point(
        data = pred,
        ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$hss_ensemble),
        size = 0.4, alpha = 0.85
      ) +
      ggplot2::scale_color_viridis_c(option = "turbo", name = "Ensemble\nHSS", limits = c(0, 1))
  }
  p <- p +
    ggplot2::labs(
      title = sprintf("Ensemble Habitat Suitability (%s)", x$method),
      subtitle = sprintf("Threshold = %.3f", x$threshold)
    ) +
    ggplot2::theme_void(
      base_size = 10,
      base_family = getOption("castSDM.font_family", "Arial")
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = getOption("castSDM.font_family", "Arial")),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.position = "right",
      legend.key.width = ggplot2::unit(0.5, "cm"),
      legend.key.height = ggplot2::unit(1.5, "cm")
    )

  p <- .add_china_dashline(p, basemap)
  p <- .add_china_outline(p, basemap) + .coord_for_basemap(basemap)
  .add_china_south_sea_inset(p, basemap)
}


#' Plot Future Climate Projection
#'
#' @param x A `cast_project` object.
#' @param scenario Character. Which scenario to plot. Default is the first.
#' @param basemap Character. `"world"`, `"china"`, or `"none"`.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
plot.cast_project <- function(x, scenario = NULL, basemap = "world", ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("sf", "for geographic mapping")

  scenario <- scenario %||% names(x$future)[1]
  if (!scenario %in% names(x$future)) {
    cli::cli_abort("Scenario {.val {scenario}} not found. Available: {.val {names(x$future)}}.")
  }

  fut <- x$future[[scenario]]
  change <- x$changes[[scenario]]
  if (!all(c("lon", "lat", "change") %in% names(change))) {
    cli::cli_abort("Change data must contain lon, lat, change columns.")
  }

  change_colors <- c(
    gain = "#27AE60", loss = "#C0392B", stable_present = "#3498DB",
    stable_absent = "grey85"
  )

  p <- ggplot2::ggplot()
  if (basemap != "none") {
    basemap_sf <- load_basemap(basemap)
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf,
        fill = "#f4f6f7", color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  if (nrow(change) > 2e5) {
    p <- p +
      ggplot2::geom_raster(
        data = change,
        ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$change)
      ) +
      ggplot2::scale_fill_manual(values = change_colors, name = "Range\nChange")
  } else {
    p <- p +
      ggplot2::geom_point(
        data = change,
        ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$change),
        size = 0.4, alpha = 0.8
      ) +
      ggplot2::scale_color_manual(values = change_colors, name = "Range\nChange")
  }
  p <- p +
    ggplot2::labs(
      title = sprintf("Range Change: %s", scenario),
      subtitle = {
        s <- x$stats[x$stats$scenario == scenario, ]
        if (nrow(s) > 0)
          sprintf("Gain=%d | Loss=%d | Stable=%d | Shift=%.1f km",
                  s$n_gain[1], s$n_loss[1], s$n_stable_present[1], s$centroid_shift_km[1])
        else ""
      }
    ) +
    ggplot2::theme_void(
      base_size = 10,
      base_family = getOption("castSDM.font_family", "Arial")
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = getOption("castSDM.font_family", "Arial")),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.position = "right"
    )

  p <- .add_china_dashline(p, basemap)
  p <- .add_china_outline(p, basemap) + .coord_for_basemap(basemap)
  .add_china_south_sea_inset(p, basemap)
}


# ========================================================================
# Internal helpers
# ========================================================================

#' castSDM Publication Theme
#' @keywords internal
#' @noRd
theme_cast <- function(base_size = 11) {
  ggplot2::theme_minimal(
    base_size = base_size,
    base_family = getOption("castSDM.font_family", "Arial")
  ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size),
      plot.subtitle = ggplot2::element_text(hjust = 0, color = "grey40", size = base_size - 2),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      strip.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )
}

#' Load Basemap Shapefile from Package
#'
#' @param type `"world"`, `"china"`, or `"dashline"`.
#' @return An `sf` object or `NULL`.
#' @keywords internal
#' @noRd
load_basemap <- function(type = "world") {
  if (!requireNamespace("sf", quietly = TRUE)) return(NULL)

  shp_name <- switch(type,
    world = "countries.shp",
    china = "china.shp",
    dashline = "dashline.shp",
    NULL
  )
  if (is.null(shp_name)) return(NULL)

  shp_path <- system.file("basemap", shp_name, package = "castSDM")
  if (shp_path == "") return(NULL)

  sf::sf_use_s2(FALSE)
  basemap <- tryCatch(
    sf::st_read(shp_path, quiet = TRUE),
    error = function(e) NULL
  )

  if (!is.null(basemap)) {
    if (is.na(sf::st_crs(basemap)) || sf::st_crs(basemap)$epsg != 4326L) {
      basemap <- sf::st_transform(basemap, 4326)
    }
    basemap <- sf::st_make_valid(basemap)
  }
  basemap
}

#' Get standard bounds for China maps
#' @keywords internal
#' @noRd
.china_main_bounds <- function() {
  china <- load_basemap("china")
  if (is.null(china)) return(NULL)
  as.numeric(sf::st_bbox(china))
}

#' Use a fixed, non-distorted extent for China maps
#' @keywords internal
#' @noRd
.coord_for_basemap <- function(basemap) {
  if (!identical(basemap, "china")) return(ggplot2::coord_sf(expand = FALSE))
  bounds <- .china_main_bounds()
  if (is.null(bounds)) return(ggplot2::coord_sf(expand = FALSE))
  ggplot2::coord_sf(
    xlim = bounds[c(1, 3)],
    ylim = bounds[c(2, 4)],
    expand = FALSE,
    datum = NA
  )
}

#' Redraw the China outline above prediction layers
#' @keywords internal
#' @noRd
.add_china_outline <- function(plot, basemap) {
  if (!identical(basemap, "china")) return(plot)
  china <- load_basemap("china")
  if (is.null(china)) return(plot)
  plot + ggplot2::geom_sf(data = china, fill = NA, colour = "#4E5963", linewidth = 0.30)
}

#' Redraw the South China Sea dashed line on the China main map
#' @keywords internal
#' @noRd
.add_china_dashline <- function(plot, basemap) {
  if (!identical(basemap, "china")) return(plot)
  dashline <- load_basemap("dashline")
  if (is.null(dashline)) return(plot)
  plot + ggplot2::geom_sf(
    data = dashline, fill = NA, colour = "#4E5963", linewidth = 0.28, linetype = "solid"
  )
}

#' Draw the South China Sea dashed-line inset inside a China map
#'
#' @keywords internal
#' @noRd
.add_china_south_sea_inset <- function(plot, basemap) {
  if (!identical(basemap, "china") ||
      !requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) return(plot)

  china <- load_basemap("china")
  dashline <- load_basemap("dashline")
  bounds <- .china_main_bounds()
  if (is.null(china) || is.null(dashline) || is.null(bounds)) return(plot)

  inset <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = china, fill = "#F7F8F8", colour = "#4E5963", linewidth = 0.20
    ) +
    ggplot2::geom_sf(
      data = dashline, fill = NA, colour = "#4E5963", linewidth = 0.22, linetype = "solid"
    ) +
    ggplot2::coord_sf(xlim = c(105, 126), ylim = c(2, 26), expand = FALSE, datum = NA) +
    ggplot2::theme_void(base_family = getOption("castSDM.font_family", "Arial")) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_rect(fill = NA, colour = "#4E5963", linewidth = 0.35),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  xmin <- bounds[3] - 10.6
  xmax <- bounds[3] - 0.7
  ymin <- bounds[2] + 0.7
  ymax <- bounds[2] + 10.6

  plot + ggplot2::annotation_custom(
    grob = ggplot2::ggplotGrob(inset),
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
  )
}
