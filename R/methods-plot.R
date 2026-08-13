# Plot Methods ----------------------------------------------------------------

#' Plot Variable Selection
#'
#' Bar chart showing predictor score and selection status.
#'
#' @param x A `cast_select` object (from [cast_select()]).
#' @param var_labels Optional named character vector for display labels.
#' @param top Integer. Show only the `top` highest-scoring predictors (union
#'   with the retained set). Default `20`; use `NULL` for all tested predictors.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_select <- function(x, var_labels = NULL, top = 20L, ...) {
  check_suggested("ggplot2", "for plotting")

  scr <- x$scores
  scr$is_selected <- scr$variable %in% x$selected

  imp_candidates <- c("freq", "abs_statistic", "cpi", "rf_importance",
                      "combined_score")
  imp_col <- imp_candidates[imp_candidates %in% names(scr)][1]
  if (is.na(imp_col) || !length(imp_col)) {
    scr$importance_plot <- as.numeric(scr$is_selected)
    imp_col <- "importance_plot"
  }
  scr[[imp_col]] <- suppressWarnings(as.numeric(scr[[imp_col]]))

  n_total <- nrow(scr)
  n_sel <- sum(scr$is_selected)

  # The conditional screen tests a candidate subset, so untested predictors
  # carry no score. Plotting all of them crushes the axis into an unreadable
  # band; keep only the evaluated predictors, then the top-scoring slice.
  tested <- scr[is.finite(scr[[imp_col]]), , drop = FALSE]
  n_tested <- nrow(tested)
  tested <- tested[order(-tested[[imp_col]]), , drop = FALSE]

  n_top <- if (is.null(top) || !is.finite(top)) n_tested else as.integer(top)
  keep <- union(utils::head(tested$variable, n_top),
                tested$variable[tested$is_selected])
  scr <- tested[tested$variable %in% keep, , drop = FALSE]
  scr <- scr[order(-scr[[imp_col]]), , drop = FALSE]

  raw_disp <- .cast_var_label(scr$variable)
  if (!is.null(var_labels)) {
    hit <- scr$variable %in% names(var_labels)
    raw_disp[hit] <- var_labels[scr$variable[hit]]
  }
  scr$display <- factor(raw_disp, levels = rev(raw_disp))

  pal <- .cast_category_palette()
  scr$category <- factor(.cast_var_category(scr$variable), levels = names(pal))

  sub_txt <- sprintf(
    "%d / %d predictors retained \u00b7 %d tested by %s",
    n_sel, n_total, n_tested, toupper(x$method %||% "CPI")
  )

  x_lab <- if (identical(imp_col, "abs_statistic")) {
    if (identical(x$method, "dml")) "|DML statistic|" else "|CPI statistic|"
  } else if (identical(imp_col, "cpi")) {
    "Conditional predictive impact"
  } else if (identical(imp_col, "freq")) {
    "Fold selection frequency"
  } else if (identical(imp_col, "importance_plot")) {
    "Selection indicator"
  } else {
    "Predictor score"
  }

  p <- ggplot2::ggplot(scr, ggplot2::aes(
    x = .data[[imp_col]], y = .data$display,
    fill = .data$category, alpha = .data$is_selected
  )) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_manual(
      values = pal, name = "Predictor class", drop = TRUE
    ) +
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.32),
      breaks = c("TRUE", "FALSE"),
      labels = c("retained", "screened out"),
      name = "Selection"
    ) +
    ggplot2::guides(
      alpha = ggplot2::guide_legend(override.aes = list(fill = "grey35"))
    ) +
    ggplot2::labs(
      title = "Variable screening",
      subtitle = sub_txt,
      x = x_lab,
      y = ""
    ) +
    theme_cast(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
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
        fill = .cast_map_bg(), color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  p <- p +
    ggplot2::geom_point(
      data = pred,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data[[hss_col]]),
      size = 0.4, alpha = 0.85
    ) +
    .cast_suitability_scale("colour", name = "HSS") +
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
  p <- .add_china_outline(p, basemap) + .coord_for_map(basemap)
  .add_china_south_sea_inset(
    p, basemap, data = pred, value_var = hss_col, bg_fill = .cast_map_bg(),
    scale = .cast_suitability_scale("colour", guide = "none")
  )
}


#' Plot Evaluation Metrics Comparison
#'
#' Cleveland-style dot plot comparing AUC, TSS, and CBI across fitted models:
#' one row per model, one dodged dot per metric, values labelled directly.
#'
#' @param x A `cast_eval` object.
#' @param metrics Character vector. Which metrics to show. Default
#'   `c("auc","tss","cbi")`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_eval <- function(x, metrics = c("auc", "tss", "cbi"), ...) {
  check_suggested("ggplot2", "for plotting")

  m <- x$metrics
  metric_colors <- c(
    AUC = "#0072B2", TSS = "#D55E00", CBI = "#009E73",
    LOGLOSS = "#CC79A7"
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
  long <- long[!is.na(long$value), , drop = FALSE]
  long$metric <- factor(long$metric, levels = toupper(sub("_mean$", "", present)))
  long$model  <- factor(long$model, levels = rev(m$model))

  # Numeric y positions with manual per-metric offsets keep stems horizontal
  n_metric <- nlevels(long$metric)
  long$y_pos <- as.numeric(long$model) +
    (as.integer(long$metric) - (n_metric + 1) / 2) * 0.18

  ggplot2::ggplot(long, ggplot2::aes(
    x = .data$value, y = .data$y_pos, color = .data$metric
  )) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$value, yend = .data$y_pos),
      linewidth = 0.5, alpha = 0.45
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$value)),
      hjust = -0.3, size = 2.6,
      family = getOption("castSDM.font_family", "Arial"),
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = metric_colors, name = NULL) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(levels(long$model)), labels = levels(long$model),
      expand = ggplot2::expansion(mult = c(0.12, 0.12))
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = "Model Performance Comparison", subtitle = src,
      x = "Score", y = NULL
    ) +
    theme_cast() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank()
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
      subtitle = sprintf(
        "%d-fold spatial (%s) block CV \u2014 folds are geographically disjoint",
        x$k, x$block_method
      ),
      x = "Spatial fold", y = toupper(mcol)
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
        data = bm, fill = .cast_map_bg(),
        color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  # Extend the palette for k > 10 folds instead of yielding NA colours.
  fold_pal <- if (x$k <= length(fold_colors)) {
    fold_colors[seq_len(x$k)]
  } else {
    grDevices::colorRampPalette(fold_colors)(x$k)
  }

  p_map <- p_map +
    ggplot2::geom_point(
      data = map_df,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$fold),
      size = 0.8, alpha = 0.7
    ) +
    ggplot2::scale_color_manual(values = fold_pal, name = "Fold") +
    ggplot2::labs(
      title = "Spatial block folds",
      subtitle = "Contiguous geographic blocks, not random points"
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 8),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  p_map <- .add_china_dashline(p_map, basemap)
  p_map <- .add_china_outline(p_map, basemap) + .coord_for_map(basemap)
  p_map <- .add_china_south_sea_inset(
    p_map, basemap, data = map_df, value_var = "fold", bg_fill = .cast_map_bg(),
    scale = ggplot2::scale_color_manual(values = fold_pal, guide = "none")
  )

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
#' @param top Integer passed to [plot.cast_select()]. Default `20`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_result <- function(x, var_labels = NULL, top = 20L, ...) {
  check_suggested("ggplot2", "for plotting")
  plot.cast_select(x$screen, var_labels = var_labels, top = top)
}


#' Plot Screening-Method Comparison
#'
#' Tile matrix of predictor (rows) by screening method (columns); a coloured
#' tile marks retention, coloured by ecological class. It contrasts the
#' conditional castSDM screen with the conventional associational baselines
#' from [cast_screen_comparison()].
#'
#' @param x A `cast_screen_comparison` object.
#' @param var_labels Optional named character vector for display labels.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_screen_comparison <- function(x, var_labels = NULL, ...) {
  check_suggested("ggplot2", "for plotting")
  mem <- x$membership
  methods <- x$methods
  if (!nrow(mem)) cli::cli_abort("No predictor was retained by any method.")

  ord <- order(-as.integer(mem[[x$cpi_method]]), -mem$n_methods)
  mem <- mem[ord, , drop = FALSE]

  raw_disp <- .cast_var_label(mem$variable)
  if (!is.null(var_labels)) {
    hit <- mem$variable %in% names(var_labels)
    raw_disp[hit] <- var_labels[mem$variable[hit]]
  }
  disp <- factor(raw_disp, levels = rev(raw_disp))
  pal <- .cast_category_palette()
  cats <- as.character(.cast_var_category(mem$variable))

  method_labels <- c(
    correlation = "Correlation\nfilter", vif = "Stepwise\nVIF",
    univariate = "Univariate\nscreen", rf = "RF\nimportance",
    cpi = "CPI\n(castSDM)"
  )

  long <- do.call(rbind, lapply(methods, function(m) {
    data.frame(
      display = disp, category = cats, method = m,
      retained = as.logical(mem[[m]]), stringsAsFactors = FALSE
    )
  }))
  long$fill_cat <- ifelse(long$retained, long$category, NA_character_)
  long$fill_cat <- factor(long$fill_cat, levels = names(pal))
  long$method <- factor(long$method, levels = methods,
                        labels = method_labels[methods])

  cnt <- vapply(methods, function(m) sum(mem[[m]]), integer(1))
  names(cnt) <- methods
  sub_txt <- sprintf(
    "Retained \u2014 correlation %d \u00b7 VIF %d \u00b7 univariate %d \u00b7 RF %d \u00b7 CPI %d",
    cnt["correlation"], cnt["vif"], cnt["univariate"], cnt["rf"], cnt["cpi"]
  )

  ggplot2::ggplot(long, ggplot2::aes(
    x = .data$method, y = .data$display, fill = .data$fill_cat
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = pal, name = "Predictor class", na.value = "grey95", drop = TRUE
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(
      title = "Conditional vs conventional predictor screening",
      subtitle = sub_txt, x = NULL, y = NULL,
      caption = paste(
        "Associational filters keep collinear bystanders;",
        "CPI retains conditionally predictive drivers."
      )
    ) +
    theme_cast(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
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
        fill = .cast_map_bg(), color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  large_grid <- nrow(pred) > 2e5
  if (large_grid) {
    p <- p +
      ggplot2::geom_raster(
        data = pred,
        ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$hss_ensemble)
      ) +
      .cast_suitability_scale("fill", name = "Ensemble\nHSS")
  } else {
    p <- p +
      ggplot2::geom_point(
        data = pred,
        ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$hss_ensemble),
        size = 0.4, alpha = 0.85
      ) +
      .cast_suitability_scale("colour", name = "Ensemble\nHSS")
  }
  p <- p +
    ggplot2::labs(
      title = sprintf("Ensemble Habitat Suitability (%s)", x$method),
      subtitle = if (is.finite(x$threshold %||% NA_real_)) sprintf("Threshold = %.3f", x$threshold)
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
  p <- .add_china_outline(p, basemap) + .coord_for_map(basemap)
  .add_china_south_sea_inset(
    p, basemap, data = pred, value_var = "hss_ensemble", raster = large_grid,
    bg_fill = .cast_map_bg(),
    scale = if (large_grid) {
      .cast_suitability_scale("fill", guide = "none")
    } else {
      .cast_suitability_scale("colour", guide = "none")
    }
  )
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

  change_bg <- "#E6E9ED"
  change_colors <- c(
    gain = "#2E9E5B", loss = "#CB4335", stable_present = "#2E86C1",
    stable_absent = change_bg
  )

  p <- ggplot2::ggplot()
  if (basemap != "none") {
    basemap_sf <- load_basemap(basemap)
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf,
        fill = change_bg, color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  large_grid <- nrow(change) > 2e5
  if (large_grid) {
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
  p <- .add_china_outline(p, basemap) + .coord_for_map(basemap)
  .add_china_south_sea_inset(
    p, basemap, data = change, value_var = "change", raster = large_grid,
    bg_fill = change_bg,
    scale = if (large_grid) {
      ggplot2::scale_fill_manual(values = change_colors, guide = "none")
    } else {
      ggplot2::scale_color_manual(values = change_colors, guide = "none")
    }
  )
}


#' Plot Causal Effects / Conditional Importance
#'
#' Forest plot of each predictor's conditional contribution. For a CPI screen
#' this is the (non-negative) conditional predictive impact; for a DML screen it
#' is the orthogonalized partial-linear effect on occurrence per one standard
#' deviation, with confidence intervals. Predictors passing FDR control are
#' highlighted.
#'
#' @param x A `cast_effect` object (from [cast_effect()]).
#' @param var_labels Optional named character vector for display labels.
#' @param top Optional integer. Show only the `top` largest-magnitude effects.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_effect <- function(x, var_labels = NULL, top = NULL, ...) {
  check_suggested("ggplot2", "for plotting")
  eff <- x$effects
  if (!is.null(top) && is.finite(top)) {
    eff <- utils::head(eff[order(-abs(eff$estimate)), , drop = FALSE], as.integer(top))
  }
  eff$sig <- ifelse(eff$selected, "significant", "not significant")
  if (!is.null(var_labels)) {
    eff$display <- ifelse(eff$variable %in% names(var_labels),
                          var_labels[eff$variable], eff$variable)
  } else {
    eff$display <- eff$variable
  }
  eff <- eff[order(eff$estimate), ]
  eff$display <- factor(eff$display, levels = eff$display)

  sig_colors <- c(significant = "#B2182B", `not significant` = "grey70")
  n_sig <- sum(eff$selected, na.rm = TRUE)
  is_cpi <- identical(x$diagnostics$measure, "cpi")

  if (is_cpi) {
    plot_title <- "Conditional predictive impact"
    plot_subtitle <- sprintf(
      "CPI (log-loss knockoff) | %d/%d significant (FDR < %.2g) | %d%% CI",
      n_sig, nrow(eff), x$alpha, round(100 * x$conf_level)
    )
    x_lab <- "Conditional predictive impact"
  } else {
    plot_title <- "Causal effect on occurrence"
    plot_subtitle <- sprintf(
      "Double machine learning | %d/%d significant (FDR < %.2g) | %d%% CI",
      n_sig, nrow(eff), x$alpha, round(100 * x$conf_level)
    )
    x_lab <- "Partial-linear effect (per +1 SD)"
  }

  ggplot2::ggplot(eff, ggplot2::aes(
    x = .data$estimate, y = .data$display, color = .data$sig
  )) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey50", linewidth = 0.4) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high),
      height = 0.25, linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_color_manual(values = sig_colors, name = NULL) +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = x_lab, y = ""
    ) +
    theme_cast(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_line(color = "grey93", linewidth = 0.3)
    )
}


#' Plot Counterfactual What-If Map
#'
#' Diverging map of the per-cell change in habitat suitability under a
#' single-predictor intervention on the current climate.
#'
#' @param x A `cast_counterfactual` object (from [cast_counterfactual()]).
#' @param basemap Character. `"world"`, `"china"`, or `"none"`.
#' @param var_label Optional display label for the intervened predictor.
#' @param title Optional plot title.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_counterfactual <- function(x, basemap = "world", var_label = NULL,
                                     title = NULL, ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("sf", "for geographic mapping")

  pred <- x$predictions
  lab <- var_label %||% x$variable
  title <- title %||% sprintf("What-if: %s + %g %s", lab, x$shift, x$shift_type)
  lim <- max(abs(pred$delta_hss), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1e-6

  p <- ggplot2::ggplot()
  if (basemap != "none") {
    basemap_sf <- load_basemap(basemap)
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf, fill = "grey95", color = "#bdc3c7", linewidth = 0.2
      )
    }
  }
  large_grid <- nrow(pred) > 2e5
  fill_scale <- ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0,
    limits = c(-lim, lim), name = "\u0394 HSS"
  )
  color_scale <- ggplot2::scale_color_gradient2(
    low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0,
    limits = c(-lim, lim), name = "\u0394 HSS"
  )
  if (large_grid) {
    p <- p + ggplot2::geom_raster(
      data = pred,
      ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$delta_hss)
    ) + fill_scale
  } else {
    p <- p + ggplot2::geom_point(
      data = pred,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$delta_hss),
      size = 0.4, alpha = 0.85
    ) + color_scale
  }
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Suitability gain in %d%% of cells (mean \u0394 = %.3f)",
                         round(100 * x$summary$frac_positive),
                         x$summary$mean_delta)
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
  p <- .add_china_outline(p, basemap) + .coord_for_map(basemap)
  .add_china_south_sea_inset(
    p, basemap, data = pred, value_var = "delta_hss", raster = large_grid,
    bg_fill = "grey95",
    scale = if (large_grid) {
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey95",
        high = "#B2182B", midpoint = 0, limits = c(-lim, lim), guide = "none")
    } else {
      ggplot2::scale_color_gradient2(low = "#2166AC", mid = "grey95",
        high = "#B2182B", midpoint = 0, limits = c(-lim, lim), guide = "none")
    }
  )
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

  # Toggle spherical geometry off for the read only; restore the user's
  # global s2 setting on exit so plotting never mutates session state.
  prev_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(prev_s2), add = TRUE)
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

#' Full national extent for China plots
#'
#' China maps always show the complete, standard national frame; a
#' range-restricted species is never cropped to its coloured cells. Seamless
#' colour across the country is achieved instead by filling the whole basemap
#' with each map's null/background colour (see the per-map `fill`), so cells
#' outside the accessible area continue the surface rather than cutting it.
#' @keywords internal
#' @noRd
.coord_for_map <- function(basemap) {
  if (!identical(basemap, "china")) return(ggplot2::coord_sf(expand = FALSE))
  bounds <- .china_main_bounds()
  if (is.null(bounds)) return(ggplot2::coord_sf(expand = FALSE))
  ggplot2::coord_sf(xlim = bounds[c(1, 3)], ylim = bounds[c(2, 4)],
                    expand = FALSE, datum = NA)
}

#' Shared neutral background fill for suitability maps
#'
#' Equal to the low anchor of the suitability ramp, so near-zero-suitability
#' cells inside the accessible area are indistinguishable from the unmodelled
#' country outside it and the map reads as one continuous surface.
#' @keywords internal
#' @noRd
.cast_map_bg <- function() "#EDF2F4"

#' Sequential suitability colour ramp (pale grey -> deep teal-blue)
#' @keywords internal
#' @noRd
.cast_suitability_scale <- function(aesthetic = c("fill", "colour"),
                                    name = "HSS", guide = "colourbar") {
  aesthetic <- match.arg(aesthetic)
  cols <- c("#EDF2F4", "#9FD3DB", "#4FA3BE", "#256E92", "#0B3D5C")
  if (aesthetic == "fill") {
    ggplot2::scale_fill_gradientn(colours = cols, limits = c(0, 1),
                                  name = name, guide = guide)
  } else {
    ggplot2::scale_colour_gradientn(colours = cols, limits = c(0, 1),
                                    name = name, guide = guide)
  }
}

#' Group a raw predictor name into a broad ecological theme
#'
#' Keyword rules over the CHECO26 naming scheme so the variable-selection
#' figure can colour predictors by class without a user-supplied dictionary.
#' @keywords internal
#' @noRd
.cast_var_category <- function(vars) {
  v <- tolower(vars)
  cat <- rep("Other", length(v))
  is_set <- rep(FALSE, length(v))
  assign_cat <- function(pattern, label) {
    hit <- !is_set & grepl(pattern, v)
    cat[hit] <<- label
    is_set[hit] <<- TRUE
  }
  assign_cat("landcover|land_cover|^lu_|_lu_|cropland|urban|forest|igbp|water_prop", "Land use / cover")
  assign_cat("^anthro|landscan|popn|population|human|hfp|footprint|nightlight", "Anthropogenic")
  assign_cat("^soil|_soil|clay|sand|silt|_ph|ocdens|bulkdens|_oc_|cec", "Soil")
  assign_cat("^hydro|river|stream|flow|discharge|topowet|_cwd|wetland|water", "Hydrology")
  assign_cat("elev|slope|aspect|^tri$|terrain|roughness|topo", "Terrain")
  assign_cat("^vege|ndvi|lai|gpp|npp|canopy", "Vegetation")
  assign_cat("bio[0-9]|temp|precip|prec_|aridity|thornthwaite|etccdi|clim|frost|gdd|season", "Climate")
  cat
}

#' Fixed palette for predictor classes
#' @keywords internal
#' @noRd
.cast_category_palette <- function() {
  c(
    "Climate"          = "#4C72B0",
    "Land use / cover" = "#55A868",
    "Soil"             = "#8C6D31",
    "Hydrology"        = "#4BB6C6",
    "Terrain"          = "#C44E52",
    "Vegetation"       = "#8172B3",
    "Anthropogenic"    = "#CCB974",
    "Other"            = "grey70"
  )
}

#' Compact human-readable label for a raw CHECO26 predictor name
#' @keywords internal
#' @noRd
.cast_var_label <- function(vars) {
  lab <- vars
  lab <- gsub("_r[0-9]+cell$", "", lab)
  repl <- c(
    "aridityindexthornthwaite" = "aridity index",
    "maxtempcoldest" = "max temp coldest month",
    "mintempwarmest" = "min temp warmest month",
    "etccdi_cwd" = "consecutive wet days",
    "topowet" = "topographic wetness",
    "landcover_igbp" = "land cover (IGBP)",
    "lu_water_prop" = "surface water proportion",
    "lu_urban" = "urban proportion",
    "lu_cropland" = "cropland proportion"
  )
  for (nm in names(repl)) lab <- gsub(nm, repl[[nm]], lab, fixed = TRUE)
  lab <- gsub("^lu_", "", lab)
  lab <- gsub("^anthro_", "human ", lab)
  lab <- gsub("_", " ", lab)
  lab <- trimws(gsub("\\s+", " ", lab))
  substr(lab, 1, 1) <- toupper(substr(lab, 1, 1))
  lab
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
#' When `data` and `value_var` are given, the inset renders the same data
#' layer as the main map (cropped to the inset extent) so colours match.
#'
#' @keywords internal
#' @noRd
.add_china_south_sea_inset <- function(plot, basemap, data = NULL,
                                       value_var = NULL, scale = NULL,
                                       raster = FALSE, bg_fill = "#F7F8F8") {
  if (!identical(basemap, "china") ||
      !requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) return(plot)

  china <- load_basemap("china")
  dashline <- load_basemap("dashline")
  bounds <- .china_main_bounds()
  if (is.null(china) || is.null(dashline) || is.null(bounds)) return(plot)

  inset <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = china, fill = bg_fill, colour = NA)

  if (!is.null(data) && !is.null(value_var) &&
      all(c("lon", "lat", value_var) %in% names(data))) {
    sub <- data[data$lon >= 104 & data$lon <= 127 &
                  data$lat >= 1 & data$lat <= 27, , drop = FALSE]
    if (nrow(sub) > 0) {
      if (isTRUE(raster)) {
        inset <- inset + ggplot2::geom_raster(
          data = sub,
          ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data[[value_var]])
        )
      } else {
        inset <- inset + ggplot2::geom_point(
          data = sub,
          ggplot2::aes(x = .data$lon, y = .data$lat, color = .data[[value_var]]),
          size = 0.2, alpha = 0.85
        )
      }
      if (!is.null(scale)) inset <- inset + scale
    }
  }

  inset <- inset +
    ggplot2::geom_sf(data = china, fill = NA, colour = "#4E5963", linewidth = 0.20) +
    ggplot2::geom_sf(
      data = dashline, fill = NA, colour = "#4E5963", linewidth = 0.22, linetype = "solid"
    ) +
    ggplot2::coord_sf(xlim = c(105, 126), ylim = c(2, 26), expand = FALSE, datum = NA) +
    ggplot2::theme_void(base_family = getOption("castSDM.font_family", "Arial")) +
    ggplot2::theme(
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
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
