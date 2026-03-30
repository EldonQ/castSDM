# Plot Methods ----------------------------------------------------------------

#' Plot a CAST DAG as a Network Graph
#'
#' Visualizes the learned causal DAG using ggraph. Nodes are colored by
#' causal role, sized by degree, and edges weighted by bootstrap strength.
#'
#' @param x A `cast_dag` object.
#' @param roles Optional [cast_roles] object for node coloring.
#' @param screen Optional [cast_screen] object to mark selected variables.
#' @param var_labels Optional named character vector mapping variable names
#'   to display labels (e.g., `c(bio02 = "Diurnal Temp Range")`).
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_dag <- function(x, roles = NULL, screen = NULL,
                          var_labels = NULL, ...) {
  check_suggested("ggraph", "for DAG plotting")
  check_suggested("igraph", "for DAG plotting")
  check_suggested("ggplot2", "for plotting")

  edges <- x$edges
  if (nrow(edges) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "No edges in DAG", size = 5,
                          fontface = "italic", color = "grey50")
    )
  }

  g <- igraph::graph_from_data_frame(
    edges[, c("from", "to"), drop = FALSE],
    directed = TRUE,
    vertices = data.frame(
      name = unique(c(edges$from, edges$to)),
      stringsAsFactors = FALSE
    )
  )

  node_names <- igraph::V(g)$name
  igraph::E(g)$strength <- edges$strength
  igraph::V(g)$deg <- igraph::degree(g, mode = "all")

  # Causal roles
  if (!is.null(roles)) {
    rmap <- stats::setNames(roles$roles$role, roles$roles$variable)
    igraph::V(g)$role <- ifelse(
      node_names %in% names(rmap), rmap[node_names], "Unscreened"
    )
  } else {
    igraph::V(g)$role <- "Unscreened"
  }

  # CAST selection
  if (!is.null(screen)) {
    igraph::V(g)$selected <- ifelse(
      node_names %in% screen$selected, "Selected", "Dropped"
    )
  } else {
    igraph::V(g)$selected <- "Selected"
  }

  # Labels
  label_fn <- function(v) {
    if (!is.null(var_labels) && v %in% names(var_labels)) {
      var_labels[v]
    } else {
      v
    }
  }
  igraph::V(g)$label <- vapply(node_names, label_fn, character(1))

  role_colors <- c(
    Root = "#2C3E50", Mediator = "#E67E22",
    Terminal = "#27AE60", Unscreened = "grey75"
  )

  set.seed(42)
  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_arc(
      ggplot2::aes(alpha = ggplot2::after_stat(index)),
      arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed"),
      end_cap = ggraph::circle(4, "mm"),
      strength = 0.15, color = "grey40"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(color = .data$role, size = .data$deg,
                   shape = .data$selected),
      alpha = 0.9
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$label),
      repel = TRUE, size = 3.2, fontface = "bold", max.overlaps = 20
    ) +
    ggplot2::scale_color_manual(values = role_colors, name = "Causal\nRole") +
    ggplot2::scale_size_continuous(range = c(3, 9), name = "Degree") +
    ggplot2::scale_shape_manual(
      values = c(Selected = 16, Dropped = 1), name = "CAST"
    ) +
    ggraph::scale_edge_alpha_continuous(
      range = c(0.3, 0.9), name = "Bootstrap\nStrength"
    ) +
    ggplot2::labs(
      title = "Causal DAG",
      subtitle = sprintf(
        "%d edges | strength >= %.2f | R = %d",
        nrow(edges), x$strength_threshold, x$boot_R
      )
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 11),
      plot.subtitle = ggplot2::element_text(
        hjust = 0, color = "grey40", size = 8.5
      ),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.position = "right"
    )
  p
}


#' Plot ATE Forest Plot
#'
#' Displays Average Treatment Effects with 95% confidence intervals.
#' Significant effects are highlighted in red.
#'
#' @param x A `cast_ate` object.
#' @param var_labels Optional named character vector for display labels.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_ate <- function(x, var_labels = NULL, ...) {
  check_suggested("ggplot2", "for plotting")

  est <- x$estimates
  est$ci_lower <- est$coef - 1.96 * est$se
  est$ci_upper <- est$coef + 1.96 * est$se
  est$sig_label <- ifelse(
    est$significant, "Significant (p < 0.05)", "Not significant"
  )

  # Labels
  if (!is.null(var_labels)) {
    est$display <- ifelse(
      est$variable %in% names(var_labels),
      var_labels[est$variable], est$variable
    )
  } else {
    est$display <- est$variable
  }
  est$display <- factor(
    est$display,
    levels = rev(est$display[order(abs(est$coef))])
  )

  n_sig <- sum(est$significant, na.rm = TRUE)
  n_total <- nrow(est)

  ggplot2::ggplot(est, ggplot2::aes(
    x = .data$coef, y = .data$display, color = .data$sig_label
  )) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey60", linewidth = 0.6
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      height = 0.3, linewidth = 0.5, alpha = 0.7
    ) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_color_manual(
      values = c(
        "Significant (p < 0.05)" = "#C0392B",
        "Not significant" = "grey60"
      ),
      name = ""
    ) +
    ggplot2::labs(
      title = "Causal Effects (ATE)",
      subtitle = sprintf(
        "DML %d-fold cross-fitting | %d/%d significant",
        x$K, n_sig, n_total
      ),
      x = "Average Treatment Effect", y = ""
    ) +
    theme_cast() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_line(
        color = "grey93", linewidth = 0.3
      )
    )
}


#' Plot CAST Screening Score Decomposition
#'
#' Stacked bar chart showing the DAG, ATE, and RF importance components
#' of the composite screening score for each variable.
#'
#' @param x A `cast_screen` object.
#' @param var_labels Optional named character vector for display labels.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_screen <- function(x, var_labels = NULL, ...) {
  check_suggested("ggplot2", "for plotting")

  scr <- x$scores
  w <- x$weights

  scr$comp_dag <- scr$score_dag * w["w_dag"]
  scr$comp_ate <- scr$score_ate * w["w_ate"]
  scr$comp_imp <- scr$score_imp * w["w_imp"]
  scr$is_selected <- scr$variable %in% x$selected

  scr <- scr[order(-scr$score_total), ]
  scr$rank <- seq_len(nrow(scr))

  # Labels
  if (!is.null(var_labels)) {
    scr$display <- ifelse(
      scr$variable %in% names(var_labels),
      var_labels[scr$variable], scr$variable
    )
  } else {
    scr$display <- scr$variable
  }
  scr$display <- factor(scr$display, levels = scr$display)

  # Pivot to long
  long <- data.frame(
    display = rep(scr$display, 3),
    component = rep(
      c("DAG topology", "ATE causal", "RF importance"),
      each = nrow(scr)
    ),
    score = c(scr$comp_dag, scr$comp_ate, scr$comp_imp),
    is_selected = rep(scr$is_selected, 3),
    stringsAsFactors = FALSE
  )
  long$component <- factor(long$component, levels = c(
    "DAG topology", "ATE causal", "RF importance"
  ))

  ggplot2::ggplot(long, ggplot2::aes(
    x = .data$display, y = .data$score,
    fill = .data$component, alpha = .data$is_selected
  )) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(
      values = c(
        "DAG topology" = "#3498DB",
        "ATE causal" = "#E74C3C",
        "RF importance" = "#2ECC71"
      ),
      name = ""
    ) +
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 0.95, "FALSE" = 0.35), guide = "none"
    ) +
    ggplot2::labs(
      title = "CAST Screening Scores",
      subtitle = sprintf(
        "Adaptive weights: w_DAG=%.2f  w_ATE=%.2f  w_RF=%.2f",
        w["w_dag"], w["w_ate"], w["w_imp"]
      ),
      x = "", y = "Composite score"
    ) +
    theme_cast() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1, size = 8
      ),
      legend.position = "bottom"
    )
}


#' Plot Spatial Habitat Suitability Map
#'
#' Displays predicted habitat suitability scores on a geographic map with
#' an optional country/province boundary basemap.
#'
#' @param x A `cast_predict` object.
#' @param model Character. Which model's HSS to plot. Default is the first
#'   model.
#' @param basemap Character. Basemap type: `"world"`, `"china"`, or `"none"`.
#'   Default `"world"`.
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

  # Basemap layer
  if (basemap != "none") {
    basemap_sf <- load_basemap(basemap)
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf,
        fill = "#f4f6f7", color = "#bdc3c7", linewidth = 0.2
      )
    }
  }

  # HSS points
  p <- p +
    ggplot2::geom_point(
      data = pred,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data[[hss_col]]),
      size = 0.4, alpha = 0.85
    ) +
    ggplot2::scale_color_viridis_c(
      option = "turbo", name = "HSS", limits = c(0, 1)
    ) +
    ggplot2::labs(title = title) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "right",
      legend.key.width = ggplot2::unit(0.5, "cm"),
      legend.key.height = ggplot2::unit(1.5, "cm")
    )

  # China dashline overlay
  if (basemap == "china") {
    dash_sf <- load_basemap("dashline")
    if (!is.null(dash_sf)) {
      p <- p + ggplot2::geom_sf(
        data = dash_sf, fill = NA, color = "#bdc3c7", linewidth = 0.3
      )
    }
  }

  p
}


#' Plot Spatial CATE Map
#'
#' Displays spatially varying Conditional Average Treatment Effects on a
#' geographic map with a diverging color scale (blue = negative, white =
#' zero, red = positive effect on species presence).
#'
#' @param x A `cast_cate` object.
#' @param variable Character. Which treatment variable to plot. Default is
#'   the first variable. Use `x$variables` to see available options.
#' @param species Character. Optional species name for the title.
#' @param basemap Character. Basemap type: `"world"`, `"china"`, or `"none"`.
#' @param var_labels Optional named character vector mapping variable names
#'   to display labels.
#' @param point_size Numeric. Point size. Default `0.5`.
#' @param legend_position Character. `"bottom"` (horizontal colorbar, like
#'   the reference figures) or `"right"`. Default `"bottom"`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_cate <- function(x, variable = NULL, species = NULL,
                           basemap = "world", var_labels = NULL,
                           point_size = 0.5,
                           legend_position = "bottom", ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("sf", "for geographic mapping")

  variable <- variable %||% x$variables[1]
  df <- x$effects[x$effects$variable == variable, ]

  if (nrow(df) == 0) {
    cli::cli_abort(c(
      "Variable {.val {variable}} not found in CATE effects.",
      i = "Available variables: {.val {x$variables}}."
    ))
  }

  if (!all(c("lon", "lat") %in% names(df))) {
    cli::cli_abort("CATE effects must contain lon/lat for spatial plotting.")
  }

  # Display label for variable
  var_display <- if (!is.null(var_labels) &&
                     variable %in% names(var_labels)) {
    var_labels[variable]
  } else {
    variable
  }

  # Symmetric limits at 98th percentile (like fig6 reference)
  cate_lim <- stats::quantile(abs(df$cate), 0.98, na.rm = TRUE)
  if (cate_lim < 1e-10) cate_lim <- max(abs(df$cate), na.rm = TRUE)

  # Title
  sp_display <- if (!is.null(species)) {
    gsub("_", " ", species)
  } else {
    NULL
  }
  main_title <- sprintf("Spatial CATE: %s", var_display)

  sub_parts <- sprintf(
    "Causal forest (%d trees) | n = %d cells",
    x$n_trees, nrow(df)
  )
  if (!is.null(sp_display)) {
    sub_parts <- paste0(sp_display, " | ", sub_parts)
  }

  p <- ggplot2::ggplot()

  # Basemap
  if (basemap != "none") {
    basemap_sf <- load_basemap(basemap)
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf,
        fill = "#f4f6f7", color = "#bdc3c7", linewidth = 0.2
      )
    }
  }

  # CATE points with diverging RdBu scale
  p <- p +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$cate),
      size = point_size, alpha = 0.85
    ) +
    ggplot2::scale_color_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      limits = c(-cate_lim, cate_lim),
      oob = scales::squish,
      name = "CATE"
    ) +
    ggplot2::labs(title = main_title, subtitle = sub_parts) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5, size = 14
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, color = "grey40", size = 9, face = "italic"
      ),
      legend.position = legend_position
    )

  # Legend styling based on position
  if (legend_position == "bottom") {
    p <- p + ggplot2::guides(
      color = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(8, "cm"),
        barheight = ggplot2::unit(0.4, "cm"),
        title.position = "top", title.hjust = 0.5
      )
    )
  } else {
    p <- p + ggplot2::guides(
      color = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(0.5, "cm"),
        barheight = ggplot2::unit(4, "cm")
      )
    )
  }

  # China dashline overlay
  if (basemap == "china") {
    dash_sf <- load_basemap("dashline")
    if (!is.null(dash_sf)) {
      p <- p + ggplot2::geom_sf(
        data = dash_sf, fill = NA, color = "#bdc3c7", linewidth = 0.3
      )
    }
  }

  p
}


#' Plot Evaluation Metrics Comparison
#'
#' Bar chart comparing AUC and TSS across fitted models.
#'
#' @param x A `cast_eval` object.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_eval <- function(x, ...) {
  check_suggested("ggplot2", "for plotting")

  m <- x$metrics
  model_colors <- c(
    cast = "#E64B35", mlp_ate = "#F39B7F", mlp = "#91D1C2",
    rf = "#4DBBD5", brt = "#3C5488", maxent = "#B09C85"
  )

  m$model <- factor(m$model, levels = names(model_colors))

  ggplot2::ggplot(m, ggplot2::aes(
    x = .data$model, y = .data$auc_mean, fill = .data$model
  )) +
    ggplot2::geom_col(width = 0.6, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$auc_mean)),
      vjust = -0.5, size = 3.5, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = model_colors, guide = "none") +
    ggplot2::labs(
      title = "Model Performance Comparison",
      subtitle = "AUC on test set",
      x = "", y = "AUC"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_cast()
}


#' Plot CAST Pipeline Result (Multi-Panel)
#'
#' Combines DAG, ATE, and screening plots into a 3-panel figure using
#' patchwork.
#'
#' @param x A `cast_result` object.
#' @param var_labels Optional named character vector for display labels.
#' @param ... Ignored.
#'
#' @return A `patchwork` combined plot.
#' @export
plot.cast_result <- function(x, var_labels = NULL, ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("patchwork", "for multi-panel layout")

  p_dag <- plot.cast_dag(x$dag, roles = x$roles, screen = x$screen,
                         var_labels = var_labels)
  p_ate <- plot.cast_ate(x$ate, var_labels = var_labels)
  p_scr <- plot.cast_screen(x$screen, var_labels = var_labels)

  combined <- p_dag | p_ate | p_scr
  combined + patchwork::plot_layout(widths = c(1.2, 1, 1))
}


# ========================================================================
# Internal helpers
# ========================================================================

#' CAST Publication Theme
#' @keywords internal
#' @noRd
theme_cast <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0, size = base_size
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0, color = "grey40", size = base_size - 2
      ),
      legend.background = ggplot2::element_rect(fill = "white", color = NA)
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
    if (is.na(sf::st_crs(basemap)) ||
        sf::st_crs(basemap)$epsg != 4326L) {
      basemap <- sf::st_transform(basemap, 4326)
    }
    basemap <- sf::st_make_valid(basemap)
  }
  basemap
}
