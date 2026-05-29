# Plot Methods ----------------------------------------------------------------

#' Plot a castSDM DAG as a Network Graph
#'
#' Visualizes the learned causal DAG using ggraph. Nodes are colored by
#' Markov Blanket role (parent/child/co_parent/predictive) when a
#' `cast_screen` object is provided, sized by degree, and edges weighted
#' by bootstrap strength. The response node (presence) is highlighted
#' when present.
#'
#' @param x A `cast_dag` object.
#' @param screen Optional `cast_screen` object from [cast_select()] to mark
#'   selected variables and color nodes by MB role.
#' @param var_labels Optional named character vector mapping variable names
#'   to display labels (e.g., `c(bio02 = "Diurnal Temp Range")`).
#' @param species Optional character string. Species name for the title
#'   (underscores are converted to spaces and italicized).
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_dag <- function(x, screen = NULL,
                          var_labels = NULL, species = NULL, ...) {
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
                          family = getOption("castSDM.font_family", "Arial"),
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
  igraph::E(g)$edge_strength <- edges$strength
  igraph::V(g)$deg <- igraph::degree(g, mode = "all")

  # Causal roles from screen (MB-based: parent/child/co_parent/predictive)
  if (!is.null(screen) && !is.null(screen$roles) && nrow(screen$roles) > 0) {
    rmap <- stats::setNames(screen$roles$role, screen$roles$variable)
    igraph::V(g)$role <- ifelse(
      node_names %in% names(rmap), rmap[node_names], "Unscreened"
    )
  } else {
    igraph::V(g)$role <- "Unscreened"
  }

  # Mark response node
  resp_node <- x$response_node
  if (!is.null(resp_node) && resp_node %in% node_names) {
    idx <- which(node_names == resp_node)
    igraph::V(g)$role[idx] <- "Response"
  }

  # Selection status
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
    parent = "#2C3E50", child = "#E67E22",
    co_parent = "#8E44AD", predictive = "#27AE60",
    Response = "#C0392B", Unscreened = "grey75"
  )

  # Title and subtitle
  n_nodes <- length(unique(c(edges$from, edges$to)))
  dag_density <- if (n_nodes > 1) {
    nrow(edges) / (n_nodes * (n_nodes - 1))
  } else {
    0
  }

  if (!is.null(species)) {
    sp_display <- gsub("_", " ", species)
    main_title <- sprintf("Causal DAG -- %s", sp_display)
  } else {
    main_title <- "Causal DAG"
  }

  sm <- x$structure_method %||% "pc"
  br <- x$boot_R
  boot_txt <- if (isTRUE(!is.na(br))) {
    sprintf("bootstrap R = %d", br)
  } else {
    "bootstrap n/a"
  }
  sub_title <- sprintf(
    "%d edges | density = %.3f | %s | method = %s",
    nrow(edges), dag_density, boot_txt, sm
  )

  set.seed(42)
  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_arc(
      ggplot2::aes(alpha = edge_strength),
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
      repel = TRUE, size = 3.2,
      family = getOption("castSDM.font_family", "Arial"),
      fontface = "bold", max.overlaps = 20
    ) +
    ggplot2::scale_color_manual(values = role_colors, name = "Causal\nRole") +
    ggplot2::scale_size_continuous(range = c(3, 9), name = "Degree") +
    ggplot2::scale_shape_manual(
      values = c(Selected = 16, Dropped = 1), name = "Selection"
    ) +
    ggraph::scale_edge_alpha_continuous(
      range = c(0.3, 0.9), name = "Edge\nStrength"
    ) +
    ggplot2::labs(title = main_title, subtitle = sub_title) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 11),
      plot.subtitle = ggplot2::element_text(
        hjust = 0, color = "grey40", size = 8.5
      ),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.position = "right"
    )
  p
}


#' Plot Variable Selection (Markov Blanket + RF Importance)
#'
#' Bar chart showing MB membership (in_mb flag) and RF permutation importance
#' for each variable. Selected variables are highlighted; roles are color-coded.
#'
#' @param x A `cast_screen` object (from [cast_select()]).
#' @param var_labels Optional named character vector for display labels.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_screen <- function(x, var_labels = NULL, ...) {
  check_suggested("ggplot2", "for plotting")

  scr <- x$scores
  scr$is_selected <- scr$variable %in% x$selected

  # Merge role info
  if (!is.null(x$roles) && nrow(x$roles) > 0) {
    scr <- merge(scr, x$roles[, c("variable", "role"), drop = FALSE],
                 by = "variable", all.x = TRUE)
    scr$role[is.na(scr$role)] <- "not selected"
  } else {
    scr$role <- ifelse(scr$is_selected, "selected", "not selected")
  }

  # Sort by importance (new cast_select uses `importance`; keep legacy fallbacks)
  imp_candidates <- c("importance", "rf_importance", "score_total", "imp_norm")
  imp_col <- imp_candidates[imp_candidates %in% names(scr)][1]
  if (is.na(imp_col) || !length(imp_col)) {
    scr$importance_plot <- as.numeric(scr$is_selected)
    imp_col <- "importance_plot"
  }
  scr[[imp_col]] <- suppressWarnings(as.numeric(scr[[imp_col]]))
  scr[[imp_col]][!is.finite(scr[[imp_col]])] <- 0
  scr <- scr[order(-scr[[imp_col]]), ]

  # Labels
  if (!is.null(var_labels)) {
    scr$display <- ifelse(
      scr$variable %in% names(var_labels),
      var_labels[scr$variable], scr$variable
    )
  } else {
    scr$display <- scr$variable
  }
  scr$display <- factor(scr$display, levels = rev(scr$display))

  role_colors <- c(
    parent = "#2C3E50", child = "#E67E22",
    co_parent = "#8E44AD", predictive = "#27AE60",
    selected = "#3498DB", `not selected` = "grey70"
  )

  n_sel <- sum(scr$is_selected)
  n_mb <- if ("in_mb" %in% names(scr)) sum(scr$in_mb, na.rm = TRUE) else NA
  sub_txt <- sprintf(
    "%d / %d variables selected",
    n_sel, nrow(scr)
  )
  if (!is.na(n_mb)) {
    sub_txt <- paste0(sub_txt, sprintf(" | %d in Markov Blanket", n_mb))
  }

  p <- ggplot2::ggplot(scr, ggplot2::aes(
    x = .data[[imp_col]], y = .data$display,
    fill = .data$role, alpha = .data$is_selected
  )) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = role_colors, name = "Role") +
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 0.95, "FALSE" = 0.35), guide = "none"
    ) +
    ggplot2::labs(
      title = "Variable Selection (MB + RF Importance)",
      subtitle = sub_txt,
      x = if (identical(imp_col, "importance_plot")) "Selection indicator" else "RF Permutation Importance",
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
#' Multi-panel bar chart comparing AUC, TSS, and CBI across fitted models.
#'
#' @param x A `cast_eval` object.
#' @param metrics Character vector. Which metrics to show. Default
#'   `c("auc","tss","cbi")`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (faceted).
#' @export
plot.cast_eval <- function(x,
                           metrics = c("auc", "tss", "cbi"),
                           ...) {
  check_suggested("ggplot2", "for plotting")

  m <- x$metrics
  model_colors <- c(
    rf = "#4DBBD5", brt = "#3C5488", maxent = "#B09C85",
    gam = "#00A087", esm = "#E64B35"
  )

  src <- if (isTRUE(x$cv_source)) "Spatial CV" else "Hold-out test"

  # Build long format: metric x model
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
  long$metric <- factor(long$metric,
                        levels = toupper(sub("_mean$", "", present)))
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
    ggplot2::labs(
      title    = "Model Performance Comparison",
      subtitle = src,
      x = "", y = "Score"
    ) +
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
#' @param basemap Character. `"world"`, `"china"`, or `"none"`. Default
#'   `"world"`.
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
    "#E41A1C", "#377EB8", "#4DAF4A",
    "#984EA3", "#FF7F00", "#A65628",
    "#F781BF", "#999999", "#66C2A5", "#FC8D62"
  )

  # -- Right panel: per-fold metric ------------------------------------------
  fd   <- x$fold_metrics
  mcol <- metric
  if (!mcol %in% names(fd)) {
    cli::cli_warn(
      "Metric '{metric}' not in fold_metrics. Using 'auc'."
    )
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
      title    = sprintf("Per-fold %s", toupper(mcol)),
      subtitle = sprintf(
        "%d-fold spatial %s CV", x$k, x$block_method
      ),
      x = "Fold", y = toupper(mcol)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_cast()

  # -- Left panel: fold map (optional) ---------------------------------------
  if (is.null(lon) || is.null(lat)) {
    return(p_metric)
  }

  map_df <- data.frame(
    lon  = lon,
    lat  = lat,
    fold = factor(x$folds)
  )

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
      ggplot2::aes(x = .data$lon, y = .data$lat,
                   color = .data$fold),
      size = 0.8, alpha = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = fold_colors[seq_len(x$k)],
      name = "Fold"
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = "Spatial fold assignment") +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5
      ),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_map + p_metric + patchwork::plot_layout(widths = c(1.4, 1))
  } else {
    p_metric
  }
}


#' Plot castSDM Pipeline Result (Multi-Panel)
#'
#' Combines DAG and variable selection plots into a 2-panel figure using
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

  p_dag <- plot.cast_dag(x$dag, screen = x$screen,
                         var_labels = var_labels)
  p_scr <- plot.cast_screen(x$screen, var_labels = var_labels)

  combined <- p_dag | p_scr
  combined + patchwork::plot_layout(widths = c(1.2, 1))
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

  p <- p +
    ggplot2::geom_point(
      data = pred,
      ggplot2::aes(x = .data$lon, y = .data$lat,
                   color = .data$hss_ensemble),
      size = 0.4, alpha = 0.85
    ) +
    ggplot2::scale_color_viridis_c(
      option = "turbo", name = "Ensemble\nHSS", limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = sprintf("Ensemble Habitat Suitability (%s)", x$method),
      subtitle = sprintf("Threshold = %.3f", x$threshold)
    ) +
    ggplot2::coord_sf(expand = FALSE) +
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


#' Plot Future Climate Projection
#'
#' @param x A `cast_project` object.
#' @param scenario Character. Which scenario to plot. Default is the first.
#' @param basemap Character. `"world"`, `"china"`, or `"none"`.
#' @param ... Ignored.
#' @return A `ggplot` object or patchwork composite.
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

  # Change map: gain / loss / stable
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

  p <- p +
    ggplot2::geom_point(
      data = change,
      ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$change),
      size = 0.4, alpha = 0.8
    ) +
    ggplot2::scale_color_manual(values = change_colors, name = "Range\nChange") +
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
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.position = "right"
    )

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
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0, size = base_size
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0, color = "grey40", size = base_size - 2
      ),
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
    if (is.na(sf::st_crs(basemap)) ||
        sf::st_crs(basemap)$epsg != 4326L) {
      basemap <- sf::st_transform(basemap, 4326)
    }
    basemap <- sf::st_make_valid(basemap)
  }
  basemap
}
