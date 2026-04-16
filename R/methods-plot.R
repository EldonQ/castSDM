# Plot Methods ----------------------------------------------------------------

#' Plot a CAST DAG as a Network Graph
#'
#' Visualizes the learned causal DAG using ggraph. Nodes are colored by
#' causal role, sized by degree, and edges weighted by bootstrap strength.
#' Follows the single-species showcase style (Fig 2) with edge alpha mapped
#' to actual bootstrap strength.
#'
#' @param x A `cast_dag` object.
#' @param roles Optional [cast_roles] object for node coloring.
#' @param screen Optional [cast_screen] object to mark selected variables.
#' @param var_labels Optional named character vector mapping variable names
#'   to display labels (e.g., `c(bio02 = "Diurnal Temp Range")`).
#' @param species Optional character string. Species name for the title
#'   (underscores are converted to spaces and italicized).
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_dag <- function(x, roles = NULL, screen = NULL,
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

  # Title and subtitle
  n_nodes <- length(unique(c(edges$from, edges$to)))
  dag_density <- if (n_nodes > 1) {
    nrow(edges) / (n_nodes * (n_nodes - 1))
  } else {
    0
  }

  if (!is.null(species)) {
    sp_display <- gsub("_", " ", species)
    main_title <- sprintf("Causal DAG \u2014 %s", sp_display)
  } else {
    main_title <- "Causal DAG"
  }

  sm <- x$structure_method %||% "bootstrap_hc"
  br <- x$boot_R
  boot_txt <- if (isTRUE(!is.na(br))) {
    sprintf("HC bootstrap R = %d", br)
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
    ggplot2::labs(title = main_title, subtitle = sub_title) +
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
#' @param type Character. `"ate"` (default) for coefficient forest plot,
#'   `"evalue"` for E-value sensitivity plot (requires [cast_evalue()] to
#'   have been called first).
#' @return A `ggplot` object.
#' @export
plot.cast_ate <- function(x, var_labels = NULL, type = c("ate", "evalue"),
                          ...) {
  check_suggested("ggplot2", "for plotting")
  type <- match.arg(type)

  if (type == "evalue") return(plot_ate_evalue(x, var_labels))

  est <- x$estimates
  est$ci_lower <- est$coef - 1.96 * est$se
  est$ci_upper <- est$coef + 1.96 * est$se
  p_method <- x$p_adjust %||% "bonferroni"
  sig_label_text <- paste0("Significant (", p_method, " p < ", x$alpha, ")")
  nonsig_label_text <- "Not significant"
  est$sig_label <- ifelse(
    est$significant, sig_label_text, nonsig_label_text
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
      values = stats::setNames(
        c("#C0392B", "grey60"),
        c(sig_label_text, nonsig_label_text)
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
#' @param hss_predict Optional [cast_predict] object. When provided together
#'   with `hss_model`, CATE values are set to `NA` wherever the chosen model's
#'   habitat suitability (`HSS_<model>`) is below `hss_threshold` (spatial
#'   masking by suitability).
#' @param hss_model Character. Which model's HSS column to use for masking.
#'   Default `"cast"`.
#' @param hss_threshold Numeric in `[0, 1]`. Default `0.1`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_cate <- function(x, variable = NULL, species = NULL,
                           basemap = "world", var_labels = NULL,
                           point_size = 0.5,
                           legend_position = "bottom",
                           hss_predict = NULL,
                           hss_model = "cast",
                           hss_threshold = 0.1,
                           ...) {
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

  if (!is.null(hss_predict)) {
    df <- .cate_mask_by_hss(df, hss_predict, hss_model, hss_threshold)
  }

  # Display label for variable
  var_display <- if (!is.null(var_labels) &&
                     variable %in% names(var_labels)) {
    var_labels[variable]
  } else {
    variable
  }

  # Symmetric limits at 98th percentile (like fig6 reference)
  finite_cate <- df$cate[is.finite(df$cate)]
  if (length(finite_cate) == 0) {
    cli::cli_warn("No finite CATE values after masking for {.val {variable}}. Returning empty plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = paste("CATE:", variable, "(no data after HSS mask)")))
  }
  cate_lim <- as.numeric(stats::quantile(abs(finite_cate), 0.98, na.rm = TRUE))
  if (!is.finite(cate_lim) || cate_lim < 1e-10) {
    cate_lim <- max(abs(finite_cate), na.rm = TRUE)
  }
  if (!is.finite(cate_lim) || cate_lim < 1e-10) cate_lim <- 0.05

  # Title
  sp_display <- if (!is.null(species)) {
    gsub("_", " ", species)
  } else {
    NULL
  }
  main_title <- sprintf("Spatial CATE: %s", var_display)

  n_vis <- sum(!is.na(df$cate))
  sub_parts <- sprintf(
    "Causal forest (%d trees) | n = %d cells displayed",
    x$n_trees, n_vis
  )
  if (!is.null(hss_predict)) {
    sub_parts <- paste0(
      sub_parts,
      sprintf(" | HSS_%s >= %.2f", hss_model, hss_threshold)
    )
  }
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
      size = point_size, alpha = 0.85, na.rm = TRUE
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
#' Multi-panel bar chart comparing AUC, TSS, CBI, SEDI, Kappa, and PRAUC
#' across fitted models.
#'
#' @param x A `cast_eval` object.
#' @param metrics Character vector. Which metrics to show. Default shows all
#'   six: `c("auc","tss","cbi","sedi","kappa","prauc")`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (faceted).
#' @export
plot.cast_eval <- function(x,
                           metrics = c("auc", "tss", "cbi",
                                       "sedi", "kappa", "prauc"),
                           ...) {
  check_suggested("ggplot2", "for plotting")

  m <- x$metrics
  model_colors <- c(
    cast = "#E64B35", mlp_ate = "#F39B7F", mlp = "#91D1C2",
    rf = "#4DBBD5", brt = "#3C5488", maxent = "#B09C85"
  )

  src <- if (isTRUE(x$cv_source)) "Spatial CV" else "Hold-out test"

  # Build long format: metric x model
  mean_cols <- paste0(metrics, "_mean")
  present   <- intersect(mean_cols, names(m))
  if (length(present) == 0) {
    # Fallback: old-style object with only auc_mean / tss_mean
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
      vjust = -0.4, size = 2.8, fontface = "bold"
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
    cast = "#E64B35", mlp_ate = "#F39B7F", mlp = "#91D1C2",
    rf = "#4DBBD5", brt = "#3C5488", maxent = "#B09C85"
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
      )
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_map + p_metric + patchwork::plot_layout(widths = c(1.4, 1))
  } else {
    p_metric
  }
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

#' Mask CATE rows by HSS threshold (match lon/lat to cast_predict grid)
#' @keywords internal
#' @noRd
.cate_mask_by_hss <- function(df, hss_predict, hss_model, hss_threshold) {
  if (!inherits(hss_predict, "cast_predict")) {
    cli::cli_abort("{.arg hss_predict} must be a {.cls cast_predict} object.")
  }
  pred <- hss_predict$predictions
  hss_col <- paste0("HSS_", hss_model)
  if (!hss_col %in% names(pred)) {
    cli::cli_abort(c(
      "Model {.val {hss_model}} not found for HSS masking.",
      i = "Available: {.val {hss_predict$models}}."
    ))
  }
  key_df <- paste0(signif(df$lon, 8), "_", signif(df$lat, 8))
  key_pr <- paste0(signif(pred$lon, 8), "_", signif(pred$lat, 8))
  mi <- match(key_df, key_pr)
  hss_vals <- pred[[hss_col]][mi]
  low <- is.na(hss_vals) | (as.numeric(hss_vals) < hss_threshold)
  df$cate[low] <- NA_real_
  df
}


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


# -- E-value sensitivity plot (internal) --------------------------------------
#' @noRd
plot_ate_evalue <- function(x, var_labels = NULL) {
  est <- x$estimates
  if (!"evalue_point" %in% names(est)) {
    cli::cli_abort("No E-value columns found. Run {.fn cast_evalue} first.")
  }

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
    levels = rev(est$display[order(est$evalue_point)])
  )

  ggplot2::ggplot(est, ggplot2::aes(
    x = .data$evalue_point, y = .data$display
  )) +
    ggplot2::geom_vline(
      xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.6
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$evalue_ci, xend = .data$evalue_point,
        y = .data$display, yend = .data$display
      ),
      linewidth = 0.8, color = "grey50"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$significant), size = 3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$evalue_ci), shape = 1, size = 2.5,
      color = "grey40"
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#C0392B", "FALSE" = "grey60"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = ""
    ) +
    ggplot2::labs(
      title = "E-value Sensitivity Analysis",
      subtitle = "Filled = point estimate; Open = CI bound",
      x = "E-value (minimum confounding strength)", y = ""
    ) +
    theme_cast() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_line(
        color = "grey93", linewidth = 0.3
      )
    )
}
