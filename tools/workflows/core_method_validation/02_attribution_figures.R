# Nature-size quantitative figures; source data are independent dataset metrics.
validation_theme <- function() ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = .3),
    axis.ticks = ggplot2::element_line(linewidth = .3),
    axis.text = ggplot2::element_text(size = 6), strip.text = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.caption = ggplot2::element_text(size = 6, hjust = 0),
    plot.tag = ggplot2::element_text(face = "bold", size = 8), legend.position = "bottom")
export_validation <- function(p, stem, height_mm = 135) {
  width <- 183 / 25.4; height <- height_mm / 25.4
  p <- p & validation_theme()
  svglite::svglite(paste0(stem, ".svg"), width = width, height = height)
  print(p); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "Arial")
  print(p); grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width, height = height,
    units = "in", res = 600, compression = "lzw")
  print(p); grDevices::dev.off()
  devices <- list(png = ragg::agg_png)
  for (ext in names(devices)) {
    args <- list(filename = paste0(stem, ".", ext), plot = p,
      width = width, height = height, units = "in", dpi = 600,
      bg = "white", device = devices[[ext]], limitsize = FALSE)
    if (ext == "pdf") args$family <- "Arial"
    if (ext == "tiff") args$compression <- "lzw"
    do.call(ggplot2::ggsave, args)
  }
  ggplot2::ggsave(paste0(stem, "_preview.png"), p, width = width, height = height,
    units = "in", dpi = 300, bg = "white", device = ragg::agg_png)
}
make_attribution_figures <- function(output) {
  for (p in c("ggplot2", "patchwork", "ragg", "svglite")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing ", p)
  library(patchwork)
  metrics <- read.csv(file.path(output, "advantage_metrics.csv"))
  reliability <- read.csv(file.path(output, "reliability_metrics.csv"))
  shards <- read.csv(file.path(output, "scenario_runs.csv"))$scenario_run
  fig <- file.path(output, "figures"); dir.create(fig, showWarnings = FALSE)
  scen_levels <- c("linear", "nonlinear", "joint", "null", "correlated", "extrapolation")
  labels <- c(linear = "Monotonic", nonlinear = "Nonlinear", joint = "Joint response",
    null = "Null driver", correlated = "Correlated", extrapolation = "Outside range")
  colours <- c(local_scenario = "#2166AC", global_response = "#888888", mean_shift = "#B35806")
  method_labels <- c(local_scenario = "Local future values", global_response = "Global PDP contrast", mean_shift = "Uniform mean shift")
  d <- metrics[metrics$domain == "all" & metrics$method %in% names(colours), ]
  d$scenario <- factor(d$scenario, levels = scen_levels)
  d$engine <- factor(d$engine, levels = c("rf", "gam"), labels = c("Random forest", "Additive GAM"))
  d$method <- factor(d$method, levels = names(colours))
  utils::write.csv(d, file.path(fig, "figure2_source.csv"), row.names = FALSE)
  # Paired error differences retain failed advantages rather than hiding them.
  key <- function(z) paste(z$scenario, z$replicate, z$engine)
  local <- d[d$method == "local_scenario", ]
  differences <- do.call(rbind, lapply(c("global_response", "mean_shift"), function(m) {
    other <- d[d$method == m, ]; idx <- match(key(local), key(other))
    data.frame(local[, c("scenario", "replicate", "engine")], baseline = m,
      improvement = other$rmse[idx] - local$rmse)
  }))
  utils::write.csv(differences, file.path(fig, "paired_advantage_source.csv"), row.names = FALSE)
  p1 <- ggplot2::ggplot(d, ggplot2::aes(.data$scenario, .data$rmse, colour = .data$method)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .15, dodge.width = .65, seed = 42),
      alpha = .3, size = .6) +
    ggplot2::stat_summary(fun = mean, geom = "point", position = ggplot2::position_dodge(.65), size = 1.8) +
    ggplot2::facet_wrap(~engine, ncol = 1) + ggplot2::scale_colour_manual(values = colours, labels = method_labels, name = NULL) +
    ggplot2::scale_x_discrete(labels = labels) + ggplot2::labs(x = NULL, y = "Spatial effect RMSE",
      title = "Accuracy of environmental effect maps")
  p2 <- ggplot2::ggplot(differences, ggplot2::aes(.data$scenario, .data$improvement, colour = .data$baseline)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = .3, linetype = 2) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .15, dodge.width = .5, seed = 42),
      alpha = .4, size = .6) +
    ggplot2::stat_summary(fun = mean, geom = "point", position = ggplot2::position_dodge(.5), size = 1.8) +
    ggplot2::facet_wrap(~engine, ncol = 1) +
    ggplot2::scale_colour_manual(values = colours[c("global_response", "mean_shift")],
      labels = method_labels[c("global_response", "mean_shift")], name = "Comparison baseline") +
    ggplot2::scale_x_discrete(labels = labels) + ggplot2::labs(x = NULL,
      y = "Baseline RMSE - local RMSE", title = "Paired advantage; negative values show failure")
  p <- p1 / p2 + patchwork::plot_annotation(tag_levels = "a",
    caption = "Small points: independent simulated datasets; large points: arithmetic means.\nSame fitted model and driver target for all baselines. No inference from grid-cell replication.")
  export_validation(p, file.path(fig, "figure2_method_advantage"), 180)
  # First replicate spatial illustration fixed by protocol.
  maps <- list()
  for (s in c("nonlinear", "joint", "extrapolation")) {
    candidates <- file.path(shards, "maps", paste0(s, "_gam.csv"))
    file <- candidates[file.exists(candidates)][1]
    md <- read.csv(file)
    for (method in c("truth", "local_scenario", "global_response", "error")) {
      maps[[length(maps) + 1L]] <- data.frame(md[, c("lon", "lat")], scenario = s,
        method, value = md[[method]])
    }
  }
  map_data <- do.call(rbind, maps)
  map_data$scenario <- factor(map_data$scenario, levels = c("nonlinear", "joint", "extrapolation"))
  map_data$method <- factor(map_data$method, levels = c("truth", "local_scenario", "global_response", "error"),
    labels = c("Known response", "Local future values", "Global PDP", "Local - truth"))
  lim <- max(abs(map_data$value), na.rm = TRUE)
  utils::write.csv(map_data, file.path(fig, "figure1_source.csv"), row.names = FALSE)
  p <- ggplot2::ggplot(map_data, ggplot2::aes(.data$lon, .data$lat, fill = .data$value)) +
    ggplot2::geom_tile(width = .05, height = .05) +
    ggplot2::facet_grid(scenario ~ method, labeller = ggplot2::labeller(scenario = labels)) +
    ggplot2::scale_fill_gradient2(low = "#B35806", mid = "#F7F7F7", high = "#2166AC",
      midpoint = 0, limits = c(-lim, lim), name = "Suitability\nchange") +
    ggplot2::coord_equal() + ggplot2::labs(x = "Synthetic spatial coordinate", y = NULL,
      title = "Spatial attribution reveals both accuracy and model limitations",
      caption = "Prespecified first dataset; additive GAM. Same colour scale across all panels.\nJoint response and extrapolation are included to expose limitations; these are simulations, not species observations.")
  export_validation(p, file.path(fig, "figure1_known_response_maps"), 142)
  reliability$scenario <- factor(reliability$scenario, levels = scen_levels)
  reliability$type <- factor(reliability$type, levels = c("resampling", "forest_seed", "model_family"),
    labels = c("New observations", "New forest seed", "RF versus GAM"))
  utils::write.csv(reliability, file.path(fig, "figure3_source.csv"), row.names = FALSE)
  p1 <- ggplot2::ggplot(reliability, ggplot2::aes(.data$scenario, .data$map_correlation, colour = .data$type)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 42), size = .7, alpha = .4) +
    ggplot2::stat_summary(fun = mean, geom = "point", position = ggplot2::position_dodge(.6), size = 1.8) +
    ggplot2::scale_colour_manual(values = c("#2166AC", "#4A9A8B", "#B35806"), name = NULL) +
    ggplot2::scale_x_discrete(labels = labels) + ggplot2::labs(x = NULL, y = "Effect-map correlation",
      title = "Repeatability and model dependence")
  p2 <- ggplot2::ggplot(reliability, ggplot2::aes(.data$scenario, .data$map_rmse, colour = .data$type)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 42), size = .7, alpha = .4) +
    ggplot2::stat_summary(fun = mean, geom = "point", position = ggplot2::position_dodge(.6), size = 1.8) +
    ggplot2::scale_colour_manual(values = c("#2166AC", "#4A9A8B", "#B35806"), name = NULL) +
    ggplot2::scale_x_discrete(labels = labels) + ggplot2::labs(x = NULL, y = "Repeat effect-map RMSE")
  p <- p1 / p2 + patchwork::plot_annotation(tag_levels = "a",
    caption = "Each point is one independent dataset pair. Correlation alone can hide magnitude errors.\nModel-family comparison measures disagreement; a high correlation does not establish ecological causality.")
  export_validation(p, file.path(fig, "figure3_reliability"), 140)
  # Range flags are a marginal diagnostic; show their errors, not just coverage.
  support <- metrics[metrics$method == "local_scenario" &
    metrics$domain %in% c("within_range", "outside_range") & metrics$n > 0, ]
  support$scenario <- factor(support$scenario, levels = scen_levels)
  support$engine <- factor(support$engine, levels = c("rf", "gam"),
    labels = c("Random forest", "Additive GAM"))
  support$domain <- factor(support$domain, levels = c("within_range", "outside_range"),
    labels = c("Inside marginal range", "Outside marginal range"))
  utils::write.csv(support, file.path(fig, "figure4_source.csv"), row.names = FALSE)
  p <- ggplot2::ggplot(support,
    ggplot2::aes(.data$scenario, .data$rmse, colour = .data$domain)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(
      jitter.width = .15, dodge.width = .5, seed = 42), size = .7, alpha = .4) +
    ggplot2::stat_summary(fun = mean, geom = "point",
      position = ggplot2::position_dodge(.5), size = 1.8) +
    ggplot2::facet_wrap(~engine, ncol = 1) +
    ggplot2::scale_colour_manual(values = c("#2166AC", "#B35806"), name = NULL) +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::labs(x = NULL, y = "Spatial effect RMSE",
      title = "Training-range flags and effect-map error",
      caption = paste0("Small points: independent datasets; large points: means. Empty domains have no points.\n",
        "Range: each predictor's training 1st-99th percentiles in both compared environments.\n",
        "Marginal range flags do not detect unsupported joint combinations or guarantee accuracy."))
  export_validation(p, file.path(fig, "figure4_support_error"), 125)
  writeLines(c("All quantitative panels trace to source CSVs. All datasets and failure scenarios retained.",
    "Arial; 183 mm width; 600 dpi PNG/TIFF; PDF/SVG editable labels. R-only exports.",
    "Simulation datasets are independent units, pixels are not replicates. Means shown without confidence intervals.",
    "Written exports await final-size visual QA; all map values and scales are unmodified."), file.path(fig, "figure_QA.txt"))
  invisible(fig)
}
