# AI-assisted restoration; scope and verification: inst/AI_ASSISTANCE.md.
# Reuses the established package cartography and the original paper layouts.
.cast_figure_theme <- function() {
  theme_cast(base_size = 7) + ggplot2::theme(
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6.5),
    plot.caption = ggplot2::element_text(size = 6, hjust = 0),
    axis.text = ggplot2::element_text(size = 6),
    axis.title = ggplot2::element_text(size = 7),
    legend.text = ggplot2::element_text(size = 6),
    legend.title = ggplot2::element_text(size = 6.5),
    strip.text = ggplot2::element_text(size = 7),
    plot.tag = ggplot2::element_text(size = 9, face = "bold"),
    legend.position = "bottom")
}

.cast_figure_wrap <- function(x, width = 140) {
  if (is.null(x)) return(NULL)
  ratio <- getOption("castSDM.figure_width_mm", 183)/183
  paste(strwrap(x, max(25, floor(width*ratio))), collapse = "\n")
}

.cast_figure_annotation <- function(title, caption = NULL) {
  patchwork::plot_annotation(title = .cast_figure_wrap(title, 85),
    caption = .cast_figure_wrap(caption),
    tag_levels = "a", theme = .cast_figure_theme())
}

.cast_figure_variable_labels <- function(x) vapply(.cast_var_label(x), function(label)
  paste(strwrap(label, 24), collapse = "\n"), character(1))

.cast_figure_compact <- function(p, discrete = FALSE) {
  p + ggplot2::labs(caption = NULL) + ggplot2::guides(fill = if (discrete)
    ggplot2::guide_legend(ncol = 1, keyheight = grid::unit(3, "mm"), keywidth = grid::unit(3, "mm")) else
    ggplot2::guide_colourbar(barheight = grid::unit(22, "mm"), barwidth = grid::unit(3, "mm"))) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 7.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 6.5),
      legend.text = ggplot2::element_text(size = 6),
      plot.margin = ggplot2::margin(3, 4, 3, 3))
}

.cast_figure_bound <- function(paths, layers, signed = FALSE) {
  values <- unlist(lapply(unique(paths), function(path) {
    r <- terra::rast(path); present <- intersect(layers, names(r))
    if (!length(present)) return(numeric())
    x <- terra::values(r[[present]], mat = FALSE); x[is.finite(x)]
  }), use.names = FALSE)
  if (!length(values)) return(if (signed) c(-.01, .01) else c(0, 1))
  upper <- max(.01, unname(stats::quantile(if (signed) abs(values) else values, .995, names = FALSE)))
  if (signed) c(-upper, upper) else c(0, upper)
}

.cast_figure_map <- function(path, layer, title, basemap, kind = "continuous",
                             limits = NULL, legend = NULL, resolution = .02) {
  title <- paste(strwrap(title, width = 32), collapse = "\n")
  cache <- getOption("castSDM.figure_map_cache")
  info <- file.info(path)
  key <- .cast_digest(list(path, info$size, as.numeric(info$mtime), layer, basemap,
    kind, limits, legend, resolution, getOption("castSDM.font_family")))
  if (is.environment(cache) && identical(cache$key, key))
    return(cache$plot + ggplot2::labs(title = title))
  stack <- terra::rast(path)
  r <- stack[[layer]]
  if (grepl("^context_|^(effect|retained)_agreement_", layer)) {
    variable <- sub("^context_(direction|gap)_|^(effect|retained)_agreement_", "", layer)
    flags <- paste0(c("support_effect_", "support_retained_"), variable)
    if (all(flags %in% names(stack)))
      r <- terra::mask(r, terra::ifel(stack[[flags[1L]]] == 1 & stack[[flags[2L]]] == 1, 1, NA))
  }
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d) <- c("lon", "lat", "value")
  if (!nrow(d)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
    label = "Unassessed: no finite eligible cells", size = 2.5) +
    ggplot2::labs(title = title) + ggplot2::theme_void())
  categorical <- kind %in% c("range", "support", "novelty", "conflict", "context")
  subtitle <- NULL
  if (categorical) {
    levels <- if (kind == "range") c(-1, 0, 1, 2) else if (kind == "context") c(-1, 0, 1, 2, 3) else c(0, 1)
    labels <- switch(kind, context = c("Consistent decrease", "No response", "Consistent increase", "Opposite signs", "One background zero"), range = c("Loss", "Stable absent", "Gain", "Stable present"),
      support = c("Outside support", "Inside support"), novelty = c("Within envelope", "Novel (MESS < 0)"),
      conflict = c("No opposite signs", "Positive/negative conflict"))
    colors <- switch(kind, context = c("#B35806", "#EDF2F4", "#2166AC", "#7B3294", "#999999"), range = c("#CB4335", "#F2F2F2", "#2E9E5B", "#2E86C1"),
      support = c("#F4A582", "#0B3D5C"), novelty = c("#EDF2F4", "#0B3D5C"),
      conflict = c("#EDF2F4", "#B35806"))
    d$value <- factor(d$value, levels = levels, labels = labels)
    scale <- function(guide) ggplot2::scale_fill_manual(values = stats::setNames(colors, labels),
      drop = FALSE, na.translate = FALSE, na.value = NA, name = legend,
      guide = if (guide == "none") "none" else "legend")
  } else {
    signed <- identical(kind, "signed")
    if (is.null(limits)) limits <- if (signed) {
      hi <- max(.01, unname(stats::quantile(abs(d$value), .995, names = FALSE))); c(-hi, hi)
    } else c(0, max(.01, unname(stats::quantile(d$value, .995, names = FALSE))))
    clipped <- sum(d$value < limits[1] | d$value > limits[2])
    subtitle <- if (clipped) sprintf("Display limits [%s, %s]; %s cells colour-clipped", signif(limits[1], 3),
      signif(limits[2], 3), format(clipped, big.mark = ",")) else NULL
    scale <- if (signed) function(guide) ggplot2::scale_fill_gradient2(
      low = "#B35806", mid = "#F2F2F2", high = "#2166AC", midpoint = 0,
      limits = limits, oob = scales::squish, na.value = NA, name = legend,
      guide = if (guide == "none") "none" else "colourbar") else function(guide)
        ggplot2::scale_fill_gradientn(colours = c("#EDF2F4", "#9FD3DB", "#4FA3BE", "#256E92", "#0B3D5C"),
          limits = limits, oob = scales::squish, na.value = NA, name = legend,
          guide = if (guide == "none") "none" else "colourbar")
  }
  if (!is.null(subtitle)) subtitle <- paste(strwrap(subtitle, width = 34), collapse = "\n")
  p <- .cast_map_surface(d, "value", scale, title, subtitle = subtitle,
    legend_name = legend, cell_size = resolution, basemap = basemap, raster = TRUE)
  p <- .cast_figure_compact(p, categorical)
  # Keep a single panel in memory: shared current/common-support atlases reuse
  # its raster, basemap and inset, while each panel keeps its own title.
  if (is.environment(cache)) { cache$key <- key; cache$plot <- p }
  p
}

.cast_figure_scenarios <- function(d) {
  s <- as.character(d$scenario)
  d$ssp <- ifelse(grepl("ssp[0-9]+", s), toupper(sub(".*(ssp[0-9]+).*", "\\1", s)), s)
  d$period <- ifelse(grepl("[0-9]{4}_[0-9]{4}$", s), sub(".*([0-9]{4})_([0-9]{4})$", "\\1-\\2", s), s)
  d$period <- factor(d$period, levels = sort(unique(d$period)))
  d
}

.cast_figure_scenario_colors <- function(groups) {
  colors <- c(SSP126 = "#009E73", SSP245 = "#E69F00", SSP370 = "#D55E00", SSP585 = "#CC79A7")
  other <- setdiff(unique(as.character(groups)), names(colors))
  if (length(other)) colors <- c(colors, stats::setNames(grDevices::hcl.colors(length(other), "Dark 3"), other))
  colors
}

.cast_figure_scenario_label <- function(scenario) {
  if (!grepl("_ssp[0-9]+_[0-9]{4}_[0-9]{4}$", scenario)) return(gsub("_", " ", scenario))
  parts <- regmatches(scenario, regexec("^(.*)_(ssp[0-9]+)_([0-9]{4})_([0-9]{4})$", scenario))[[1L]]
  paste(toupper(gsub("_", "-", sub("^joint_", "", parts[2L]))),
    toupper(parts[3L]), paste(parts[4L], parts[5L], sep = "-"))
}

.cast_figure_trend <- function(data, components = c("climate_main", "lulc_main", "joint_change"),
                               y = "pct_of_current", title = "Thresholded range change", domain = NULL,
                               group_labels = c("First group", "Second group")) {
  components <- intersect(components, unique(data$component))
  d <- data[data$component %in% components, , drop = FALSE]
  if (!is.null(domain) && "domain" %in% names(d)) d <- d[d$domain %in% domain, , drop = FALSE]
  d <- .cast_figure_scenarios(d)
  labels <- c(climate_main = paste(group_labels[1L], "baseline effect"),
    lulc_main = paste(group_labels[2L], "baseline effect"), joint_change = "Future change", interaction = "Interaction",
    climate_conditional = paste(group_labels[1L], "conditional effect"),
    lulc_conditional = paste(group_labels[2L], "conditional effect"),
    climate_average = paste(group_labels[1L], "average effect"),
    lulc_average = paste(group_labels[2L], "average effect"))
  d$component_label <- factor(unname(labels[d$component]), levels = unname(labels[components]))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$period, .data[[y]], colour = .data$ssp, group = .data$ssp)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = .3, colour = "grey65") +
    (if (nlevels(d$period) > 1L) ggplot2::geom_line(linewidth = .5)) + ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~component_label, ncol = if (length(components) <= 3L) 1 else 2) +
    ggplot2::scale_colour_manual(values = .cast_figure_scenario_colors(d$ssp)) +
    ggplot2::labs(title = title, x = NULL, y = if (y == "pct_of_current") "Net range change (%)" else y,
      colour = NULL) + .cast_figure_theme()
  p
}

.cast_figure_performance <- function(data) {
  d <- data$curves
  roc <- ggplot2::ggplot(d, ggplot2::aes(.data$fpr, .data$recall,
    colour = .data$model, linetype = .data$dataset)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = .3) +
    ggplot2::geom_line(linewidth = .45) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(title = "ROC", x = "False positive rate", y = "True positive rate", colour = NULL, linetype = NULL)
  pr <- ggplot2::ggplot(d, ggplot2::aes(.data$recall, .data$precision,
    colour = .data$model, linetype = .data$dataset)) +
    ggplot2::geom_hline(data = unique(d[, c("dataset", "sampled_prevalence")]),
      ggplot2::aes(yintercept = .data$sampled_prevalence, linetype = .data$dataset),
      colour = "grey65", linewidth = .3, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = .45) + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(title = "Precision-recall", x = "Recall", y = "Precision", colour = NULL, linetype = NULL)
  b <- data$reliability
  b <- b[b$n > 0 & b$dataset == "Hold-out", , drop = FALSE]
  calibration <- ggplot2::ggplot(b, ggplot2::aes(.data$mean_score,
    .data$sampled_positive_fraction, colour = .data$model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = .3) +
    ggplot2::geom_line(linewidth = .45) + ggplot2::geom_point(ggplot2::aes(size = .data$n)) +
    ggplot2::scale_colour_brewer(palette = "Dark2") + ggplot2::guides(size = "none", colour = "none") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = "Hold-out reliability\n(sampled responses)", x = "Mean predicted score",
      y = "Sampled positive fraction", colour = NULL)
  caption <- "Grey PR lines mark sampled positive fractions (AP baselines). OOF curves are descriptive; hold-out labels do not choose weights or thresholds. Reliability concerns sampled labels, not ecological occurrence probabilities."
  if ("legacy_clamp_mismatch" %in% data$diagnostics$notes$code) caption <- paste(caption,
    "Legacy validation was unclamped while maps were clamped; historical metrics are preserved.")
  patchwork::wrap_plots(roc, pr, calibration, ncol = 3, guides = "collect") +
    .cast_figure_annotation("Model performance and score interpretation",
      caption) & .cast_figure_theme()
}

.cast_figure_registry <- function(data, basemap) {
  cmp <- data$comparison; report <- data$scenario
  group_labels <- names(cmp$groups)
  if (length(group_labels) != 2L) group_labels <- c("First group", "Second group")
  group_labels <- vapply(group_labels, function(x) if (tolower(x) == "lulc") "LULC" else
    tools::toTitleCase(gsub("_", " ", x)), character(1))
  path <- data$displays[[report]]; att <- data$attribution
  suitability_limits <- local({
    value <- NULL
    function() {
      if (is.null(value)) value <<- if (length(data$displays))
        .cast_figure_bound(data$displays, c("current", "joint")) else c(0, 1)
      value
    }
  })
  paired <- identical(cmp$design, "paired")
  effects <- if (paired) "joint_change" else .cast_factorial_names()
  names_short <- c(current = "Current suitability", joint = "Joint future suitability",
    climate = paste(group_labels[1L], "only suitability"), lulc = paste(group_labels[2L], "only suitability"),
    climate_main = paste(group_labels[1L], "baseline effect"), lulc_main = paste(group_labels[2L], "baseline effect"),
    joint_change = "Joint suitability change", interaction = "Two-group interaction",
    climate_conditional = paste(group_labels[1L], "conditional effect"), lulc_conditional = paste(group_labels[2L], "conditional effect"),
    climate_average = paste(group_labels[1L], "average effect"), lulc_average = paste(group_labels[2L], "average effect"))
  if (paired) names_short[c("joint", "joint_change")] <- c("Future suitability", "Future suitability change")
  registry <- list()
  add <- function(name, build, height, available = TRUE, reason = "Required saved product is unavailable.") {
    registry[[name]] <<- list(build = build, height_mm = height, available = available, reason = reason)
  }
  map <- function(layer, title = unname(names_short[layer]), kind = "continuous", limits = NULL,
                  source = path, legend = "Suitability") .cast_figure_map(source, layer,
    title, basemap, kind, limits, legend, data$display_resolution)
  add("01_model_validation", function() {
    info <- data$folds
    cv <- data$cv
    if ("replicate" %in% names(cv$fold_metrics)) {
      cv$fold_metrics <- cv$fold_metrics[cv$fold_metrics$replicate == 1L, ]
      if (!is.null(cv$repeat_status)) cv$fold_status <- cv$repeat_status$status[cv$repeat_status$replicate == 1L]
    }
    p <- plot(cv, lon = info$lon, lat = info$lat, basemap = basemap, metric = "auc")
    p & .cast_figure_theme()
  }, 100, !is.null(data$cv))
  add("02_habitat_change", function() {
    plots <- list(map("current", limits = suitability_limits()), map("joint", limits = suitability_limits()),
      map("joint_change", kind = "signed", legend = "Future - current"),
      .cast_figure_trend(cmp$binary, group_labels = group_labels))
    patchwork::wrap_plots(plots, ncol = 2) + .cast_figure_annotation(
      "Predicted environmental suitability change",
      paste("Reporting scenario:", report, ". Shared suitability scale. Range trends use the frozen training OOF threshold; continuous maps retain sub-threshold values. Maps are display-only reprojections."))
  }, 155, !is.null(path) && !is.null(cmp$binary))
  add("03_environmental_driver", function() {
    only <- att$summary[att$summary$contrast == "only", ]
    changing <- only$driver[is.finite(only$mean_abs) & only$mean_abs > 0]
    top <- if (length(changing)) changing[which.max(only$mean_abs[match(changing, only$driver)])] else att$drivers[1L]
    limits <- .cast_figure_bound(data$attribution_display,
      paste0(c("effect_", "retained_"), top), signed = TRUE)
    bars <- function(contrast) {
      hi <- max(.01, ceiling(max(abs(att$summary$mean_signed[
        att$summary$contrast %in% c("only", "retained")]), na.rm = TRUE)*100)/100)
      plot(att, contrast = contrast) + .cast_figure_theme() +
        ggplot2::labs(title = if (contrast == "only") "Driver-only future effects" else "Retained-driver contributions",
          caption = NULL) +
        ggplot2::scale_x_continuous(limits = c(-hi, hi), breaks = pretty(c(-hi, hi), 3),
          labels = scales::label_number(accuracy = .001)) +
        ggplot2::scale_y_discrete(labels = function(x) vapply(x, function(label)
          paste(strwrap(label, 27), collapse = "\n"), character(1)))
    }
    maps <- list(map(paste0("effect_", top), paste("Driver-only:", .cast_var_label(top)),
      "signed", limits, source = data$attribution_display, legend = "Suitability change"),
      map(paste0("retained_", top), paste("Retained-driver:", .cast_var_label(top)),
        "signed", limits, source = data$attribution_display, legend = "Suitability change"))
    plots <- c(maps, list(patchwork::wrap_plots(bars("only"), bars("retained"), ncol = 2)))
    patchwork::wrap_plots(plots, ncol = 1, heights = c(1, 1, 1.15)) + .cast_figure_annotation(
      "Environmental explanation under frozen fitted models",
      "Driver-only and retained-driver contrasts are non-additive predictive explanations. Static predictors remain fixed. These contrasts do not identify ecological causal effects.")
  }, 220, !is.null(att) && length(att$drivers) > 0L && !is.null(data$attribution_display),
    if (!is.null(att) && !length(att$drivers)) "No fitted predictor changes in the supplied scenario; per-factor scenario maps are inapplicable." else
      "Driver attribution needs raw scenario predictors or a supplied saved cast_attribution result.")
  add("04_support_and_disagreement", function() {
    patchwork::wrap_plots(map("support", "Training marginal support", "support", legend = "Training range"),
      map("sd_joint_change", "Paired cross-model disagreement", limits = .cast_figure_bound(path, "sd_joint_change"),
        legend = "SD of change"), ncol = 2) + .cast_figure_annotation(
        "Where scenario interpretation needs caution",
        "Support is the marginal training-quantile envelope on raw inputs. SD is calculated after paired within-model changes; it is not a standard error or confidence interval. Uncoloured basemap areas are geographic context, not model predictions.")
  }, 115, !is.null(path))
  if (!paired) add("05_factorial_attribution", function() {
    plots <- list(.cast_figure_trend(cmp$binary, c("climate_main", "lulc_main", "joint_change", "interaction"),
      title = "Factorial range contrasts", group_labels = group_labels), map("climate_main", kind = "signed", legend = "Suitability change"),
      map("lulc_main", kind = "signed", legend = "Suitability change"),
      map("interaction", kind = "signed", legend = "Suitability change"))
    patchwork::wrap_plots(plots, ncol = 2) + .cast_figure_annotation(
      paste("Paired scenario contrasts:", paste(group_labels, collapse = " and ")),
      "Baseline effects compare each partial scenario with current. Interaction = joint - first-only - second-only + current. Model-conditional contrasts, not identified ecological causation. Inspect average and conditional effects in the supplementary atlases.")
  }, 150, !is.null(path) && !is.null(cmp$binary))
  add("06_response_curves", function() {
    p <- plot(data$responses$ale %||% data$responses$pdp, ensemble = FALSE) +
      ggplot2::facet_wrap(~variable, scales = "free", ncol = 3,
        labeller = ggplot2::as_labeller(.cast_figure_variable_labels)) +
      ggplot2::labs(title = if (!is.null(data$responses$ale)) "Accumulated local effects" else "Partial-dependence response curves",
        subtitle = sprintf("All %d fitted predictors; deterministic %d-row training reference", length(data$variables), data$pdp_reference_n),
        caption = .cast_figure_wrap(if (!is.null(data$responses$ale)) "Published ALE (iml): observed local changes are accumulated and centred. Curves describe fitted score responses, not causal effects." else "Individual engine curves are unweighted. Correlated predictors can create unsupported combinations. Model-implied profiles are not causal effects."))
    p + .cast_figure_theme()
  }, max(110, ceiling(length(data$variables)/3)*39), !is.null(data$responses$ale %||% data$responses$pdp))
  add("07_model_performance", function() .cast_figure_performance(data), 100,
    !is.null(data$curves) && !is.null(data$reliability))
  add("08_range_change_class", function() {
    layers <- names(terra::rast(path)); use <- intersect(c("range_class", "range_class_supported",
      "range_class_common_support"), layers)
    titles <- c(range_class = "Full domain", range_class_supported = "Scenario-supported domain",
      range_class_common_support = "All-scenario common support")
    plots <- lapply(use, function(layer) map(layer, titles[[layer]], "range", legend = "Range change"))
    patchwork::wrap_plots(plots, ncol = length(plots)) + .cast_figure_annotation(
      "Thresholded range loss, gain and persistence",
      sprintf("Classification uses native-grid current/future scores >= %.6g before any display reprojection. Supported-domain panels exclude unsupported cells; the full domain retains them.", cmp$threshold))
  }, 105, !is.null(path) && "range_class" %in% names(terra::rast(path)))
  add("09_extrapolation_mess", function() {
    patchwork::wrap_plots(map("novel_current", "Current environmental novelty", "novelty", legend = "MESS < 0", source = data$mess_display),
      map("novel_future", "Future environmental novelty", "novelty", legend = "MESS < 0", source = data$mess_display), ncol = 2) +
      .cast_figure_annotation("Environmental novelty beyond the training envelope",
        "MESS is calculated on original native-grid predictors. MESS < 0 indicates at least one out-of-range predictor. This is distinct from marginal quantile support and multivariate AOA.")
  }, 110, !is.null(data$mess_display),
    "Raw predictor stacks or a saved MESS diagnostic product are required.")
  add("10_future_suitability", function() {
    plots <- lapply(names(data$displays), function(s) map("joint", .cast_figure_scenario_label(s),
      limits = suitability_limits(), source = data$displays[[s]]))
    patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
      .cast_figure_annotation("Future suitability across every requested scenario",
        "All panels use the same colour limits and frozen model weights. Uncoloured basemap areas provide geographic context. Source scenario inputs and time matching are recorded by the workflow.")
  }, max(90, ceiling(length(data$displays)/3)*65), length(data$displays) > 0)
  add("11_variable_selection", function() {
    a <- plot(data$screen, top = 30L) + .cast_figure_theme() +
      ggplot2::scale_y_discrete(labels = .cast_figure_variable_labels)
    a <- a + ggplot2::labs(subtitle = paste(strwrap(a$labels$subtitle, 40), collapse = "\n"), caption = NULL)
    freq <- data$cv$selection_freq
    b <- ggplot2::ggplot(freq, ggplot2::aes(.data$freq,
      stats::reorder(.data$variable, .data$freq))) + ggplot2::geom_col(fill = "#256E92", width = .72) +
      ggplot2::labs(title = "Spatial-fold selection stability", x = "Selection frequency", y = NULL) +
      ggplot2::scale_y_discrete(labels = .cast_figure_variable_labels) +
      .cast_figure_theme()
    patchwork::wrap_plots(a, b, ncol = 2) + .cast_figure_annotation(
      "Variable screening and selection stability", "Overview: the top 30 ranks plus every retained predictor, and all fold-selected predictors. All scored predictors appear in the paginated S_selection catalogue; unscored exclusions remain in the source table. Screening rank is not causal importance.")
  }, 190,
    !is.null(data$screen) && !is.null(data$cv$selection_freq))
  add("12_threshold_sensitivity", function() {
    d <- data$thresholds[data$thresholds$component == "joint_change", ]; d <- .cast_figure_scenarios(d)
    d$threshold_label <- factor(signif(d$threshold, 4))
    ggplot2::ggplot(d, ggplot2::aes(.data$period, .data$pct_of_current, colour = .data$ssp,
      group = interaction(.data$ssp, .data$threshold_label), linetype = .data$threshold_label)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = .3, colour = "grey65") +
      (if (nlevels(d$period) > 1L) ggplot2::geom_line(linewidth = .45)) + ggplot2::geom_point(size = 1.2) +
      ggplot2::facet_wrap(~domain, ncol = 1) +
      ggplot2::scale_colour_manual(values = .cast_figure_scenario_colors(d$ssp)) +
      ggplot2::labs(title = "Threshold and support sensitivity", x = NULL, y = "Net range change (%)",
        colour = NULL, linetype = "Threshold", caption = .cast_figure_wrap("Each threshold/domain has its own explicitly recorded native-grid denominator. No threshold is selected from these future projections.")) + .cast_figure_theme()
  }, 145, !is.null(data$thresholds))
  add("13_model_weights_and_scores", function() {
    d <- data.frame(model = names(cmp$weights), weight = as.numeric(cmp$weights))
    a <- ggplot2::ggplot(d, ggplot2::aes(.data$model, .data$weight)) +
      ggplot2::geom_col(fill = "#256E92", width = .65) + ggplot2::labs(title = "Frozen model weights", x = NULL, y = "Weight")
    e <- data$diagnostics$holdout
    b <- ggplot2::ggplot(e, ggplot2::aes(stats::reorder(.data$model, .data$brier), .data$brier)) +
      ggplot2::geom_col(fill = "#4FA3BE", width = .65) + ggplot2::coord_flip() +
      ggplot2::labs(title = "Hold-out sampled-response Brier score", x = NULL, y = "Brier score")
    patchwork::wrap_plots(a, b, ncol = 2) + .cast_figure_annotation(
      "Ensemble policy and sampled-response loss", "Weights are frozen from training CV. The constant baseline uses training prevalence. Brier score concerns sampled labels, not calibrated ecological occurrence probability.") & .cast_figure_theme()
  }, 100, length(cmp$weights) > 0 && !is.null(data$diagnostics$holdout))
  add("14_fixed_response_profiles", function() plot(data$responses$fixed, ensemble = FALSE) +
    ggplot2::facet_wrap(~variable, scales = "free", ncol = 3,
      labeller = ggplot2::as_labeller(.cast_figure_variable_labels)) +
    ggplot2::labs(title = "Fixed response profiles", caption = .cast_figure_wrap("Each variable varies across its reference range; all others are fixed at training medians. These are not partial-dependence averages or causal effects.")) +
    .cast_figure_theme(), max(110, ceiling(length(data$variables)/3)*39), !is.null(data$responses$fixed))
  add("15_area_of_applicability", function() {
    a <- data$cv$aoa
    caption <- "The dashed line is the recorded training-side AOA threshold. OOF diagnostics refer to repeat 1 and do not by themselves explain a weak spatial fold."
    if (length(data$training_version) == 1L && !is.na(data$training_version) &&
        utils::compareVersion(data$training_version, "0.27.0") < 0)
      caption <- paste(caption, "Stored legacy AOA; version 0.27 corrected nearest-neighbour distances. Recompute before scientific interpretation.")
    d <- rbind(data.frame(di = a$di_train, set = "Training"), data.frame(di = a$di_oof, set = "Training OOF"))
    d <- d[is.finite(d$di), ]
    ggplot2::ggplot(d, ggplot2::aes(.data$di, colour = .data$set)) +
      ggplot2::stat_ecdf(linewidth = .5) + ggplot2::geom_vline(xintercept = a$threshold, linetype = 2, colour = "#B35806") +
      ggplot2::labs(title = "Area-of-applicability calibration", x = "Dissimilarity index", y = "Empirical cumulative fraction",
        colour = NULL, caption = .cast_figure_wrap(caption)) + .cast_figure_theme()
  }, 95, !is.null(data$cv$aoa) && is.finite(data$cv$aoa$threshold %||% NA_real_))
  if (!is.null(data$screen) && is.data.frame(data$screen$scores)) {
    scored <- plot(data$screen, top = NULL)$data
    pages <- split(seq_len(nrow(scored)), ceiling(seq_len(nrow(scored))/45L))
    for (page in seq_along(pages)) local({
      number <- page; rows <- pages[[number]]
      add(sprintf("S_selection_%02d", number), function() {
        screen <- data$screen
        screen$scores <- screen$scores[match(scored$variable[rows], screen$scores$variable), , drop = FALSE]
        plot(screen, top = NULL) + .cast_figure_theme() +
          ggplot2::scale_y_discrete(labels = .cast_figure_variable_labels) +
          ggplot2::labs(title = sprintf("Complete screening ranks: page %d of %d", number, length(pages)),
            subtitle = sprintf("Ranks %d-%d of %d scored predictors", min(rows), max(rows), nrow(scored)),
            caption = .cast_figure_wrap("Every scored predictor is retained across these pages. Scores describe fitted selector components (historical screens retain their original score definitions); original values and qualification exclusions remain in screening.csv. Display ordering does not determine the selected set or identify causal importance."))
      }, max(100, length(rows)*3.5+45))
    })
  }
  if (length(data$displays)) {
    layers <- names(terra::rast(path))
    atlas_layers <- setdiff(layers, c("current", "joint",
      "novel_current", "novel_future", "mess_current", "mess_future"))
    for (layer in atlas_layers) local({
      v <- layer
      available <- data$displays[vapply(data$displays, function(p) v %in% names(terra::rast(p)), logical(1))]
      if (identical(cmp$domain, "common") && (startsWith(v, "current_") || v == "common_support"))
        available[] <- available[1L]
      kind <- if (startsWith(v, "range_class")) "range" else if (v %in% effects || startsWith(v, "change_")) "signed" else if (v %in% c("support", "common_support")) "support" else
        if (v == "sign_conflict_joint") "conflict" else "continuous"
      add(paste0("S_", v), function() {
        limits <- if (kind %in% c("continuous", "signed")) .cast_figure_bound(available, v, kind == "signed") else NULL
        plots <- lapply(names(available), function(s) map(v, .cast_figure_scenario_label(s), kind,
          limits, source = available[[s]], legend = gsub("_", " ", v)))
        patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
          .cast_figure_annotation(gsub("_", " ", v),
            "Every requested scenario is retained; shared display limits. Native values and domain denominators are preserved in source tables. SD and sign agreement describe model disagreement, not statistical confidence.")
      }, max(90, ceiling(length(available)/3)*65),
      available = !startsWith(v, "sd_") || length(cmp$weights) > 1L,
      reason = "Cross-model SD is unassessed with fewer than two positive-weight models.")
    })
  }
  if (!is.null(data$attribution_display)) for (driver in att$drivers) local({
    v <- driver
    add(paste0("D_", v), function() {
      layers <- paste0(c("effect_", "retained_"), v)
      limits <- .cast_figure_bound(data$attribution_display, layers, TRUE)
      plots <- lapply(seq_along(layers), function(i) map(layers[i],
        paste(c("Driver-only:", "Retained-driver:")[i], .cast_var_label(v)), "signed", limits,
        source = data$attribution_display, legend = "Suitability change"))
      patchwork::wrap_plots(plots, ncol = 2) + .cast_figure_annotation(.cast_var_label(v),
        "Frozen-model single-driver scenario contrasts. No refit or probability recalibration; static inputs remain current.")
    }, 105)
    needed <- paste0(c("context_direction_", "context_gap_", "effect_agreement_", "retained_agreement_"), v)
    add(paste0("E_", v), function() {
      gap <- .cast_figure_bound(data$attribution_display, needed[2L], TRUE)
      plots <- list(
        map(needed[1L], "Response direction in both backgrounds", "context",
          source = data$attribution_display, legend = "Background consistency"),
        map(needed[2L], "Difference between background responses", "signed", gap,
          source = data$attribution_display, legend = "Response difference"),
        map(needed[3L], "Engine agreement: current background", limits = c(0, 1),
          source = data$attribution_display, legend = "Agreement fraction"),
        map(needed[4L], "Engine agreement: future background", limits = c(0, 1),
          source = data$attribution_display, legend = "Agreement fraction"))
      patchwork::wrap_plots(plots, ncol = 2) + .cast_figure_annotation(
        paste("Environmental response consistency:", .cast_var_label(v)),
        "The same factor is compared in current and future backgrounds. Both-background range support is required for display. Sign reversals indicate background dependence; engine agreement is descriptive, not confidence or causal identification. Counts and denominators are in response_consistency.csv.")
    }, 165, all(needed %in% names(terra::rast(data$attribution_display))),
      "Background consistency requires attribution diagnostics = TRUE.")
  })
  registry
}
