# Four result figures: validation, habitat change, environmental driver, support.
options(castSDM.font_family = "Arial")
paper_theme <- function() {
  cast_theme <- tryCatch(getFromNamespace("theme_cast", "castSDM"),
                         error = function(e) NULL)
  if (is.function(cast_theme)) return(cast_theme(base_size = 7))
  ggplot2::theme(text = ggplot2::element_text(family = "Arial", size = 7),
    plot.title = ggplot2::element_text(family = "Arial", size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(family = "Arial", size = 6),
    plot.caption = ggplot2::element_text(family = "Arial", size = 6, hjust = 0),
    axis.text = ggplot2::element_text(family = "Arial", size = 6),
    axis.title = ggplot2::element_text(family = "Arial", size = 7),
    legend.text = ggplot2::element_text(family = "Arial", size = 6),
    legend.title = ggplot2::element_text(family = "Arial", size = 6),
    strip.text = ggplot2::element_text(family = "Arial", size = 7),
    plot.tag = ggplot2::element_text(family = "Arial", size = 9, face = "bold"),
    legend.position = "bottom")
}

figure_paths <- function(path) {
  paste0(tools::file_path_sans_ext(path), c(".png", ".pdf", ".svg", ".tiff", "_preview.png"))
}

save_plot <- function(plot, path, width = 9, height = 6, apply_theme = TRUE) {
  paths <- figure_paths(path)
  height <- 183 / 25.4 * height / width
  width <- 183 / 25.4
  wrap_label <- function(x, n) {
    if (is.null(x)) return(NULL)
    paste(strwrap(gsub("_", " ", x, fixed = TRUE), width = n), collapse = "\n")
  }
  if (inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
    plot <- plot + ggplot2::labs(title = wrap_label(plot$labels$title, 70),
      subtitle = wrap_label(plot$labels$subtitle, 95),
      caption = wrap_label(plot$labels$caption, 130))
  } else if (inherits(plot, "patchwork")) {
    annotation <- plot$patches$annotation
    if (!is.null(annotation)) {
      for (name in c("title", "subtitle", "caption")) {
        n <- if (name == "caption") 130 else 70
        annotation[[name]] <- wrap_label(annotation[[name]], n)
      }
      plot$patches$annotation <- annotation
    }
  }
  if (apply_theme) {
    plot <- if (inherits(plot, "patchwork")) plot & paper_theme() else plot + paper_theme()
  }
  devices <- list(ragg::agg_png, function(...) grDevices::cairo_pdf(...),
                  function(...) svglite::svglite(...), ragg::agg_tiff, ragg::agg_png)
  rows <- lapply(seq_along(paths), function(i) {
    notes <- character()
    status <- tryCatch(withCallingHandlers({
      args <- list(filename = paths[i], plot = plot, width = width, height = height,
                   units = "in", dpi = if (i == 5L) 300 else 600,
                   device = devices[[i]], limitsize = FALSE, bg = "white")
      if (i == 2L) args$family <- "Arial"
      if (i == 4L) args$compression <- "lzw"
      do.call(ggplot2::ggsave, args)
      if (!file.exists(paths[i]) || file.info(paths[i])$size == 0) stop("Empty export")
      "written_not_visually_checked"
    }, warning = function(w) {
      notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) {
      notes <<- c(notes, conditionMessage(e)); "failed"
    })
    data.frame(path = paths[i], status = status, reason = paste(unique(notes), collapse = "; "),
               width_mm = 183, height_mm = height * 25.4)
  })
  do.call(rbind, rows)
}

compact_panel <- function(p) {
  p + ggplot2::labs(caption = NULL) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(
      barheight = grid::unit(22, "mm"), barwidth = grid::unit(3, "mm"))) + ggplot2::theme(
    plot.title = ggplot2::element_text(size = 7.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6),
    legend.title = ggplot2::element_text(size = 6.5),
    legend.text = ggplot2::element_text(size = 6),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(3, "mm"),
    plot.margin = ggplot2::margin(3, 4, 3, 3))
}

make_species_figures <- function(cfg, computation_source = NULL, output_dir = NULL) {
  out_dir <- file.path(cfg$output_root, cfg$species_name)
  bundle_file <- file.path(out_dir, "castSDM_result.rds")
  bundle <- readRDS(bundle_file)
  if (is.null(computation_source)) {
    require_current_bundle(bundle, cfg)
  } else {
    if (!file.exists(computation_source)) stop("Missing original computation receipt")
    receipt <- readRDS(computation_source)
    expected <- bundle$computation_provenance %||% bundle$provenance
    if (is.null(expected) || !identical(receipt, expected)) {
      stop("Computation receipt does not match the saved result bundle")
    }
    if (is.null(receipt$source_hashes) && is.null(receipt$hashes) &&
        is.null(receipt$source) && is.null(receipt$signature)) {
      stop("Original computation receipt contains no source identity")
    }
  }
  fig_dir <- output_dir %||% file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  result <- bundle$result
  scenario <- cfg$report_scenario
  if (is.null(scenario) || !scenario %in% names(bundle$attribution)) {
    stop("report_scenario must explicitly name a computed scenario")
  }
  at <- bundle$attribution[[scenario]]
  projection <- bundle$projection
  scenario_title <- gsub("_", " ", scenario)
  table_dir <- file.path(out_dir, "tables")
  inputs <- c(bundle_file, at$raster_file, projection$current$hss_path,
    file.path(out_dir, "rasters", paste0(scenario, "_change_class.tif")),
    file.path(table_dir, c("training_data.csv", "cv_fold_metrics.csv", "projection_stats.csv")))
  if (any(!file.exists(inputs))) stop("Required result sources are missing")
  train <- utils::read.csv(file.path(table_dir, "training_data.csv"))
  # Coordinate values can be native projected x/y despite their legacy names.
  if (any(abs(train$lon) > 180 | abs(train$lat) > 90)) {
    crs <- terra::crs(terra::rast(projection$current$hss_path))
    pts <- sf::st_as_sf(train, coords = c("lon", "lat"), crs = crs)
    xy <- sf::st_coordinates(sf::st_transform(pts, 4326))
    train$lon <- xy[, 1]; train$lat <- xy[, 2]
  }
  status <- list()
  write_figure <- function(name, plot, height) {
    status[[name]] <<- save_plot(plot, file.path(fig_dir, paste0(name, ".png")),
                                 width = 183, height = height, apply_theme = FALSE)
    invisible(plot)
  }
  annotation <- function(title, caption) patchwork::plot_annotation(
    title = title, caption = caption, tag_levels = "a", theme = paper_theme())

  fold_plot <- plot(result$cv, lon = train$lon, lat = train$lat,
                    basemap = "china", metric = "auc")
  write_figure("01_model_validation", fold_plot & paper_theme(), 95)

  r <- terra::rast(at$raster_file)
  display_max <- unname(stats::quantile(c(
    terra::spatSample(r[["current"]], size = 50000, method = "regular", values = TRUE,
                      na.rm = TRUE)[, 1],
    terra::spatSample(r[["future"]], size = 50000, method = "regular", values = TRUE,
                      na.rm = TRUE)[, 1]), .995, na.rm = TRUE))
  display_max <- max(.01, display_max)
  current <- cast_plot_effect_map(at$raster_file, "current", "Current suitability",
    "Suitability", limits = c(0, display_max), breaks = pretty(c(0, display_max), 3))
  future <- cast_plot_effect_map(at$raster_file, "future", "Future suitability",
    "Suitability", limits = c(0, display_max), breaks = pretty(c(0, display_max), 3))
  display_change <- max(.01, unname(stats::quantile(abs(
    terra::spatSample(r[["total_change"]], size = 50000, method = "regular",
                      values = TRUE, na.rm = TRUE)[, 1]), .995, na.rm = TRUE)))
  change <- cast_plot_effect_map(at$raster_file, "total_change", "Suitability change",
    "Future - current", signed = TRUE, limits = c(-display_change, display_change))
  stats <- projection$stats
  parse_scenario <- function(s) {
    found <- regmatches(s, regexec("^(.*)_(ssp[0-9]+)_([0-9]+)_([0-9]+)$", s))[[1]]
    if (length(found) != 5L) stop("Unrecognised scenario name: ", s)
    c(climate_model = found[2], ssp = toupper(found[3]),
      period = paste(found[4], found[5], sep = "-"))
  }
  parsed <- as.data.frame(t(vapply(stats$scenario, parse_scenario, character(3))))
  stats <- cbind(stats, parsed)
  stats$period <- factor(stats$period, levels = unique(stats$period))
  area <- ggplot2::ggplot(stats, ggplot2::aes(period, pct_change,
    colour = ssp, group = ssp)) + ggplot2::geom_hline(yintercept = 0,
      linewidth = .3, colour = "grey65") + ggplot2::geom_line(linewidth = .5) +
    ggplot2::geom_point(size = 1.5) + ggplot2::facet_wrap(~climate_model, ncol = 1) +
    ggplot2::scale_colour_manual(values = c(SSP126 = "#009E73",
      SSP245 = "#E69F00", SSP585 = "#CC79A7")) +
    ggplot2::labs(title = "Thresholded range change", x = NULL,
      y = "Net change (%)", colour = NULL) + paper_theme() +
      ggplot2::theme(legend.position = "bottom") + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))
  habitat <- patchwork::wrap_plots(lapply(list(current, future, change, area), compact_panel),
    ncol = 2) + annotation("Predicted habitat change",
      paste0("Map scenario: ", scenario_title, ". Shared suitability display 0-", signif(display_max, 3),
        "; change clipped at +/-", signif(display_change, 3), " (p99.5).\n",
        "Range trends use the recorded threshold; continuous maps retain sub-threshold values."))
  write_figure("02_habitat_change", habitat, 155)

  s <- at$summary
  changing <- unique(s$driver[s$contrast == "only" & is.finite(s$mean_abs) & s$mean_abs > 0])
  if (!length(changing)) stop("The reporting scenario contains no non-zero driver response")
  contribution <- plot(at, contrast = "retained") +
    ggplot2::labs(title = "Future-value driver comparisons", caption = NULL) + paper_theme()
  short_labels <- c(bio05 = "Warmest-month temperature",
    vege_ch2019_ch_f20000 = "Vegetation height (41 x 41 km window)",
    lulc_glu1km_2020_water_prop_f1000 = "Surface-water cover (3 x 3 km window)",
    hydro_h90m_spi_max = "Stream-power index",
    dist_h90m_stream_dist_up_farth_max = "Upstream distance",
    dist_grip4_road = "Distance to road",
    hydro_h90m_stream_diff_up_farth = "Upstream elevation difference",
    anthro_foothii_2017_hii = "Human influence",
    hydro_h90m_outlet_diff_dw_basin = "Outlet elevation difference")
  label_key <- stats::setNames(short_labels, castSDM:::.cast_var_label(names(short_labels)))
  contribution <- contribution + ggplot2::scale_y_discrete(labels = function(x) {
    known <- x %in% names(label_key); x[known] <- label_key[x[known]]; x
  })
  # A single changing predictor has identical total/only/retained contrasts.
  # Display its map once, and state the resulting limitation explicitly.
  only <- s[s$contrast == "only", , drop = FALSE]
  top <- changing[which.max(only$mean_abs[match(changing, only$driver)])]
  bound <- max(.01, unname(stats::quantile(abs(terra::spatSample(
    r[[paste0("retained_", top)]], size = 50000, method = "regular",
    values = TRUE, na.rm = TRUE)[, 1]), .995, na.rm = TRUE)))
  driver <- cast_plot_effect_map(at$raster_file, paste0("retained_", top),
    paste("Scenario response:", castSDM:::.cast_var_label(top)), "Suitability change",
    signed = TRUE, limits = c(-bound, bound))
  note <- if (length(changing) == 1L) {
    "Only one predictor changes: total, driver-only and retained-driver contrasts coincide."
  } else "Individual driver contrasts are non-additive; the remaining joint response is reported in the data."
  driver_plot <- patchwork::wrap_plots(compact_panel(driver), compact_panel(contribution),
    ncol = 2, widths = c(1.1, 1.3)) + annotation("Environmental explanation of the scenario",
      paste0(note, "\nMap display clipped symmetrically at +/-", signif(bound, 3), " (p99.5); raw values retained."))
  write_figure("03_environmental_driver", driver_plot, 100)

  coverage <- cast_plot_effect_map(at$raster_file, "support_total",
    "Training-range coverage", "Training range", limits = c(0, 1), breaks = c(0, 1), method = "near")
  bound_sd <- max(.01, unname(stats::quantile(terra::spatSample(r[["total_sd"]],
    size = 50000, method = "regular", values = TRUE, na.rm = TRUE)[, 1], .995, na.rm = TRUE)))
  disagreement <- cast_plot_effect_map(at$raster_file, "total_sd",
    "Cross-model disagreement", "SD of change", limits = c(0, bound_sd)) +
    ggplot2::labs(subtitle = paste0("Display upper limit ", signif(bound_sd, 3), " (p99.5)"))
  coverage <- compact_panel(coverage) + ggplot2::guides(fill = ggplot2::guide_legend())
  reliable <- patchwork::wrap_plots(coverage, compact_panel(disagreement),
    ncol = 2) + annotation("Where the scenario explanation needs caution",
      "Range coverage checks marginal training quantiles. Model SD measures disagreement among engines.\nNeither establishes sampling uncertainty or ecological causal effects.")
  write_figure("04_support_and_disagreement", reliable, 100)

  for (name in c("cv_metrics", "cv_fold_metrics", "cv_fold_status")) {
    f <- file.path(table_dir, paste0(name, ".csv"))
    if (file.exists(f)) file.copy(f, file.path(fig_dir, paste0(name, ".csv")), overwrite = TRUE)
  }
  utils::write.csv(stats, file.path(fig_dir, "habitat_change.csv"), row.names = FALSE)
  utils::write.csv(s, file.path(fig_dir, "environmental_driver.csv"), row.names = FALSE)
  utils::write.csv(at$weights, file.path(fig_dir, "model_weights.csv"))
  exports <- do.call(rbind, status)
  exports$figure <- rep(names(status), lengths(lapply(status, function(x) x$path)))
  source_data <- data.frame(path = inputs, status = "present",
                             md5 = unname(tools::md5sum(inputs)))
  utils::write.csv(exports, file.path(fig_dir, "figure_status.csv"), row.names = FALSE)
  utils::write.csv(source_data, file.path(fig_dir, "source_data.csv"), row.names = FALSE)
  saveRDS(list(figure_status = exports, source_data = source_data,
    plotting_source = workflow_source_signature(cfg),
    computation_source = if (is.null(computation_source)) bundle$provenance else receipt,
    visual_QA = "not performed"), file.path(fig_dir, "figure_manifest.rds"))
  if (is.null(computation_source)) file.copy(file.path(fig_dir, "figure_manifest.rds"),
    file.path(out_dir, "figure_manifest.rds"), overwrite = TRUE)
  chosen <- stats[stats$scenario == scenario, , drop = FALSE]
  total_row <- s[s$contrast == "total", , drop = FALSE]
  cv_good <- sum(result$cv$fold_status == "evaluated")
  writeLines(c("# \u5148\u8bfb\u8fd9\u91cc\uff1a\u73af\u5883\u60c5\u666f\u5982\u4f55\u6539\u53d8\u9884\u6d4b\u9002\u5b9c\u5ea6", "",
    sprintf("\u672c\u6b21\u4e3b\u56fe\u60c5\u666f\uff1a%s\u3002\u9608\u503c\u5316\u9002\u5b9c\u9762\u79ef\u7684\u51c0\u53d8\u5316\u4e3a %+.2f%%\uff1b\u8fd9\u662f\u6a21\u578b\u60c5\u666f\u9884\u6d4b\uff0c\u4e0d\u662f\u5df2\u53d1\u751f\u7684\u5206\u5e03\u53d8\u5316\u3002",
      scenario, chosen$pct_change[1]),
    sprintf("\u7a7a\u95f4\u4ea4\u53c9\u9a8c\u8bc1\u5b9e\u9645\u8bc4\u4f30 %d/%d \u6298\uff1b\u65e0\u6548\u6298\u53ca\u9010\u6298\u6307\u6807\u89c1\u56fe 1 \u548c cv_fold_status.csv\u3002",
      cv_good, result$cv$k),
    sprintf("\u8be5\u60c5\u666f\u7684\u975e\u96f6\u53d8\u5316\u54cd\u5e94\u6d89\u53ca %d \u4e2a\u56e0\u5b50\uff1a%s\u3002\u8bad\u7ec3\u8fb9\u9645\u8303\u56f4\u8986\u76d6\u7387\u4e3a %.1f%%\u3002",
      length(changing), paste(changing, collapse = ", "), 100 * total_row$support[1]), "",
    "\u6309\u79d1\u5b66\u95ee\u9898\u987a\u5e8f\u770b\u56db\u5f20\u56fe\uff1a",
    "1. 01_model_validation\uff1a\u7a7a\u95f4\u5206\u7ec4\u4e0e\u9010\u6298\u9884\u6d4b\u8868\u73b0\uff1b\u4e0d\u8981\u53ea\u770b\u4e00\u4e2a\u5e73\u5747 AUC\u3002",
    "2. 02_habitat_change\uff1a\u5f53\u524d\u4e0e\u672a\u6765\u9002\u5b9c\u5ea6\u3001\u7a7a\u95f4\u53d8\u5316\u3001\u5168\u90e8 18 \u4e2a\u60c5\u666f\u7684\u9762\u79ef\u51c0\u53d8\u5316\u3002",
    "3. 03_environmental_driver\uff1a\u5b9e\u9645\u672a\u6765\u503c\u53d8\u5316\u5982\u4f55\u5f71\u54cd\u6a21\u578b\u54cd\u5e94\uff1b\u5b8c\u6574\u53d8\u91cf\u540d\u4e0e\u6307\u6807\u5728 environmental_driver.csv\u3002",
    "4. 04_support_and_disagreement\uff1a\u8bad\u7ec3\u8303\u56f4\u8986\u76d6\u4e0e\u6a21\u578b\u95f4\u5206\u6b67\uff1b\u4e24\u8005\u5747\u4e0d\u662f\u7f6e\u4fe1\u533a\u95f4\u3002", "",
    note,
    "\u5f53\u524d Schizopygopsis \u4f7f\u7528\u5c5e\u7ea7\u8bb0\u5f55\uff1b\u5168\u56fd\u9646\u5730\u80cc\u666f\u548c\u51fa\u73b0\u70b9\u7f13\u51b2\u533a\u5c1a\u4e0d\u80fd\u4ee3\u8868\u9c7c\u7c7b\u53ef\u8fbe\u6cb3\u6d41\u6816\u606f\u5730\u3002",
    "\u53ea\u6709\u4e00\u4e2a\u73af\u5883\u56e0\u5b50\u53d8\u5316\u65f6\uff0c\u603b\u53d8\u5316\u3001\u5355\u56e0\u5b50\u53d8\u5316\u548c\u4fdd\u7559\u56e0\u5b50\u5bf9\u6bd4\u76f8\u540c\uff1b\u5176\u4f59\u96f6\u503c\u8868\u793a\u6ca1\u6709\u63d0\u4f9b\u8be5\u56e0\u5b50\u7684\u672a\u6765\u53d8\u5316\uff0c\u4e0d\u8868\u793a\u8be5\u56e0\u5b50\u6ca1\u6709\u751f\u6001\u4f5c\u7528\u3002",
    "\u56fe\u7531\u5305\u5185\u7ed8\u56fe\u51fd\u6570\u751f\u6210\u3002\u91cd\u7ed8\u4fdd\u7559\u539f\u59cb\u8ba1\u7b97\u6765\u6e90\uff1b\u91cd\u7ed8\u5386\u53f2\u6570\u503c\u7ed3\u679c\u4e0d\u7b49\u4e8e\u91cd\u65b0\u5b8c\u6210\u62df\u5408\u5b9e\u9a8c\u3002",
    "\u5bfc\u51fa\u6210\u529f\u4e0e\u89c6\u89c9\u68c0\u67e5\u5206\u522b\u8bb0\u5f55\u3002"), file.path(fig_dir, "READ_FIRST.md"), useBytes = TRUE)
  if (any(exports$status == "failed")) stop("Figure export failed; inspect figure_status.csv")
  invisible(fig_dir)
}

