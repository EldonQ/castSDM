# AI-assisted restoration; scope and verification: inst/AI_ASSISTANCE.md.
.cast_figure_pdf_backend <- function(prefer = getOption("castSDM.figure_pdf_device", "auto")) {
  if (identical(prefer, "pdf")) return("pdf")
  if (identical(prefer, "cairo")) return("cairo")
  if (!isTRUE(capabilities("cairo"))) return("pdf")
  path <- tempfile(fileext = ".pdf")
  device <- grDevices::dev.cur()
  on.exit({
    if (grDevices::dev.cur() != device) grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  ready <- tryCatch(suppressWarnings({
    grDevices::cairo_pdf(path, family = "Arial")
    grid::grid.text("castSDM", gp = grid::gpar(fontfamily = "Arial"))
    grDevices::dev.off()
    TRUE
  }), error = function(e) FALSE)
  if (ready) "cairo" else "pdf"
}

# Construct the layout once on the first export device. Reusing the grob
# avoids rebuilding every large raster panel for the matching PDF.
.cast_figure_export <- function(plot, path, width_mm, height_mm, dpi, grob = NULL) {
  if (endsWith(path, ".png")) {
    ragg::agg_png(path, width = width_mm, height = height_mm, units = "mm",
      res = dpi, background = "white")
  } else if (endsWith(path, ".pdf")) {
    if (identical(getOption("castSDM.figure_pdf_device", "cairo"), "pdf")) {
      grDevices::pdf(path, width = width_mm/25.4, height = height_mm/25.4,
        family = "Helvetica", bg = "white", useDingbats = FALSE)
    } else {
      grDevices::cairo_pdf(path, width = width_mm/25.4, height = height_mm/25.4,
        family = "Arial", bg = "white")
    }
  } else if (endsWith(path, ".svg")) {
    svglite::svglite(path, width = width_mm/25.4, height = height_mm/25.4,
      bg = "white", system_fonts = list(sans = "Arial"))
  } else {
    ragg::agg_tiff(path, width = width_mm, height = height_mm, units = "mm",
      res = dpi, compression = "lzw", background = "white")
  }
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  if (is.null(grob)) grob <- if (inherits(plot, "patchwork"))
    patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)
  grid::grid.newpage()
  grid::grid.draw(grob)
  grob
}

#' Export the Complete Saved-Result Figure Collection
#'
#' Restores the numbered validation, habitat, driver, support, factorial,
#' response, performance, range, novelty and future-scenario figures, plus
#' selection, threshold, weight, AOA and all available contrast/driver atlases.
#' Uses the package cartography, including China outlines and South China Sea
#' inset when `basemap = "china"`. Numerical diagnostic products are cached
#' separately from PNG/PDF rendering; fitted models are never refitted.
#' Unchanged installed-package exports are reused after checking their data,
#' code, rendering settings and file hashes. Missing or damaged exports are
#' regenerated individually.
#'
#' @param object A `cast_result` from [cast_run()].
#' @param comparison A `cast_comparison` from [cast_compare()].
#' @param output_dir Figure directory. Defaults to a sibling `figures` directory
#'   beside the training output, or a temporary directory for in-memory runs.
#' @param scenario Reporting scenario. Defaults to the last requested future.
#' @param current,future Optional raw predictor rasters for native-grid MESS
#'   and driver maps. `future` may be a single raster or named scenario list.
#'   Omitted static predictors inherit `current`.
#' @param attribution Optional saved `cast_attribution` for the reporting
#'   scenario, avoiding repeat driver predictions. It must use the same frozen
#'   models/weights and extrapolation policy.
#' @param basemap `"world"`, `"china"` or `"none"`.
#' @param include `"all"` (default), exact figure names, or numbered prefixes
#'   such as `c("01", "02", "07")`. All available raster layers and all
#'   requested futures are included in the default supplementary atlases.
#' @param render_only Require an existing matching numerical figure-data
#'   bundle. No fitting, diagnostic prediction or response calculation occurs.
#' @param response_grid Response-curve grid size, default `50`. Partial
#'   dependence uses [cast_response_curves()]' declared deterministic reference
#'   cap of 2,000 rows; every fitted predictor and engine is retained.
#' @param display_resolution WGS84 display grid resolution in degrees, default
#'   `0.02`. Original scientific rasters and native-grid range counts are
#'   unchanged. Discrete layers use nearest-neighbour display resampling.
#' @param width_mm,dpi Export width and PNG resolution, defaults `183` and
#'   `600`; matching PDF exports retain vector boundaries/text and embed dense
#'   gridded data as rasters. Display colour clipping is reported on maps.
#'   Cairo/Arial is preferred; systems without a working Cairo text device use
#'   standard PDF/Helvetica. The selected backend/font is recorded.
#' @param resources Resource budget for diagnostic prediction; see [cast_compare()].
#' @param verbose Print diagnostic and figure progress.
#' @param formats Export formats. Defaults to `c("png", "pdf")` to preserve the
#'   established collection. Optional `"svg"` (requires `svglite`) and `"tiff"`
#'   reuse the same layout; no duplicate diagnostic calculations are needed.
#' @return A `cast_figures` list with `status`, `requested`, `tables`, `data_file`, output
#'   directory and the full figure inventory. Missing inputs are reported as
#'   `unavailable` with a reason; export failures abort after recording status.
#' @details
#' The first complete export may compute expensive driver/PDP products from
#' frozen fits. Subsequent complete exports reuse verified products. Use
#' `render_only = TRUE` for theme/format changes. Native-grid classifications
#' use the declared threshold before reprojection. Marginal support, MESS and
#' AOA remain distinct; paired model SD/sign agreement are not confidence.
#' Hold-out reliability refers to sampled presence/background responses.
#' @examples
#' # cast_figures(run, comparison, "results/figures", basemap = "china",
#' #   current = current_stack, future = future_stack)
#' # cast_figures(run, comparison, "results/figures", render_only = TRUE)
#' @export
cast_figures <- function(object, comparison, output_dir = NULL,
    scenario = utils::tail(names(comparison$results), 1L), current = NULL,
    future = NULL, attribution = NULL, basemap = c("world", "china", "none"),
    include = "all", render_only = FALSE, response_grid = 50L,
    display_resolution = .02, width_mm = 183, dpi = 600,
    resources = list(), verbose = TRUE, formats = c("png", "pdf")) {
  if (!inherits(object, "cast_result") || !inherits(comparison, "cast_comparison"))
    cli::cli_abort("Supply a cast_result and its cast_comparison.")
  if (!is.character(scenario) || length(scenario) != 1L || !scenario %in% names(comparison$results))
    cli::cli_abort("scenario must name a saved comparison scenario.")
  basemap <- match.arg(basemap)
  for (pkg in c("ggplot2", "patchwork", "sf", "terra", "scales", "ragg"))
    check_suggested(pkg, "for the complete figure collection")
  if (!is.character(formats) || !length(formats) || anyNA(formats) ||
      any(!formats %in% c("png", "pdf", "svg", "tiff")))
    cli::cli_abort("formats must contain png, pdf, svg or tiff.")
  formats <- unique(formats)
  if ("svg" %in% formats) check_suggested("svglite", "for editable SVG figures")
  width_mm <- .cast_check_scalar(width_mm, "width_mm", lower = 40)
  dpi <- .cast_check_scalar(dpi, "dpi", lower = 72, integer = TRUE)
  response_grid <- .cast_check_scalar(response_grid, "response_grid", lower = 2, integer = TRUE)
  display_resolution <- .cast_check_scalar(display_resolution, "display_resolution", lower = .0001)
  if (!is.null(current) && !inherits(current, "SpatRaster") && !is.character(current))
    cli::cli_abort("current must be a raw predictor SpatRaster or raster file path(s).")
  if (is.list(future) && !inherits(future, "SpatRaster")) {
    if (!scenario %in% names(future)) cli::cli_abort("The future list must contain the reporting scenario.")
    value <- future[[scenario]]
  } else value <- future
  if (!is.null(value) && !inherits(value, "SpatRaster") && !is.character(value))
    cli::cli_abort("future must contain raw predictor SpatRaster objects or raster file path(s).")
  if (!is.null(attribution) && (!inherits(attribution, "cast_attribution") ||
      !identical(attribution$weights, comparison$weights)))
    cli::cli_abort("attribution must be a cast_attribution with the comparison's frozen model weights.")
  if (is.null(output_dir)) output_dir <- if (!is.null(object$output_dir))
    file.path(dirname(object$output_dir), "figures") else tempfile("cast-figures-")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  lock <- .cast_lock(file.path(output_dir, ".figures.lock"))
  on.exit(filelock::unlock(lock), add = TRUE)
  pdf_device <- if ("pdf" %in% formats) .cast_figure_pdf_backend() else "cairo"
  font_family <- if (pdf_device == "pdf") "Helvetica" else "Arial"
  old_options <- options(castSDM.font_family = font_family, castSDM.figure_width_mm = width_mm,
    castSDM.figure_pdf_device = pdf_device)
  on.exit(options(old_options), add = TRUE)
  data_dir <- file.path(output_dir, ".data")
  data <- .cast_figure_data(object, comparison, scenario, data_dir,
    current, future, attribution, render_only, response_grid,
    display_resolution, resources, verbose)
  registry <- .cast_figure_registry(data, basemap)
  if (!is.character(include) || !length(include) || anyNA(include))
    cli::cli_abort("include must be 'all', figure names or numbered prefixes.")
  selected <- if (identical(include, "all")) names(registry) else names(registry)[
    vapply(names(registry), function(name) any(vapply(include, function(i)
      identical(i, name) || startsWith(name, paste0(i, "_")), logical(1))), logical(1))]
  if (!length(selected)) cli::cli_abort("No requested figure names match the inventory.")
  table_dir <- file.path(output_dir, "source_data")
  dir.create(table_dir, showWarnings = FALSE)
  tables <- list(screening = data$screen$scores, selection_frequency = data$cv$selection_freq,
    spatial_fold_metrics = data$cv$fold_metrics, spatial_fold_coordinates = data$folds,
    holdout_metrics = data$diagnostics$holdout, oof_score_views = data$diagnostics$score_views_oof,
    metric_definitions = data$diagnostics$metric_definitions, diagnostic_notes = data$diagnostics$notes,
    roc_pr_curves = data$curves, reliability = data$reliability,
    response_curves = as.data.frame(data$responses$pdp),
    fixed_profiles = as.data.frame(data$responses$fixed),
    factorial_summary = comparison$summary, factorial_binary = comparison$binary,
    factorial_binary_supported = comparison$binary_supported,
    common_support_summary = comparison$summary_common_support,
    common_support_binary = comparison$binary_common_support,
    model_agreement = comparison$agreement, threshold_sensitivity = data$thresholds,
    environmental_driver = data$attribution$summary,
    model_weights = data.frame(model = names(comparison$weights), weight = as.numeric(comparison$weights)),
    figure_provenance = data.frame(key = c("training_software_version", "figure_software_version",
      "report_scenario", "native_binary_threshold", "display_resolution_degrees", "display_resampling",
      "response_grid_size", "pdp_reference_rows", "full_training_reference_rows"),
      value = c(data$training_version, data$source_version, scenario,
        comparison$threshold %||% NA_real_, display_resolution,
        "bilinear continuous; nearest categorical; native scientific rasters unchanged",
        response_grid, data$pdp_reference_n, data$response_reference_n)),
    unavailable_products = data$notes)
  table_rows <- list()
  for (name in names(tables)) if (is.data.frame(tables[[name]])) {
    path <- file.path(table_dir, paste0(name, ".csv"))
    utils::write.csv(tables[[name]], path, row.names = FALSE)
    table_rows[[name]] <- data.frame(path = file.path("source_data", basename(path)),
      rows = nrow(tables[[name]]), columns = ncol(tables[[name]]), md5 = unname(tools::md5sum(path)))
  }
  manifest_path <- file.path(output_dir, "figure_manifest.rds")
  previous_manifest <- if (file.exists(manifest_path))
    tryCatch(readRDS(manifest_path), error = function(e) NULL) else NULL
  data_md5 <- unname(tools::md5sum(file.path(data_dir, "figure_data.rds")))
  code_file <- file.path(getNamespaceInfo(asNamespace("castSDM"), "path"), "R", "castSDM.rdb")
  # Installed code is immutable for a run. Development namespaces always draw.
  render_seed <- .cast_digest(list(data_md5, basemap, width_mm, dpi, font_family, pdf_device, as.character(getRversion()),
    as.character(utils::packageVersion("castSDM")),
    if (file.exists(code_file)) unname(tools::md5sum(code_file)),
    vapply(c("ggplot2", "patchwork", "sf", "terra", "ragg", "svglite"), function(p)
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "missing", character(1))))
  previous_status <- previous_manifest$status
  status <- list()
  for (name in selected) {
    spec <- registry[[name]]
    paths <- file.path(output_dir, paste0(name, ".", formats))
    if (!spec$available) {
      status[[name]] <- data.frame(figure = name, path = basename(paths), status = "unavailable",
        reason = spec$reason, warnings = "", width_mm = width_mm, height_mm = spec$height_mm,
        cache_hit = FALSE, render_key = "", md5 = NA_character_)
      if (verbose) cli::cli_inform("{name}: unavailable ({spec$reason})")
      next
    }
    keys <- vapply(paths, function(path) .cast_digest(list(render_seed, name, spec$height_mm,
      tools::file_ext(path))), character(1))
    reusable <- rep(FALSE, length(paths))
    previous_rows <- rep(NA_integer_, length(paths))
    if (file.exists(code_file) && is.data.frame(previous_status) &&
        all(c("render_key", "md5") %in% names(previous_status))) {
      previous_rows <- match(basename(paths), previous_status$path)
      for (i in which(!is.na(previous_rows) & file.exists(paths))) {
        old <- previous_status[previous_rows[i], ]
        reusable[i] <- identical(old$status, "written_not_visually_checked") &&
          identical(old$render_key, unname(keys[i])) &&
          identical(old$md5, unname(tools::md5sum(paths[i])))
      }
    }
    if (all(reusable)) {
      rows <- previous_status[previous_rows, , drop = FALSE]; rows$cache_hit <- TRUE
      status[[name]] <- rows
      if (verbose) cli::cli_inform("{name}: verified {length(paths)} figure exports.")
      next
    }
    if (verbose) cli::cli_inform("Rendering {name} ({width_mm} mm; {dpi} dpi; {paste(toupper(formats), collapse = ' + ')}).")
    plot <- tryCatch(spec$build(), error = identity)
    if (inherits(plot, "error")) {
      status[[name]] <- data.frame(figure = name, path = basename(paths), status = "failed",
        reason = conditionMessage(plot), warnings = "", width_mm = width_mm, height_mm = spec$height_mm,
        cache_hit = FALSE, render_key = "", md5 = NA_character_)
      next
    }
    grob <- NULL
    export_rows <- vector("list", length(paths))
    for (i in seq_along(paths)) {
      if (reusable[i]) {
        export_rows[[i]] <- previous_status[previous_rows[i], , drop = FALSE]
        export_rows[[i]]$cache_hit <- TRUE
        next
      }
      notes <- character()
      outcome <- tryCatch(withCallingHandlers({
        grob <- .cast_figure_export(plot, paths[i], width_mm, spec$height_mm, dpi, grob)
        if (!file.exists(paths[i]) || file.info(paths[i])$size == 0) stop("Empty figure export.")
        NULL
      }, warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = identity)
      export_rows[[i]] <- data.frame(figure = name, path = basename(paths[i]),
        status = if (inherits(outcome, "error")) "failed" else "written_not_visually_checked",
        reason = if (inherits(outcome, "error")) conditionMessage(outcome) else "",
        warnings = paste(unique(notes), collapse = "; "), width_mm = width_mm, height_mm = spec$height_mm,
        cache_hit = FALSE, render_key = unname(keys[i]),
        md5 = if (!inherits(outcome, "error")) unname(tools::md5sum(paths[i])) else NA_character_)
    }
    status[[name]] <- do.call(rbind, export_rows)
    rm(plot, grob); invisible(gc(FALSE))
  }
  status <- do.call(rbind, status); rownames(status) <- NULL
  # A subset re-render must not erase the export record for the other figures.
  # Restrict retained rows to the same numerical bundle and current inventory.
  if (identical(previous_manifest$data_md5, data_md5) &&
      is.data.frame(previous_manifest$status) &&
      identical(names(previous_manifest$status), names(status))) {
    retained <- previous_manifest$status[!previous_manifest$status$path %in% status$path &
      previous_manifest$status$figure %in% names(registry), , drop = FALSE]
    missing <- retained$status == "written_not_visually_checked" &
      !file.exists(file.path(output_dir, retained$path))
    retained$status[missing] <- "missing"
    retained$reason[missing] <- "Previously exported file is missing; include this figure in a re-render."
    status <- rbind(retained, status)
    status <- status[order(match(status$figure, names(registry)), status$path), , drop = FALSE]
    rownames(status) <- NULL
  }
  table_manifest <- if (length(table_rows)) do.call(rbind, table_rows) else data.frame()
  utils::write.csv(status, file.path(output_dir, "figure_status.csv"), row.names = FALSE)
  utils::write.csv(table_manifest, file.path(output_dir, "table_manifest.csv"), row.names = FALSE)
  inventory <- do.call(rbind, lapply(names(registry), function(name) data.frame(
    figure = name, available = registry[[name]]$available, reason = if (registry[[name]]$available) "" else registry[[name]]$reason)))
  utils::write.csv(inventory, file.path(output_dir, "figure_inventory.csv"), row.names = FALSE)
  result <- structure(list(status = status, requested = selected, tables = table_manifest, inventory = inventory,
    output_dir = output_dir, data_file = file.path(data_dir, "figure_data.rds"),
    data_md5 = data_md5,
    scenario = scenario, basemap = basemap, dpi = dpi, width_mm = width_mm, formats = formats,
    font_family = font_family, pdf_device = pdf_device,
    display_resolution = display_resolution, render_only = render_only,
    interpretation = "Frozen-model diagnostics; source scientific rasters are unchanged.",
    visual_QA = "not performed automatically"), class = "cast_figures")
  .cast_atomic_rds(result, manifest_path)
  if (any(status$status == "failed")) {
    failed <- status[status$status == "failed", ]
    details <- unique(paste0(failed$path, ": ", failed$reason))
    cli::cli_abort(c("Figure exports failed: {paste(utils::head(details, 3L), collapse = '; ')}",
      i = "Inspect figure_status.csv for all reasons; other successful exports are retained."))
  }
  if (verbose) cli::cli_inform("Figure collection: {sum(status$status == 'written_not_visually_checked')} exports ready ({sum(status$cache_hit & status$figure %in% selected)} requested exports reused); {sum(status$status == 'unavailable')} unavailable outputs documented.")
  result
}

#' @export
print.cast_figures <- function(x, ...) {
  cli::cli_h1("castSDM figure collection")
  cli::cli_inform("{length(unique(x$status$figure))} figure sets; {sum(x$status$status == 'written_not_visually_checked')} files written; {sum(x$status$status == 'unavailable')} unavailable outputs recorded.")
  cli::cli_inform("Figures, source tables and export status: {.path {x$output_dir}}")
  invisible(x)
}
