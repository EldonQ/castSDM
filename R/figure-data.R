# AI-assisted restoration; scope and verification: inst/AI_ASSISTANCE.md.
# Numerical products are cached independently of themes and file formats.
.cast_figure_split <- function(object) {
  path <- file.path(object$output_dir %||% "", "stages", "01_prepare.rds")
  if (file.exists(path)) readRDS(path) else NULL
}

.cast_figure_curves <- function(predictions, dataset) {
  if (!is.data.frame(predictions) || !"obs" %in% names(predictions)) return(NULL)
  rows <- lapply(grep("^HSS_", names(predictions), value = TRUE), function(col) {
    p <- predictions[[col]]; y <- predictions$obs
    keep <- is.finite(p) & y %in% c(0, 1)
    p <- p[keep]; y <- y[keep]
    if (!length(p) || length(unique(y)) < 2L) return(NULL)
    order <- order(-p); p <- p[order]; y <- y[order]
    last <- !duplicated(p, fromLast = TRUE)
    tp <- cumsum(y == 1)[last]; fp <- cumsum(y == 0)[last]
    data.frame(model = sub("^HSS_", "", col), dataset = dataset,
      threshold = c(Inf, p[last]), fpr = c(0, fp/sum(y == 0)),
      recall = c(0, tp/sum(y == 1)), precision = c(1, tp/(tp+fp)),
      sampled_prevalence = mean(y), n = length(y))
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)) do.call(rbind, rows) else NULL
}

.cast_figure_classes <- function(result, threshold, path) {
  r <- terra::rast(result$raster_file)
  cols <- c("current", "joint", "support", intersect("common_support", names(r)))
  r <- r[[cols]]
  terra::readStart(r); on.exit(terra::readStop(r), add = TRUE)
  layers <- c("range_class", "range_class_supported",
    if ("common_support" %in% cols) "range_class_common_support")
  writer <- .cast_raster_writer(r, path, layers, datatype = "INT2S", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  for (first in seq.int(1L, terra::nrow(r), by = 256L)) {
    nr <- min(256L, terra::nrow(r)-first+1L)
    x <- terra::readValues(r, row = first, nrows = nr, mat = TRUE)
    ok <- is.finite(x[, 1L]) & is.finite(x[, 2L])
    a <- x[, 1L] >= threshold; b <- x[, 2L] >= threshold
    value <- ifelse(a & b, 2, ifelse(a, -1, ifelse(b, 1, 0)))
    value[!ok] <- NA_real_
    z <- cbind(value, ifelse(x[, 3L] == 1 & is.finite(x[, 3L]), value, NA_real_))
    if (ncol(x) == 4L) z <- cbind(z,
      ifelse(x[, 4L] == 1 & is.finite(x[, 4L]), value, NA_real_))
    .cast_writer_block(writer, z, first, nr)
  }
  .cast_writer_finish(writer)
  path
}

.cast_figure_mess <- function(fit, current, future, domain, path, block_size = 250000L) {
  r <- terra::rast(domain)[["current"]]
  writer <- .cast_raster_writer(r, path,
    c("mess_current", "mess_future", "novel_current", "novel_future"),
    datatype = "FLT8S", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  rows <- max(1L, as.integer(block_size/terra::ncol(r)))
  for (first in seq.int(1L, terra::nrow(r), by = rows)) {
    nr <- min(rows, terra::nrow(r)-first+1L)
    eligible <- is.finite(terra::values(r, row = first, nrows = nr, mat = FALSE))
    z <- matrix(NA_real_, length(eligible), 4L)
    if (any(eligible)) {
      a <- .cast_raster_block(current, first, nr)
      b <- .cast_raster_block(future, first, nr)
      ok <- eligible & .cast_finite_rows(a) & .cast_finite_rows(b)
      if (any(ok)) {
        ma <- .cast_mess(fit$scaling$reference, a[ok, , drop = FALSE])
        mb <- .cast_mess(fit$scaling$reference, b[ok, , drop = FALSE])
        z[ok, ] <- cbind(ma, mb, as.numeric(ma < 0), as.numeric(mb < 0))
      }
    }
    .cast_writer_block(writer, z, first, nr)
  }
  .cast_writer_finish(writer)
  path
}

.cast_figure_display <- function(paths, output, resolution = .02) {
  stack <- do.call(c, lapply(paths, terra::rast))
  discrete <- grepl("^support|common_support|^range_class|^novel_|^sign_", names(stack))
  # The saved model rasters are never resampled or overwritten. These files
  # are a declared display-only grid, shared by the PNG/PDF panels.
  template <- terra::project(stack[[1L]], "EPSG:4326", res = resolution, method = "bilinear")
  continuous <- if (any(!discrete)) terra::project(stack[[which(!discrete)]], template, method = "bilinear")
  categorical <- if (any(discrete)) terra::project(stack[[which(discrete)]], template, method = "near")
  display <- if (is.null(continuous)) categorical else if (is.null(categorical)) continuous else c(continuous, categorical)
  display <- display[[names(stack)]]
  terra::writeRaster(display, output, overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES", "INTERLEAVE=BAND")))
  output
}

.cast_figure_data <- function(object, comparison, scenario, directory,
    current = NULL, future = NULL, attribution = NULL, render_only = FALSE,
    response_grid = 50L, display_resolution = .02, resources = list(), verbose = TRUE,
    products = c("maps", "responses", "attribution", "mess"), response_method = "ale") {
  fit <- object$fit_full %||% object$fit
  .cast_check_fit(fit)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  bundle_path <- file.path(directory, "figure_data.rds")
  provenance <- comparison$provenance
  provenance$session <- NULL
  base_key <- .cast_digest(list(schema = "cast_figure_data/3", .cast_fit_identity(fit), provenance,
    object[c("screen", "cv", "eval")], object$protocol$clamp,
    comparison$weights, comparison$threshold, comparison$threshold_sensitivity,
    scenario, response_grid, display_resolution, response_method,
    .cast_code_identity(cast_response_curves), .cast_code_identity(.cast_mess),
    .cast_code_identity(.cast_figure_classes), .cast_code_identity(cast_diagnose),
    .cast_code_identity(.cast_diagnostic_metrics), .cast_code_identity(compute_average_precision)))
  old <- if (file.exists(bundle_path)) tryCatch(readRDS(bundle_path), error = function(e) NULL) else NULL
  valid_old <- !is.null(old) && identical(old$key, base_key) &&
    !is.null(.cast_cache_read(bundle_path, base_key))
  verify_raw <- !render_only && any(products %in% c("attribution", "mess")) &&
    !is.null(current) && !is.null(future)
  products <- union(products, if (valid_old) old$value$products)
  if (valid_old && is.null(attribution) && !verify_raw && all(products %in% old$value$products)) {
    if (verbose) cli::cli_inform("Verified figure data: no fitting or diagnostic prediction required.")
    return(old$value)
  }
  if (valid_old && is.null(attribution) && !verify_raw && "attribution" %in% old$value$products)
    attribution <- old$value$attribution
  if (render_only) cli::cli_abort("A matching figure-data cache is missing. Run once with render_only = FALSE to compute frozen-model diagnostic products.")
  budget <- .cast_resources(resources)
  previous <- options(castSDM.predict_threads = budget$cores)
  on.exit(options(previous), add = TRUE)
  split <- .cast_figure_split(object)
  ev <- object$eval
  if (is.null(ev$predictions) && !is.null(split$test)) {
    if (verbose) cli::cli_inform("Preparing legacy hold-out scores from the saved fit (no refit).")
    ev <- cast_evaluate(object$fit, split$test, thresholds = object$cv$thresholds,
      clamp = object$eval$clamp %||% FALSE)
  }
  diagnostic_object <- object; diagnostic_object$eval <- ev
  diagnostic <- cast_diagnose(diagnostic_object, comparison)
  oof <- object$cv$oof; heldout <- ev$predictions
  if (length(comparison$weights)) for (kind in c("oof", "heldout")) {
    d <- if (kind == "oof") oof else heldout
    if (!is.data.frame(d)) next
    cols <- paste0("HSS_", names(comparison$weights))
    if (all(cols %in% names(d))) {
      x <- as.matrix(d[, cols, drop = FALSE]); pred <- as.numeric(x %*% comparison$weights)
      pred[!.cast_finite_rows(x)] <- NA_real_
      d$HSS_ensemble <- pred
    }
    if (kind == "oof") oof <- d else heldout <- d
  }
  curves <- rbind(.cast_figure_curves(oof, "Training OOF"),
                   .cast_figure_curves(heldout, "Hold-out"))
  reliability <- list()
  for (kind in c("Training OOF", "Hold-out")) {
    d <- if (kind == "Training OOF") oof else heldout
    if (!is.data.frame(d)) next
    for (col in grep("^HSS_", names(d), value = TRUE)) {
      bins <- .cast_diagnostic_bins(d[[col]], d$obs, sub("^HSS_", "", col), "declared")
      bins$dataset <- kind; reliability[[length(reliability)+1L]] <- bins
    }
  }
  cached <- function(name, key, compute, artifact = character()) {
    path <- file.path(directory, paste0(name, ".rds"))
    value <- .cast_cache_read(path, key)
    if (!is.null(value)) return(value)
    value <- force(compute)
    .cast_cache_write(value, path, key, artifact)
    value
  }
  responses <- NULL
  if ("responses" %in% products) {
    response_key <- .cast_digest(list(.cast_fit_identity(fit), response_grid,
      object$protocol$clamp, .cast_code_identity(cast_response_curves)))
    response_key <- .cast_digest(list(response_key, response_method))
    if (verbose) cli::cli_inform("Response curves: fixed profiles and {toupper(response_method)}; checking verified cache.")
    responses <- cached(paste0("responses_", response_method), response_key, {
      value <- list(fixed = cast_response_curves(fit, type = "fixed", grid_size = response_grid,
        clamp = object$protocol$clamp %||% TRUE))
      value[[response_method]] <- cast_response_curves(fit, type = response_method,
        grid_size = response_grid, clamp = object$protocol$clamp %||% TRUE)
      value
    })
  }
  outputs <- character(); displays <- list(); classes <- list(); notes <- list()
  note <- function(product, reason) notes[[length(notes)+1L]] <<-
    data.frame(product = product, reason = reason)
  if (!is.null(current)) {
    if (is.character(current)) current <- terra::rast(current)
    current <- current[[fit$env_vars]]
  }
  if (!is.null(future)) {
    if (is.list(future) && !inherits(future, "SpatRaster")) future <- future[[scenario]]
    if (is.character(future)) future <- terra::rast(future)
    for (v in setdiff(fit$env_vars, names(future))) future <- c(future, current[[v]])
    future <- future[[fit$env_vars]]
  }
  report <- comparison$results[[scenario]]
  if ("attribution" %in% products && is.null(attribution) && !is.null(current) && !is.null(future)) {
    if (verbose) cli::cli_inform("Preparing reporting-scenario driver maps; reuse requires matching models, inputs, mask and settings.")
    domain <- if (!is.null(report$raster_file)) terra::rast(report$raster_file)[["current"]] else NULL
    if (!is.null(domain)) domain <- terra::ifel(is.finite(domain), 1, NA)
    attribution <- cast_attribution(fit, current, future,
      drivers = unlist(comparison$groups, use.names = FALSE),
      weights = comparison$weights, mask = domain,
      clamp = object$protocol$clamp %||% FALSE, support_probs = comparison$support_probs %||% c(.01, .99),
      resources = resources, block_size = 250000L,
      overwrite = TRUE,
      filename = if (!is.null(domain)) file.path(directory, "driver_attribution.tif"))
  }
  if (is.null(attribution)) note("03_environmental_driver", "Current/future predictor inputs or an existing cast_attribution object are required.")
  mess <- NULL
  if ("mess" %in% products && !is.null(current) && !is.null(future) && !is.null(report$raster_file)) {
    path <- file.path(directory, "mess_native.tif")
    memo <- new.env(parent = emptyenv())
    key <- .cast_digest(list(.cast_fit_identity(fit), .cast_data_identity(current, memo),
      .cast_data_identity(future, memo), provenance, scenario, .cast_code_identity(.cast_figure_mess)))
    mess <- cached("mess", key, .cast_figure_mess(fit, current, future, report$raster_file, path), path)
    outputs <- c(outputs, mess)
  } else note("09_extrapolation_mess", "Raw current/future predictor stacks are required for native-grid MESS; support flags cannot substitute for MESS.")
  for (s in if ("maps" %in% products) names(comparison$results) else character()) {
    result <- comparison$results[[s]]
    if (is.null(result$raster_file)) next
    if (verbose) cli::cli_inform("Preparing display grid: {s} (original scientific raster unchanged).")
    input <- result$raster_file
    if (!is.null(comparison$threshold)) {
      path <- file.path(directory, paste0("range_", s, ".tif"))
      key <- .cast_digest(list(unname(tools::md5sum(input)), comparison$threshold, .cast_code_identity(.cast_figure_classes)))
      classes[[s]] <- cached(paste0("range_", s), key,
        .cast_figure_classes(result, comparison$threshold, path), path)
      outputs <- c(outputs, classes[[s]]); input <- c(input, classes[[s]])
    }
    path <- file.path(directory, paste0("display_", s, ".tif"))
    key <- .cast_digest(list(unname(tools::md5sum(input)),
      display_resolution, .cast_code_identity(.cast_figure_display)))
    displays[[s]] <- cached(paste0("display_", s), key,
      .cast_figure_display(input, path, display_resolution), path)
    outputs <- c(outputs, displays[[s]])
    invisible(gc(FALSE))
  }
  mess_display <- NULL
  if (!is.null(mess)) {
    path <- file.path(directory, "display_mess.tif")
    key <- .cast_digest(list(unname(tools::md5sum(mess)), display_resolution,
      .cast_code_identity(.cast_figure_display)))
    mess_display <- cached("display_mess", key,
      .cast_figure_display(mess, path, display_resolution), path)
    outputs <- c(outputs, mess_display)
  }
  attribution_display <- NULL
  if (!is.null(attribution$raster_file)) {
    outputs <- c(outputs, attribution$raster_file)
    path <- file.path(directory, "display_drivers.tif")
    key <- .cast_digest(list(unname(tools::md5sum(attribution$raster_file)),
      display_resolution, .cast_code_identity(.cast_figure_display)))
    attribution_display <- cached("display_drivers", key,
      .cast_figure_display(attribution$raster_file, path, display_resolution), path)
    outputs <- c(outputs, attribution_display)
  }
  info <- object$cv$fold_info
  if (is.null(info) && !is.null(split$train) && length(object$cv$folds) == nrow(split$train))
    info <- data.frame(fold = object$cv$folds, lon = split$train$lon, lat = split$train$lat,
      obs = split$train$presence)
  crs <- object$protocol$coords_crs %||% if (!is.null(report$raster_file)) terra::crs(terra::rast(report$raster_file)) else NULL
  if (!is.null(info) && any(abs(info$lon) > 180 | abs(info$lat) > 90)) {
    if (is.null(crs) || !nzchar(crs)) cli::cli_abort("Projected fold coordinates require a raster CRS for geographic maps.")
    p <- sf::st_as_sf(info, coords = c("lon", "lat"), crs = crs)
    xy <- sf::st_coordinates(sf::st_transform(p, 4326)); info$lon <- xy[, 1L]; info$lat <- xy[, 2L]
  }
  thresholds <- summary(comparison, thresholds = unique(c(comparison$threshold, .3, .5, .7)))
  origin_path <- file.path(object$output_dir %||% "", "manifest.rds")
  origin <- if (file.exists(origin_path)) readRDS(origin_path) else NULL
  completed <- products
  if (is.null(attribution)) completed <- setdiff(completed, "attribution")
  if (is.null(mess)) completed <- setdiff(completed, "mess")
  if (!length(displays)) completed <- setdiff(completed, "maps")
  data <- list(schema = "cast_figure_data/3", analysis_key = base_key, products = completed, screen = object$screen, cv = object$cv,
    evaluation = ev, diagnostics = diagnostic, folds = info, responses = responses,
    curves = curves, reliability = if (length(reliability)) do.call(rbind, reliability),
    comparison = comparison, attribution = attribution, attribution_display = attribution_display,
    displays = displays, mess_display = mess_display, range_files = classes, thresholds = thresholds,
    scenario = scenario, variables = fit$env_vars, model_names = names(fit$models),
    response_reference_n = nrow(fit$scaling$reference), pdp_reference_n = min(2000L, nrow(fit$scaling$reference)),
    notes = if (length(notes)) do.call(rbind, notes), display_resolution = display_resolution,
    source_version = as.character(utils::packageVersion("castSDM")),
    training_version = origin$cast_version %||% NA_character_)
  # Receipt file verification is restricted to sibling products, by design.
  sibling <- outputs[file.exists(outputs) & dirname(outputs) == directory]
  .cast_cache_write(data, bundle_path, base_key, unique(sibling))
  data
}

# Public presets and individual names share the same dependency resolver.
.cast_figure_choose <- function(include, inventory) {
  if (!is.character(include) || !length(include) || anyNA(include))
    cli::cli_abort("include must name figure presets or individual figures.")
  unique(unlist(lapply(include, function(x) {
    if (x == "all") return(inventory)
    pattern <- switch(x,
      basic = "^(01_|07_|11_|12_|13_|15_|S_selection_)",
      maps = "^(02_|04_|05_|08_|10_|S_(?!selection_))",
      responses = "^(06_|14_)", drivers = "^(03_|D_)", novelty = "^09_", NULL)
    if (!is.null(pattern)) return(grep(pattern, inventory, value = TRUE, perl = TRUE))
    inventory[inventory == x | startsWith(inventory, paste0(x, "_"))]
  }), use.names = FALSE))
}

.cast_figure_products <- function(include) {
  if (!is.character(include) || !length(include) || anyNA(include))
    cli::cli_abort("include must name figure presets or individual figures.")
  if ("all" %in% include) return(c("maps", "responses", "attribution", "mess"))
  out <- character()
  for (x in include) {
    if (x %in% c("maps", "novelty") || grepl("^(02|04|05|08|09|10)(_|$)", x) ||
        (startsWith(x, "S_") && !startsWith(x, "S_selection")) || x == "S") out <- c(out, "maps")
    if (x == "responses" || grepl("^(06|14)(_|$)", x)) out <- c(out, "responses")
    if (x == "drivers" || grepl("^03(_|$)|^D_", x)) out <- c(out, "attribution")
    if (x == "novelty" || grepl("^09(_|$)", x)) out <- c(out, "mess")
  }
  unique(out)
}
