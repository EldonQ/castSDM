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

.cast_figure_mess <- function(fit, current, future, domain, path, block_size = 250000L,
                              prepared = NULL) {
  r <- terra::rast(domain)[["current"]]
  writer <- .cast_raster_writer(r, path,
    c("mess_current", "mess_future", "novel_current", "novel_future"),
    datatype = "FLT8S", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  reference <- .cast_mess_reference(fit$scaling$reference)
  shared <- names(reference)[vapply(names(reference), function(v)
    .cast_shared_raster_layers(current, future, v), logical(1))]
  changing <- setdiff(names(reference), shared)
  blocks <- prepared$blocks
  if (!is.null(blocks) && blocks$rows*terra::ncol(r) > block_size) blocks <- NULL
  rows <- if (!is.null(blocks)) blocks$rows else max(1L, as.integer(block_size/terra::ncol(r)))
  for (first in seq.int(1L, terra::nrow(r), by = rows)) {
    nr <- min(rows, terra::nrow(r)-first+1L)
    eligible <- is.finite(.cast_raster_block(r, first, nr)[[1L]])
    z <- matrix(NA_real_, length(eligible), 4L)
    if (any(eligible)) {
      saved <- if (!is.null(blocks)) .cast_compare_current_block(blocks, first, nr) else
        list(indices = which(eligible), values = .cast_raster_block(current, first, nr)[eligible, , drop = FALSE])
      use <- eligible[saved$indices]; positions <- saved$indices[use]
      a <- saved$values[use, , drop = FALSE]; b <- a
      if (length(changing) && length(positions)) b[changing] <-
        .cast_raster_cells(future[[changing]], first, nr, positions)
      ok <- .cast_finite_rows(a) & .cast_finite_rows(b)
      if (any(ok)) {
        evaluate <- function(vars, x) if (length(vars))
          .cast_mess_apply(reference[vars], x[, vars, drop = FALSE]) else rep(Inf, nrow(x))
        static <- evaluate(shared, a[ok, , drop = FALSE])
        ma <- pmin(static, evaluate(changing, a[ok, , drop = FALSE]))
        mb <- pmin(static, evaluate(changing, b[ok, , drop = FALSE]))
        z[positions[ok], ] <- cbind(ma, mb, as.numeric(ma < 0), as.numeric(mb < 0))
      }
    }
    .cast_writer_block(writer, z, first, nr)
  }
  .cast_writer_finish(writer)
  path
}

.cast_figure_display <- function(paths, output, resolution = .02, layers = NULL,
    references = NULL, shared_key = NULL) {
  stack <- if (inherits(paths, "SpatRaster")) paths else do.call(c, lapply(paths, terra::rast))
  if (!is.null(layers)) stack <- stack[[layers]]
  if (is.null(references)) references <- stats::setNames("full_domain", names(stack)[1L])
  cache <- getOption("castSDM.figure_display_cache")
  template_key <- paste0("template_", .cast_digest(list(terra::crs(stack),
    as.vector(terra::ext(stack)), terra::res(stack), resolution)))
  template <- if (is.environment(cache)) cache[[template_key]]
  if (is.null(template)) {
    template <- terra::project(terra::rast(stack[[1L]]), "EPSG:4326", res = resolution)
    if (is.environment(cache)) cache[[template_key]] <- template
  }
  original <- stack
  shared <- if (is.environment(cache) && !is.null(shared_key)) cache[[shared_key]]
  shared_layers <- intersect(grep("^current($|_)|^common_support$", names(stack), value = TRUE),
    if (!is.null(shared)) names(terra::rast(shared$path)) else names(stack))
  if (!is.null(shared)) stack <- stack[[setdiff(names(stack), shared_layers)]]
  discrete <- grepl("^support|common_support|^range_class|^novel_|^sign_|^context_direction_", names(stack))
  # The saved model rasters are never resampled or overwritten. These files
  # are a declared display-only grid, shared by the PNG/PDF panels.
  support <- grepl("^support|^common_support|^novel_", names(stack))
  continuous <- if (any(!discrete)) terra::project(stack[[which(!discrete)]], template, method = "average")
  categorical <- if (any(discrete & !support)) terra::project(stack[[which(discrete & !support)]], template, method = "mode")
  # Preserve mixed support explicitly. A coarse support flag must never mask
  # effects: their eligibility is applied on the native grid before averaging.
  novel <- startsWith(names(stack), "novel_")
  flags <- if (any(support & !novel)) {
    r <- stack[[which(support & !novel)]]
    low <- terra::project(r, template, method = "min")
    high <- terra::project(r, template, method = "max")
    terra::ifel(low == high, low, 2)
  }
  novelty <- if (any(novel)) terra::project(stack[[which(novel)]], template, method = "max")
  display <- do.call(c, Filter(Negate(is.null), list(continuous, categorical, flags, novelty)))
  if (!is.null(shared) && length(shared_layers)) display <- c(display, terra::rast(shared$path)[[shared_layers]])
  display <- display[[names(original)]]
  terra::writeRaster(display, output, overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES", "INTERLEAVE=BAND")))
  coverage <- if (!is.null(shared) && identical(references, shared$references)) shared$coverage else
    .cast_display_coverage(original, display, resolution, references)
  utils::write.csv(coverage, paste0(output, ".csv"), row.names = FALSE)
  if (is.environment(cache) && !is.null(shared_key) && is.null(shared))
    cache[[shared_key]] <- list(path = output, coverage = coverage, references = references)
  output
}

# Separate full-domain D_ maps from native-masked main/E_ display products.
# Masks have exactly the native scopes used by their numerical summaries.
.cast_figure_driver_display <- function(attribution, output, resolution, drivers = attribution$drivers, include_diagnostics = TRUE) {
  r <- terra::rast(attribution$raster_file)
  contrasts <- c("total_change", paste0("effect_", drivers), paste0("retained_", drivers), "remainder")
  layers <- c(paste0("effect_", drivers), if (include_diagnostics) paste0("retained_", drivers))
  n_contrasts <- length(layers)
  references <- stats::setNames("all_driver_paired_support", layers[1L])
  diagnostics <- list()
  for (v in if (include_diagnostics) drivers else character()) {
    cols <- paste0(c("context_direction_", "context_gap_", "effect_agreement_", "retained_agreement_"), v)
    if (!all(cols %in% names(r))) next
    diagnostics[[v]] <- cols
    layers <- c(layers, cols)
    references <- c(references, stats::setNames(paste0("driver_paired_support:", v), cols[1L]))
  }
  # Read the saved effects once in bounded blocks. Do not materialise the
  # national multiband stack or repeatedly mask it for each driver.
  needed <- unique(c(contrasts, "support_total", layers,
    if (length(diagnostics)) paste0("support_effect_", names(diagnostics)),
    if (length(diagnostics)) paste0("support_retained_", names(diagnostics))))
  r <- r[[needed]]
  terra::readStart(r)
  on.exit(terra::readStop(r), add = TRUE)
  path <- tempfile("cast-native-mask-", fileext = ".tif")
  on.exit(unlink(path), add = TRUE)
  writer <- .cast_raster_writer(r, path, layers, datatype = "FLT8S", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  rows <- max(1L, as.integer(100000L/terra::ncol(r)))
  for (first in seq.int(1L, terra::nrow(r), by = rows)) {
    nr <- min(rows, terra::nrow(r)-first+1L)
    x <- as.data.frame(terra::readValues(r, row = first, nrows = nr, mat = TRUE), check.names = FALSE)
    common <- x$support_total == 1 & .cast_finite_rows(x[, contrasts, drop = FALSE])
    z <- as.matrix(x[, layers, drop = FALSE])
    z[which(is.na(common) | !common), seq_len(n_contrasts)] <- NA_real_
    for (v in names(diagnostics)) {
      ok <- x[[paste0("support_effect_", v)]] == 1 & x[[paste0("support_retained_", v)]] == 1
      z[which(is.na(ok) | !ok), diagnostics[[v]]] <- NA_real_
    }
    .cast_writer_block(writer, z, first, nr)
  }
  .cast_writer_finish(writer)
  .cast_figure_display(path, output, resolution, references = references)
}

.cast_figure_data <- function(object, comparison, scenario, directory,
    current = NULL, future = NULL, attribution = NULL, render_only = FALSE,
    response_grid = 50L, display_resolution = .02, resources = list(), verbose = TRUE,
    products = c("maps", "responses", "attribution", "mess"), response_method = "ale",
    response_variables = NULL, response_types = c(response_method, "fixed"),
    display_request = NULL) {
  fit <- object$fit_full %||% object$fit
  .cast_check_fit(fit)
  if (!render_only && any(products %in% c("responses", "attribution", "mess")) &&
      !identical(comparison$prediction_policy, "regular_gam/1"))
    cli::cli_abort("Recompute the comparison with the current prediction policy before preparing new response products. Do not combine earlier approximate GAM projections with regular predictions.")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  bundle_path <- file.path(directory, "figure_data.rds")
  provenance <- comparison$provenance
  provenance$session <- NULL
  if (is.null(response_variables)) response_variables <- fit$env_vars
  base_key <- .cast_digest(list(schema = "cast_figure_data/7", .cast_fit_identity(fit), provenance,
    comparison[c("groups", "design", "static", "domain", "clamp", "support_probs")],
    object[c("screen", "cv", "eval")], object$protocol$clamp,
    comparison$weights, comparison$threshold, comparison$threshold_sensitivity,
    scenario, response_grid, display_resolution, response_method,
    .cast_code_identity(cast_response_curves), .cast_code_identity(.cast_response_predictor), .cast_code_identity(.cast_mess),
    .cast_code_identity(.cast_figure_classes), .cast_code_identity(cast_diagnose),
    .cast_code_identity(.cast_response_validation_summary),
    .cast_code_identity(.cast_diagnostic_metrics), .cast_code_identity(compute_average_precision),
    .cast_code_identity(.cast_figure_display), .cast_code_identity(.cast_figure_driver_display),
    .cast_code_identity(.cast_display_coverage)))
  old <- .cast_cache_read(bundle_path, base_key)
  valid_old <- !is.null(old)
  required_products <- products
  if (!length(unlist(comparison$groups, use.names = FALSE)))
    required_products <- setdiff(required_products, "attribution")
  verify_raw <- !render_only && any(products %in% c("attribution", "mess")) &&
    !is.null(current) && !is.null(future)
  curves_ready <- !"responses" %in% products || (valid_old &&
    all(vapply(response_types, function(type) !is.null(old$responses[[type]]) &&
      all(response_variables %in% old$responses[[type]]$variable), logical(1))))
  if (valid_old && is.null(attribution) && !verify_raw && all(required_products %in% old$products) && curves_ready &&
      identical(old$display_request, display_request)) {
    if (verbose) cli::cli_inform("Verified figure data: no fitting or diagnostic prediction required.")
    return(old)
  }
  if (!"attribution" %in% products) attribution <- if (valid_old) old$attribution else NULL
  if (valid_old && is.null(attribution) && !verify_raw && "attribution" %in% products &&
      "attribution" %in% old$products)
    attribution <- old$attribution
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
    if (!is.null(value)) {
      if (verbose) cli::cli_inform("Verified numerical cache: {name}.")
      return(value)
    }
    value <- force(compute)
    .cast_cache_write(value, path, key, artifact)
    value
  }
  responses <- if (valid_old) old$responses else NULL
  if ("responses" %in% products && length(response_variables)) {
    response_key <- .cast_digest(list(.cast_fit_identity(fit), response_grid,
      response_variables, object$protocol$clamp, .cast_code_identity(cast_response_curves),
      .cast_code_identity(.cast_response_predictor)))
    if (verbose) cli::cli_inform("Response curves: {length(response_variables)} requested predictors; {nrow(fit$scaling$reference)} reference rows; {paste(toupper(response_types), collapse = ', ')}.")
    for (type in response_types) if (is.null(responses[[type]]) ||
        !all(response_variables %in% responses[[type]]$variable))
      responses[[type]] <- cached(paste0("responses_", type), .cast_digest(list(response_key, type)),
        cast_response_curves(fit, variables = response_variables, type = type,
          grid_size = response_grid, clamp = object$protocol$clamp %||% TRUE))
  }
  outputs <- if (valid_old) old$artifact_files else character()
  displays <- if (valid_old) old$displays else list()
  classes <- if (valid_old) old$range_files else list()
  notes <- list()
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
  prepared <- if (!is.null(current) && !is.null(future) &&
      any(products %in% c("attribution", "mess"))) .cast_compare_reuse(fit, current, future,
        comparison, comparison$weights, object$protocol$clamp %||% FALSE) else NULL
  if ("attribution" %in% products && length(unlist(comparison$groups, use.names = FALSE)) &&
      is.null(attribution) && !is.null(current) && !is.null(future)) {
    if (verbose) cli::cli_inform("Checking reporting-scenario response maps; a verified receipt skips numerical recomputation.")
    domain <- if (!is.null(report$raster_file)) terra::rast(report$raster_file)[["current"]] else NULL
    if (!is.null(domain)) domain <- terra::ifel(is.finite(domain), 1, NA)
    attribution <- .cast_attribution(fit, current, future,
      drivers = as.character(unlist(comparison$groups, use.names = FALSE)),
      weights = comparison$weights, mask = domain,
      clamp = object$protocol$clamp %||% FALSE, support_probs = comparison$support_probs %||% c(.01, .99),
      resources = resources, block_size = 250000L,
      overwrite = TRUE, diagnostics = TRUE, prepared = prepared,
      filename = if (!is.null(domain)) file.path(directory, "driver_attribution.tif"))
  }
  if (!is.null(attribution$raster_file) && !is.null(current) && !is.null(future)) {
    contract <- list(fit = .cast_fit_identity(fit), scenario = scenario,
      weights = comparison$weights, clamp = comparison$clamp,
      support_probs = comparison$support_probs,
      current = prepared$current_key %||% .cast_data_identity(current),
      future = prepared$future_key %||% .cast_data_identity(future),
      comparison = unname(tools::md5sum(report$raster_file)))
    .cast_cache_write(list(contract = contract, attribution = attribution),
      file.path(directory, "driver_catalog.rds"), .cast_digest(contract), attribution$raster_file)
    outputs <- c(outputs, file.path(directory, "driver_catalog.rds"))
  }
  if ("attribution" %in% products && is.null(attribution)) note("03_environmental_driver",
    if (!length(unlist(comparison$groups, use.names = FALSE))) "No fitted predictor changes; scenario driver maps are inapplicable." else
      "Current/future predictor inputs or an existing cast_attribution object are required.")
  mess <- if (valid_old) old$mess_native else NULL
  if ("mess" %in% products && !is.null(current) && !is.null(future) && !is.null(report$raster_file)) {
    path <- file.path(directory, "mess_native.tif")
    memo <- new.env(parent = emptyenv())
    key <- .cast_digest(list(.cast_fit_identity(fit), prepared$current_key %||% .cast_data_identity(current, memo),
      prepared$future_key %||% .cast_data_identity(future, memo), provenance, scenario,
      .cast_code_identity(.cast_figure_mess), .cast_code_identity(.cast_mess_apply),
      .cast_code_identity(.cast_mess_reference)))
    mess <- cached("mess", key, .cast_figure_mess(fit, current, future, report$raster_file, path, prepared = prepared), path)
    outputs <- c(outputs, mess)
  } else if ("mess" %in% products) note("09_extrapolation_mess", "Raw current/future predictor stacks are required for native-grid MESS; support flags cannot substitute for MESS.")
  shared_display <- if (identical(comparison$domain, "common")) paste0("shared_",
    .cast_digest(list(provenance, comparison$weights, comparison$clamp,
      comparison$support_probs, display_resolution)))
  for (s in if ("maps" %in% products) names(comparison$results) else character()) {
    result <- comparison$results[[s]]
    artifact_id <- result$alias_of %||% s
    if (is.null(result$raster_file)) next
    if (verbose) cli::cli_inform("Display grid {s}: checking cache (scientific raster unchanged).")
    input <- result$raster_file
    requested <- if (s == scenario) display_request$report_layers else display_request$other_layers
    if (!is.null(requested) && !length(requested)) next
    if (!is.null(comparison$threshold) && (is.null(display_request) || isTRUE(display_request$range))) {
      path <- file.path(directory, paste0("range_", artifact_id, ".tif"))
      key <- .cast_digest(list(unname(tools::md5sum(input)), comparison$threshold, .cast_code_identity(.cast_figure_classes)))
      classes[[s]] <- cached(paste0("range_", artifact_id), key,
        .cast_figure_classes(result, comparison$threshold, path), path)
      outputs <- c(outputs, classes[[s]]); input <- c(input, classes[[s]])
    }
    path <- file.path(directory, paste0("display_", artifact_id, ".tif"))
    key <- .cast_digest(list(unname(tools::md5sum(input)),
      display_resolution, requested, .cast_code_identity(.cast_figure_display), .cast_code_identity(.cast_display_coverage)))
    displays[[s]] <- cached(paste0("display_", artifact_id), key,
      .cast_figure_display(input, path, display_resolution, layers = requested, shared_key = shared_display),
      c(path, paste0(path, ".csv")))
    outputs <- c(outputs, displays[[s]])
    invisible(gc(FALSE))
  }
  mess_display <- if (valid_old) old$mess_display else NULL
  if ("mess" %in% products && !is.null(mess)) {
    path <- file.path(directory, "display_mess.tif")
    key <- .cast_digest(list(unname(tools::md5sum(mess)), display_resolution,
      .cast_code_identity(.cast_figure_display), .cast_code_identity(.cast_display_coverage)))
    mess_display <- cached("display_mess", key,
      .cast_figure_display(mess, path, display_resolution), c(path, paste0(path, ".csv")))
    outputs <- c(outputs, mess_display)
  }
  attribution_display <- if (valid_old) old$attribution_display else NULL
  attribution_supported_display <- if (valid_old) old$attribution_supported_display else NULL
  focus_driver <- NULL
  if ("attribution" %in% products && !is.null(attribution$raster_file)) {
    outputs <- c(outputs, attribution$raster_file)
    drivers <- attribution$drivers
    focus_only <- identical(display_request$driver_mode, "focus")
    if (length(drivers)) {
      z <- attribution$common_summary %||% attribution$summary
      z <- z[z$contrast == "only", , drop = FALSE]
      focus_driver <- if (any(is.finite(z$mean_abs) & z$mean_abs > 0))
        z$driver[which.max(replace(z$mean_abs, !is.finite(z$mean_abs), -Inf))] else drivers[1L]
      if (focus_only) drivers <- focus_driver
    }
    layers <- c(paste0("effect_", drivers), if (!focus_only) paste0("retained_", drivers))
    path <- file.path(directory, "display_drivers.tif")
    key <- .cast_digest(list(unname(tools::md5sum(attribution$raster_file)),
      layers, display_resolution, .cast_code_identity(.cast_figure_display), .cast_code_identity(.cast_display_coverage)))
    attribution_display <- cached("display_drivers", key,
      .cast_figure_display(attribution$raster_file, path, display_resolution, layers = layers),
      c(path, paste0(path, ".csv")))
    outputs <- c(outputs, attribution_display)
    if (!is.null(attribution$common_summary)) {
      path <- file.path(directory, "display_drivers_supported.tif")
      key <- .cast_digest(list(key, .cast_code_identity(.cast_figure_driver_display)))
      attribution_supported_display <- cached("display_drivers_supported", key,
        .cast_figure_driver_display(attribution, path, display_resolution, drivers = drivers,
          include_diagnostics = !focus_only), c(path, paste0(path, ".csv")))
      outputs <- c(outputs, attribution_supported_display)
      if (verbose) {
        coverage <- utils::read.csv(paste0(path, ".csv"))[1L, ]
        cli::cli_inform("Native-supported driver display: {coverage$native_valid_cells} cells; {round(100*coverage$native_centres_represented_fraction, 2)}% of native centres represented. Support masks precede display aggregation.")
      }
    }
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
  if (is.null(responses)) completed <- setdiff(completed, "responses")
  if (!length(displays)) completed <- setdiff(completed, "maps")
  if (valid_old) completed <- union(old$products, completed)
  coverage_files <- paste0(c(unlist(displays, use.names = FALSE), attribution_display,
    attribution_supported_display, mess_display), ".csv")
  coverage_files <- unique(coverage_files[file.exists(coverage_files)])
  display_coverage <- if (length(coverage_files)) do.call(rbind, lapply(coverage_files, function(p) {
    value <- utils::read.csv(p); value$product <- basename(sub("[.]csv$", "", p)); value
  }))
  outputs <- c(outputs, coverage_files)
  data <- list(schema = "cast_figure_data/7", display_coverage = display_coverage, analysis_key = base_key, display_request = display_request, products = completed, screen = object$screen, cv = object$cv,
    evaluation = ev, diagnostics = diagnostic, folds = info, responses = responses,
    curves = curves, reliability = if (length(reliability)) do.call(rbind, reliability),
    comparison = comparison, attribution = attribution, focus_driver = focus_driver %||% old$focus_driver,
    attribution_display = attribution_display,
    attribution_supported_display = attribution_supported_display,
    displays = displays, mess_native = mess, mess_display = mess_display, range_files = classes, thresholds = thresholds,
    scenario = scenario, variables = fit$env_vars, model_names = names(fit$models),
    response_reference_n = nrow(fit$scaling$reference),
    notes = if (length(notes)) do.call(rbind, notes), display_resolution = display_resolution,
    source_version = as.character(utils::packageVersion("castSDM")),
    training_version = origin$cast_version %||% NA_character_)
  # Receipt file verification is restricted to sibling products, by design.
  sibling <- outputs[file.exists(outputs) & dirname(outputs) == directory]
  data$artifact_files <- unique(sibling)
  .cast_cache_write(data, bundle_path, base_key, unique(sibling))
  data
}

# Public presets and individual names share the same dependency resolver.
.cast_figure_choose <- function(include, inventory) {
  if (!is.character(include) || !length(include) || anyNA(include))
    cli::cli_abort("include must name figure presets or individual figures.")
  unique(unlist(lapply(include, function(x) {
    if (x == "core") return(intersect(inventory, c("01_model_validation", "02_habitat_change",
      "03_environmental_driver", "04_support_and_disagreement")))
    if (x == "all") return(inventory)
    pattern <- switch(x,
      basic = "^(01_|07_|11_|12_|13_|15_|S_selection_|S_stability_)",
      maps = "^(02_|04_|05_|08_|10_|S_(?!selection_|stability_))",
      responses = "^(06_|14_|R_)", drivers = "^(03_|D_|E_)", novelty = "^09_", NULL)
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
    if (x == "core") out <- c(out, "maps", "attribution")
    if (x %in% c("maps", "novelty") || grepl("^(02|04|05|08|09|10)(_|$)", x) ||
        (startsWith(x, "S_") && !grepl("^S_(selection|stability)", x)) || x == "S") out <- c(out, "maps")
    if (x == "responses" || startsWith(x, "R_") || grepl("^(06|14)(_|$)", x)) out <- c(out, "responses")
    if (x == "drivers" || grepl("^03(_|$)|^[DE]_", x)) out <- c(out, "attribution")
    if (x == "novelty" || grepl("^09(_|$)", x)) out <- c(out, "mess")
  }
  unique(out)
}

# Figure selection controls computation; it never changes fitted predictors.
.cast_figure_response_request <- function(include, comparison, variables, method) {
  complete <- any(include %in% c("all", "responses") | startsWith(include, "R_"))
  vars <- if (complete) variables else intersect(unlist(comparison$groups, use.names = FALSE), variables)
  both <- any(include %in% c("all", "responses"))
  types <- c(if (both || any(grepl("^06(_|$)|^R_(ale|pdp)", include))) method,
    if (both || any(grepl("^14(_|$)|^R_fixed", include))) "fixed")
  list(variables = vars, types = unique(types))
}

# Core maps use the reporting scenario. Full atlases remain explicit.
.cast_figure_display_request <- function(include) {
  if (any(include %in% c("all", "maps", "S", "drivers")) ||
      any(grepl("^[SDE]_", include))) return(NULL)
  has <- function(id) any(grepl(paste0("^", id, "(_|$)"), include))
  report <- character()
  if ("core" %in% include || has("02")) report <- c(report, "current", "joint", "joint_change")
  if ("core" %in% include || has("04")) report <- c(report, "support", "sd_joint_change")
  if (has("05")) report <- c(report, "climate_main", "lulc_main", "interaction")
  if (has("08")) report <- c(report, "range_class", "range_class_supported", "range_class_common_support")
  future <- if (has("10")) "joint" else character()
  list(report_layers = unique(c(report, future)), other_layers = future,
    range = has("08"), driver_mode = "focus")
}
