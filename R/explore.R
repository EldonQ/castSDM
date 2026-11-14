# AI-assisted development (2026-11-11); see inst/AI_ASSISTANCE.md.
# File paths and ordinary parameters cross process boundaries; rasters do not.

.cast_explore_project <- function(directory) {
  if (!is.character(directory) || length(directory) != 1L || is.na(directory))
    cli::cli_abort("Supply a completed single-species case directory.")
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  path <- file.path(directory, "case.rds")
  if (!file.exists(path)) cli::cli_abort("No case.rds in {.path {directory}}. Run the configured workflow first.")
  if (!file.exists(file.path(directory, "completion.rds")))
    cli::cli_abort("This case is incomplete or is being updated. Finish the research workflow first.")
  case <- readRDS(path)
  if (!identical(case$schema, "cast_single_species/1"))
    cli::cli_abort("Unsupported case schema.")
  run <- readRDS(file.path(directory, case$training))
  comparison <- readRDS(file.path(directory, case$comparison))
  fit <- run$fit_full %||% run$fit
  .cast_check_fit(fit)
  pool <- utils::read.csv(file.path(directory, "training_pool.csv"), stringsAsFactors = FALSE)
  selected <- pool[match(fit$env_vars, pool$layer), , drop = FALSE]
  if (anyNA(selected$path)) cli::cli_abort("The case inventory lacks fitted predictors.")
  metadata <- data.frame(layer = fit$env_vars, product = selected$product,
    definition = NA_character_, units = NA_character_, source = NA_character_, notes = NA_character_, label = fit$env_vars)
  dictionaries <- c(file.path(directory, "predictor_metadata.csv"),
    file.path(case$cfg$example_root, "config", "predictor_metadata.csv"))
  dictionaries <- dictionaries[file.exists(dictionaries)]
  if (length(dictionaries)) {
    dictionary <- utils::read.csv(dictionaries[1L], stringsAsFactors = FALSE, encoding = "UTF-8")
    required <- c("layer", "product", "definition", "units", "source", "notes")
    if (!all(required %in% names(dictionary)) ||
        anyDuplicated(paste(dictionary$product, dictionary$layer)))
      cli::cli_abort("Predictor metadata needs unique product/layer keys and definition, units, source, notes.")
    matched <- match(paste(metadata$product, metadata$layer), paste(dictionary$product, dictionary$layer))
    found <- !is.na(matched)
    metadata[found, required] <- dictionary[matched[found], required]
    if ("label" %in% names(dictionary)) {
      labelled <- found & !is.na(dictionary$label[matched]) & nzchar(dictionary$label[matched])
      metadata$label[labelled] <- dictionary$label[matched[labelled]]
    }
  }
  scenarios <- utils::read.csv(file.path(directory, "scenario_table.csv"), stringsAsFactors = FALSE)
  comparison$results <- lapply(comparison$results, function(x) {
    if (!is.null(x$artifact)) x$raster_file <- file.path(directory, "comparison", "rasters", x$artifact)
    x
  })
  list(directory = directory, case = case, run = run, fit = fit,
    comparison = comparison, pool = selected, scenarios = scenarios, metadata = metadata)
}

.cast_explore_stacks <- function(project, scenario) {
  paths <- project$pool$path
  vars <- project$fit$env_vars
  tab <- project$scenarios[project$scenarios$scenario == scenario, , drop = FALSE]
  if (!nrow(tab)) cli::cli_abort("No predictor mapping for scenario {.val {scenario}}.")
  future <- paths
  for (i in seq_len(nrow(tab))) {
    k <- match(tab$layer[i], vars)
    if (!is.na(k)) future[k] <- tab$path[i]
  }
  missing <- unique(c(paths, future)[!file.exists(c(paths, future))])
  if (length(missing)) cli::cli_abort("Missing environmental sources: {.path {missing}}.")
  a <- terra::rast(paths); b <- terra::rast(future)
  names(a) <- names(b) <- vars
  terra::compareGeom(a, b, lyrs = FALSE, stopOnError = TRUE)
  list(current = a, future = b)
}

.cast_explore_selection <- function(project, variables) {
  vars <- project$fit$env_vars
  if (is.null(variables)) return(vars)
  if (!is.character(variables) || !length(variables) || anyNA(variables) ||
      anyDuplicated(variables) || !all(variables %in% vars))
    cli::cli_abort("{.arg variables} must name unique fitted predictors.")
  variables
}

.cast_explore_sample <- function(project, a, b, n) {
  info <- project$run$cv$fold_info
  if (is.null(info) || !all(c("lon", "lat", "fold") %in% names(info)))
    cli::cli_abort("Response groups require the saved spatial CV fold_info (lon, lat, fold).")
  # Evenly spaced row indices are a deterministic subset of existing records,
  # not new observations or manufactured absence labels.
  at <- unique(as.integer(round(seq(1L, nrow(info), length.out = min(n, nrow(info))))))
  xy <- as.matrix(info[at, c("lon", "lat"), drop = FALSE])
  # The case workflow records these coordinates in the environmental grid CRS.
  current <- terra::extract(a, xy)
  future <- terra::extract(b, xy)
  list(current = current, future = future, xy = xy, fold = info$fold[at], row = at)
}

.cast_explore_group_map <- function(groups, predictors, support, filename, block_size = 50000L) {
  writer <- .cast_raster_writer(predictors, filename, "response_group", datatype = "INT4S", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  rows <- max(1L, as.integer(block_size / terra::ncol(predictors)))
  for (first in seq.int(1L, terra::nrow(predictors), by = rows)) {
    count <- min(rows, terra::nrow(predictors) - first + 1L)
    x <- .cast_raster_block(predictors, first, count)
    ok <- .cast_finite_rows(x) & .cast_mask_rows(support, first, count, nrow(x))
    z <- rep(NA_integer_, nrow(x))
    if (any(ok)) z[ok] <- stats::predict(groups, x[ok, , drop = FALSE])
    .cast_writer_block(writer, matrix(z, ncol = 1L), first, count)
  }
  .cast_writer_finish(writer)
  filename
}

.cast_explore_preview <- function(path, output, pixels = 250000L) {
  r <- if (inherits(path, "SpatRaster")) path else terra::rast(path)
  coarse <- terra::spatSample(r, size = min(pixels, terra::ncell(r)),
    method = "regular", as.raster = TRUE)
  categorical <- grepl("^support|^response_group$|^context_direction_|^(effect|retained)_agreement_", names(coarse))
  template <- leaflet::projectRasterForLeaflet(coarse[[1L]], method = "bilinear")
  out <- lapply(seq_len(terra::nlyr(coarse)), function(i)
    terra::project(coarse[[i]], template, method = if (categorical[i]) "near" else "bilinear"))
  display <- do.call(c, out); names(display) <- names(r)
  terra::writeRaster(display, output, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES")))
  output
}

.cast_explore_pair <- function(project, scenario, variables, shifts, context) {
  stacks <- .cast_explore_stacks(project, scenario)
  # Choosing a factor changes the explanation, not the imported joint scenario.
  if (is.null(shifts)) return(stacks)
  a <- b <- if (context == "current") stacks$current else stacks$future
  for (v in variables) {
    a[[v]] <- stacks$current[[v]]
    b[[v]] <- if (is.null(shifts)) stacks$future[[v]] else stacks$current[[v]]
  }
  list(current = a, future = b)
}

#' Explore Environmental Responses in a Completed Case
#'
#' Reuses frozen models, recorded spatial folds and actual environmental
#' rasters from the single-species workflow. Computes paired scenario or
#' explicitly hypothetical finite-change responses and descriptive groups.
#'
#' @param case_dir Directory containing the completed workflow's `case.rds`,
#'   training result, current inventory and scenario table.
#' @param scenario Saved scenario identifier. Default the case reporting scenario.
#' @param variables Fitted predictors to explain. Scenario mode preserves the
#'   complete imported current/future pair; selecting one factor does not
#'   freeze the other scenario changes. `NULL` uses all predictors.
#'   Related derived predictors should be changed together using a consistent
#'   scenario. A single derived-variable change is a sensitivity analysis.
#' @param shifts Optional named numeric vector of raw-unit predictor shifts.
#'   When supplied, replaces the imported future change for those variables
#'   with these hypothetical changes. All other selected variables must be
#'   omitted: `variables` is inferred from `names(shifts)`.
#' @param context For hypothetical shifts only, `"current"` keeps other predictors current. `"future"` keeps
#'   other predictors at the selected future values in both members of the pair.
#' @param sample_n Maximum number of saved CV locations used for descriptive
#'   grouping. Default 5000, selected deterministically from real saved rows.
#' @param grouping Optional descriptive tree grouping, default `FALSE`.
#' @param maxdepth Maximum response tree depth. See [cast_response_groups()].
#' @param maps Compute native-grid response/group rasters. `FALSE` only computes
#'   the real-location response table and grouping; it never predicts a full grid.
#' @param focus Optional selected driver for detailed maps. `NULL` automatically
#'   uses the largest mean absolute driver-only response on the shared support
#'   domain. Only this driver's two background responses, model sign agreement
#'   and background diagnostics are mapped, avoiding full maps for every driver.
#' @param output_dir Separate directory for exploration artifacts. Default
#'   `interactive/` within the case. Research fits and existing figures are not overwritten.
#' @param resources Resource budget, as in [cast_attribution()].
#' @param resume Reuse content-verified exploration artifacts.
#' @return A list of artifact paths, frozen settings, groups and interpretation.
#'   Existing spatial folds tune a descriptive response tree, not an ecological
#'   causal effect. Raster previews are separate from native-grid outputs.
#' @export
cast_explore <- function(case_dir, scenario = NULL, variables = NULL, shifts = NULL,
                         context = c("current", "future"), sample_n = 5000L,
                         maxdepth = 3L, maps = TRUE, output_dir = NULL,
                         resources = list(), resume = TRUE, focus = NULL, grouping = FALSE) {
  check_suggested("terra", "for saved-case exploration")
  if (!is.logical(grouping) || length(grouping) != 1L || is.na(grouping)) cli::cli_abort("grouping must be TRUE or FALSE.")
  if (grouping) check_suggested("rpart", "for descriptive response grouping")
  if (!is.logical(maps) || length(maps) != 1L || is.na(maps))
    cli::cli_abort("{.arg maps} must be TRUE or FALSE.")
  if (maps) check_suggested("leaflet", "for interactive display rasters")
  context <- match.arg(context)
  .cast_check_scalar(sample_n, "sample_n", lower = 4, integer = TRUE)
  .cast_check_scalar(maxdepth, "maxdepth", lower = 1, upper = 30, integer = TRUE)
  case_lock <- .cast_lock(file.path(case_dir, ".case.lock"))
  on.exit(filelock::unlock(case_lock), add = TRUE)
  project <- .cast_explore_project(case_dir)
  scenario <- scenario %||% project$case$report_scenario
  if (length(scenario) != 1L || is.na(scenario) || !scenario %in% names(project$comparison$results))
    cli::cli_abort("Select one saved scenario.")
  if (!is.null(shifts)) {
    if (!is.numeric(shifts) || !length(shifts) || is.null(names(shifts)) ||
        anyNA(names(shifts)) || anyDuplicated(names(shifts)) || any(!is.finite(shifts)) ||
        any(shifts == 0) || !all(names(shifts) %in% project$fit$env_vars))
      cli::cli_abort("{.arg shifts} must name fitted predictors with finite non-zero raw-unit changes.")
    if (!is.null(variables) && !setequal(variables, names(shifts)))
      cli::cli_abort("With shifts, variables must equal names(shifts), or be NULL.")
    variables <- names(shifts)
  }
  variables <- .cast_explore_selection(project, variables)
  requested_focus <- focus
  if (!is.null(focus) && (!is.character(focus) || length(focus) != 1L ||
      is.na(focus) || !focus %in% variables))
    cli::cli_abort("{.arg focus} must name one of the compared variables.")
  pair <- .cast_explore_pair(project, scenario, variables, shifts, context)
  a <- pair$current; b <- pair$future
  # Hash in the worker, with a shared memo so repeated/static files are read once.
  memo <- new.env(parent = emptyenv())
  invisible(.cast_resources(resources))
  settings <- list(scenario = scenario, variables = variables, shifts = shifts,
    context = context, sample_n = sample_n, maxdepth = maxdepth, maps = maps, focus = focus, grouping = grouping,
    weights = project$comparison$weights,
    clamp = project$run$protocol$clamp %||% FALSE,
    support_probs = project$comparison$support_probs,
    group_population = "Saved training-CV presence/background locations; not an area-weighted ecological population",
    interpretation = "Paired fitted-suitability response; not ecological causal effects or occurrence probability change")
  mask_file <- file.path(project$case$preparation, "analysis_domain_mask.tif")
  if (!file.exists(mask_file)) cli::cli_abort("Missing recorded analysis domain: {.path {mask_file}}.")
  mask <- terra::rast(mask_file)
  key <- .cast_digest(list(schema = "cast_explore/1", settings = settings,
    fit = .cast_fit_identity(project$fit),
    current = .cast_data_identity(a, memo), future = .cast_data_identity(b, memo),
    mask = .cast_data_identity(mask, memo), folds = project$run$cv$fold_info,
    metadata = project$metadata,
    grouping = .cast_code_identity(cast_response_groups),
    membership = .cast_code_identity(.cast_response_nodes),
    evidence = .cast_code_identity(.cast_response_evidence),
    direction = .cast_code_identity(.cast_response_direction),
    compute = .cast_code_identity(cast_explore), attribution = .cast_code_identity(cast_attribution),
    sample = .cast_code_identity(.cast_explore_sample), pair = .cast_code_identity(.cast_explore_pair),
    group_map = .cast_code_identity(.cast_explore_group_map),
    preview = .cast_code_identity(.cast_explore_preview),
    rpart = as.character(utils::packageVersion("rpart"))))
  root <- output_dir %||% file.path(project$directory, "interactive")
  directory <- file.path(root, key)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  receipt <- file.path(directory, "exploration.rds")
  lock <- .cast_lock(paste0(receipt, ".lock"))
  on.exit(filelock::unlock(lock), add = TRUE)
  if (resume) {
    cached <- .cast_cache_read(receipt, key)
    if (!is.null(cached)) {
      for (field in c("raster_file", "group_file", "preview", "points_file"))
        if (!is.null(cached[[field]])) cached[[field]] <- file.path(directory, basename(cached[[field]]))
      cached$directory <- directory; cached$case_dir <- project$directory
      return(cached)
    }
  }
  cli::cli_inform("Reading real saved CV locations and aligned scenario values.")
  sample <- .cast_explore_sample(project, a, b, sample_n)
  response <- cast_attribution(project$fit, sample$current, sample$future,
    drivers = variables, weights = settings$weights, shifts = shifts,
    support_probs = settings$support_probs, clamp = settings$clamp,
    resources = resources, verbose = FALSE, diagnostics = TRUE)
  actual_future <- sample$future
  for (v in names(shifts)) actual_future[[v]] <- actual_future[[v]] + shifts[[v]]
  evidence <- .cast_response_evidence(response$data, variables, sample$current, actual_future)
  ranked <- evidence$summary$driver[!is.na(evidence$summary$rank) & evidence$summary$rank == 1L]
  if (is.null(focus) && length(ranked) == 1L) focus <- ranked
  settings$focus <- focus
  groups <- NULL; group_note <- "Optional descriptive grouping was not requested."
  if (grouping) {
    cli::cli_inform("Selecting optional descriptive response groups using saved spatial folds.")
    supported <- response$data$support_total == 1
    groups <- tryCatch(cast_response_groups(sample$current, response$data$total_change,
      sample$fold, supported = supported, maxdepth = maxdepth), error = function(e) e)
    group_note <- if (inherits(groups, "error")) conditionMessage(groups) else NULL
    if (!is.null(group_note)) groups <- NULL
  }
  points <- data.frame(source_cv_row = sample$row, x = sample$xy[, 1L], y = sample$xy[, 2L],
    fold = sample$fold, response$data, evidence$locations,
    group = if (is.null(groups)) NA_integer_ else groups$group)
  utils::write.csv(points, file.path(directory, "responses.csv"), row.names = FALSE)
  utils::write.csv(response$summary, file.path(directory, "summary.csv"), row.names = FALSE)
  utils::write.csv(evidence$summary, file.path(directory, "drivers_common.csv"), row.names = FALSE)
  utils::write.csv(evidence$population, file.path(directory, "support_population.csv"), row.names = FALSE)
  utils::write.csv(evidence$drivers, file.path(directory, "driver_ids.csv"), row.names = FALSE)
  utils::write.csv(project$metadata, file.path(directory, "predictor_metadata.csv"), row.names = FALSE)
  if (!is.null(groups)) {
    utils::write.csv(groups$summary, file.path(directory, "groups.csv"), row.names = FALSE)
    utils::write.csv(groups$cv, file.path(directory, "group_tuning.csv"), row.names = FALSE)
    utils::write.csv(groups$fold_losses, file.path(directory, "group_fold_losses.csv"), row.names = FALSE)
  }
  raw_file <- preview <- group_file <- NULL
  if (maps) {
    cli::cli_inform("Computing native-grid paired responses in bounded blocks.")
    raw_file <- file.path(directory, "response.tif")
    cast_attribution(project$fit, a, b, drivers = focus %||% character(), weights = settings$weights,
      shifts = shifts, support_probs = settings$support_probs, clamp = settings$clamp,
      mask = mask, filename = raw_file, overwrite = TRUE, resources = resources,
      resume = resume, diagnostics = TRUE)
    if (!is.null(groups)) group_file <- .cast_explore_group_map(groups, a,
      terra::rast(raw_file)[["support_total"]], file.path(directory, "groups.tif"))
    cli::cli_inform("Preparing a separate interactive preview (not the analysis grid).")
    preview_source <- terra::rast(raw_file)
    if (!is.null(group_file)) preview_source <- c(preview_source, terra::rast(group_file))
    # The preview helper accepts raster objects without crossing process boundaries.
    preview <- .cast_explore_preview(preview_source, file.path(directory, "preview.tif"))
  }
  result <- list(schema = "cast_explore/1", directory = directory, key = key,
    case_dir = project$directory, settings = settings, summary = response$summary,
    evidence = evidence,
    groups = groups, group_note = group_note, raster_file = raw_file,
    group_file = group_file, preview = preview,
    points_file = file.path(directory, "responses.csv"),
    package_version = as.character(utils::packageVersion("castSDM")))
  .cast_atomic_rds(settings, file.path(directory, "settings.rds"))
  replay <- c("# Uses the original real case; hypothetical shifts remain explicitly recorded.",
    "library(castSDM)", "args <- ", utils::capture.output(dput(list(case_dir = project$directory,
      scenario = scenario, variables = variables, shifts = shifts, context = context,
      sample_n = sample_n, maxdepth = maxdepth, maps = maps, output_dir = root, grouping = grouping,
      resources = resources, focus = requested_focus))), "result <- do.call(cast_explore, args)")
  writeLines(replay, file.path(directory, "reproduce.R"))
  files <- list.files(directory, full.names = TRUE)
  files <- files[!grepl("\\.lock$|exploration\\.rds$", files)]
  .cast_cache_write(result, receipt, key, files)
  result
}

.cast_app_view <- function(project, scenario) {
  source <- project$comparison$results[[scenario]]$raster_file
  if (is.null(source) || !file.exists(source)) cli::cli_abort("Saved scenario raster is missing.")
  r <- terra::rast(source)
  r <- r[[c("current", "joint", "joint_change", "support", "sd_joint_change")]]
  names(r) <- c("current", "future", "total_change", "support_total", "total_sd")
  key <- .cast_digest(list(scenario = scenario, source = unname(tools::md5sum(source)),
    fit = .cast_fit_identity(project$fit), weights = project$comparison$weights,
    clamp = project$run$protocol$clamp %||% FALSE, support = project$comparison$support_probs,
    .cast_code_identity(.cast_explore_preview), .cast_code_identity(.cast_app_view),
    terra = as.character(utils::packageVersion("terra")),
    leaflet = as.character(utils::packageVersion("leaflet")), "cast_app_view/1"))
  directory <- file.path(project$directory, "interactive", paste0("view_", key))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  lock <- .cast_lock(file.path(directory, ".lock"))
  on.exit(filelock::unlock(lock), add = TRUE)
  receipt <- file.path(directory, "view.rds")
  cached <- .cast_cache_read(receipt, key)
  if (!is.null(cached)) {
    cached$directory <- directory; cached$case_dir <- project$directory
    cached$preview <- file.path(directory, "preview.tif")
    cached$raster_file <- source
    return(cached)
  }
  preview <- .cast_explore_preview(r, file.path(directory, "preview.tif"))
  result <- list(schema = "cast_app_view/1", directory = directory,
    case_dir = project$directory, preview = preview, raster_file = source,
    settings = list(scenario = scenario, context = "current", variables = project$fit$env_vars,
      shifts = NULL, weights = project$comparison$weights,
      clamp = project$run$protocol$clamp %||% FALSE,
      interpretation = "Saved joint scenario; paired fitted suitability"),
    groups = NULL, group_note = "Run response grouping to calculate descriptive subgroups.")
  .cast_cache_write(result, receipt, key, preview)
  result
}

.cast_app_task <- function(request) {
  if (identical(request$operation, "setup")) {
    request$operation <- NULL
    return(do.call(cast_workflow, request))
  }
  if (identical(request$operation, "view")) {
    lock <- .cast_lock(file.path(request$case_dir, ".case.lock"))
    on.exit(filelock::unlock(lock), add = TRUE)
    return(.cast_app_view(.cast_explore_project(request$case_dir), request$scenario))
  }
  request$operation <- NULL
  do.call(cast_explore, request)
}
