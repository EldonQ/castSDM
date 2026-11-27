# AI-assisted development (2026-11-11); see inst/AI_ASSISTANCE.md.
# File paths and ordinary parameters cross process boundaries; rasters do not.

.cast_completed_case <- function(directory) {
  if (!is.character(directory) || length(directory) != 1L || is.na(directory))
    cli::cli_abort("Supply a completed single-species case directory.")
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  path <- file.path(directory, "case.rds")
  if (!file.exists(path)) cli::cli_abort("No case.rds in {.path {directory}}. Run the configured workflow first.")
  marker <- file.path(directory, "completion.rds")
  if (!file.exists(marker))
    cli::cli_abort("This case is incomplete or is being updated. Finish the research workflow first.")
  case <- readRDS(path)
  if (!identical(case$schema, "cast_single_species/1"))
    cli::cli_abort("Unsupported case schema.")
  completed <- readRDS(marker)
  if (!is.null(case$analysis_key) || !is.null(completed$analysis_key)) {
    claim <- readRDS(file.path(directory, "run_claim.rds"))
    if (is.null(case$analysis_key) || !identical(completed$analysis_key, case$analysis_key) ||
        !identical(claim$key, case$analysis_key) ||
        !identical(completed$case_md5, unname(tools::md5sum(path))))
      cli::cli_abort("The saved case and completion record belong to different analysis inputs.")
  }
  case
}

.cast_explore_project <- function(directory) {
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  case <- .cast_completed_case(directory)
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
  comparison$output_dir <- file.path(directory, "comparison")
  list(directory = directory, case = case, run = run, fit = fit,
    comparison = comparison, pool = selected, scenarios = scenarios, metadata = metadata)
}

.cast_explore_stacks <- function(project, scenario, memo = new.env(parent = emptyenv())) {
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
  proven <- project$comparison$provenance
  if (!identical(.cast_data_identity(a, memo), proven$current) ||
      !identical(.cast_data_identity(b, memo), proven$future[[scenario]]))
    cli::cli_abort("Environmental inputs have changed since the saved comparison. Use a completed case with matching files to compute a response.")
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

.cast_paired_extract <- function(current, future, xy) {
  before <- terra::extract(current, xy)
  if (!identical(names(current), names(future)))
    return(list(current = before, future = terra::extract(future, xy)))
  vars <- names(current)
  changed <- vars[!vapply(vars, function(v)
    .cast_shared_raster_layers(current, future, v), logical(1))]
  after <- before
  if (length(changed)) after[changed] <- terra::extract(future[[changed]], xy)
  list(current = before, future = after)
}

.cast_explore_sample <- function(project, a, b, n) {
  info <- project$run$cv$fold_info
  if (is.null(info) || !all(c("lon", "lat", "fold") %in% names(info)))
    cli::cli_abort("Response groups require the saved spatial CV fold_info (lon, lat, fold).")
  # Evenly spaced row indices are a deterministic subset of existing records,
  # not new observations or manufactured absence labels.
  at <- if (is.null(n)) seq_len(nrow(info)) else
    unique(as.integer(round(seq(1L, nrow(info), length.out = min(n, nrow(info))))))
  xy <- as.matrix(info[at, c("lon", "lat"), drop = FALSE])
  # The case workflow records these coordinates in the environmental grid CRS.
  c(.cast_paired_extract(a, b, xy), list(xy = xy, fold = info$fold[at], row = at))
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
  # A geometric template, not regular-cell sampling: narrow valid features
  # contribute to their overlapping display cells instead of being missed.
  template <- terra::project(terra::rast(r[[1L]]), "EPSG:3857")
  ratio <- max(1, sqrt(terra::ncell(template)/pixels))
  template <- terra::rast(nrows = max(1, floor(terra::nrow(template)/ratio)),
    ncols = max(1, floor(terra::ncol(template)/ratio)),
    ext = terra::ext(template), crs = terra::crs(template))
  fields <- names(r)
  flags <- grepl("^support_", fields)
  categorical <- grepl("^response_group$|^context_direction_|^sign_conflict", fields)
  eligible <- fields[!flags & fields != "response_group"]
  support_for <- function(v) {
    if (grepl("^context_(gap|direction)_", v)) {
      driver <- sub("^context_(gap|direction)_", "", v)
      return(paste0(c("support_effect_", "support_retained_"), driver))
    }
    for (prefix in c("effect_", "retained_")) if (startsWith(v, prefix)) {
      driver <- sub("^(effect|retained)_(agreement_|sd_)?", "", v)
      return(paste0("support_", prefix, driver))
    }
    "support_total"
  }
  masks <- stats::setNames(lapply(eligible, support_for), eligible)
  masks <- masks[vapply(masks, function(x) all(x %in% fields), logical(1))]
  source <- tempfile("cast-preview-native-", fileext = ".tif")
  on.exit(unlink(source), add = TRUE)
  layers <- c(fields, if (length(masks)) paste0("supported_", names(masks)))
  writer <- .cast_raster_writer(r, source, layers, datatype = "FLT4S", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  rows <- max(1L, as.integer(250000L/terra::ncol(r)))
  terra::readStart(r); on.exit(terra::readStop(r), add = TRUE)
  for (first in seq.int(1L, terra::nrow(r), by = rows)) {
    count <- min(rows, terra::nrow(r)-first+1L)
    x <- terra::readValues(r, row = first, nrows = count, mat = TRUE)
    z <- matrix(NA_real_, nrow(x), length(layers), dimnames = list(NULL, layers))
    z[, fields] <- x
    for (v in names(masks)) {
      support <- x[, masks[[v]], drop = FALSE]
      keep <- rowSums(!is.finite(support) | support != 1) == 0L
      z[keep, paste0("supported_", v)] <- x[keep, v]
    }
    .cast_writer_block(writer, z, first, count)
  }
  .cast_writer_finish(writer)
  input <- terra::rast(source)
  discrete <- c(fields[categorical], if (any(fields[categorical] %in% names(masks)))
    paste0("supported_", intersect(fields[categorical], names(masks))))
  continuous <- setdiff(layers, c(fields[flags], discrete))
  parts <- list()
  if (length(continuous)) parts <- c(parts, list(terra::project(input[[continuous]], template, method = "average")))
  if (length(discrete)) parts <- c(parts, list(terra::project(input[[discrete]], template, method = "mode")))
  if (any(flags)) {
    low <- terra::project(input[[fields[flags]]], template, method = "min")
    high <- terra::project(input[[fields[flags]]], template, method = "max")
    mixed <- terra::ifel(low == high, low, 2); names(mixed) <- fields[flags]
    parts <- c(parts, list(mixed))
  }
  display <- do.call(c, parts)[[layers]]
  terra::writeRaster(display, output, overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES")))
  output
}

.cast_explore_pair <- function(project, scenario, variables, shifts, context,
                               memo = new.env(parent = emptyenv())) {
  stacks <- .cast_explore_stacks(project, scenario, memo)
  # Choosing a factor changes the explanation, not the imported joint scenario.
  if (is.null(shifts)) return(stacks)
  a <- b <- if (context == "current") stacks$current else stacks$future
  for (v in variables) {
    a[[v]] <- stacks$current[[v]]
    b[[v]] <- stacks$current[[v]]
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
#'   freeze the other scenario changes. `NULL` uses the fitted scenario-changing predictors.
#'   Related derived predictors should be changed together using a consistent
#'   scenario. A single derived-variable change is a sensitivity analysis.
#' @param shifts Optional named numeric vector of raw-unit predictor shifts.
#'   When supplied, replaces the imported future change for those variables
#'   with these hypothetical changes. All other selected variables must be
#'   omitted: `variables` is inferred from `names(shifts)`.
#' @param context For hypothetical shifts only, `"current"` keeps other predictors current. `"future"` keeps
#'   other predictors at the selected future values in both members of the pair.
#' @param sample_n Optional maximum number of saved training-CV locations used
#'   for response summaries, driver ranking and optional grouping. Default
#'   `NULL` retains every saved location. A supplied limit selects deterministic
#'   existing rows; these summaries are not area-weighted raster summaries.
#' @param grouping Optional descriptive tree grouping, default `FALSE`.
#' @param maxdepth Maximum response tree depth. See [cast_response_groups()].
#' @param maps Compute native-grid response/group rasters. `FALSE` only computes
#'   the real-location response table and grouping; it never predicts a full grid.
#'   Matching saved publication response maps remain available even with FALSE;
#'   native common-support summaries and location summaries are kept separate.
#' @param focus Optional selected driver for detailed maps. `NULL` automatically
#'   uses the largest mean absolute driver-only response on the shared support
#'   subset of the saved training-CV locations. Only this driver's two background responses, model sign agreement
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
                         context = c("current", "future"), sample_n = NULL,
                         maxdepth = 3L, maps = TRUE, output_dir = NULL,
                         resources = list(), resume = TRUE, focus = NULL, grouping = FALSE) {
  if (!is.environment(getOption("castSDM.fit_identity_cache"))) {
    identity_options <- options(castSDM.fit_identity_cache = new.env(parent = emptyenv()))
    on.exit(options(identity_options), add = TRUE)
  }
  check_suggested("terra", "for saved-case exploration")
  if (!is.logical(grouping) || length(grouping) != 1L || is.na(grouping)) cli::cli_abort("grouping must be TRUE or FALSE.")
  if (grouping) check_suggested("rpart", "for descriptive response grouping")
  if (!is.logical(maps) || length(maps) != 1L || is.na(maps))
    cli::cli_abort("{.arg maps} must be TRUE or FALSE.")
  if (maps) check_suggested("leaflet", "for interactive display rasters")
  context <- match.arg(context)
  if (!is.null(sample_n)) .cast_check_scalar(sample_n, "sample_n", lower = 4, integer = TRUE)
  if (grouping) .cast_check_scalar(maxdepth, "maxdepth", lower = 1, upper = 30, integer = TRUE)
  if (!grouping) maxdepth <- NULL
  case_lock <- .cast_lock(file.path(case_dir, ".case.lock"))
  on.exit(filelock::unlock(case_lock), add = TRUE)
  project <- .cast_explore_project(case_dir)
  if (!identical(project$comparison$prediction_policy, "regular_gam/1"))
    cli::cli_abort("This saved case used an earlier prediction approximation. Run the current workflow in its new analysis directory before computing new responses. Saved maps remain available for browsing.")
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
  if (is.null(variables) && is.null(shifts)) {
    variables <- as.character(unlist(project$comparison$groups, use.names = FALSE))
    if (!length(variables)) cli::cli_abort("No fitted predictor changes in this scenario. Choose an explicit hypothetical shift to inspect sensitivity.")
  }
  variables <- .cast_explore_selection(project, variables)
  requested_focus <- focus
  if (!is.null(focus) && (!is.character(focus) || length(focus) != 1L ||
      is.na(focus) || !focus %in% variables))
    cli::cli_abort("{.arg focus} must name one of the compared variables.")
  # Hash in the worker, with a shared memo so repeated/static files are read once.
  memo <- new.env(parent = emptyenv())
  pair <- .cast_explore_pair(project, scenario, variables, shifts, context, memo)
  a <- pair$current; b <- pair$future
  invisible(.cast_resources(resources))
  settings <- list(scenario = scenario, variables = variables, shifts = shifts,
    context = context, sample_n = sample_n, maxdepth = maxdepth, maps = maps, focus = focus, grouping = grouping,
    weights = project$comparison$weights,
    clamp = project$run$protocol$clamp %||% FALSE,
    support_probs = project$comparison$support_probs,
    group_population = "Saved training-CV presence/background locations; not an area-weighted ecological population",
    ranking_population = "Range-supported saved training-CV locations; not all native raster cells",
    interpretation = "Paired fitted-suitability response; not ecological causal effects or occurrence probability change")
  mask_file <- file.path(project$case$preparation, "analysis_domain_mask.tif")
  if (!file.exists(mask_file)) cli::cli_abort("Missing recorded analysis domain: {.path {mask_file}}.")
  mask <- terra::rast(mask_file)
  mask_key <- .cast_data_identity(mask, memo)
  recorded_mask <- project$comparison$provenance$mask
  if (!is.null(recorded_mask) && !identical(mask_key, recorded_mask))
    cli::cli_abort("The analysis mask has changed since the saved comparison. Use the recorded mask or complete a new comparison.")
  saved <- if (is.null(shifts)) .cast_explore_saved(project, scenario, a, b, memo) else NULL
  if (!is.null(saved) && !all(variables %in% saved$drivers)) saved <- NULL
  key <- .cast_digest(list(schema = "cast_explore/2", settings = settings,
    saved = if (!is.null(saved)) unname(tools::md5sum(saved$raster_file)),
    fit = .cast_fit_identity(project$fit),
    current = .cast_data_identity(a, memo), future = .cast_data_identity(b, memo),
    mask = mask_key, folds = project$run$cv$fold_info,
    metadata = project$metadata,
    grouping = .cast_code_identity(cast_response_groups),
    membership = .cast_code_identity(.cast_response_nodes),
    evidence = .cast_code_identity(.cast_response_evidence),
    direction = .cast_code_identity(.cast_response_direction),
    compute = .cast_code_identity(cast_explore), attribution = .cast_code_identity(cast_attribution),
    attribution_execution = .cast_code_identity(.cast_attribution),
    sample = .cast_code_identity(.cast_explore_sample), extraction = .cast_code_identity(.cast_paired_extract),
    pair = .cast_code_identity(.cast_explore_pair),
    group_map = .cast_code_identity(.cast_explore_group_map),
    preview = .cast_code_identity(.cast_explore_preview),
    rpart = if (grouping) as.character(utils::packageVersion("rpart"))))
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
      if (!is.null(saved)) cached$raster_file <- saved$raster_file
      cached$directory <- directory; cached$case_dir <- project$directory
      return(cached)
    }
  }
  cli::cli_inform("Reading saved training-CV locations and aligned scenario values; rankings describe these locations, not habitat area.")
  sample <- .cast_explore_sample(project, a, b, sample_n)
  response <- if (!is.null(saved)) {
    value <- saved
    value$data <- terra::extract(terra::rast(saved$raster_file), sample$xy)
    value
  } else cast_attribution(project$fit, sample$current, sample$future,
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
  utils::write.csv(response$model_summary, file.path(directory, "response_by_model.csv"), row.names = FALSE)
  utils::write.csv(evidence$summary, file.path(directory, "drivers_common.csv"), row.names = FALSE)
  utils::write.csv(evidence$population, file.path(directory, "support_population.csv"), row.names = FALSE)
  utils::write.csv(evidence$drivers, file.path(directory, "driver_ids.csv"), row.names = FALSE)
  utils::write.csv(project$metadata, file.path(directory, "predictor_metadata.csv"), row.names = FALSE)
  if (!is.null(groups)) {
    utils::write.csv(groups$summary, file.path(directory, "groups.csv"), row.names = FALSE)
    utils::write.csv(groups$cv, file.path(directory, "group_tuning.csv"), row.names = FALSE)
    utils::write.csv(groups$fold_losses, file.path(directory, "group_fold_losses.csv"), row.names = FALSE)
  }
  mapped <- saved
  raw_file <- saved$raster_file
  preview <- group_file <- NULL
  if (maps) {
    if (is.null(saved)) {
      cli::cli_inform("Computing native-grid paired responses in bounded blocks.")
      raw_file <- file.path(directory, "response.tif")
      prepared <- if (is.null(shifts)) .cast_compare_reuse(project$fit, a, b,
        project$comparison, settings$weights, settings$clamp, memo)
      mapped <- .cast_attribution(project$fit, a, b, drivers = focus %||% character(), weights = settings$weights,
        shifts = shifts, support_probs = settings$support_probs, clamp = settings$clamp,
        mask = mask, filename = raw_file, overwrite = TRUE, resources = resources,
        resume = resume, diagnostics = TRUE, prepared = prepared)
    }
    if (!is.null(groups)) group_file <- .cast_explore_group_map(groups, a,
      terra::rast(raw_file)[["support_total"]], file.path(directory, "groups.tif"))
    cli::cli_inform("Preparing a separate interactive preview (not the analysis grid).")
    preview_source <- terra::rast(raw_file)
    wanted <- c("current", "future", "total_change", "support_total", "total_sd",
      if (!is.null(focus)) paste0(c("effect_", "retained_", "support_effect_", "support_retained_",
        "context_direction_", "context_gap_", "effect_agreement_", "retained_agreement_"), focus))
    preview_source <- preview_source[[intersect(wanted, names(preview_source))]]
    if (is.null(shifts)) {
      models <- terra::rast(project$comparison$results[[scenario]]$raster_file)
      detail <- grep("^change_|^sign_conflict_joint$", names(models), value = TRUE)
      if (length(detail)) preview_source <- c(preview_source, models[[detail]])
    }
    if (!is.null(group_file)) preview_source <- c(preview_source, terra::rast(group_file))
    # The preview helper accepts raster objects without crossing process boundaries.
    preview <- .cast_explore_preview(preview_source, file.path(directory, "preview.tif"))
  }
  if (!is.null(mapped)) {
    utils::write.csv(mapped$common_summary, file.path(directory,"native_common_summary.csv"), row.names=FALSE)
    if (!is.null(mapped$model_summary_common)) utils::write.csv(mapped$model_summary_common,
      file.path(directory,"native_response_by_model.csv"), row.names=FALSE)
  }
  if (!is.null(preview)) settings$response_limits <- .cast_app_response_limits(project, raw_file)
  result <- list(schema = "cast_explore/1", directory = directory, key = key,
    case_dir = project$directory, settings = settings, summary = response$summary,
    model_summary = response$model_summary,
    model_raster_file = if (is.null(shifts)) project$comparison$results[[scenario]]$raster_file,
    native_summary = mapped$common_summary,
    model_summary_population = if (!is.null(mapped$model_summary_common)) "Native grid cells; common paired support" else
      if (!is.null(saved)) "Native grid cells; legacy per-driver support" else settings$group_population,
    model_summary_common = mapped$model_summary_common,
    raster_hashes = list(raster_file = if (!is.null(raw_file)) unname(tools::md5sum(raw_file)),
      model_raster_file = if (is.null(shifts)) unname(tools::md5sum(project$comparison$results[[scenario]]$raster_file))),
    reused_response_maps = !is.null(saved), evidence = evidence,
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
  detail <- grep("^change_|^sign_conflict_joint$", names(r), value = TRUE)
  r <- r[[c("current", "joint", "joint_change", "support", "sd_joint_change", detail)]]
  names(r) <- c("current", "future", "total_change", "support_total", "total_sd", detail)
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
    cached$raster_file <- source; cached$model_raster_file <- source
    return(cached)
  }
  preview <- .cast_explore_preview(r, file.path(directory, "preview.tif"))
  result <- list(schema = "cast_app_view/1", directory = directory,
    case_dir = project$directory, preview = preview, raster_file = source, model_raster_file = source,
    raster_hashes = list(raster_file = unname(tools::md5sum(source)), model_raster_file = unname(tools::md5sum(source))),
    settings = list(scenario = scenario, context = "current", variables = project$fit$env_vars,
      shifts = NULL, weights = project$comparison$weights,
      response_limits = .cast_app_response_limits(project),
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

# A small, content-verified catalog connects publication maps to exploration.
# The numerical data and extrapolation policy must match, not just filenames.
.cast_explore_saved <- function(project, scenario, current, future, memo) {
  path <- file.path(project$directory, "figures", ".data", "driver_catalog.rds")
  if (!file.exists(path)) return(NULL)
  envelope <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(envelope$key)) return(NULL)
  value <- .cast_cache_read(path, envelope$key)
  if (is.null(value)) return(NULL)
  expected <- list(fit = .cast_fit_identity(project$fit), scenario = scenario,
    weights = project$comparison$weights, clamp = project$comparison$clamp,
    support_probs = project$comparison$support_probs,
    current = .cast_data_identity(current, memo), future = .cast_data_identity(future, memo),
    comparison = unname(tools::md5sum(project$comparison$results[[scenario]]$raster_file)))
  if (!identical(value$contract, expected) || !identical(envelope$key, .cast_digest(expected))) return(NULL)
  value$attribution$raster_file <- file.path(dirname(path), basename(value$attribution$raster_file))
  cli::cli_inform("Verified saved environmental responses: reuse native maps and extract location values; no response prediction required.")
  value$attribution
}

# One observed magnitude scale across every saved scenario and available driver
# map. Hypothetical maps may expand it; no values are clipped or overwritten.
.cast_app_response_limits <- function(project, extra = NULL) {
  catalog <- file.path(project$directory, "figures", ".data", "driver_attribution.tif")
  paths <- unique(c(catalog, vapply(project$comparison$results, `[[`, character(1), "raster_file"), extra))
  high <- 0
  for (p in paths[file.exists(paths)]) {
    r <- terra::rast(p)
    use <- grep("^change_|^joint_change$|^total_change$|^effect_|^retained_|^context_gap_", names(r), value = TRUE)
    use <- use[!grepl("agreement|_sd_", use)]
    if (!length(use)) next
    extrema <- terra::minmax(r[[use]], compute = TRUE)
    high <- max(high, abs(extrema[is.finite(extrema)]))
  }
  if (!is.finite(high) || high == 0) high <- .Machine$double.eps
  c(-high, high)
}
