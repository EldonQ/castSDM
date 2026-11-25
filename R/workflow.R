.cast_workflow_run_id <- function() {
  parts <- strsplit(as.character(utils::packageVersion("castSDM")), ".", fixed = TRUE)[[1L]]
  # Patch releases resume the same user workspace; artifact keys check inputs.
  parts <- utils::head(parts, 2L)
  paste0("cast", paste(parts, collapse = ""), "_chelsa_ensmean")
}

# One selection contract for script entry points and Shiny. Configuration
# survives R restarts; a stale session override must not switch the data set.
.cast_workflow_selection <- function(path, allow_unconfigured = FALSE) {
  file <- file.path(path, "workflow_settings.rds")
  settings <- if (file.exists(file)) tryCatch(readRDS(file), error = function(e)
    cli::cli_abort("Cannot read workflow settings {.path {file}}: {conditionMessage(e)}")) else NULL
  valid <- function(x) is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl("^[A-Za-z][A-Za-z0-9_-]*$", x)
  if (!is.null(settings) && (!is.list(settings) || !valid(settings$species)))
    cli::cli_abort("The saved workflow settings need one explicit species identifier. Configure the workspace with cast_workflow().")
  override <- trimws(Sys.getenv("CAST_SPECIES", ""))
  if (nzchar(override) && !valid(override)) cli::cli_abort("Invalid CAST_SPECIES identifier.")
  if (!is.null(settings) && nzchar(override) && !identical(override, settings$species))
    cli::cli_abort("CAST_SPECIES ({override}) conflicts with the saved species ({settings$species}). Unset CAST_SPECIES or reconfigure this workspace with cast_workflow(); no other species will be run.")
  species <- settings$species %||% if (nzchar(override)) override else NULL
  if (is.null(species) && !allow_unconfigured)
    cli::cli_abort("No species is configured. Use cast_workflow(path, data_dir, occurrence_file, species, run = FALSE) or the Shiny setup form. There is no default fish or squirrel analysis.")
  list(species = species, settings = settings,
    source = if (!is.null(settings)) "saved project settings" else if (!is.null(species)) "session override (not persistent)" else "unconfigured")
}

#' Initialise, Check or Run the CHELSA Workflow
#'
#' Creates a workspace from the installed, versioned workflow. Checks the real
#' input inventory without fitting; execution is opt-in. Current nonclimate
#' predictors remain fixed in all future projections.
#' @param path Writable workflow workspace, separate from the environmental library.
#' @param data_dir Environmental library containing `covariates/reg` and `mask_china.tif`.
#' @param occurrence_file Real occurrence CSV/PSV file supported by the data adapter.
#' @param species Explicit taxon identifier. Required; no species is inferred
#'   from a filename or a previous case. Input format and domain are separate
#'   settings and are never inferred from the taxon name.
#' @param settings Named list: `domain` (`national`, `freshwater_river` or
#'   `custom_mask`), `domain_mask` (path to an aligned study-domain raster;
#'   positive cells are active),
#'   `input_type` (`gbif` or `nsdm_occ`), `lulc_product` (`glu` or `hlu`),
#'   `models`, `cv_k`, `block_method` (`blockcv` or `grid`), `block_size`
#'   (metres for `blockcv`; NULL estimates the range), `holdout_block_size`
#'   (optional distinct outer scale), `background_points`, `thin_km` (distance thinning in km),
#'   `predictor_min_coverage` (actual-domain coverage), `vif_threshold`, `predictor_missing`, and `resources` (cores, workers, memory_mb).
#' @param run Execute preparation, fitting, projections and figures after checks.
#'   Default `FALSE` only indexes and checks physical inputs.
#' @return Workspace, entry, case directory and checked settings, invisibly.
#' @details Future files must be CHELSA v2.1 ensemble-mean climate for SSP126,
#'   SSP370 and SSP585 in 2011-2040, 2041-2070 and 2071-2100, relative to
#'   1981-2010. Predictions at mean climate do not estimate GCM uncertainty.
#'   No occurrence labels or environmental measurements are manufactured.
#'   Existing helper scripts must match this installed version; settings and
#'   outputs are preserved. Source inputs are read only.
#' @export
cast_workflow <- function(path, data_dir, occurrence_file, species = NULL,
                          settings = list(), run = FALSE) {
  if (!is.character(species) || length(species) != 1L || is.na(species) || !grepl("^[A-Za-z][A-Za-z0-9_-]*$", species))
    cli::cli_abort("Use a taxon identifier with letters, digits, underscores or hyphens.")
  for (p in c("terra", "data.table", "sf", "blockCV", "ranger", "maxnet", "mgcv",
              "iml", "pROC", "FNN", "ggplot2", "patchwork", "scales", "ragg"))
    check_suggested(p, "for the CHELSA workflow")
  paths <- list(path = path, data_dir = data_dir, occurrence_file = occurrence_file)
  if (any(!vapply(paths, function(x) is.character(x) && length(x) == 1L &&
      !is.na(x) && nzchar(trimws(x)), logical(1)))) cli::cli_abort("Supply workspace, data and occurrence paths.")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  occurrence_file <- normalizePath(occurrence_file, winslash = "/", mustWork = TRUE)
  if (startsWith(tolower(paste0(path, "/")), tolower(paste0(data_dir, "/"))))
    cli::cli_abort("Keep the workspace outside the environmental data library.")
  allowed <- c("domain", "domain_mask", "input_type", "lulc_product", "models", "cv_k", "block_method",
    "block_size", "holdout_block_size", "background_points", "thin_km", "predictor_min_coverage", "vif_threshold", "predictor_missing", "resources")
  if (!is.list(settings) || (length(settings) && (is.null(names(settings)) ||
      any(!names(settings) %in% allowed) || anyDuplicated(names(settings)))))
    cli::cli_abort("Unknown or duplicated workflow settings.")
  defaults <- list(domain = "national",
    input_type = "gbif", lulc_product = "glu", models = c("rf", "maxent", "gam"),
    cv_k = 5L, background_points = 10000L, thin_km = 1, predictor_min_coverage = .95, vif_threshold = 10, predictor_missing = "complete_cases", resources = list(cores = 4L, workers = 1L, memory_mb = 4096))
  settings <- utils::modifyList(defaults, settings)
  for (n in c("domain", "input_type", "lulc_product", "block_method")) {
    if (n == "block_method" && is.null(settings[[n]])) next
    choices <- switch(n, domain = c("national", "freshwater_river", "custom_mask"),
      input_type = c("gbif", "nsdm_occ"), lulc_product = c("glu", "hlu"), block_method = c("blockcv", "grid"))
    if (length(settings[[n]]) != 1L || !settings[[n]] %in% choices) cli::cli_abort("Invalid workflow {.val {n}} setting.")
  }
  if (identical(settings$domain, "custom_mask")) {
    if (!is.character(settings$domain_mask) || length(settings$domain_mask) != 1L ||
        is.na(settings$domain_mask) || !file.exists(settings$domain_mask))
      cli::cli_abort("custom_mask requires an existing domain_mask raster.")
    settings$domain_mask <- normalizePath(settings$domain_mask, winslash = "/", mustWork = TRUE)
    mask <- terra::rast(settings$domain_mask)
    reference <- terra::rast(file.path(data_dir, "mask_china.tif"))
    if (terra::nlyr(mask) != 1L || !isTRUE(terra::compareGeom(mask, reference, stopOnError = FALSE)))
      cli::cli_abort("domain_mask must be one raster layer aligned with the environmental library; no resampling is applied automatically.")
  }
  # Reuse package validation rather than maintain a second numerical contract.
  cast_protocol(models = settings$models, cv_k = settings$cv_k,
    block_method = settings$block_method %||% "blockcv", block_size = settings$block_size,
    holdout_block_size = settings$holdout_block_size,
    num_threads = settings$resources$cores, select_vif_threshold = settings$vif_threshold)
  if (length(settings$predictor_missing) != 1L || is.na(settings$predictor_missing) || !settings$predictor_missing %in% c("drop_predictors", "complete_cases")) cli::cli_abort("Invalid predictor_missing policy.")
  settings$background_points <- .cast_check_scalar(settings$background_points, "background_points", lower = 1, integer = TRUE)
  settings$thin_km <- .cast_check_scalar(settings$thin_km, "thin_km", lower = 0)
  settings$predictor_min_coverage <- .cast_check_scalar(settings$predictor_min_coverage, "predictor_min_coverage", lower = 0, upper = 1, min_inclusive = FALSE)
  .cast_resources(settings$resources)
  settings <- c(list(species = species, data_dir = data_dir, occurrence_file = occurrence_file), settings)
  template <- system.file("workflows", "single_species", package = "castSDM", mustWork = TRUE)
  scripts <- list.files(template, full.names = TRUE)
  scripts <- scripts[!file.info(scripts)$isdir]
  existing <- file.path(path, basename(scripts)); found <- file.exists(existing)
  if (any(found) && !identical(unname(tools::md5sum(existing[found])), unname(tools::md5sum(scripts[found]))))
    cli::cli_abort("Workspace helper scripts differ from this package. Use a new workspace or update its scripts from the installed workflow; keep settings and outputs.")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (any(!file.copy(scripts[!found], path))) cli::cli_abort("Cannot initialise the workflow workspace.")
  lock <- .cast_lock(file.path(path, ".workflow.lock"))
  on.exit(filelock::unlock(lock), add = TRUE)
  .cast_atomic_rds(settings, file.path(path, "workflow_settings.rds"))
  variables <- c("CAST_EXAMPLE_ROOT", "CAST_DATA_ROOT", "CAST_SPECIES", "CAST_RUN_ID", "CAST_OUTPUT_ROOT", "CAST_DRY_RUN")
  old <- Sys.getenv(variables, unset = NA_character_)
  on.exit({ Sys.unsetenv(variables); present <- !is.na(old)
    if (any(present)) do.call(Sys.setenv, as.list(stats::setNames(old[present], variables[present]))) }, add = TRUE)
  run_id <- Sys.getenv("CAST_RUN_ID", .cast_workflow_run_id())
  Sys.setenv(CAST_EXAMPLE_ROOT = path, CAST_DATA_ROOT = data_dir, CAST_SPECIES = species,
    CAST_RUN_ID = run_id, CAST_OUTPUT_ROOT = file.path(path, "output"),
    CAST_DRY_RUN = if (isTRUE(run)) "0" else "1")
  entry <- file.path(path, "00_RUN_ALL.R")
  sys.source(entry, envir = new.env(parent = globalenv()), keep.source = FALSE)
  invisible(list(path = path, entry = entry, case_dir = file.path(path, "output", run_id, species),
    checked = TRUE, run = isTRUE(run), settings = settings))
}
