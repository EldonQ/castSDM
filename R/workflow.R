#' Initialise, Check or Run the CHELSA Workflow
#'
#' Creates a workspace from the installed, versioned workflow. Checks the real
#' input inventory without fitting; execution is opt-in. Current nonclimate
#' predictors remain fixed in all future projections.
#' @param path Writable workflow workspace, separate from the environmental library.
#' @param data_dir Environmental library containing `covariates/reg` and `mask_china.tif`.
#' @param occurrence_file Real occurrence CSV/PSV file supported by the data adapter.
#' @param species Taxon identifier. The Schizothoracinae preset is a pooled subfamily.
#' @param settings Named list: `domain` (`national` or `freshwater_river`),
#'   `input_type` (`gbif` or `nsdm_occ`), `lulc_product` (`glu` or `hlu`),
#'   `models`, `cv_k`, `background_points`, and `resources` (cores, workers, memory_mb).
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
cast_workflow <- function(path, data_dir, occurrence_file, species = "Schizothoracinae",
                          settings = list(), run = FALSE) {
  for (p in c("terra", "data.table", "sf", "blockCV", "ranger", "maxnet", "mgcv",
              "glmnet", "iml", "pROC", "FNN", "ggplot2", "patchwork", "scales", "ragg"))
    check_suggested(p, "for the CHELSA workflow")
  paths <- list(path = path, data_dir = data_dir, occurrence_file = occurrence_file)
  if (any(!vapply(paths, function(x) is.character(x) && length(x) == 1L &&
      !is.na(x) && nzchar(trimws(x)), logical(1)))) cli::cli_abort("Supply workspace, data and occurrence paths.")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  occurrence_file <- normalizePath(occurrence_file, winslash = "/", mustWork = TRUE)
  if (startsWith(tolower(paste0(path, "/")), tolower(paste0(data_dir, "/"))))
    cli::cli_abort("Keep the workspace outside the environmental data library.")
  if (length(species) != 1L || is.na(species) || !grepl("^[A-Za-z][A-Za-z0-9_-]*$", species))
    cli::cli_abort("Use a taxon identifier with letters, digits, underscores or hyphens.")
  allowed <- c("domain", "input_type", "lulc_product", "models", "cv_k", "background_points", "resources")
  if (!is.list(settings) || (length(settings) && (is.null(names(settings)) ||
      any(!names(settings) %in% allowed) || anyDuplicated(names(settings)))))
    cli::cli_abort("Unknown or duplicated workflow settings.")
  defaults <- list(domain = if (species == "Schizothoracinae") "freshwater_river" else "national",
    input_type = "gbif", lulc_product = "glu", models = c("rf", "maxent", "gam"),
    cv_k = 5L, background_points = 10000L, resources = list(cores = 4L, workers = 1L, memory_mb = 4096))
  settings <- utils::modifyList(defaults, settings)
  for (n in c("domain", "input_type", "lulc_product")) {
    choices <- switch(n, domain = c("national", "freshwater_river"),
      input_type = c("gbif", "nsdm_occ"), lulc_product = c("glu", "hlu"))
    if (length(settings[[n]]) != 1L || !settings[[n]] %in% choices) cli::cli_abort("Invalid workflow {.val {n}} setting.")
  }
  # Reuse package validation rather than maintain a second numerical contract.
  cast_protocol(models = settings$models, cv_k = settings$cv_k,
    num_threads = settings$resources$cores)
  settings$background_points <- .cast_check_scalar(settings$background_points, "background_points", lower = 1, integer = TRUE)
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
  Sys.setenv(CAST_EXAMPLE_ROOT = path, CAST_DATA_ROOT = data_dir, CAST_SPECIES = species,
    CAST_RUN_ID = "cast033_chelsa_ensmean", CAST_OUTPUT_ROOT = file.path(path, "output"),
    CAST_DRY_RUN = if (isTRUE(run)) "0" else "1")
  entry <- file.path(path, "00_RUN_ALL.R")
  sys.source(entry, envir = new.env(parent = globalenv()), keep.source = FALSE)
  invisible(list(path = path, entry = entry, case_dir = file.path(path, "output", "cast033_chelsa_ensmean", species),
    checked = TRUE, run = isTRUE(run), settings = settings))
}
