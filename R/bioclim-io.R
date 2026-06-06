#' Load CHNECO26 Bioclim Raster Stack
#'
#' Locates and loads bioclimatic raster layers from a CHNECO26-style directory
#' structure where each variable lives in its own subdirectory. Supports two
#' discovery modes: metadata CSV (preferred) or recursive directory scan.
#'
#' @section CHNECO26 Directory Layout:
#' ```
#' bioclim_root/
#'   CHNECO26_datalayers_details_bioclim.csv   (metadata)
#'   bioclim/cnclim1km/present/1991_2020/
#'     bio01/bioclim_cnclim1km_present_1991_2020_bio01.tif
#'     bio02/bioclim_cnclim1km_present_1991_2020_bio02.tif
#'     ...
#'   bioclim/cnclim1km_ensmean/future/
#'     2021_2040/ssp126/bio01/...tif
#'     ...
#' ```
#'
#' @param bioclim_root Character. Root directory of the bioclim dataset
#'   (e.g., `"E:/CHECO26/Bioclim/CHNECO26_Bioclim"`).
#' @param variables Character vector. Variable names to load
#'   (e.g., `sprintf("bio%02d", 1:19)`). Default loads standard 19 bioclim.
#' @param period Character. `"present"` or `"future"`. Default `"present"`.
#' @param sub_period Character. Time period (e.g., `"1991_2020"`,
#'   `"2021_2040"`). Default `"1991_2020"`.
#' @param scenario Character or `NULL`. SSP scenario for future periods
#'   (e.g., `"ssp126"`). Required when `period = "future"`.
#' @param model Character or `NULL`. Climate model (e.g., `"ensmean"`,
#'   `"access_cm2"`). For future periods. Default `"ensmean"`.
#' @param metadata_csv Character or `NULL`. Path to the CHNECO26 metadata
#'   CSV. If `NULL`, attempts to find it at
#'   `<bioclim_root>/CHNECO26_datalayers_details_bioclim.csv`.
#'   If the CSV doesn't exist, falls back to recursive directory scanning.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A named `terra::SpatRaster` stack with one layer per variable.
#'   Layer names match the `variables` argument.
#'
#' @details
#' **Metadata CSV mode** (preferred): Reads the CSV, filters by period,
#' sub_period, scenario, model, and variable, then constructs absolute
#' paths from relative_path + bioclim_root.
#'
#' **Recursive scan mode** (fallback): Searches recursively for TIF files
#' matching the variable names, excluding `.tmp.tif` files.
#'
#' @export
cast_load_bioclim <- function(bioclim_root,
                              variables   = sprintf("bio%02d", 1:19),
                              period      = "present",
                              sub_period  = "1991_2020",
                              scenario    = NULL,
                              model       = "ensmean",
                              metadata_csv = NULL,
                              verbose     = TRUE) {

  # ── Validate inputs ──
  if (!dir.exists(bioclim_root)) {
    cli::cli_abort("Bioclim root directory not found: {.path {bioclim_root}}")
  }
  period <- match.arg(period, c("present", "future"))
  if (period == "future" && is.null(scenario)) {
    cli::cli_abort("{.arg scenario} is required for future period (e.g., 'ssp126').")
  }

  # ── Try metadata CSV mode first ──
  if (is.null(metadata_csv)) {
    metadata_csv <- file.path(bioclim_root, "CHNECO26_datalayers_details_bioclim.csv")
  }

  if (file.exists(metadata_csv)) {
    paths <- .load_bioclim_from_metadata(
      bioclim_root, metadata_csv, variables,
      period, sub_period, scenario, model, verbose
    )
  } else {
    if (verbose) {
      cli::cli_alert_info("Metadata CSV not found, using recursive scan")
    }
    paths <- .load_bioclim_recursive(
      bioclim_root, variables, period, sub_period,
      scenario, model, verbose
    )
  }

  # ── Load as SpatRaster stack ──
  if (!requireNamespace("terra", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg terra} is required for raster loading.")
  }

  stack <- terra::rast(paths)
  names(stack) <- variables

  if (verbose) {
    cli::cli_alert_success(
      "Loaded {length(variables)} bioclim layers ({period}, {sub_period})"
    )
    cli::cli_alert_info(
      "Resolution: {paste(round(terra::res(stack), 4), collapse = ' x ')}"
    )
  }

  stack
}


#' Load Bioclim Paths from Metadata CSV
#' @noRd
.load_bioclim_from_metadata <- function(bioclim_root, metadata_csv, variables,
                                         period, sub_period, scenario, model,
                                         verbose) {
  meta <- utils::read.csv(metadata_csv, stringsAsFactors = FALSE)

  # Normalise NA strings
  if ("scenario" %in% names(meta)) {
    meta$scenario <- ifelse(
      is.na(meta$scenario) | meta$scenario == "NA",
      NA_character_, meta$scenario
    )
  }
  if ("model" %in% names(meta)) {
    meta$model <- ifelse(
      is.na(meta$model) | meta$model == "NA",
      NA_character_, meta$model
    )
  }

  # Filter rows
  if (period == "present") {
    rows <- meta[
      meta$period == "present" &
        meta$sub_period == sub_period &
        meta$variable %in% variables &
        is.na(meta$scenario) &
        is.na(meta$model),
      ,
      drop = FALSE
    ]
  } else {
    # For future: match model name in dataset field (e.g., cnclim1km_ensmean)
    model_pattern <- if (!is.null(model)) model else ""
    rows <- meta[
      meta$period == "future" &
        meta$sub_period == sub_period &
        meta$scenario == scenario &
        meta$variable %in% variables,
      ,
      drop = FALSE
    ]
    # Filter by model if specified
    if (!is.null(model) && nzchar(model)) {
      rows <- rows[grepl(model, rows$dataset, fixed = TRUE), , drop = FALSE]
    }
  }

  # Deduplicate and order
  rows <- rows[!duplicated(rows$variable), , drop = FALSE]
  rows <- rows[order(match(rows$variable, variables)), , drop = FALSE]

  # Check completeness
  found_vars <- rows$variable
  missing <- setdiff(variables, found_vars)
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Missing bioclim variables in metadata:",
      "x" = "Not found: {paste(missing, collapse = ', ')}",
      "i" = "Period: {period}, sub_period: {sub_period}",
      "i" = "Scenario: {scenario %||% 'NA'}, model: {model %||% 'NA'}"
    ))
  }

  # Build absolute paths
  paths <- file.path(
    normalizePath(bioclim_root, winslash = "/", mustWork = FALSE),
    gsub("\\\\", "/", rows$relative_path)
  )

  # Verify existence
  exists_ok <- file.exists(paths)
  if (!all(exists_ok)) {
    bad <- paths[!exists_ok]
    cli::cli_abort(c(
      "Bioclim raster file(s) not found on disk:",
      "x" = "{paste(head(bad, 3), collapse = '\\n')}"
    ))
  }

  if (verbose) {
    cli::cli_alert_info("Located {length(paths)} layers via metadata CSV")
  }

  paths
}


#' Load Bioclim Paths via Recursive Directory Scan
#' @noRd
.load_bioclim_recursive <- function(bioclim_root, variables, period,
                                     sub_period, scenario, model, verbose) {
  # Build expected directory path
  if (period == "present") {
    search_dir <- file.path(bioclim_root, "bioclim", "cnclim1km",
                            "present", sub_period)
  } else {
    model_dir <- if (!is.null(model) && nzchar(model)) {
      paste0("cnclim1km_", model)
    } else {
      "cnclim1km_ensmean"
    }
    search_dir <- file.path(bioclim_root, "bioclim", model_dir,
                            "future", sub_period, scenario)
  }

  if (!dir.exists(search_dir)) {
    cli::cli_abort("Bioclim directory not found: {.path {search_dir}}")
  }

  # Find TIF files recursively
  all_tifs <- list.files(search_dir, pattern = "\\.tif$",
                         recursive = TRUE, full.names = TRUE)
  # Exclude .tmp.tif files
  all_tifs <- all_tifs[!grepl("\\.tmp\\.tif$", all_tifs)]

  # Match each variable
  paths <- character(length(variables))
  names(paths) <- variables

  for (i in seq_along(variables)) {
    v <- variables[i]
    # Match pattern: variable name in filename, e.g., *_bio01.tif
    pattern <- paste0("_", v, "\\.tif$")
    matches <- grep(pattern, all_tifs, value = TRUE)
    if (length(matches) == 0) {
      # Try broader match: /bio01/ directory
      pattern2 <- paste0("/", v, "/")
      matches <- grep(pattern2, all_tifs, value = TRUE)
    }
    if (length(matches) == 0) {
      cli::cli_abort("Could not find TIF for variable {.val {v}} in {.path {search_dir}}")
    }
    paths[i] <- matches[1]  # Take first match
  }

  if (verbose) {
    cli::cli_alert_info("Located {length(paths)} layers via recursive scan")
  }

  unname(paths)
}


#' Build Future Bioclim Raster List for Projection
#'
#' Convenience wrapper that loads multiple future scenario raster stacks
#' into a named list suitable for [cast_project_raster()].
#'
#' @param bioclim_root Character. Root bioclim directory.
#' @param variables Character vector. Variables to load.
#' @param scenarios Character vector. SSP scenarios (e.g., `c("ssp126", "ssp245", "ssp585")`).
#' @param periods Character vector. Future periods (e.g., `c("2021_2040", "2041_2070", "2071_2100")`).
#' @param model Character. Climate model. Default `"ensmean"`.
#' @param metadata_csv Character or `NULL`. Metadata CSV path.
#' @param verbose Logical.
#'
#' @return Named list of `terra::SpatRaster` stacks. Names follow the pattern
#'   `"{model}_{scenario}_{period}"`.
#'
#' @export
cast_load_future_bioclim <- function(bioclim_root,
                                     variables = sprintf("bio%02d", 1:19),
                                     scenarios = c("ssp126", "ssp245", "ssp585"),
                                     periods   = c("2021_2040", "2041_2070", "2071_2100"),
                                     model     = "ensmean",
                                     metadata_csv = NULL,
                                     verbose   = TRUE) {

  future_rasters <- list()

  for (ssp in scenarios) {
    for (prd in periods) {
      scen_name <- paste(model, ssp, prd, sep = "_")

      stk <- tryCatch(
        cast_load_bioclim(
          bioclim_root = bioclim_root,
          variables    = variables,
          period       = "future",
          sub_period   = prd,
          scenario     = ssp,
          model        = model,
          metadata_csv = metadata_csv,
          verbose      = FALSE
        ),
        error = function(e) {
          if (verbose) {
            cli::cli_alert_warning(
              "Skipping {scen_name}: {conditionMessage(e)}"
            )
          }
          NULL
        }
      )

      if (!is.null(stk)) {
        future_rasters[[scen_name]] <- stk
        if (verbose) {
          cli::cli_alert_success("Loaded future: {scen_name}")
        }
      }
    }
  }

  if (length(future_rasters) == 0) {
    cli::cli_warn("No future bioclim scenarios found.")
  } else if (verbose) {
    cli::cli_alert_success(
      "Total: {length(future_rasters)} future scenario stacks loaded"
    )
  }

  future_rasters
}
