#' Download CMIP6 Bioclimatic Data for Future Projections
#'
#' Wraps \code{geodata::cmip6_world()} to download and organize CMIP6
#' downscaled bioclimatic variables for multiple GCM x SSP x time-period
#' combinations. Optionally computes a multi-GCM ensemble mean per scenario.
#'
#' Downloaded rasters are cached in `path` and reused on subsequent calls.
#'
#' @param gcms Character vector. CMIP6 General Circulation Models. Default
#'   includes five models spanning low-to-high climate sensitivity:
#'   `"ACCESS-CM2"`, `"CMCC-ESM2"`, `"MIROC6"`, `"MRI-ESM2-0"`,
#'   `"IPSL-CM6A-LR"`.
#' @param ssps Character vector. SSP scenarios (number only). Default
#'   `c("245", "585")`. Options: `"126"`, `"245"`, `"370"`, `"585"`.
#' @param periods Character vector. Future time periods. Default
#'   `c("2041-2060", "2061-2080")`. Options: `"2021-2040"`, `"2041-2060"`,
#'   `"2061-2080"`.
#' @param var Character. Variable set: `"bioc"` (bioclimatic, default),
#'   `"tmin"`, `"tmax"`, `"prec"`.
#' @param res Numeric. Spatial resolution in arc-minutes. Default `2.5`.
#'   Options: `2.5`, `5`, `10`.
#' @param ensemble Logical. Compute multi-GCM ensemble mean per SSP x period.
#'   Default `TRUE`.
#' @param path Character. Download/cache directory. Default `"data/cmip6"`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A named list of `SpatRaster` objects. Names follow the pattern
#'   `"ssp{ssp}_{period}"` for ensemble means, or
#'   `"ssp{ssp}_{period}_{gcm}"` for individual GCM rasters.
#'
#' @details
#' ## Recommended GCM selection
#'
#' Following IPCC AR6 guidance, using 5+ GCMs spanning the range of climate
#' sensitivity reduces projection uncertainty. The defaults represent:
#'
#' - **ACCESS-CM2** (Australia, mid-high sensitivity)
#' - **CMCC-ESM2** (Italy, medium sensitivity)
#' - **MIROC6** (Japan, low sensitivity)
#' - **MRI-ESM2-0** (Japan, low sensitivity)
#' - **IPSL-CM6A-LR** (France, high sensitivity)
#'
#' ## SSP scenarios
#'
#' - **SSP1-2.6**: Sustainability path, ~1.5C warming
#' - **SSP2-4.5**: Middle of the road
#' - **SSP3-7.0**: Regional rivalry
#' - **SSP5-8.5**: Fossil-fueled development
#'
#' @seealso [cast_prepare_future_env()], [cast_project()]
#'
#' @export
cast_download_cmip6 <- function(
    gcms     = c("ACCESS-CM2", "CMCC-ESM2", "MIROC6",
                 "MRI-ESM2-0", "IPSL-CM6A-LR"),
    ssps     = c("245", "585"),
    periods  = c("2041-2060", "2061-2080"),
    var      = "bioc",
    res      = 2.5,
    ensemble = TRUE,
    path     = "data/cmip6",
    verbose  = TRUE) {
  check_suggested("geodata", "for CMIP6 climate data download")
  check_suggested("terra", "for raster operations")

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  cmip6_world <- utils::getFromNamespace("cmip6_world", "geodata")
  result <- list()

  for (ssp in ssps) {
    for (period in periods) {
      scenario_name <- sprintf("ssp%s_%s", ssp, gsub("-", "_", period))
      gcm_rasters <- list()

      for (gcm in gcms) {
        label <- sprintf("ssp%s_%s_%s", ssp, gsub("-", "_", period), gcm)
        if (verbose) {
          cli::cli_inform("Downloading {.val {gcm}} SSP{ssp} {period}...")
        }

        r <- tryCatch(
          cmip6_world(
            model = gcm, ssp = ssp, time = period,
            var = var, res = res, path = path
          ),
          error = function(e) {
            cli::cli_warn(
              "Failed to download {gcm} SSP{ssp} {period}: {e$message}"
            )
            NULL
          }
        )

        if (!is.null(r)) {
          result[[label]] <- r
          gcm_rasters[[gcm]] <- r
        }
      }

      # Compute multi-GCM ensemble mean
      if (ensemble && length(gcm_rasters) >= 2L) {
        if (verbose) {
          cli::cli_inform(
            "Computing ensemble mean for {scenario_name} ({length(gcm_rasters)} GCMs)..."
          )
        }
        ens <- tryCatch({
          stk <- terra::rast(gcm_rasters)
          # Average across GCMs for each bioclim layer
          n_layers <- terra::nlyr(gcm_rasters[[1]])
          layer_names <- names(gcm_rasters[[1]])
          ens_layers <- list()
          for (i in seq_len(n_layers)) {
            idx <- seq(i, terra::nlyr(stk), by = n_layers)
            ens_layers[[i]] <- terra::mean(stk[[idx]], na.rm = TRUE)
          }
          ens_r <- terra::rast(ens_layers)
          names(ens_r) <- layer_names

          # Save ensemble as GeoTIFF
          ens_dir <- file.path(path, "ensemble")
          dir.create(ens_dir, recursive = TRUE, showWarnings = FALSE)
          ens_path <- file.path(ens_dir, paste0(scenario_name, ".tif"))
          terra::writeRaster(ens_r, ens_path, overwrite = TRUE)
          if (verbose) cli::cli_inform("  Saved: {ens_path}")
          ens_r
        }, error = function(e) {
          cli::cli_warn("Ensemble computation failed: {e$message}")
          NULL
        })

        if (!is.null(ens)) {
          result[[scenario_name]] <- ens
        }
      } else if (length(gcm_rasters) == 1L) {
        # Only one GCM: use it directly as the scenario raster
        result[[scenario_name]] <- gcm_rasters[[1]]
      }
    }
  }

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Downloaded {length(result)} raster(s) to {.path {path}}."
    ))
  }

  result
}


#' Prepare Future Environment Data for Projection
#'
#' Extracts environmental values at species distribution grid points from
#' raster data, joins with static (non-climate) variables, and renames
#' columns to match training data variable names.
#'
#' @param rasters A named list of `SpatRaster` objects (from
#'   [cast_download_cmip6()] or user-supplied), or a single `SpatRaster`.
#'   When a named list, each element becomes one future scenario.
#' @param coords A `data.frame` with `lon` and `lat` columns defining the
#'   prediction grid (same grid used for current predictions).
#' @param static_vars Optional `data.frame` with the same number of rows as
#'   `coords`, containing non-climate variables (elevation, land cover, etc.)
#'   that are assumed constant across future scenarios. Default `NULL`.
#' @param var_mapping Optional named character vector mapping raster layer
#'   names to training data variable names. For example:
#'   `c("wc2.1_2.5m_bioc_01" = "bio01", "wc2.1_2.5m_bioc_02" = "bio02")`.
#'   When `NULL` (default), attempts automatic matching by extracting
#'   trailing bioclim numbers (e.g., `"_01"` -> `"bio01"`).
#' @param save_dir Optional path. When provided, saves each future
#'   environment data.frame as a CSV file. Default `NULL`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A named list of `data.frame`s, each containing `lon`, `lat`, and
#'   environmental variables ready for [cast_project()].
#'
#' @details
#' ## Variable name mapping
#'
#' CMIP6 raster layers from `geodata` have long names like
#' `"wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2041-2060_01"`. The `var_mapping`
#' parameter or the auto-detect logic strips these to match your training
#' data column names (e.g., `"bio01"`).
#'
#' ## Static variables
#'
#' Topographic variables (elevation, slope, terrain ruggedness) and land
#' cover are typically assumed constant in future projections. Pass them
#' via `static_vars` — they are joined by row position (must match
#' `coords` row count).
#'
#' @seealso [cast_download_cmip6()], [cast_project()]
#'
#' @export
cast_prepare_future_env <- function(
    rasters,
    coords,
    static_vars = NULL,
    var_mapping = NULL,
    save_dir    = NULL,
    verbose     = TRUE) {
  check_suggested("terra", "for raster extraction")

  if (!is.data.frame(coords) || !all(c("lon", "lat") %in% names(coords))) {
    cli::cli_abort("{.arg coords} must be a data.frame with lon and lat columns.")
  }

  # Handle single raster
  if (inherits(rasters, "SpatRaster")) {
    rasters <- list(scenario = rasters)
  }
  if (!is.list(rasters) || length(rasters) == 0) {
    cli::cli_abort("{.arg rasters} must be a non-empty named list of SpatRaster objects.")
  }

  # Filter to only SpatRaster entries (skip individual GCM rasters if ensemble exists)
  raster_list <- Filter(function(x) inherits(x, "SpatRaster"), rasters)
  if (length(raster_list) == 0) {
    cli::cli_abort("No valid SpatRaster objects found in {.arg rasters}.")
  }

  pts <- terra::vect(
    as.matrix(coords[, c("lon", "lat")]),
    type = "points",
    crs = "EPSG:4326"
  )

  if (!is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }

  result <- list()

  for (name in names(raster_list)) {
    r <- raster_list[[name]]
    if (verbose) cli::cli_inform("Extracting {.val {name}}: {terra::nlyr(r)} layers...")

    vals <- terra::extract(r, pts)
    if ("ID" %in% names(vals)) vals$ID <- NULL

    # Apply variable name mapping
    if (!is.null(var_mapping)) {
      old_names <- names(vals)
      for (i in seq_along(old_names)) {
        if (old_names[i] %in% names(var_mapping)) {
          names(vals)[i] <- var_mapping[old_names[i]]
        }
      }
    } else {
      # Auto-detect: extract bioclim numbers from layer names
      names(vals) <- .auto_rename_bioclim(names(vals))
    }

    # Build output data.frame
    out <- data.frame(lon = coords$lon, lat = coords$lat)
    out <- cbind(out, vals)

    # Join static variables
    if (!is.null(static_vars)) {
      if (nrow(static_vars) != nrow(coords)) {
        cli::cli_abort(
          "{.arg static_vars} must have the same number of rows as {.arg coords} ({nrow(coords)})."
        )
      }
      # Avoid duplicating lon/lat or climate vars
      static_cols <- setdiff(names(static_vars), c("lon", "lat", names(vals)))
      if (length(static_cols) > 0) {
        out <- cbind(out, static_vars[, static_cols, drop = FALSE])
      }
    }

    result[[name]] <- out

    if (!is.null(save_dir)) {
      csv_path <- file.path(save_dir, paste0(name, ".csv"))
      utils::write.csv(out, csv_path, row.names = FALSE)
      if (verbose) cli::cli_inform("  Saved: {csv_path}")
    }
  }

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Prepared {length(result)} future environment data.frame(s)."
    ))
  }

  result
}


#' Auto-rename bioclim layers to standard short names
#'
#' Converts long geodata layer names like
#' `"wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2041-2060_01"` to `"bio01"`.
#'
#' @param x Character vector of layer names.
#' @return Character vector of renamed layer names.
#' @keywords internal
#' @noRd
.auto_rename_bioclim <- function(x) {
  vapply(x, function(name) {
    # Try to extract trailing number (bioclim index)
    m <- regmatches(name, regexpr("_(\\d{1,2})$", name))
    if (length(m) == 1 && nzchar(m)) {
      num <- sub("^_", "", m)
      return(sprintf("bio%s", num))
    }
    # Try pattern like "bioc_01" or "bio_01"
    m2 <- regmatches(name, regexpr("bio[c]?_(\\d{1,2})", name))
    if (length(m2) == 1 && nzchar(m2)) {
      num <- sub("^bio[c]?_", "", m2)
      return(sprintf("bio%s", num))
    }
    # No match — keep original
    name
  }, character(1), USE.NAMES = FALSE)
}
