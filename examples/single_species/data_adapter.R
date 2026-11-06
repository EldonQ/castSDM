# ==============================================================================
# Data adapter for the castSDM single-species workflow.
#
# The adapter owns the data-library layout; the model and scenario resolver
# only consume its outputs:
#
#   * nsdm_load_current_stack(): a terra::SpatRaster with one layer per
#     eligible current-pool variable, named by the pool `layer` id, with the
#     layer table attached as attr(., "nsdm_layers").
#   * nsdm_build_scenario_stacks(): one SpatRaster per resolved scenario, with
#     layers ordered and named exactly as `selected_vars`; climate/LULC layers
#     come from the scenario table, static layers are reused from the current
#     stack.
#   * wf_analysis_domain(): one binary analysis domain (national land mask,
#     river-network buffer with a local occurrence envelope, or a custom mask)
#     used by background sampling, projection, attribution and maps.
#
# Data source: the national 1 km Albers environmental library. Current and
# future layer metadata is written by 00_INVENTORY.R; no variable names are
# stored here.
#
# CRS policy (HARD CONSTRAINT on this machine):
#   In-process PROJ calls segfault here; raw occurrence points use external cs2cs.
# ==============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

nsdm_default_data_root <- function() Sys.getenv("CAST_DATA_ROOT", "E:/CHECO26/N-SDM/data")

nsdm_project_occurrences <- function(lonlat, target_crs) {
  exe <- Sys.which("cs2cs")
  if (!nzchar(exe)) stop("PROJ cs2cs is required to transform WGS84 occurrence points.")
  input <- tempfile(fileext = ".txt")
  output <- tempfile(fileext = ".txt")
  errors <- tempfile(fileext = ".log")
  on.exit(unlink(c(input, output, errors)), add = TRUE)
  utils::write.table(lonlat, input, row.names = FALSE, col.names = FALSE, quote = FALSE)
  args <- c("-f", "%.8f", "+proj=longlat", "+datum=WGS84", "+to",
            strsplit(target_crs, " ", fixed = TRUE)[[1]], shQuote(input))
  status <- system2(exe, args, stdout = output, stderr = errors)
  if (status != 0L) stop("Coordinate conversion failed: ", paste(readLines(errors), collapse = " "))
  xy <- as.matrix(utils::read.table(output)[, 1:2, drop = FALSE])
  if (nrow(xy) != nrow(lonlat) || !all(is.finite(xy))) {
    stop("Coordinate conversion did not return finite coordinates for every record.")
  }
  colnames(xy) <- c("lon", "lat")
  xy
}

nsdm_default_pool_csv <- function() {
  file.path(Sys.getenv("CAST_EXAMPLE_ROOT", getwd()), "config", "pool_current_eligible_glu.csv")
}

# ---- Current stack -----------------------------------------------------------

# Read the generated current-pool manifest. Returns a data.frame with columns
# group/layer/path/is_climate (plus any extra QC columns from the CSV).
# Paths in the generated manifest are authoritative. Rebuild it after moving data.
nsdm_pool_layers <- function(pool_csv = nsdm_default_pool_csv()) {
  if (!file.exists(pool_csv)) stop("Pool layer table not found: ", pool_csv)
  pool <- utils::read.csv(pool_csv, stringsAsFactors = FALSE)
  need <- c("group", "layer", "path", "is_climate")
  if (!all(need %in% names(pool))) {
    stop("Pool table ", pool_csv, " lacks column(s): ",
         paste(setdiff(need, names(pool)), collapse = ", "))
  }
  if ("status" %in% names(pool)) pool <- pool[pool$status %in% c("kept", NA), , drop = FALSE]
  if (anyDuplicated(pool$layer)) stop("Duplicate layer IDs in pool manifest: ", pool_csv)
  pool$orig_layer <- pool$layer
  pool$layer <- nsdm_sanitize_layer_id(pool$layer)
  missing_paths <- pool$path[!file.exists(pool$path)]
  if (length(missing_paths)) {
    stop("Missing pool raster(s):\n  ", paste(missing_paths, collapse = "\n  "))
  }
  rownames(pool) <- NULL
  pool
}

# Layer ids are the single naming choke point of this adapter: the current
# stack, model_data columns, selection/attribution output and future stacks all
# take their names from the pool `layer` column. Hydro layer stems contain
# hyphens (e.g. "h90m_channel-grad-dw-seg"), which break every formula-interface
# fitter downstream -- ranger aborts with "Illegal column names in formula
# interface", and gbm/maxnet term parsing misreads "a-b" as subtraction.
# Sanitize once here so the whole pipeline sees syntactic ids; `orig_layer`
# preserves the file-stem id for provenance and future-path resolution.
nsdm_sanitize_layer_id <- function(x) {
  y <- gsub("[^0-9A-Za-z]+", "_", x)
  y <- sub("^_+", "", y)
  y <- ifelse(grepl("^[0-9]", y), paste0("x", y), y)
  if (anyDuplicated(y)) {
    stop("Layer-id sanitization collided: ",
         paste(unique(x[duplicated(y)]), collapse = ", "))
  }
  y
}

# Metadata-only grid check: every layer must align with the reference mask.
# Never reprojects (banned here); aborts and lists the misaligned layers.
nsdm_check_geom <- function(paths, reference, label = "layer") {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to load environmental rasters.")
  }
  if (!file.exists(reference)) stop("Reference mask not found: ", reference)
  ref <- terra::rast(reference)
  bad <- character()
  for (p in paths) {
    ok <- isTRUE(tryCatch(
      terra::compareGeom(ref, terra::rast(p), stopOnError = FALSE, messages = FALSE),
      error = function(e) FALSE))
    if (!ok) bad <- c(bad, p)
  }
  if (length(bad)) {
    stop(label, "(s) not aligned with ", basename(reference),
         " (reprojection is disabled on this machine; fix the grid upstream,",
         " e.g. gdalwarp):\n  ", paste(bad, collapse = "\n  "))
  }
  invisible(TRUE)
}

# Load the current environmental pool as one SpatRaster (named by
# pool `layer`). Mirrors checo26_load_current_stack(mode = "all_modeling",
# include_lulc = TRUE): the pool already contains its 5 static lulc
# proportions, so there is no separate LULC merge.
nsdm_load_current_stack <- function(varnames = NULL,
                                    pool_csv = nsdm_default_pool_csv(),
                                    data_root = nsdm_default_data_root(),
                                    check_geom = TRUE,
                                    verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to load environmental rasters.")
  }
  pool <- nsdm_pool_layers(pool_csv)
  if (!is.null(varnames)) {
    missing_vars <- setdiff(varnames, pool$layer)
    if (length(missing_vars)) {
      stop("Requested variables missing from the pool: ",
           paste(missing_vars, collapse = ", "))
    }
    pool <- pool[match(varnames, pool$layer), , drop = FALSE]
  }
  if (isTRUE(check_geom)) {
    nsdm_check_geom(pool$path, file.path(data_root, "mask_china.tif"),
                    label = "pool layer")
  }
  if (isTRUE(verbose)) message("Current environmental stack: ", nrow(pool), " layers")
  r <- terra::rast(pool$path)
  names(r) <- pool$layer
  attr(r, "nsdm_layers") <- pool
  r
}

# ---- Future stacks (resolved scenario table only) ------------------------------

# ---- Analysis domain (background = projection = attribution = maps) ---------

# One binary domain is used by background sampling and projection. This is the
# contract that prevents the old "buffer-trained, national-projected" mismatch.
nsdm_domain_study_area <- function(mask, occ, method) {
  n_cells <- terra::global(mask, "sum", na.rm = TRUE)[[1]][1]
  if (!is.finite(n_cells) || n_cells <= 0L) {
    stop("Analysis domain has no positive cells.")
  }
  structure(
    list(mask = mask, boundary = NULL, n_cells = as.integer(n_cells),
         method = method, buffer_km = NA_real_, n_occurrences = nrow(occ)),
    class = "cast_study_area")
}

nsdm_national_domain <- function(template, national_mask_file, verbose = TRUE) {
  if (!file.exists(national_mask_file)) {
    stop("National land mask not found: ", national_mask_file)
  }
  m <- terra::rast(national_mask_file)
  if (terra::nlyr(m) > 1L) m <- m[[1L]]
  if (!isTRUE(tryCatch(terra::compareGeom(template, m, stopOnError = FALSE,
                                           messages = FALSE),
                       error = function(e) FALSE))) {
    stop("National land mask is not aligned with the current raster grid: ",
         national_mask_file)
  }
  m <- terra::ifel(m > 0, 1L, NA)
  n_cells <- terra::global(m, "sum", na.rm = TRUE)[[1]][1]
  if (verbose) {
    message("Analysis domain: national land mask (", format(n_cells, big.mark = ","),
            " cells)")
  }
  m
}

nsdm_freshwater_river_domain <- function(occ, river_distance_raster, template,
                                         national_mask_file,
                                         river_buffer_m = 100,
                                         occurrence_buffer_km = 10,
                                         min_cells = 1000L, verbose = TRUE) {
  if (!file.exists(river_distance_raster)) {
    stop("Freshwater river-distance raster not found: ", river_distance_raster)
  }
  if (!is.numeric(river_buffer_m) || length(river_buffer_m) != 1L ||
      !is.finite(river_buffer_m) || river_buffer_m < 0) {
    stop("freshwater_river_buffer_m must be one non-negative number.")
  }
  if (!is.numeric(occurrence_buffer_km) || length(occurrence_buffer_km) != 1L ||
      !is.finite(occurrence_buffer_km) || occurrence_buffer_km < 0) {
    stop("freshwater_occurrence_buffer_km must be one non-negative number.")
  }
  d <- terra::rast(river_distance_raster)
  if (terra::nlyr(d) > 1L) d <- d[[1L]]
  if (!isTRUE(tryCatch(terra::compareGeom(template, d, stopOnError = FALSE,
                                           messages = FALSE),
                       error = function(e) FALSE))) {
    stop("Freshwater river-distance raster is not aligned with the current raster grid: ",
         river_distance_raster)
  }
  # River-network skeleton: cells within the stream-distance buffer.
  river <- d <= river_buffer_m
  # Coarse literature coordinates often miss the 1 km stream cell. A local
  # occurrence envelope keeps those records in the domain without reverting to
  # a whole catchment. The river buffer remains the freshwater core.
  pts <- terra::vect(as.matrix(occ[, c("lon", "lat"), drop = FALSE]),
                     type = "points", crs = terra::crs(template))
  local <- terra::rasterize(
    terra::buffer(pts, width = occurrence_buffer_km * 1000),
    template, field = 1L, background = NA)
  local <- terra::ifel(!is.na(local), 1L, NA)
  m <- terra::ifel(river | !is.na(local), 1L, NA)
  if (!is.null(national_mask_file) && file.exists(national_mask_file)) {
    land <- terra::rast(national_mask_file)
    if (terra::nlyr(land) > 1L) land <- land[[1L]]
    if (!isTRUE(tryCatch(terra::compareGeom(template, land, stopOnError = FALSE,
                                             messages = FALSE),
                         error = function(e) FALSE))) {
      stop("National land mask is not aligned with the freshwater river raster.")
    }
    m <- terra::mask(m, land, maskvalues = c(0, NA))
  }
  n_cells <- terra::global(m, "sum", na.rm = TRUE)[[1]][1]
  if (!is.finite(n_cells) || n_cells < min_cells) {
    stop("Freshwater river domain has only ", n_cells,
         " cells; minimum is ", min_cells, ".")
  }
  if (verbose) {
    message("Analysis domain: freshwater river-network buffer (",
            format(river_buffer_m, big.mark = ","), " m) + ",
            occurrence_buffer_km, " km occurrence locality envelope; ",
            format(n_cells, big.mark = ","), " cells")
  }
  m
}

nsdm_custom_domain <- function(mask_file, template, verbose = TRUE) {
  if (!file.exists(mask_file)) stop("analysis_domain_mask_file not found: ", mask_file)
  m <- terra::rast(mask_file)
  if (terra::nlyr(m) > 1L) m <- m[[1L]]
  if (!isTRUE(tryCatch(terra::compareGeom(template, m, stopOnError = FALSE,
                                           messages = FALSE),
                       error = function(e) FALSE))) {
    stop("analysis_domain_mask_file is not aligned with the current raster grid: ",
         mask_file)
  }
  m <- terra::ifel(m > 0, 1L, NA)
  n_cells <- terra::global(m, "sum", na.rm = TRUE)[[1]][1]
  if (!is.finite(n_cells) || n_cells <= 0L) {
    stop("analysis_domain_mask_file has no positive cells after masking: ", mask_file)
  }
  if (verbose) {
    message("Analysis domain: custom mask (", format(n_cells, big.mark = ","),
            " cells)")
  }
  m
}

# ---- Workflow dispatchers (called by 02_pipeline.R) ---------------------------

wf_load_current_stack <- function(cfg, varnames = NULL, verbose = TRUE) {
  nsdm_load_current_stack(varnames = varnames,
                          pool_csv = cfg$nsdm_pool_csv,
                          data_root = cfg$nsdm_data_root,
                          verbose = verbose)
}

wf_analysis_domain <- function(cfg, occ, raster_template, verbose = TRUE) {
  mode <- cfg$domain_mode %||% "national"
  national_mask <- file.path(cfg$nsdm_data_root, "mask_china.tif")
  if (identical(mode, "national")) {
    mask <- nsdm_national_domain(raster_template, national_mask, verbose = verbose)
    source <- national_mask
  } else if (identical(mode, "freshwater_river")) {
    if (is.null(cfg$freshwater_river_distance_raster) ||
        !nzchar(cfg$freshwater_river_distance_raster)) {
      stop("domain_mode = 'freshwater_river' requires ",
           "freshwater_river_distance_raster.")
    }
    mask <- nsdm_freshwater_river_domain(
      occ, cfg$freshwater_river_distance_raster, raster_template,
      national_mask_file = national_mask,
      river_buffer_m = cfg$freshwater_river_buffer_m %||% 100,
      occurrence_buffer_km = cfg$freshwater_occurrence_buffer_km %||% 10,
      verbose = verbose)
    source <- cfg$freshwater_river_distance_raster
  } else if (identical(mode, "custom_mask")) {
    if (is.null(cfg$analysis_domain_mask_file) ||
        !nzchar(cfg$analysis_domain_mask_file)) {
      stop("domain_mode = 'custom_mask' requires analysis_domain_mask_file.")
    }
    mask <- nsdm_custom_domain(cfg$analysis_domain_mask_file, raster_template,
                               verbose = verbose)
    source <- cfg$analysis_domain_mask_file
  } else {
    stop("Unknown domain_mode: ", mode,
         " (use 'national', 'freshwater_river' or 'custom_mask').")
  }
  list(mask = mask, mode = mode, source = source)
}

# Build lazy future stacks from an exact mapping table. Unlisted static
# layers retain their current sources; cast_compare constructs partial corners.
nsdm_build_scenario_stacks <- function(current_stack, selected_vars,
                                       scenario_table_csv, verbose = TRUE) {
  required <- c("scenario", "family", "climate_model", "ssp", "period",
                "climate_year", "lulc_year", "layer", "path", "source")
  tab <- utils::read.csv(scenario_table_csv, stringsAsFactors = FALSE)
  if (!all(required %in% names(tab))) {
    stop("scenario table is missing columns: ",
         paste(setdiff(required, names(tab)), collapse = ", "))
  }
  missing_sel <- setdiff(selected_vars, names(current_stack))
  if (length(missing_sel)) {
    stop("Selected variables are absent from current stack: ",
         paste(missing_sel, collapse = ", "))
  }
  bad_layers <- setdiff(unique(tab$layer), selected_vars)
  if (length(bad_layers)) {
    stop("Scenario table maps layers outside the selected model set: ",
         paste(bad_layers, collapse = ", "))
  }
  bad_paths <- unique(tab$path[!file.exists(tab$path)])
  if (length(bad_paths)) {
    stop("Missing scenario raster(s):\n  ", paste(bad_paths, collapse = "\n  "))
  }
  dup <- duplicated(tab[, c("scenario", "layer")])
  if (any(dup)) {
    stop("Duplicate scenario/layer rows: ",
         paste(unique(paste(tab$scenario[dup], tab$layer[dup], sep = "/")),
               collapse = ", "))
  }
  scenarios <- unique(tab$scenario)
  if (!length(scenarios)) stop("Scenario table has no scenarios.")
  out <- lapply(scenarios, function(scen) {
    combined <- current_stack[[selected_vars]]
    source_tab <- data.frame(
      scenario = scen, layer = selected_vars, source = "current",
      path = vapply(selected_vars, function(v)
        as.character(terra::sources(combined[[v]]))[[1L]], character(1)),
      stringsAsFactors = FALSE)
    rows <- tab[tab$scenario == scen, , drop = FALSE]
    meta_cols <- c("family", "climate_model", "ssp", "period",
                   "climate_year", "lulc_year")
    source_tab <- cbind(
      rows[rep(1L, nrow(source_tab)), meta_cols, drop = FALSE],
      source_tab)
    for (i in seq_len(nrow(rows))) {
      v <- rows$layer[[i]]
      r <- terra::rast(rows$path[[i]])
      if (!terra::compareGeom(r, combined[[v]], stopOnError = FALSE)) {
        stop("Scenario raster does not match current grid: ", rows$path[[i]])
      }
      names(r) <- v
      combined[[v]] <- r
      j <- match(v, source_tab$layer)
      source_tab$source[j] <- rows$source[[i]]
      source_tab$path[j] <- rows$path[[i]]
    }
    combined <- combined[[selected_vars]]
    names(combined) <- selected_vars
    attr(combined, "future_layer_source") <- source_tab
    combined
  })
  names(out) <- scenarios
  source_manifest <- do.call(rbind, lapply(out, attr, "future_layer_source"))
  out <- lapply(out, function(x) { attr(x, "future_layer_source") <- NULL; x })
  names(out) <- scenarios
  attr(out, "future_layer_source") <- source_manifest
  if (isTRUE(verbose)) {
    fam <- unique(tab[, c("scenario", "family")])$family
    cat("Scenario stacks:", length(out), " (",
        paste(names(table(fam)), table(fam), sep = "=", collapse = ", "),
        ")\n", sep = "")
  }
  out
}
