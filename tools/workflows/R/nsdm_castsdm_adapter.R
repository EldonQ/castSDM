# ==============================================================================
# N-SDM data adapter for the castSDM single-species workflow.
#
# Drop-in replacement for the CHECO26 Final helpers
# (R/checo26_castsdm_helpers.R, deprecated with E:/CHECO26/Final). The loader
# and future-stack builder below return EXACTLY the structures 02_pipeline.R
# consumes from the checo26_* functions:
#
#   * nsdm_load_current_stack()  ~ checo26_load_current_stack():
#       a terra::SpatRaster with one layer per pool variable, named by the pool
#       `layer` id, with the layer table attached as attr(., "nsdm_layers")
#       (checo26 used attr "checo26_layers"; nothing downstream reads it).
#   * nsdm_build_future_stacks_with_static()  ~
#       checo26_build_future_stacks_with_static():
#       a named list "<model>_<ssp>_<period>" -> SpatRaster whose layers are
#       ordered and named exactly as `selected_vars`; climate layers come from
#       the requested scenario, all other layers are reused from the current
#       stack (static).
#   * nsdm_study_area()  ~ cast_study_area():
#       an S3 "cast_study_area" object (mask/boundary/n_cells/method/
#       buffer_km/n_occurrences), accepted by cast_background().
#
# Data source: E:/CHECO26/N-SDM/data (1 km Albers China grid, 4861 x 4062).
# Layer manifest: single_species/pool_layers.csv -- the 68
# QC'd layers written by v2_common.R::build_candidate_pool() (all
# status == "kept", all on the mask_china.tif grid).
#
# CRS policy (HARD CONSTRAINT on this machine):
#   In-process PROJ calls segfault here; raw occurrence points use external cs2cs.
# ==============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

nsdm_default_data_root <- function() "E:/CHECO26/N-SDM/data"

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
  "E:/Package/castSDM_workflows/single_species/pool_layers.csv"
}

# N-SDM future bioclim is stored per scenario YEAR (2030/2055/2085). The
# workflow keeps the CHECO26-era sub-period labels in scenario names so every
# downstream artefact (raster prefixes, projection_stats.csv, 03_figures.R
# globs) stays byte-compatible with the old pipeline; each label is the
# CHECO26 sub-period bin the N-SDM year falls into:
#   2030 -> "2021_2040"   2055 -> "2041_2070"   2085 -> "2071_2100"
.nsdm_period_dir <- c("2021_2040" = "2030",
                      "2041_2070" = "2055",
                      "2071_2100" = "2085")

# GCM/ensemble product id used in directory and file names:
# "bcc_csm2_mr" -> "cnclim1km-bcc-csm2-mr", "ensmean" -> "cnclim1km-ensmean".
.nsdm_product_id <- function(model) paste0("cnclim1km-", gsub("_", "-", model))

# ---- Current stack -----------------------------------------------------------

# Read the 68-layer pool manifest. Returns a data.frame with columns
# group/layer/path/is_climate (plus any extra QC columns from the CSV).
# `data_root` rewrites the default path prefix if the data live elsewhere.
nsdm_pool_layers <- function(pool_csv = nsdm_default_pool_csv(),
                             data_root = nsdm_default_data_root()) {
  if (!file.exists(pool_csv)) stop("N-SDM pool layer table not found: ", pool_csv)
  pool <- utils::read.csv(pool_csv, stringsAsFactors = FALSE)
  need <- c("group", "layer", "path", "is_climate")
  if (!all(need %in% names(pool))) {
    stop("Pool table ", pool_csv, " lacks column(s): ",
         paste(setdiff(need, names(pool)), collapse = ", "))
  }
  if ("status" %in% names(pool)) pool <- pool[pool$status %in% c("kept", NA), , drop = FALSE]
  pool <- pool[!duplicated(pool$layer), , drop = FALSE]
  pool$orig_layer <- pool$layer
  pool$layer <- nsdm_sanitize_layer_id(pool$layer)
  default_root <- nsdm_default_data_root()
  if (!identical(data_root, default_root)) {
    pool$path <- sub(default_root, data_root, pool$path, fixed = TRUE)
  }
  missing_paths <- pool$path[!file.exists(pool$path)]
  if (length(missing_paths)) {
    stop("Missing N-SDM pool raster(s):\n  ", paste(missing_paths, collapse = "\n  "))
  }
  rownames(pool) <- NULL
  pool
}

# Layer ids are the single naming choke point of this adapter: the current
# stack, model_data columns, selection/attribution output and future stacks all
# take their names from the pool `layer` column. N-SDM hydro stems contain
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
    stop("Package 'terra' is required to load N-SDM rasters.")
  }
  if (!file.exists(reference)) stop("N-SDM reference mask not found: ", reference)
  ref <- terra::rast(reference)
  bad <- character()
  for (p in paths) {
    ok <- isTRUE(tryCatch(
      terra::compareGeom(ref, terra::rast(p), stopOnError = FALSE, messages = FALSE),
      error = function(e) FALSE))
    if (!ok) bad <- c(bad, p)
  }
  if (length(bad)) {
    stop("N-SDM ", label, "(s) not aligned with ", basename(reference),
         " (reprojection is disabled on this machine; fix the grid upstream,",
         " e.g. gdalwarp):\n  ", paste(bad, collapse = "\n  "))
  }
  invisible(TRUE)
}

# Load the current N-SDM modelling pool as one SpatRaster (68 layers, named by
# pool `layer`). Mirrors checo26_load_current_stack(mode = "all_modeling",
# include_lulc = TRUE): the N-SDM pool already contains its 5 static lulc
# proportions, so there is no separate LULC merge.
nsdm_load_current_stack <- function(varnames = NULL,
                                    pool_csv = nsdm_default_pool_csv(),
                                    data_root = nsdm_default_data_root(),
                                    check_geom = TRUE,
                                    verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to load N-SDM rasters.")
  }
  pool <- nsdm_pool_layers(pool_csv, data_root)
  if (!is.null(varnames)) {
    missing_vars <- setdiff(varnames, pool$layer)
    if (length(missing_vars)) {
      stop("Requested variables missing from N-SDM pool: ",
           paste(missing_vars, collapse = ", "))
    }
    pool <- pool[match(varnames, pool$layer), , drop = FALSE]
  }
  if (isTRUE(check_geom)) {
    nsdm_check_geom(pool$path, file.path(data_root, "mask_china.tif"),
                    label = "pool layer")
  }
  if (isTRUE(verbose)) message("N-SDM current stack: ", nrow(pool), " layers")
  r <- terra::rast(pool$path)
  names(r) <- pool$layer
  attr(r, "nsdm_layers") <- pool
  r
}

# ---- Future stacks -----------------------------------------------------------

# Locate one future bioclim layer:
#   <root>/covariates/glo/bioclim/<product>/scenario/<year>/<ssp>/<var>/
#     glo_bioclim_<product>_scenario_<year>_<ssp>_<var>.tif
nsdm_future_bioclim_path <- function(varname, model, ssp, period,
                                     data_root = nsdm_default_data_root()) {
  year <- unname(.nsdm_period_dir[period])
  if (is.na(year)) {
    stop("Unknown period '", period, "'. Known: ",
         paste(names(.nsdm_period_dir), collapse = ", "))
  }
  product <- .nsdm_product_id(model)
  file.path(data_root, "covariates", "glo", "bioclim", product,
            "scenario", year, ssp, varname,
            sprintf("glo_bioclim_%s_scenario_%s_%s_%s.tif",
                    product, year, ssp, varname))
}

# Build one future SpatRaster per "<model>_<ssp>_<period>": scenario-specific
# bioclim for the selected climate layers (pool is_climate == TRUE: bio01-23
# minus bio07, prec, sunp) + the CURRENT layers for everything else (topo,
# soil, lulc, vege, anthro, hydro, dist are static -- N-SDM provides no future
# scenarios for them). Layer order/names follow `selected_vars` exactly, as
# checo26_build_future_stacks_with_static() did.
nsdm_build_future_stacks_with_static <- function(
    current_stack,
    selected_vars = names(current_stack),
    future_models = c("bcc_csm2_mr", "ensmean"),
    ssp_scenarios = c("ssp126", "ssp245", "ssp585"),
    periods = c("2021_2040", "2041_2070", "2071_2100"),
    pool_csv = nsdm_default_pool_csv(),
    data_root = nsdm_default_data_root(),
    check_geom = TRUE,
    verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to load N-SDM rasters.")
  }
  selected_vars <- intersect(selected_vars, names(current_stack))
  pool <- nsdm_pool_layers(pool_csv, data_root)
  climate_use <- intersect(pool$layer[pool$is_climate], selected_vars)
  # Future bioclim directories/files are keyed by the ORIGINAL variable id;
  # `climate_use` carries the sanitized layer id (identical for bio01-23, but
  # never assume it stays that way for future hyphenated climate products).
  climate_orig <- if ("orig_layer" %in% names(pool)) {
    pool$orig_layer[match(climate_use, pool$layer)]
  } else {
    climate_use
  }
  static_use <- setdiff(selected_vars, climate_use)
  static_stack <- if (length(static_use)) current_stack[[static_use]] else NULL

  scenario_names <- unlist(lapply(future_models, function(model) {
    as.vector(outer(ssp_scenarios, periods, function(ssp, period) {
      paste(model, ssp, period, sep = "_")
    }))
  }), use.names = FALSE)

  # Resolve and verify every future climate file once, before building.
  fut_paths <- list()
  if (length(climate_use)) {
    for (model in future_models) for (ssp in ssp_scenarios) for (period in periods) {
      scen <- paste(model, ssp, period, sep = "_")
      paths <- vapply(climate_orig, nsdm_future_bioclim_path, character(1),
                      model = model, ssp = ssp, period = period,
                      data_root = data_root)
      missing_paths <- paths[!file.exists(paths)]
      if (length(missing_paths)) {
        stop("Missing N-SDM future bioclim raster(s) for ", scen, ":\n  ",
             paste(missing_paths, collapse = "\n  "))
      }
      fut_paths[[scen]] <- paths
    }
    if (isTRUE(check_geom)) {
      nsdm_check_geom(unique(unlist(fut_paths, use.names = FALSE)),
                      file.path(data_root, "mask_china.tif"),
                      label = "future bioclim layer")
    }
  }

  out <- lapply(scenario_names, function(scen) {
    climate_part <- if (length(climate_use)) {
      r <- terra::rast(fut_paths[[scen]])
      names(r) <- climate_use
      r
    }
    combined <- if (is.null(climate_part)) {
      static_stack
    } else if (is.null(static_stack)) {
      climate_part
    } else {
      c(climate_part, static_stack)
    }
    combined <- combined[[selected_vars]]
    names(combined) <- selected_vars
    combined
  })
  names(out) <- scenario_names
  if (isTRUE(verbose)) {
    message("N-SDM future stacks: ", length(out), " scenarios (",
            length(climate_use), " climate + ", length(static_use),
            " static layers each)")
  }
  out
}

# ---- Study area (zero-reprojection cast_study_area equivalent) ---------------

# Accessible-area M for N-SDM data. Returns a genuine "cast_study_area" S3
# object, but built with GEOS-only sf calls: the occurrence coordinates are
# ALREADY on the template's Albers grid, so no st_transform / terra::project
# is ever needed (both segfault on this machine; cast_study_area() itself
# would attempt one identity st_transform). `occ` is the data.frame produced
# by read_occurrences() for input_type "nsdm_occ", i.e. lon/lat hold Albers
# x/y in metres.
nsdm_study_area <- function(occ, raster_template,
                            method = c("buffer", "bbox_buffer", "full"),
                            buffer_km = 300, min_cells = 1000L,
                            verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
  method <- match.arg(method)
  crs_albers <- sf::st_crs(terra::crs(raster_template))
  n_occ <- nrow(occ)
  if (method != "full" && n_occ < 3L) {
    stop("Need at least 3 occurrence points (got ", n_occ, ").")
  }

  build_mask <- function(dist_km) {
    if (method == "full") {
      m <- !is.na(raster_template[[1]])
      m[m == 0] <- NA
      return(m)
    }
    occ_sf <- sf::st_as_sf(occ, coords = c("lon", "lat"), crs = crs_albers)
    base <- if (method == "buffer") {
      sf::st_convex_hull(sf::st_union(occ_sf))
    } else {
      sf::st_as_sfc(sf::st_bbox(occ_sf))
    }
    # Projected CRS: GEOS buffer directly in metres, no equal-area detour.
    boundary <- sf::st_buffer(base, dist = dist_km * 1000)
    m <- terra::mask(raster_template[[1]], terra::vect(boundary))
    terra::classify(m, cbind(-Inf, Inf, 1))
  }

  used_buffer <- if (method == "full") NA_real_ else buffer_km
  mask_r <- build_mask(buffer_km)
  n_cells <- sum(terra::values(mask_r, mat = FALSE) == 1, na.rm = TRUE)
  while (method != "full" && n_cells < min_cells && used_buffer < buffer_km * 3) {
    used_buffer <- used_buffer * 1.5
    if (verbose) {
      message("Study area has only ", n_cells, " cells (< ", min_cells,
              "); expanding buffer to ", round(used_buffer), " km.")
    }
    mask_r <- build_mask(used_buffer)
    n_cells <- sum(terra::values(mask_r, mat = FALSE) == 1, na.rm = TRUE)
  }
  if (n_cells < min_cells) {
    warning("Study area has only ", n_cells, " valid cells (< min_cells = ",
            min_cells, "); proceeding anyway.")
  }

  boundary_sf <- if (method == "full") {
    e <- terra::ext(raster_template)
    coords <- matrix(c(e[1], e[3], e[2], e[3], e[2], e[4], e[1], e[4], e[1], e[3]),
                     ncol = 2, byrow = TRUE)
    sf::st_sfc(sf::st_polygon(list(coords)), crs = crs_albers)
  } else {
    sf::st_buffer(
      if (method == "buffer") {
        sf::st_convex_hull(sf::st_union(
          sf::st_as_sf(occ, coords = c("lon", "lat"), crs = crs_albers)))
      } else {
        sf::st_as_sfc(sf::st_bbox(
          sf::st_as_sf(occ, coords = c("lon", "lat"), crs = crs_albers)))
      },
      dist = used_buffer * 1000)
  }

  if (verbose) {
    message("Study area defined: ", method, " method | buffer: ",
            ifelse(is.na(used_buffer), "N/A", paste0(round(used_buffer), " km")),
            " | valid cells: ", format(n_cells, big.mark = ","),
            " | occurrences: ", n_occ)
  }
  structure(
    list(mask = mask_r, boundary = boundary_sf, n_cells = n_cells,
         method = method, buffer_km = used_buffer, n_occurrences = n_occ),
    class = "cast_study_area"
  )
}

# ---- Workflow dispatchers (called by 02_pipeline.R) ---------------------------

wf_load_current_stack <- function(cfg, verbose = TRUE) {
  nsdm_load_current_stack(pool_csv = cfg$nsdm_pool_csv,
                          data_root = cfg$nsdm_data_root,
                          verbose = verbose)
}

wf_build_future_stacks <- function(cfg, current_stack, selected_vars,
                                   verbose = TRUE) {
  nsdm_build_future_stacks_with_static(
    current_stack, selected_vars = selected_vars,
    future_models = cfg$future_models, ssp_scenarios = cfg$ssp_scenarios,
    periods = cfg$periods, pool_csv = cfg$nsdm_pool_csv,
    data_root = cfg$nsdm_data_root, verbose = verbose)
}

wf_study_area <- function(cfg, occ, raster_template, verbose = TRUE) {
  nsdm_study_area(occ, raster_template, method = cfg$study_area_method,
                  buffer_km = cfg$study_area_buffer_km, verbose = verbose)
}
