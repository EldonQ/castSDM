# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# Case-specific occurrence QC, domain and background preparation.
# Numerical fitting, validation and scenario comparison belong to castSDM.
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

case_file_manifest <- function(stacks) {
  files <- if (inherits(stacks, "SpatRaster")) unique(terra::sources(stacks)) else unique(stacks)
  if (!length(files) || any(!file.exists(files))) stop("Missing file-backed data input.")
  data.frame(path = normalizePath(files, winslash = "/"), size = file.info(files)$size,
    sha256 = vapply(files, function(p) digest::digest(file = p, algo = "sha256"), character(1)),
    stringsAsFactors = FALSE)
}

case_preparation_source <- function(cfg) {
  paths <- file.path(cfg$example_root, c("prepare_data.R", "data_adapter.R"))
  packages <- c("terra", "sf", "data.table", "FNN")
  list(files = unname(tools::md5sum(paths)),
    code = castSDM:::.cast_code_dependencies(c("cast_thin", "cast_background")),
    versions = vapply(packages, function(p) if (requireNamespace(p, quietly = TRUE))
      as.character(utils::packageVersion(p)) else "not_installed", character(1)))
}

ensure_workflow_tempdir <- function() {
  td <- tempdir()
  if (!dir.exists(td)) dir.create(td, recursive = TRUE)
  if (file.access(td, 2L) != 0L) stop("R temporary directory is not writable.")
}

gbif_taxon_ok <- function(raw, target) {
  target <- tolower(trimws(as.character(target)[1]))
  if (!nzchar(target)) return(rep(FALSE, nrow(raw)))
  if ("target_taxon" %in% names(raw)) {
    tt <- raw$target_taxon
    if (is.logical(tt) || is.numeric(tt)) {
      return(!is.na(tt) & as.logical(tt))
    }
    return(tolower(trimws(as.character(tt))) == target)
  }
  candidates <- intersect(c("species", "genus", "family", "subfamily"), names(raw))
  if (!length(candidates)) return(rep(FALSE, nrow(raw)))
  ok <- rep(FALSE, nrow(raw))
  for (col in candidates) {
    v <- tolower(trimws(as.character(raw[[col]])))
    ok <- ok | (!is.na(v) & (v == target | startsWith(v, paste0(target, " "))))
  }
  ok
}

case_year_checks <- function(year, cfg) {
  limits <- cfg$occurrence_year_range
  if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2L ||
      any(!is.finite(limits)) || limits[1L] > limits[2L]))
    stop("occurrence_year_range must be NULL or c(first_year, last_year).")
  known <- is.finite(year) & year > 0
  list(outside_declared_years = if (is.null(limits)) rep(FALSE, length(year)) else
    known & (year < limits[1L] | year > limits[2L]),
    unknown_year = !known & !isTRUE(cfg$occurrence_keep_unknown_year %||% TRUE))
}

case_source_metadata <- function(raw, cfg) {
  file <- cfg$occurrence_metadata_file
  if (is.null(file) || !file.exists(file) || !"datasetKey" %in% names(raw)) return(raw)
  use <- which(grepl("^SchiSOFT_24638538", raw$datasetKey))
  if (!length(use)) return(raw)
  if (!all(c("gbifID", "species", "decimalLongitude", "decimalLatitude") %in% names(raw)))
    stop("SchiSOFT metadata matching requires IDs, species names and original coordinates.")
  original <- data.table::fread(file, data.table = FALSE, colClasses = "character", na.strings = c("NA", ""))
  required <- c("occurrenceID", "scientificName", "genus", "taxonomicStatus", "longitudeX", "latitudeY", "sampledDate", "remarksOfSampledDate", "sourceType")
  if (!all(required %in% names(original)) || anyDuplicated(original$occurrenceID))
    stop("Invalid original SchiSOFT metadata: ", file)
  at <- match(as.character(raw$gbifID[use]), original$occurrenceID)
  if (anyNA(at) || any(raw$species[use] != original$scientificName[at]) ||
      any(abs(as.numeric(raw$decimalLongitude[use])-as.numeric(original$longitudeX[at])) > 1e-6) ||
      any(abs(as.numeric(raw$decimalLatitude[use])-as.numeric(original$latitudeY[at])) > 1e-6))
    stop("SchiSOFT IDs, names or coordinates do not match the original source. No metadata were inferred.")
  fields <- c(sampled_date = "sampledDate", sampling_date_remarks = "remarksOfSampledDate", source_type = "sourceType")
  for (name in names(fields)) {
    raw[[name]] <- NA_character_
    raw[[name]][use] <- original[[fields[[name]]]][at]
  }
  if (!"genus" %in% names(raw)) raw$genus <- NA_character_
  if (!"taxonRank" %in% names(raw)) raw$taxonRank <- NA_character_
  if (!"subfamily" %in% names(raw)) raw$subfamily <- NA_character_
  raw$genus[use] <- original$genus[at]
  raw$taxonRank[use] <- toupper(original$taxonomicStatus[at])
  raw$subfamily[use] <- "Schizothoracinae"
  raw$sampling_date_status <- "unverified_record_year"
  single <- !is.na(raw$sampled_date) & grepl("^[0-9]{4}([0-9]{2}){0,2}$", raw$sampled_date)
  actual <- use[is.na(raw$sampling_date_remarks[use]) & single[use]]
  raw$sampling_date_status[use] <- ifelse(is.na(raw$sampled_date[use]), "unknown_sampling_date",
    ifelse(!is.na(raw$sampling_date_remarks[use]), "source_date_substitute", "reported_sampling_interval"))
  raw$sampling_date_status[actual] <- "reported_sampling_date"
  raw$sampling_year_verified <- NA_real_
  raw$sampling_year_verified[actual] <- as.numeric(substr(raw$sampled_date[actual], 1L, 4L))
  intervals <- use[is.na(raw$sampling_date_remarks[use]) & !single[use] & !is.na(raw$sampled_date[use])]
  for (j in intervals) {
    ends <- strsplit(raw$sampled_date[j], "_", fixed = TRUE)[[1L]]
    if (length(ends) == 2L && all(grepl("^[0-9]{4}([0-9]{2}){0,2}$", ends)) &&
        length(unique(substr(ends, 1L, 4L))) == 1L) {
      raw$sampling_year_verified[j] <- as.numeric(substr(ends[1L], 1L, 4L))
      raw$sampling_date_status[j] <- "reported_interval_within_one_year"
    }
  }
  # This is a precision warning, not a claim that a coordinate was rounded.
  raw$qc_tenth_degree_grid <- abs(as.numeric(raw$decimalLongitude)*10-round(as.numeric(raw$decimalLongitude)*10)) < 1e-7 &
    abs(as.numeric(raw$decimalLatitude)*10-round(as.numeric(raw$decimalLatitude)*10)) < 1e-7
  raw$source_metadata_file <- raw$source_metadata_md5 <- NA_character_
  raw$source_metadata_file[use] <- normalizePath(file, winslash = "/")
  raw$source_metadata_md5[use] <- unname(tools::md5sum(file))
  message("Original SchiSOFT metadata matched: ", length(use), " records; reported single sampling dates: ", length(actual),
    "; additional intervals within one calendar year: ", sum(raw$sampling_date_status == "reported_interval_within_one_year"),
    ". Multi-year intervals and substitute dates are not converted to precise sampling years.")
  raw
}

read_occurrences <- function(cfg) {
  if (!length(cfg$occurrence_file) || any(!file.exists(cfg$occurrence_file))) {
    stop("Occurrence files are missing for ", cfg$species_name)
  }
  raw <- as.data.frame(data.table::rbindlist(lapply(cfg$occurrence_file, function(p) {
    d <- data.table::fread(p, data.table = FALSE,
                          colClasses = if (cfg$input_type == "gbif") "character" else NULL)
    names(d) <- sub("^\ufeff", "", names(d))
    d$source_file <- normalizePath(p, winslash = "/")
    d
  }), fill = TRUE))
  raw <- case_source_metadata(raw, cfg)
  if (cfg$input_type == "gbif") {
    required <- c("gbifID", "decimalLongitude", "decimalLatitude",
                  "countryCode", "occurrenceStatus", "hasGeospatialIssues",
                  "coordinateUncertaintyInMeters", "basisOfRecord")
    if (!all(required %in% names(raw))) stop("Missing GBIF fields: ",
      paste(setdiff(required, names(raw)), collapse = ", "))
    taxon_cols <- c("species", "genus", "family", "subfamily", "target_taxon")
    if (!any(taxon_cols %in% names(raw))) {
      stop("GBIF input lacks a taxon column; expected one of: ",
           paste(taxon_cols, collapse = ", "))
    }
    for (v in c("decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters")) {
      raw[[v]] <- as.numeric(raw[[v]])
    }
    raw$hasGeospatialIssues <- tolower(raw$hasGeospatialIssues) %in% c("true", "t", "1")
    taxon_ok <- gbif_taxon_ok(raw, cfg$species_name)
    if (!any(taxon_ok) && nrow(raw) > 0L) {
      stop("No GBIF record matches target taxon '", cfg$species_name,
           "' in columns ", paste(intersect(taxon_cols, names(raw)), collapse = ", "),
           ". For subfamily-level raw GBIF data, provide a subfamily or target_taxon ",
           "column (or preprocess the target label).")
    }
    basis_norm <- toupper(trimws(as.character(raw$basisOfRecord)))
    accepted_basis <- toupper(as.character(cfg$occurrence_basis_of_record %||%
      c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OBSERVATION",
        "PRESERVED_SPECIMEN", "MATERIAL_SAMPLE", "MATERIAL_CITATION",
        "LITERATURE", "OCCURRENCE")))
    unsupported_basis <- if (length(accepted_basis)) {
      is.na(basis_norm) | !basis_norm %in% accepted_basis
    } else {
      rep(FALSE, nrow(raw))
    }
    status_norm <- tolower(trimws(as.character(raw$occurrenceStatus)))
    country_norm <- toupper(trimws(as.character(raw$countryCode)))
    geo_fields <- intersect(c("georeferenceSources", "georeferenceRemarks", "locality"), names(raw))
    geo <- if (length(geo_fields)) {
      do.call(paste, c(as.list(raw[, geo_fields, drop = FALSE]), sep = " "))
    } else {
      rep("", nrow(raw))
    }
    raw$qc_centroid <- grepl("centroid|country cent", geo, ignore.case = TRUE)
    raw$qc_unknown_uncertainty <- is.na(raw$coordinateUncertaintyInMeters)
    valid_xy <- is.finite(raw$decimalLongitude) & is.finite(raw$decimalLatitude) &
      abs(raw$decimalLongitude) <= 180 & abs(raw$decimalLatitude) < 90
    raw$year <- if ("year" %in% names(raw)) suppressWarnings(as.numeric(raw$year)) else NA_real_
    checks <- list(
      invalid_coordinates = !valid_xy,
      wrong_taxon = !taxon_ok,
      outside_country = is.na(country_norm) | country_norm != "CN",
      not_present = is.na(status_norm) | status_norm != "present",
      unsupported_basis = unsupported_basis,
      geospatial_issue = raw$hasGeospatialIssues %in% TRUE,
      centroid_coordinates = raw$qc_centroid,
      excessive_uncertainty = !raw$qc_unknown_uncertainty &
        (raw$coordinateUncertaintyInMeters < 0 |
         raw$coordinateUncertaintyInMeters > cfg$occurrence_uncertainty_m),
      unknown_uncertainty = raw$qc_unknown_uncertainty & !cfg$occurrence_keep_unknown_uncertainty)
    checks <- c(checks, case_year_checks(raw$year, cfg))
    raw$qc_reason <- "kept"
    for (reason in names(checks)) raw$qc_reason[raw$qc_reason == "kept" & checks[[reason]]] <- reason
    mask <- terra::rast(file.path(cfg$nsdm_data_root, "mask_china.tif"))
    raw$x_albers <- raw$y_albers <- raw$grid_cell <- NA_real_
    eligible <- which(raw$qc_reason == "kept")
    if (length(eligible)) {
      xy <- nsdm_project_occurrences(raw[eligible, c("decimalLongitude", "decimalLatitude")],
                                    terra::crs(mask, proj = TRUE))
      raw$x_albers[eligible] <- xy[, 1]
      raw$y_albers[eligible] <- xy[, 2]
      raw$grid_cell[eligible] <- terra::cellFromXY(mask, xy)
      mask_values <- terra::extract(mask, xy)[[1]]
      raw$qc_reason[eligible[!is.finite(mask_values) | mask_values <= 0]] <- "outside_mask"
    }
    eligible <- which(raw$qc_reason == "kept")
    eligible <- eligible[order(raw$coordinateUncertaintyInMeters[eligible],
                              raw$gbifID[eligible], na.last = TRUE)]
    raw$qc_reason[eligible[duplicated(raw$gbifID[eligible])]] <- "duplicate_record"
    eligible <- eligible[raw$qc_reason[eligible] == "kept"]
    raw$qc_reason[eligible[duplicated(raw$grid_cell[eligible])]] <- "duplicate_cell"
    kept <- raw$qc_reason == "kept"
    fields <- intersect(c("source_file", "gbifID", "datasetKey", "species", "genus", "subfamily",
      "taxonRank", "basisOfRecord", "decimalLongitude", "decimalLatitude", "year", "eventDate",
      "coordinateUncertaintyInMeters", "georeferenceSources", "georeferenceRemarks",
      "locality", "countryCode", "occurrenceStatus", "hasGeospatialIssues",
      "sampled_date", "sampling_date_remarks", "sampling_date_status", "sampling_year_verified",
      "source_type", "source_metadata_file", "source_metadata_md5", "qc_tenth_degree_grid",
      "x_albers", "y_albers", "grid_cell", "qc_centroid", "qc_unknown_uncertainty", "qc_reason"), names(raw))
    out <- data.frame(lon = raw$x_albers[kept], lat = raw$y_albers[kept])
    attr(out, "occurrence_audit") <- raw[, fields, drop = FALSE]
    return(out)
  } else {
    lon <- grep("^lon|^longitude$", names(raw), ignore.case = TRUE, value = TRUE)[1]
    lat <- grep("^lat", names(raw), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(lon) || is.na(lat)) stop("Longitude/latitude columns not found: ", cfg$occurrence_file)
    year_col <- grep("^year$", names(raw), ignore.case = TRUE, value = TRUE)[1]
    year <- if (is.na(year_col)) rep(NA_integer_, nrow(raw)) else
      suppressWarnings(as.integer(raw[[year_col]]))
    out <- data.frame(
      lon = as.numeric(raw[[lon]]), lat = as.numeric(raw[[lat]]),
      year = year, source_file = raw$source_file, stringsAsFactors = FALSE)
    occurrence_audit <- data.frame(
      source_file = out$source_file,
      x_albers = out$lon, y_albers = out$lat,
      lon = out$lon, lat = out$lat, year = out$year,
      grid_cell = NA_integer_, qc_reason = "kept", used_in_model = FALSE,
      stringsAsFactors = FALSE)
    checks <- c(list(invalid_coordinates = !is.finite(out$lon) | !is.finite(out$lat)),
      case_year_checks(out$year, cfg))
    for (reason in names(checks)) occurrence_audit$qc_reason[
      occurrence_audit$qc_reason == "kept" & checks[[reason]]] <- reason
    eligible <- which(occurrence_audit$qc_reason == "kept")
    occurrence_audit$qc_reason[eligible[duplicated(out[eligible, c("lon", "lat")])]] <- "duplicate_coordinate"
    out <- out[occurrence_audit$qc_reason == "kept", c("lon", "lat"), drop = FALSE]
    attr(out, "occurrence_audit") <- occurrence_audit
    return(out)
  }
}

prep_fingerprint_inputs <- function(cfg, raster_manifest,
                                    target_group_file = character(0),
                                    target_group_md5 = NA_character_) {
  domain_mode <- cfg$domain_mode %||% "national"
  domain_source <- switch(
    domain_mode,
    national = file.path(cfg$nsdm_data_root, "mask_china.tif"),
    freshwater_river = cfg$freshwater_river_distance_raster %||% NA_character_,
    custom_mask = cfg$analysis_domain_mask_file %||% NA_character_,
    NA_character_)
  domain_md5 <- if (length(domain_source) &&
                    all(!is.na(domain_source)) &&
                    all(nzchar(domain_source)) &&
                    all(file.exists(domain_source))) {
    unname(tools::md5sum(domain_source))
  } else {
    NA_character_
  }
  list(
    species = cfg$species_name,
    occurrence = tools::md5sum(cfg$occurrence_file),
    occurrence_metadata = if (!is.null(cfg$occurrence_metadata_file) && file.exists(cfg$occurrence_metadata_file))
      tools::md5sum(cfg$occurrence_metadata_file),
    pool = tools::md5sum(cfg$nsdm_pool_csv),
    rasters = raster_manifest,
    target_group = list(
      enabled = isTRUE(cfg$bg_target_group),
      file = if (length(target_group_file)) as.character(target_group_file)[1] else NA_character_,
      md5 = if (length(target_group_md5)) as.character(target_group_md5)[1] else NA_character_
    ),
    analysis_domain = list(
      mode = domain_mode,
      source = domain_source,
      md5 = domain_md5
    ),
    background_mode = paste0(domain_mode, "_",
                             if (isTRUE(cfg$bg_target_group)) "target_group" else
                               cfg$bg_strategy),
    preprocessing_scope = cfg$preprocessing_scope %||% "unspecified",
    config = cfg[c("input_type", "occurrence_uncertainty_m", "occurrence_keep_unknown_uncertainty",
                   "occurrence_year_range", "occurrence_keep_unknown_year", "occurrence_year_semantics",
                   "occurrence_thin_km", "occurrence_thin_order", "occurrence_basis_of_record",
                   "preprocessing_scope",
                   "domain_mode", "freshwater_river_distance_raster",
                   "freshwater_river_buffer_m", "freshwater_occurrence_buffer_km",
                   "analysis_domain_mask_file", "domain_min_presence_coverage",
                   "bg_n", "predictor_min_coverage", "predictor_missing",
                   "bg_strategy", "bg_bin_k", "bg_sre_quantile",
                   "bg_target_group", "bg_target_group_file", "seed")],
    source = case_preparation_source(cfg)
  )
}

analysis_domain_label <- function(mode) {
  switch(
    mode,
    national = "national China land mask",
    freshwater_river = "river-containing grid cells within observed-locality scope",
    custom_mask = "custom analysis-domain mask",
    stop("Unknown domain_mode: ", mode))
}

case_preparation_inputs <- function(cfg, current) {
  # A target-group background is a different sampling design, not a fallback
  # option. Silently reverting to random background would mislabel every
  # downstream product, so an enabled-but-unusable source is a hard error.
  tg_file <- if (isTRUE(cfg$bg_target_group)) cfg$bg_target_group_file %||% character(0) else character(0)
  if (isTRUE(cfg$bg_target_group) &&
      (length(tg_file) != 1L || is.na(tg_file) || !nzchar(tg_file) || !file.exists(tg_file))) {
    stop("bg_target_group = TRUE but bg_target_group_file is missing or unreadable: ",
         if (length(tg_file)) tg_file else "<not set>")
  }
  tg_md5 <- if (length(tg_file) == 1L && file.exists(tg_file)) unname(tools::md5sum(tg_file)) else NA_character_
  prep_fingerprint_inputs(cfg, case_file_manifest(current), tg_file, tg_md5)
}

# Prepare the case data once, independently of fitted models and future inputs.
prepare_species <- function(cfg, current, prep_inputs) {
  ensure_workflow_tempdir()
  tg_file <- prep_inputs$target_group$file
  preparation_parent <- file.path(cfg$preparation_root, cfg$species_name)
  dir.create(preparation_parent, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(preparation_parent, ".preparation.lock"), timeout = 0)
  if (is.null(lock)) stop("Another process is preparing this case: ", preparation_parent)
  on.exit(filelock::unlock(lock), add = TRUE)
  prep_fp <- digest::digest(prep_inputs, algo = "sha256")
  # A new cohort/domain never replaces the preparation referenced by an old run.
  out_dir <- file.path(preparation_parent, prep_fp)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data_file <- file.path(out_dir, paste0(cfg$species_name, "_model_data.csv"))
  domain_file <- file.path(out_dir, "analysis_domain_mask.tif")
  manifest_file <- file.path(out_dir, "preparation_manifest.rds")
  old <- if (file.exists(manifest_file)) tryCatch(readRDS(manifest_file), error = function(e) NULL) else NULL
  artifacts <- c(data_file, domain_file, file.path(out_dir,
    c("occurrence_qc.csv", "occurrence_qc_counts.csv", "analysis_domain_summary.csv",
      "predictor_qualification.csv", "predictor_missing_qc.csv", "occurrence_environment_qc.csv", "background_summary.csv")))
  if (!is.null(old) && identical(old$fingerprint, prep_fp) && all(file.exists(artifacts)) &&
      identical(unname(tools::md5sum(artifacts)), old$artifacts)) {
    message("Verified preparation cache: ", cfg$species_name)
    return(invisible(out_dir))
  }
  n_layers <- terra::nlyr(current)
  occ <- read_occurrences(cfg)
  occurrence_audit <- attr(occ, "occurrence_audit")
  # The nsdm_occ branch has unique coordinates but not yet unique grid cells.
  # Collapse one record per 1 km cell before thinning so the audit and the
  # model data use the same cell universe.
  if (!is.null(occurrence_audit) && nrow(occ) > 0L) {
    cells <- terra::cellFromXY(current[[1]],
      as.matrix(occ[, c("lon", "lat")]))
    dup <- duplicated(cells)
    if (any(dup)) {
      keep_occ <- which(occurrence_audit$qc_reason == "kept")
      occurrence_audit$qc_reason[keep_occ[dup]] <- "duplicate_cell"
      occ <- occ[!dup, , drop = FALSE]
    }
  }
  if (!is.null(occurrence_audit)) {
    data.table::fwrite(occurrence_audit, file.path(out_dir, "occurrence_qc.csv"))
    data.table::fwrite(as.data.frame(table(occurrence_audit$qc_reason)),
                      file.path(out_dir, "occurrence_qc_counts.csv"))
  }
  if (nrow(occ) < 10L) stop("Fewer than 10 unique eligible occurrence cells; inspect occurrence_qc.csv.")
  # --- Occurrence spatial thinning ---------------------------------------------
  # Grid deduplication (one record per 1 km cell) is not distance thinning:
  # adjacent-cell records can sit metres apart, so spatial autocorrelation
  # survives it. Greedy minimum-distance thinning is the standard remedy
  # (Aiello-Lammens et al. 2015 spThin).
  thin_km <- cfg$occurrence_thin_km %||% 0
  if (thin_km > 0 && nrow(occ) > 1L) {
    thin_year <- if (!is.null(occurrence_audit) &&
                     identical(cfg$occurrence_thin_order, "year_desc")) {
      yr <- occurrence_audit$year[occurrence_audit$qc_reason == "kept"]
      if (length(yr) == nrow(occ)) yr else NULL
    } else NULL
    candidates <- data.frame(occ, .case_row = seq_len(nrow(occ)),
      .case_year = if (is.null(thin_year)) NA_real_ else thin_year)
    thinned <- castSDM::cast_thin(candidates, min_dist_km = thin_km, projected = TRUE,
      order = cfg$occurrence_thin_order, year_col = ".case_year", seed = cfg$seed, verbose = FALSE)
    keep_i <- thinned$.case_row
    n_before <- nrow(occ)
    occ <- occ[keep_i, , drop = FALSE]
    if (!is.null(occurrence_audit)) {
      eligible <- which(occurrence_audit$qc_reason == "kept")
      occurrence_audit$qc_reason[eligible[setdiff(seq_len(n_before), keep_i)]] <- "spatially_thinned"
    }
    message("Spatial thinning: ", n_before, " -> ", nrow(occ), " records at ", thin_km, " km.")
    if (nrow(occ) < 10L) stop("Fewer than 10 occurrence cells survive thinning; lower occurrence_thin_km.")
  }
  # One analysis domain for background, projection, attribution and maps.
  domain_info <- wf_analysis_domain(cfg, occ, current[[1]], verbose = TRUE)
  domain <- domain_info$mask
  terra::writeRaster(domain, domain_file, overwrite = TRUE,
                     wopt = list(datatype = "INT1U", gdal = "COMPRESS=LZW"))
  # Occurrence cells outside the domain cannot inform a model that is fitted
  # and projected on that domain; exclude them and record the QC reason.
  in_domain <- terra::extract(domain, as.matrix(occ[, c("lon", "lat"), drop = FALSE]))[[1L]]
  in_domain <- is.finite(in_domain) & in_domain == 1L
  occurrence_coverage_input <- mean(in_domain)
  if (!is.finite(occurrence_coverage_input) ||
      (occurrence_coverage_input < cfg$domain_min_presence_coverage &&
       !identical(cfg$domain_mode, "custom_mask"))) {
    stop(sprintf(
      paste0("Analysis domain covers only %.1f%% of occurrence cells (%d/%d); ",
             "minimum is %.1f%%. Use a broader domain mask or fix the occurrence ",
             "coordinates before rerunning."),
      100 * occurrence_coverage_input, sum(in_domain), length(in_domain),
      100 * cfg$domain_min_presence_coverage), call. = FALSE)
  }
  n_outside_domain <- sum(!in_domain)
  if (n_outside_domain > 0L) {
    if (!is.null(occurrence_audit)) {
      kept_idx <- which(occurrence_audit$qc_reason == "kept")
      occurrence_audit$qc_reason[kept_idx[!in_domain]] <- "outside_domain"
    }
    message(n_outside_domain, " occurrence cell(s) outside the analysis domain were excluded.")
    occ <- occ[in_domain, , drop = FALSE]
  }
  if (nrow(occ) < 10L) {
    stop("Fewer than 10 occurrence cells remain inside the analysis domain.")
  }
  qualification <- case_domain_qualification(current, domain, cfg$predictor_min_coverage %||% .95,
    manifest = prep_inputs$rasters, cache_dir = file.path(cfg$preparation_root, ".coverage"))
  data.table::fwrite(qualification, file.path(out_dir, "predictor_qualification.csv"))
  kept <- qualification$layer[qualification$eligible]
  if (!length(kept)) stop("No predictor meets declared coverage on the actual analysis mask.")
  current <- current[[kept]]
  message("Actual-domain coverage: ", n_layers, " candidates -> ", length(kept),
    "; excluded: ", paste(qualification$layer[!qualification$eligible], collapse = ", "))
  domain_area <- nsdm_domain_study_area(domain, occ, method = domain_info$mode)
  domain_cells_value <- as.integer(domain_area$n_cells)
  # Optional target-group locations must be supplied for a comparable survey
  # group; no all-taxa source is treated automatically as fish survey effort.
  tg_table <- NULL
  if (isTRUE(cfg$bg_target_group) && length(cfg$bg_target_group_file) &&
      file.exists(cfg$bg_target_group_file)) {
    tg_table <- {
      tg <- data.table::fread(cfg$bg_target_group_file, sep = "|",
                              select = c("X", "Y"))
      names(tg) <- c("lon", "lat")
      # Snap to the 1 km template grid and keep unique valid cells.
      tmpl <- current[[1]]
      res_x <- terra::res(tmpl)[1]; res_y <- terra::res(tmpl)[2]
      col <- floor((tg$lon - terra::xmin(tmpl)) / res_x) + 1
      row <- floor((tg$lat - terra::ymin(tmpl)) / res_y) + 1
      ncol_g <- terra::ncol(tmpl); nrow_g <- terra::nrow(tmpl)
      inside <- col >= 1 & col <= ncol_g & row >= 1 & row <= nrow_g &
        is.finite(tg$lon) & is.finite(tg$lat)
      cell <- (row[inside] - 1) * ncol_g + col[inside]
      # National land mask: keep cells that are non-NA in the stack.
      land <- !is.na(terra::values(tmpl, mat = FALSE))
      valid_cell <- sort(unique(cell[land[cell] %in% TRUE]))
      cell_xy <- data.frame(
        lon = terra::xFromCol(tmpl, (valid_cell - 1) %% ncol_g + 1),
        lat = terra::yFromRow(tmpl, (valid_cell - 1) %/% ncol_g + 1))
      message(sprintf(
        "Target-group background pool: %s records -> %d unique 1 km cells on the national mask.",
        format(nrow(tg), big.mark = ","), nrow(cell_xy)))
      cell_xy
    }
  }
  if (!is.null(tg_table)) {
    tg_cells <- terra::cellFromXY(domain, as.matrix(tg_table[, c("lon", "lat")]))
    in_domain_tg <- !is.na(tg_cells) &
      terra::values(domain, mat = FALSE)[tg_cells] == 1
    in_domain_tg[is.na(in_domain_tg)] <- FALSE
    tg_table <- tg_table[in_domain_tg, , drop = FALSE]
    if (!nrow(tg_table)) {
      stop("Target-group background has no points inside the analysis domain.")
    }
    message("Target-group background restricted to the analysis domain: ",
            nrow(tg_table), " candidate cells.")
  }
  model_data <- castSDM::cast_background(
    occ, domain_area, current, strategy = cfg$bg_strategy,
    bin_k = cfg$bg_bin_k, sre_quantile = cfg$bg_sre_quantile,
    user_table = tg_table,
    n_bg = cfg$bg_n, missing = cfg$predictor_missing %||% "complete_cases",
    cell_thin = TRUE, exclude_presence = TRUE, seed = cfg$seed, verbose = TRUE
  )
  data.table::fwrite(attr(model_data, "missing_qc"), file.path(out_dir, "predictor_missing_qc.csv"))
  location_qc <- attr(model_data, "presence_qc")
  data.table::fwrite(location_qc, file.path(out_dir, "occurrence_environment_qc.csv"))
  observed <- attr(model_data, "missing_qc")
  at <- match(qualification$layer, observed$variable)
  qualification$n_presence_checked <- observed$n_presence[at]
  qualification$n_presence_missing <- observed$n_presence_missing[at]
  qualification$presence_coverage <- observed$presence_coverage[at]
  data.table::fwrite(qualification, file.path(out_dir, "predictor_qualification.csv"))
  current <- current[[intersect(names(current), names(model_data))]]
  data.table::fwrite(data.frame(requested_background = cfg$bg_n,
    actual_background = sum(model_data$presence == 0L),
    actual_presence = sum(model_data$presence == 1L),
    ratio = sum(model_data$presence == 0L)/sum(model_data$presence == 1L),
    strategy = cfg$bg_strategy), file.path(out_dir, "background_summary.csv"))
  vars <- castSDM::get_env_vars(model_data)
  missing_file <- file.path(out_dir, "predictor_missing_qc.csv")
  if (file.exists(missing_file)) {
    missing_qc <- data.table::fread(missing_file, data.table = FALSE)
    affected <- missing_qc[missing_qc$n_presence_missing > 0, , drop = FALSE]
    if (nrow(affected)) {
      affected <- affected[order(-affected$n_presence_missing), , drop = FALSE]
      message("Occurrence missingness by layer (counts overlap): ",
        paste(paste0(utils::head(affected$variable, 3), "=", utils::head(affected$n_presence_missing, 3)), collapse = ", "),
        ". Missingness may be structural or unobserved; no values are filled.")
    }
  }
  dropped <- setdiff(setdiff(names(model_data), c("lon", "lat", "presence")), vars)
  if (length(dropped)) {
    warning(sprintf(
      "%d/%d pool layers are constant in the complete prepared sample (zero variance or all NA) and were dropped; first few: %s",
      length(dropped), terra::nlyr(current), paste(utils::head(dropped, 15), collapse = ", ")),
      call. = FALSE)
    data.table::fwrite(data.frame(layer = dropped),
                       file.path(out_dir, "zero_variance_layers.csv"))
  }
  if (!is.null(occurrence_audit)) {
    eligible <- which(occurrence_audit$qc_reason == "kept")
    idx <- match(terra::cellFromXY(current[[1]], as.matrix(occ[, c("lon", "lat")])), location_qc$cell)
    if (anyNA(idx)) stop("Occurrence QC does not match the sampled grid cells.")
    occurrence_audit$missing_predictors <- ""
    occurrence_audit$missing_predictors[eligible] <- location_qc$missing_predictors[idx]
    pres_cells <- terra::cellFromXY(current,
      as.matrix(model_data[model_data$presence == 1L, c("lon", "lat")]))
    if (!"grid_cell" %in% names(occurrence_audit) ||
        all(is.na(occurrence_audit$grid_cell))) {
      occurrence_audit$grid_cell <- terra::cellFromXY(current[[1]],
        cbind(occurrence_audit$x_albers, occurrence_audit$y_albers))
    }
    occurrence_audit$used_in_model <- FALSE
    occurrence_audit$used_in_model[eligible] <- occurrence_audit$grid_cell[eligible] %in% pres_cells
    occurrence_audit$qc_reason[eligible[!location_qc$used_in_model[idx]]] <- "missing_environment"
    data.table::fwrite(occurrence_audit, file.path(out_dir, "occurrence_qc.csv"))
    data.table::fwrite(as.data.frame(table(occurrence_audit$qc_reason)),
                      file.path(out_dir, "occurrence_qc_counts.csv"))
    stopifnot(sum(occurrence_audit$used_in_model) == sum(model_data$presence == 1L))
  }
  if (sum(model_data$presence == 1L) < 10L) {
    stop("Fewer than 10 complete occurrence cells remain after environmental extraction.")
  }
  data.table::fwrite(model_data, data_file)
  presence_xy <- as.matrix(model_data[model_data$presence == 1L, c("lon", "lat")])
  presence_in_domain <- if (nrow(presence_xy)) {
    terra::extract(domain, presence_xy)[[1]]
  } else {
    numeric(0)
  }
  presence_coverage <- if (length(presence_in_domain)) {
    mean(presence_in_domain == 1, na.rm = TRUE)
  } else {
    NA_real_
  }
  national_cells <- terra::global(terra::rast(file.path(cfg$nsdm_data_root, "mask_china.tif")) > 0, "sum", na.rm = TRUE)[[1L]]
  domain_summary <- data.frame(
    domain_mode = domain_info$mode,
    domain_label = analysis_domain_label(domain_info$mode),
    domain_source = paste(domain_info$source, collapse = ";"),
    domain_mask_file = domain_file,
    domain_mask_md5 = unname(tools::md5sum(domain_file)),
    domain_cells = domain_cells_value,
    national_mask_cells = national_cells, domain_share_national = domain_cells_value/national_cells,
    area_interpretation = "Coverage of native grid cells within the declared study domain; not species occupancy or physical habitat surface area.",
    occurrence_coverage_input = occurrence_coverage_input,
    occurrence_excluded_outside_domain = n_outside_domain,
    pool_layers = as.integer(n_layers),
    prepared_pool_columns = length(vars),
    presence_cells = sum(model_data$presence == 1L),
    background_cells = sum(model_data$presence == 0L),
    presence_in_domain = sum(presence_in_domain == 1, na.rm = TRUE),
    presence_coverage = presence_coverage,
    stringsAsFactors = FALSE)
  data.table::fwrite(domain_summary,
                     file.path(out_dir, "analysis_domain_summary.csv"))
  message("Analysis domain covers ", sprintf("%.2f%%", 100*domain_cells_value/national_cells),
    " of national-mask cells. Blank areas outside this domain were not predicted; coloured grid area is not occupied habitat area.")
  if (!is.finite(presence_coverage) ||
      presence_coverage < cfg$domain_min_presence_coverage) {
    stop(sprintf(
      paste0("Analysis domain covers only %.1f%% of presence cells (%d/%d); ",
             "minimum is %.1f%%. The supplied domain excludes occurrence ",
             "records and would invalidate fitting. Use a domain that covers ",
             "the data (larger water buffer, river-network/catchment mask, or ",
             "the accessible land area) and rerun with a new output_root."),
      100 * presence_coverage, sum(presence_in_domain == 1, na.rm = TRUE),
      length(presence_in_domain), 100 * cfg$domain_min_presence_coverage),
      call. = FALSE)
  }
  prep_inputs$dropped_layers <- dropped
  saveRDS(list(fingerprint = prep_fp, inputs = prep_inputs,
    artifacts = unname(tools::md5sum(artifacts))), manifest_file)
  invisible(out_dir)
}

# Exact coverage on the actual analysis mask; no national-coverage proxy or imputation.
case_domain_qualification <- function(current, domain, minimum, manifest = NULL, cache_dir = NULL) {
  if (length(minimum) != 1L || !is.finite(minimum) || minimum <= 0 || minimum > 1)
    stop("predictor_min_coverage must be in (0, 1].")
  mask <- terra::values(domain, mat = FALSE)
  cells <- which(is.finite(mask) & mask != 0)
  if (!length(cells)) stop("The analysis mask has no active cells.")
  path <- NULL; keys <- NULL; previous <- numeric()
  counts <- rep(NA_real_, terra::nlyr(current))
  if (!is.null(cache_dir)) {
    if (is.null(manifest)) manifest <- case_file_manifest(current)
    key <- digest::digest(list(cells, terra::crs(domain), as.vector(terra::ext(domain)),
      dim(domain), terra::res(domain), as.character(utils::packageVersion("terra")),
      castSDM:::.cast_code_identity(case_domain_qualification)), algo = "sha256")
    path <- file.path(cache_dir, paste0(key, ".rds"))
    lock <- castSDM:::.cast_lock(paste0(path, ".lock"))
    on.exit(filelock::unlock(lock), add = TRUE)
    previous <- castSDM:::.cast_cache_read(path, key) %||% numeric()
    sources <- terra::sources(current, bands = TRUE)
    ix <- match(normalizePath(sources$source, winslash = "/"), manifest$path)
    if (nrow(sources) != length(counts) || anyNA(ix)) stop("Coverage cache needs all source layer identities.")
    keys <- paste(manifest$sha256[ix], sources$bands, sep = ":")
    counts <- unname(previous[keys])
  }
  nc <- terra::ncol(current)
  # Open each source once. Repeated 100,000-cell extraction reopened a national
  # raster up to 95 times per layer. Sequential bounded windows count exactly
  # the same native cells without those reopen/extraction costs.
  rows <- max(1L, as.integer(1000000L/nc))
  windows <- lapply(seq.int(1L, terra::nrow(current), by = rows), function(first) {
    count <- min(rows, terra::nrow(current)-first+1L)
    at <- cells[cells > (first-1L)*nc & cells <= (first+count-1L)*nc] - (first-1L)*nc
    if (!length(at)) return(NULL)
    col <- (at-1L) %% nc + 1L
    left <- min(col); width <- max(col)-left+1L
    list(first=first, count=count, left=left, width=width,
      index=((at-1L) %/% nc)*width + col-left+1L)
  })
  windows <- Filter(Negate(is.null), windows)
  missing <- which(is.na(counts))
  message("Exact analysis-mask coverage: ", length(counts)-length(missing),
    " cached; ", length(missing), " layers to scan.")
  for (i in missing) {
    layer <- current[[i]]
    terra::readStart(layer)
    counts[i] <- tryCatch(sum(vapply(windows, function(w) {
      z <- terra::readValues(layer, row=w$first, nrows=w$count,
        col=w$left, ncols=w$width, mat=FALSE)
      sum(is.finite(z[w$index]))
    }, numeric(1))), finally=terra::readStop(layer))
  }
  if (!is.null(path) && length(missing)) {
    previous[keys[missing]] <- counts[missing]
    castSDM:::.cast_cache_write(previous, path, key)
  }
  coverage <- counts/length(cells)
  data.frame(layer = names(current), n_domain = length(cells), n_finite = counts,
    coverage = coverage, minimum = minimum, eligible = coverage >= minimum,
    reason = ifelse(coverage >= minimum, "eligible", "insufficient_actual_domain_coverage"))
}
