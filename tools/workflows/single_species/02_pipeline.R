required_workflow_packages <- function() {
  packages <- c(
    "devtools", "pkgload", "parallelly", "data.table", "terra", "sf", "ranger",
    "xgboost", "maxnet", "mgcv", "ggplot2", "digest",
    "future.apply", "patchwork", "svglite", "ragg", "blockCV",
    "glmnet", "RRF"   # covsel Step B embedding (elastic net + guided RF)
  )
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
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
  if (cfg$input_type == "gbif") {
    required <- c("gbifID", "decimalLongitude", "decimalLatitude", "genus",
                  "countryCode", "occurrenceStatus", "hasGeospatialIssues",
                  "coordinateUncertaintyInMeters")
    if (!all(required %in% names(raw))) stop("Missing GBIF fields: ",
      paste(setdiff(required, names(raw)), collapse = ", "))
    for (v in c("decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters")) {
      raw[[v]] <- as.numeric(raw[[v]])
    }
    raw$hasGeospatialIssues <- tolower(raw$hasGeospatialIssues) %in% c("true", "t", "1")
    geo_fields <- intersect(c("georeferenceSources", "georeferenceRemarks", "locality"), names(raw))
    geo <- do.call(paste, c(as.list(raw[, geo_fields, drop = FALSE]), sep = " "))
    raw$qc_centroid <- grepl("centroid|country cent", geo, ignore.case = TRUE)
    raw$qc_unknown_uncertainty <- is.na(raw$coordinateUncertaintyInMeters)
    valid_xy <- is.finite(raw$decimalLongitude) & is.finite(raw$decimalLatitude) &
      abs(raw$decimalLongitude) <= 180 & abs(raw$decimalLatitude) < 90
    checks <- list(
      invalid_coordinates = !valid_xy,
      wrong_taxon = is.na(raw$genus) | raw$genus != cfg$species_name,
      outside_country = is.na(raw$countryCode) | raw$countryCode != "CN",
      not_present = is.na(raw$occurrenceStatus) | raw$occurrenceStatus != "PRESENT",
      geospatial_issue = raw$hasGeospatialIssues %in% TRUE,
      centroid_coordinates = raw$qc_centroid,
      excessive_uncertainty = !raw$qc_unknown_uncertainty &
        (raw$coordinateUncertaintyInMeters < 0 |
         raw$coordinateUncertaintyInMeters > cfg$occurrence_uncertainty_m),
      unknown_uncertainty = raw$qc_unknown_uncertainty & !cfg$occurrence_keep_unknown_uncertainty)
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
    fields <- intersect(c("source_file", "gbifID", "datasetKey", "species", "genus",
      "taxonRank", "basisOfRecord", "decimalLongitude", "decimalLatitude", "year", "eventDate",
      "coordinateUncertaintyInMeters", "georeferenceSources", "georeferenceRemarks",
      "locality", "countryCode", "occurrenceStatus", "hasGeospatialIssues",
      "x_albers", "y_albers", "grid_cell", "qc_centroid", "qc_unknown_uncertainty", "qc_reason"), names(raw))
    out <- data.frame(lon = raw$x_albers[kept], lat = raw$y_albers[kept])
    attr(out, "occurrence_audit") <- raw[, fields, drop = FALSE]
    return(out)
  } else if (cfg$input_type == "nsdm_occ") {
    # N-SDM species_occ tables are QC'd + 1 km-thinned and carry BOTH WGS84
    # (lon/lat) and Albers 1 km-grid (x_albers/y_albers) coordinates. The N-SDM
    # predictor stack IS that Albers grid and reprojection is banned on this
    # machine (terra::project / sf::st_transform segfault), so the Albers
    # columns are passed through under the lon/lat names castSDM expects.
    # Every downstream consumer (terra::cellFromXY / xyFromCell in
    # cast_background, as.data.frame(xy = TRUE) grids) treats lon/lat as plain
    # raster-CRS coordinates, keeping presences, background and prediction
    # grids mutually consistent with zero reprojection.
    if (!all(c("x_albers", "y_albers") %in% names(raw))) {
      stop("x_albers/y_albers columns not found: ", cfg$occurrence_file)
    }
    out <- data.frame(lon = as.numeric(raw$x_albers), lat = as.numeric(raw$y_albers))
  } else {
    lon <- grep("^lon|^longitude$", names(raw), ignore.case = TRUE, value = TRUE)[1]
    lat <- grep("^lat", names(raw), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(lon) || is.na(lat)) stop("Longitude/latitude columns not found: ", cfg$occurrence_file)
    out <- data.frame(lon = as.numeric(raw[[lon]]), lat = as.numeric(raw[[lat]]))
  }
  unique(out[is.finite(out$lon) & is.finite(out$lat), , drop = FALSE])
}

claim_species_output <- function(cfg) {
  if (length(cfg$output_root) != 1L || is.na(cfg$output_root) || !nzchar(cfg$output_root))
    stop("Choose an isolated output_root before preparing a case.")
  out_dir <- file.path(cfg$output_root, cfg$species_name)
  claim_file <- file.path(out_dir, "run_source.rds")
  claim <- list(config = cfg, source = workflow_source_signature(cfg))
  if (file.exists(file.path(out_dir, "completion.rds")))
    stop("Completed case is immutable; choose a new output_root: ", out_dir)
  if (file.exists(claim_file)) {
    if (!identical(readRDS(claim_file), claim))
      stop("Existing case has different source or settings; choose a new output_root: ", out_dir)
    return(invisible(out_dir))
  }
  if (file.exists(out_dir) && (!dir.exists(out_dir) ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE))))
    stop("Unclaimed output is not empty; refusing to overwrite: ", out_dir)
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE))
    stop("Cannot create case output directory: ", out_dir)
  saveRDS(claim, claim_file)
  invisible(out_dir)
}

prepare_species <- function(cfg) {
  out_dir <- claim_species_output(cfg)
  data_file <- file.path(out_dir, paste0(cfg$species_name, "_model_data.csv"))
  study_file <- file.path(out_dir, "study_area.rds")
  manifest_file <- file.path(out_dir,"preparation_manifest.rds")
  current <- wf_load_current_stack(cfg, verbose = TRUE)
  prep_inputs <- list(occurrence=tools::md5sum(cfg$occurrence_file),
    pool=tools::md5sum(cfg$nsdm_pool_csv), rasters=raster_input_manifest(current),
    config=cfg[c("input_type", "occurrence_uncertainty_m", "occurrence_keep_unknown_uncertainty",
                 "study_area_method","study_area_buffer_km","bg_ratio","bg_min","bg_max",
                 # The background-sampling design is part of the prepared data.
                 # Leaving these out let bg_strategy change from "random" to
                 # "environmental" while prep_fp stayed identical, so the run
                 # resumed and reused a background sample drawn by the OLD
                 # strategy -- silently mislabelled prepared data.
                 "bg_strategy","bg_bin_k","bg_sre_quantile","seed")],
    source=workflow_source_signature(cfg)$signature)
  prep_fp <- digest::digest(prep_inputs)
  old_fp <- if (file.exists(manifest_file)) readRDS(manifest_file)$fingerprint else NULL
  if (!is.null(old_fp) && !identical(old_fp, prep_fp))
    stop("Preparation inputs changed; choose a new output_root instead of overwriting this run.")
  n_layers <- as.integer(cfg$n_stack_layers)
  if (identical(old_fp,prep_fp) && file.exists(data_file) && file.exists(study_file)) {
    header <- names(data.table::fread(data_file, nrows = 1L))
    got <- length(setdiff(header, c("lon", "lat", "presence")))
    # prepare() drops zero-variance / all-NA pool layers, so the stored table
    # legitimately has fewer columns than the pool. Requiring an exact match
    # made the resume branch unreachable for every production configuration.
    dropped <- 0L
    man <- file.path(out_dir, "preparation_manifest.rds")
    if (file.exists(man)) {
      dropped <- readRDS(man)$inputs$dropped_layers
      dropped <- if (is.null(dropped)) 0L else length(dropped)
    }
    if (identical(got, n_layers - dropped)) {
      message("[resume] prepared data: ", cfg$species_name, " (", got,
              " predictors, ", dropped, " pool layers dropped)")
      return(invisible(out_dir))
    }
  }
  occ <- read_occurrences(cfg)
  occurrence_audit <- attr(occ, "occurrence_audit")
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
  # (Aiello-Lammens et al. 2015 spThin; N-SDM nsdm.disaggregate).
  thin_km <- cfg$occurrence_thin_km %||% 0
  if (thin_km > 0 && nrow(occ) > 1L) {
    mean_nn <- function(xy) mean(vapply(seq_len(nrow(xy)), function(i) {
      d <- sqrt((xy$lon[-i] - xy$lon[i])^2 + (xy$lat[-i] - xy$lat[i])^2)
      min(d)
    }, numeric(1)), na.rm = TRUE)
    nn_before <- mean_nn(occ)
    thin_year <- if (!is.null(occurrence_audit) &&
                     identical(cfg$occurrence_thin_order, "year_desc")) {
      yr <- occurrence_audit$year[occurrence_audit$qc_reason == "kept"]
      if (length(yr) == nrow(occ)) yr else NULL
    } else NULL
    keep_i <- which(.cast_thin_points(occ$lon, occ$lat, thin_km,
                                order = cfg$occurrence_thin_order %||% "random",
                                year = thin_year, seed = cfg$seed))
    n_before <- nrow(occ)
    occ <- occ[keep_i, , drop = FALSE]
    if (!is.null(occurrence_audit)) {
      eligible <- which(occurrence_audit$qc_reason == "kept")
      occurrence_audit$qc_reason[eligible[setdiff(seq_len(n_before), keep_i)]] <- "spatially_thinned"
    }
    nn_after <- if (nrow(occ) > 1L) mean_nn(occ) else NA_real_
    message(sprintf(
      "Spatial thinning (%.1f km, %s order): %d -> %d records; mean nearest-neighbour distance %.0f -> %.0f m.",
      thin_km, cfg$occurrence_thin_order %||% "random", n_before, nrow(occ),
      nn_before, nn_after))
    if (nrow(occ) < 10L) stop("Fewer than 10 occurrence cells survive thinning; lower occurrence_thin_km.")
  }
  if (terra::nlyr(current) != n_layers) {
    stop("Expected ", n_layers, " raster layers; found ", terra::nlyr(current))
  }
  # Reporting domain stays the accessible area; the BACKGROUND domain is a
  # separate decision (cfg$background_domain). With "national" the model is
  # shown the full climate contrast, so its variable choice is made on data
  # that can actually discriminate the niche; the selector's verdict either
  # way is then a finding, not an artefact of a truncated contrast.
  study <- wf_study_area(
    cfg, occ, current[[1]], verbose = TRUE
  )
  if (identical(cfg$background_domain, "national")) {
    study_bg <- nsdm_study_area(occ, current[[1]], method = "full", verbose = TRUE)
    message("Background domain: national land mask (cfg$background_domain = 'national')")
  } else {
    study_bg <- study
  }
  # --- Target-group background (Phillips et al. 2009) ---------------------------
  # Uniform-random background assumes uniform survey effort, which literature
  # records violate (roads, towns, accessible valleys). Sampling the
  # background from a BROADER taxon group's pooled occurrences (here: all
  # GBIF records in China) makes the background inherit the same spatial
  # bias as the presence data. The 555 MB source is reduced once -- snap to
  # the 1 km grid, keep cells inside the national mask, drop duplicate cells
  # -- and cached; cast_background() then draws n_bg from it with its own
  # filters (presence-cell exclusion, NA environments).
  tg_table <- NULL
  if (isTRUE(cfg$bg_target_group) && length(cfg$bg_target_group_file) &&
      file.exists(cfg$bg_target_group_file)) {
    tg_fp <- digest::digest(list(file = cfg$bg_target_group_file,
                                 size = file.info(cfg$bg_target_group_file)$size,
                                 mtime = file.info(cfg$bg_target_group_file)$mtime,
                                 seed = cfg$seed))
    tg_table <- cache_step(file.path(cache_dir, "target_group_cells.rds"), tg_fp, {
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
    })
  }
  model_data <- cast_background(
    occ, study_bg, current, ratio = cfg$bg_ratio,
    min_bg = cfg$bg_min, max_bg = cfg$bg_max, strategy = cfg$bg_strategy,
    bin_k = cfg$bg_bin_k, sre_quantile = cfg$bg_sre_quantile,
    user_table = tg_table,
    n_bg = if (!is.null(tg_table)) cfg$bg_max else NULL,
    cell_thin = TRUE, exclude_presence = TRUE, seed = cfg$seed, verbose = TRUE
  )
  vars <- get_env_vars(model_data)
  dropped <- setdiff(setdiff(names(model_data), c("lon", "lat", "presence")), vars)
  if (length(dropped)) {
    warning(sprintf(
      "%d/%d pool layers carry no training-sample signal (zero variance or all NA) and were dropped; first few: %s",
      length(dropped), n_layers, paste(utils::head(dropped, 15), collapse = ", ")),
      call. = FALSE)
    data.table::fwrite(data.frame(layer = dropped),
                       file.path(out_dir, "zero_variance_layers.csv"))
  }
  if (!is.null(occurrence_audit)) {
    eligible <- which(occurrence_audit$qc_reason == "kept")
    env <- as.matrix(terra::extract(current, as.matrix(occ[, c("lon", "lat")])))
    missing <- !is.finite(env)
    occurrence_audit$missing_predictors <- ""
    occurrence_audit$missing_predictors[eligible] <- apply(missing, 1L, function(z) {
      paste(names(current)[z], collapse = ";")
    })
    pres_cells <- terra::cellFromXY(current,
      as.matrix(model_data[model_data$presence == 1L, c("lon", "lat")]))
    occurrence_audit$used_in_model <- FALSE
    occurrence_audit$used_in_model[eligible] <- occurrence_audit$grid_cell[eligible] %in% pres_cells
    occurrence_audit$qc_reason[eligible[rowSums(missing) > 0L]] <- "missing_environment"
    data.table::fwrite(occurrence_audit, file.path(out_dir, "occurrence_qc.csv"))
    data.table::fwrite(as.data.frame(table(occurrence_audit$qc_reason)),
                      file.path(out_dir, "occurrence_qc_counts.csv"))
    stopifnot(sum(occurrence_audit$used_in_model) == sum(model_data$presence == 1L))
  }
  if (sum(model_data$presence == 1L) < 10L) {
    stop("Fewer than 10 complete occurrence cells remain after environmental extraction.")
  }
  data.table::fwrite(model_data, data_file)
  saveRDS(study, study_file)
  prep_inputs$dropped_layers <- dropped
  prep_inputs$background_domain <- cfg$background_domain
  saveRDS(list(fingerprint=prep_fp,inputs=prep_inputs),manifest_file)
  invisible(out_dir)
}

cache_step <- function(path, fingerprint, expr) {
  if (file.exists(path)) {
    old <- readRDS(path)
    if (identical(old$fingerprint, fingerprint)) {
      message("[resume] ", basename(path))
      return(old$value)
    }
  }
  value <- force(expr)
  tmp <- paste0(path, ".tmp")
  saveRDS(list(fingerprint = fingerprint, value = value), tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Cannot commit cache: ", path)
  value
}

# Preflight: prove the loaded castSDM actually provides every entry point and
# argument this workflow calls, BEFORE any expensive stage runs.
#
# Without this, a castSDM upgrade that renames or drops a single argument
# surfaces as an R "unused argument" error deep inside a cached step, or --
# worse -- as a silently ignored call whose product no longer matches the one
# the figures were built from. The workflow hardcodes package internals
# (castSDM:::theme_cast, .cast_map_bg) and passes a long tail of named
# arguments, so it needs an explicit compatibility contract with the package.
#
# Returns a record that goes into the run provenance, so every delivered
# bundle states the exact API surface it was produced against.
assert_cast_api <- function(cfg) {
  ns <- asNamespace("castSDM")
  exported <- getNamespaceExports("castSDM")

  # Every castSDM entry point this workflow depends on, with the arguments it
  # passes by name. Add a row here when a call is added; the preflight then
  # protects it automatically.
  required <- list(
    cast_prepare          = c("train_fraction", "seed", "env_vars",
                              "block_method", "verbose"),
    cast_select           = c("method", "corcut", "algorithms", "ncov",
                              "maxncov", "keep", "ensure", "env_vars",
                              "var_id", "category", "seed", "verbose"),
    cast_fit              = c("screen", "models", "seed", "verbose"),
    cast_cv               = c("select_method", "select_args", "k", "models",
                              "block_method", "block_size", "buffer",
                              "buffer_unit", "n_repeat",
                              "tune", "prevalence_target", "threshold_method",
                              "aoa", "parallel", "seed", "verbose"),
    cast_evaluate         = c("test_data", "threshold_method", "thresholds"),
    cast_predict          = c("new_data", "models", "cv", "aoa", "clamp"),
    cast_ensemble         = c("cv", "new_data", "method", "min_score", "decay",
                              "fallback", "min_metric_folds"),
    cast_ensemble_importance = c("cv", "new_data", "n_perm", "metric", "seed"),
    cast_project_raster   = c("cv", "current_raster", "future_rasters",
                              "output_dir", "mask", "method", "min_score",
                              "fallback", "min_metric_folds", "overwrite",
                              "clamp"),
    cast_effect_table     = c("newdata", "drivers", "shift", "shift_type",
                              "support_probs", "min_support"),
    cast_effect_map       = c("current_stack", "drivers", "shift", "shift_type",
                              "filename", "overwrite", "support_probs", "min_support"),
    cast_response_curves  = c("variables", "type"),
    cast_background       = c("study_area", "raster_stack", "ratio", "min_bg",
                              "max_bg", "strategy", "bin_k", "cell_thin",
                              "exclude_presence", "seed"),
    cast_importance       = "object",
    cast_report_odmap     = c("path", "meta", "verbose"),
    get_env_vars          = "data",
    new_cast_result       = c("screen", "fit", "eval", "cv", "fit_full")
  )

  problems <- character(0)
  missing_fn <- character(0)
  for (fn in names(required)) {
    if (!exists(fn, envir = ns, inherits = FALSE)) {
      missing_fn <- c(missing_fn, fn)
      next
    }
    if (fn %in% exported) next  # formals below cover exported functions
    # internal helpers are allowed but must still exist
  }
  if (length(missing_fn)) {
    problems <- c(problems, sprintf("missing castSDM function(s): %s",
                                    paste(missing_fn, collapse = ", ")))
  }
  for (fn in names(required)) {
    if (!exists(fn, envir = ns, inherits = FALSE)) next
    fml <- names(formals(get(fn, envir = ns)))
    need <- required[[fn]]
    absent <- setdiff(need, fml)
    if (length(absent)) {
      problems <- c(problems, sprintf("%s() lacks argument(s): %s", fn,
                                      paste(absent, collapse = ", ")))
    }
  }
  # Internals the figures reach into. These are unversioned by definition, so
  # they are checked explicitly rather than assumed.
  for (internal in c("theme_cast", ".cast_map_bg")) {
    if (!exists(internal, envir = ns, inherits = FALSE)) {
      problems <- c(problems, sprintf("missing castSDM internal: %s", internal))
    }
  }

  # Version floor. A newer package may still be fine; an OLDER one cannot
  # support the calls above, and the formals check cannot see a whole feature
  # that never existed (aoa, decay, permutation importance).
  have <- utils::packageVersion("castSDM")
  # packageVersion() takes a PACKAGE NAME; the floor is a version string, so it
  # has to be parsed with package_version() instead.
  want <- tryCatch(base::package_version(as.character(cfg$min_cast_version)),
                   error = function(e) NULL)
  if (is.null(want) || length(want) != 1L) {
    problems <- c(problems, sprintf(
      "cfg$min_cast_version (%s) is not a single version string",
      format(cfg$min_cast_version)))
    # Safe sentinel that always parses: the run must fail on the problem just
    # recorded, not on an error raised while describing it.
    want <- base::package_version("0.0.0")
  }
  version_ok <- have >= want
  if (!version_ok) {
    problems <- c(problems, sprintf(
      "castSDM %s is older than the required %s (cfg$min_cast_version)",
      as.character(have), as.character(want)))
  }

  # Enum-valued settings. castSDM validates these itself, but only once the
  # expensive stack has been extracted; checking them in the preflight turns a
  # late abort into an immediate, named one.
  check_enum <- function(value, arg, allowed) {
    if (length(value) != 1L || !is.character(value) || !nzchar(value)) {
      problems <<- c(problems, sprintf("cfg$%s must be a single string", arg))
    } else if (!value %in% allowed) {
      problems <<- c(problems, sprintf("cfg$%s = %s is not one of: %s", arg,
                                       format(value),
                                       paste(allowed, collapse = ", ")))
    }
  }
  check_enum(cfg$cv_threshold_method, "cv_threshold_method",
             eval(formals(cast_cv)$threshold_method))
  check_enum(cfg$cv_block_method, "cv_block_method",
             eval(formals(cast_cv)$block_method))
  check_enum(cfg$cv_buffer_unit, "cv_buffer_unit",
             eval(formals(cast_cv)$buffer_unit))
  check_enum(cfg$bg_strategy, "bg_strategy",
             eval(formals(cast_background)$strategy))
  check_enum(cfg$ensemble_method, "ensemble_method",
             eval(formals(cast_ensemble)$method))
  check_enum(cfg$ensemble_fallback, "ensemble_fallback",
             eval(formals(cast_ensemble)$fallback))
  check_enum(cfg$select_method, "select_method",
             eval(formals(cast_select)$method))
  check_enum(cfg$ensemble_perm_metric, "ensemble_perm_metric",
             eval(formals(cast_ensemble_importance)$metric))

  # Unit coherence. This workflow passes Albers METRE coordinates, but
  # cast_cv() only understands "deg" and "km". A nonzero buffer would be read
  # as degrees or kilometres, so refuse rather than silently run a fold with
  # a 1000x-wrong exclusion band.
  if (!is.null(cfg$cv_buffer) && is.numeric(cfg$cv_buffer) &&
      length(cfg$cv_buffer) == 1L && isTRUE(cfg$cv_buffer != 0)) {
    problems <- c(problems, sprintf(paste0(
      "cfg$cv_buffer = %s is not usable: cast_cv(buffer_unit = ) accepts ",
      "only \"deg\" and \"km\", but this workflow's coordinates are Albers ",
      "METRES, so a nonzero buffer cannot be expressed correctly. Keep it 0 ",
      "or pre-buffer the training rows before calling cast_cv()."),
      format(cfg$cv_buffer)))
  }

  record <- list(
    cast_version = as.character(have),
    min_cast_version = as.character(want),
    version_ok = version_ok,
    checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    functions = sort(names(required)),
    problems = problems
  )
  if (length(problems)) {
    stop("castSDM API preflight failed; the loaded package cannot run this ",
         "workflow as written:\n  - ", paste(problems, collapse = "\n  - "),
         "\n\nResolved castSDM: ", as.character(have),
         " (required >= ", as.character(want), ")")
  }
  # One string: cli_inform() treats each element of its argument as a
  # SEPARATE message, so splitting this would print a truncated line.
  cli::cli_inform("castSDM API preflight OK (v{as.character(have)}; {length(required)} entry points checked).")
  record
}


# Reporting products that derive entirely from the saved result + audit, so they
# cost nothing to regenerate and never need a refit.
#
#   odmap_report.md            ODMAP protocol (Zurell et al. 2020). An SDM study
#                              is expected to ship one; cast_report_odmap()
#                              fills everything the pipeline knows and leaves
#                              "[report]" placeholders for the rest.
#   attribution_vs_importance  the two-stage screen's forward path
#                              (cast_importance(): admission step_added with
#                              loss_gain) joined to the single-shift
#                              hard-masked audit, so the case study shows
#                              screening and shift-response side by side.
write_species_reports <- function(cfg, result, audit) {
  out_dir <- file.path(cfg$output_root, cfg$species_name)
  table_dir <- file.path(out_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  cast_report_odmap(result, path = file.path(out_dir, "odmap_report.md"),
                    meta = list(
                      taxon = paste(cfg$species_name, cfg$taxon_level),
                      title = paste(cfg$species_name, "- castSDM", cfg$taxon_level, "case study"),
                      purpose = "Environmental drivers and spatial attribution of current/future habitat suitability change; paired scenarios and explicit training-range diagnostics.",
                      data_description = paste0(
                        "N-SDM 1 km Albers China stack, ", cfg$n_stack_layers,
                        "-layer pool; presence/background data, one record per raster cell.",
                        " See occurrence_qc.csv for raw-source exclusions, unknown location uncertainty",
                        " and records lost during environmental extraction."),
                      assumptions = paste("Presence/background scores are not absolute occurrence probabilities.",
                        "This opportunistic record collection has no identified sampling model or sufficient causal adjustment set.",
                        "Holding other predictors fixed does not establish causal identification.",
                        "Range-box masking is a range-extrapolation diagnostic, not an identification test.",
                        "National land background does not represent fish-accessible river habitat; reporting is restricted to the occurrence buffer.")),
                    verbose = TRUE)

  imp <- if (isTRUE(cfg$single_shift_diagnostics)) tryCatch(cast_importance(result), error = function(e) {
    cli::cli_warn("cast_importance() unavailable for {cfg$species_name}: {conditionMessage(e)}")
    NULL
  }) else NULL
  if (!is.null(imp) && !is.null(audit$audit)) {
    ie <- as.data.frame(imp$effects)
    # Forward path: admission step with inner-CV loss improvement; excluded
    # survivors carry missing steps and gains.
    cmp <- merge(audit$audit,
                 ie[, c("variable", "step_added", "loss_gain",
                        "selected")],
                 by.x = "driver", by.y = "variable", all.x = TRUE)
    names(cmp)[names(cmp) == "selected"] <- "kept_by_screen"
    cmp$rank_admission <- rank(cmp$step_added,
                               na.last = "keep", ties.method = "min")
    cmp$rank_effect <- rank(-cmp$mean_abs_dHSS, na.last = "keep",
                            ties.method = "min")
    cmp <- cmp[order(cmp$rank_effect), ]
    rownames(cmp) <- NULL
    data.table::fwrite(cmp, file.path(table_dir, "attribution_vs_importance.csv"))
  }
  invisible(out_dir)
}


# Lightweight cache probe: TRUE when `path` already holds a result whose stored
# fingerprint matches. Lets us skip spawning parallel workers on a clean resume.
cache_is_fresh <- function(path, fingerprint) {
  if (!file.exists(path)) return(FALSE)
  old <- tryCatch(readRDS(path), error = function(e) NULL)
  !is.null(old) && identical(old$fingerprint, fingerprint)
}

run_species <- function(cfg) {
  source_info <- workflow_source_signature(cfg)
  out_dir <- prepare_species(cfg)
  cache_dir <- file.path(out_dir, "cache_v060")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  data_file <- file.path(out_dir, paste0(cfg$species_name, "_model_data.csv"))
  d <- data.table::fread(data_file, data.table = FALSE)
  vars <- get_env_vars(d)
  if (length(vars) < as.integer(cfg$n_stack_layers)) {
    warning(sprintf(
      "%d varying predictors of %d pool layers (zero-variance layers dropped); continuing.",
      length(vars), cfg$n_stack_layers), call. = FALSE)
  }
  # --- covsel selection metadata (derived once, reused per CV fold) --------
  # var_id: map repeated measurements onto their construct. Year tokens are
  # stripped wherever they sit in the name (anthro_foothii_2001_hii and
  # anthro_foothii_2015_hii are the same construct at two dates), then the
  # focal-scale / min-max suffixes (end-anchored, from config). Step A keeps
  # one covariate per construct -- but only after temporal matching has
  # already chosen the construct's date slice.
  var_id <- vars
  var_id <- gsub("_(19|20)[0-9]{2}(?=_|$)", "", var_id, perl = TRUE)
  for (pat in (cfg$select_varid_strip %||% character(0))) {
    var_id <- sub(paste0(pat, "$"), "", var_id)
  }
  names(var_id) <- vars
  category <- rep(NA_character_, length(vars))
  names(category) <- vars
  climate_layers <- character(0)
  is_climate_flag <- rep(FALSE, length(vars))
  names(is_climate_flag) <- vars
  layer_year <- rep(NA_integer_, length(vars))
  names(layer_year) <- vars
  if (length(cfg$nsdm_pool_csv) && file.exists(cfg$nsdm_pool_csv)) {
    pool <- data.table::fread(cfg$nsdm_pool_csv, data.table = FALSE)
    pool <- pool[pool$status == "kept", , drop = FALSE]
    pool$layer <- gsub("-", "_", pool$layer)
    hit <- match(vars, pool$layer)
    # Unmatched layers (legacy columns absent from the manifest) get a generic
    # bucket rather than NA -- cast_select() rejects NA category labels.
    category <- ifelse(is.na(hit), "unassigned", pool$group[hit])
    names(category) <- vars
    climate_layers <- pool$layer[pool$is_climate %in% TRUE]
    # Climate flag from the pool manifest's own declaration (provenance, not a
    # variable name) -- feeds the generic coverage rule in cast_select().
    is_climate_flag[!is.na(hit)] <- pool$is_climate[hit[!is.na(hit)]] %in% TRUE
    # Layer year: the year directory in the library path (glo/<group>/<family>/
    # <year>/...), else a trailing _YYYY in the layer name. NA = static layer.
    py <- suppressWarnings(as.integer(sub("^.*/((19|20)[0-9]{2})/.*$", "\\1",
                                          pool$path)))
    py[is.na(py)] <- suppressWarnings(as.integer(
      sub("^.*_((19|20)[0-9]{2})(_.*)?$", "\\1", pool$layer[is.na(py)])))
    layer_year[!is.na(hit)] <- py[hit[!is.na(hit)]]
  }
  # Record-year reference: median of all occurrence years.
  ref_year <- NULL
  occ_files <- (cfg$occurrence_file %||% character(0))[file.exists(cfg$occurrence_file %||% character(0))]
  if (length(occ_files)) {
    yrs <- unlist(lapply(occ_files, function(f) {
      hdr <- names(data.table::fread(f, nrows = 0))
      ycol <- grep("^year$", hdr, ignore.case = TRUE, value = TRUE)
      if (!length(ycol)) return(NULL)
      data.table::fread(f, select = ycol[1])[[1]]
    }))
    yrs <- suppressWarnings(as.numeric(yrs))
    yrs <- yrs[is.finite(yrs) & yrs > 1500 & yrs < 2100]
    if (length(yrs)) ref_year <- stats::median(yrs)
  }
  # npts: within each construct keep only the nearest-previous slice (fallback:
  # the earliest slice when every slice postdates the reference).
  tmatch <- cfg$select_temporal_match %||% "npts"
  if (!is.null(ref_year) && tmatch %in% c("npts", "nts") &&
      any(!is.na(layer_year))) {
    keep_t <- rep(TRUE, length(vars))
    for (v in unique(var_id)) {
      ix <- which(var_id == v)
      yrs <- layer_year[ix]
      if (all(is.na(yrs))) next
      cand <- yrs[!is.na(yrs)]
      m <- if (identical(tmatch, "npts")) {
        if (any(cand <= ref_year)) max(cand[cand <= ref_year]) else min(cand)
      } else {
        cand[which.min(abs(cand - ref_year))]
      }
      keep_t[ix[!is.na(yrs) & yrs != m]] <- FALSE
    }
    n_tdrop <- sum(!keep_t)
    if (n_tdrop > 0L) {
      message(sprintf(
        "Temporal matching (%s, record-year reference %d): dropped %d of %d layers postdating or off the matched slice.",
        tmatch, as.integer(ref_year), n_tdrop, length(vars)))
      vars <- vars[keep_t]
      var_id <- var_id[keep_t]
      category <- category[keep_t]
      layer_year <- layer_year[keep_t]
      is_climate_flag <- is_climate_flag[keep_t]
    }
  } else if (is.null(ref_year)) {
    message("No occurrence years available; temporal matching skipped (all year-tagged layers stay in the pool).")
  }
  # Fingerprint only the settings that drive selection/fitting/CV. Projection-only
  # options (future_models, ssp_scenarios, periods) run after the cached steps, so
  # they must never invalidate the cached modeling steps and are dropped here;
  # "basemap" only affects figure rendering and is dropped for the same reason.
  # Reporting-only knobs are dropped for the same reason: ensemble permutation
  # importance carries its OWN fingerprint, so leaving these in would let a
  # toggle that cannot change a single predicted value throw away a good fit.
  cfg_model <- cfg[setdiff(names(cfg), c("future_models", "ssp_scenarios", "periods",
                                        "basemap",
                                        "ensemble_permutation", "ensemble_perm_n",
                                        "ensemble_perm_metric",
                                        "min_cast_version"))]
  # Package version is read dynamically so a dev pkgload::load_all() of castSDM
  # (e.g. 0.17.0) fingerprints correctly and stale caches invalidate on upgrade.
  fp <- digest::digest(list(cfg_model, file.info(data_file)[, c("size", "mtime")], vars,
                            select_method = cfg$select_method,
                            select_temporal_match = tmatch,
                            record_year_ref = if (is.null(ref_year)) NA_integer_ else as.integer(ref_year),
                            data_hash = tools::md5sum(data_file),
                            source = workflow_source_signature(cfg)$signature))

  split <- cache_step(file.path(cache_dir, "prepare.rds"), fp,
                      cast_prepare(d, train_fraction = cfg$train_fraction,
                                   seed = cfg$seed, env_vars = vars,
                                   block_method = cfg$cv_block_method,
                                   verbose = FALSE))
  # Align the selection metadata with what the training split ACTUALLY
  # supports. cast_prepare drops layers that are constant in the training
  # sample (e.g. land-cover patch statistics over a uniform plateau window):
  # declared = upper bound, split$env_vars = what can be fitted. Aligning here
  # keeps var_id/category/is_climate_flag the same length as the predictors
  # cast_select() sees.
  fit_vars <- split$env_vars
  drop_vars <- setdiff(vars, fit_vars)
  if (length(drop_vars)) {
    message(sprintf("Training-sample variance screen dropped %d layer(s): %s%s",
                    length(drop_vars), paste(drop_vars, collapse = ", "),
                    if (length(drop_vars) > 5) ", ..." else ""))
    vars <- fit_vars
    var_id <- var_id[match(fit_vars, names(var_id))]
    category <- category[match(fit_vars, names(category))]
    is_climate_flag <- is_climate_flag[match(fit_vars, names(is_climate_flag))]
    climate_layers <- intersect(climate_layers, vars)
  }
  # Selection arguments per method. "covsel" is the two-step embedded
  # selection of Adde et al. (2023); the target count follows the package
  # default ceiling(log2(n_presences)), capped by maxncov. "full" keeps every
  # varying predictor (prespecified sets only).
  select_args <- switch(cfg$select_method,
    covsel = list(corcut = cfg$select_corcut,
                  algorithms = cfg$select_algorithms,
                  ncov = if (is.null(cfg$select_ncov)) NULL else as.integer(cfg$select_ncov),
                  maxncov = as.integer(cfg$select_maxncov),
                  keep = cfg$select_keep %||% character(0),
                  # Generic coverage rule: the climate flag comes from the pool
                  # manifest (provenance), never a hardcoded variable name, so
                  # the mechanism is identical for every species while the
                  # admitted covariate is whatever the data ranks strongest.
                  ensure = if (isTRUE(cfg$select_ensure_climate) &&
                               any(is_climate_flag)) {
                    list(climate = is_climate_flag)
                  } else {
                    NULL
                  },
                  # Explicit predictor contract (as in cast_prepare): a layer
                  # constant inside one CV fold's training rows is dropped
                  # there rather than causing a length mismatch. The named
                  # var_id/category vectors are re-indexed automatically.
                  env_vars = vars,
                  var_id = var_id, category = category),
    full = list())
  screen <- cache_step(file.path(cache_dir, "screen.rds"), fp,
    do.call(cast_select, c(list(data = split$train, method = cfg$select_method,
                                seed = cfg$seed, verbose = TRUE), select_args)))
  if (!length(screen$selected)) stop("No varying predictors retained; inspect the selection result before fitting.")
  # --- Projection-degeneracy guard ------------------------------------------
  # A model whose selected covariates are all static makes every future
  # scenario identical to the present -- verified on the 20261010 run: 0%
  # change under all 18 scenarios, pixel-identical maps. When future
  # scenarios are requested, require at least one selected covariate that
  # varies between current and future stacks (the pool declares is_climate).
  # The remedy is the generic coverage rule (select_ensure_climate) or a
  # researcher-declared set in select_keep (covsel force), not a
  # species-specific code change.
  climate_layers <- intersect(climate_layers, vars)
  if (length(cfg$future_models %||% character(0)) > 0 && length(climate_layers) &&
      !any(screen$selected %in% climate_layers)) {
    stop(sprintf(
      paste0("None of the %d selected covariates varies between current and future stacks;",
             " every future scenario would be identical to the present (0%% change by construction).",
             " Keep select_ensure_climate = TRUE (generic coverage rule) or declare the",
             " bioclimatic set in cfg$select_keep (covsel force), or clear cfg$future_models."),
      length(screen$selected)), call. = FALSE)
  }
  fit <- cache_step(file.path(cache_dir, "fit.rds"), fp,
    cast_fit(split$train, screen = screen, models = cfg$models,
             rf_ntree = cfg$rf_ntree, rf_mtry = cfg$rf_mtry,
             brt_n_trees = cfg$brt_n_trees, brt_depth = cfg$brt_depth,
             brt_shrinkage = cfg$brt_shrinkage,
             maxent_classes = cfg$maxent_classes, maxent_regmult = cfg$maxent_regmult,
             num_threads = cfg$model_num_threads, tune = cfg$cv_tune,
             tune_folds = cfg$tune_folds,
             prevalence_target = if (is.na(cfg$prevalence_target)) NULL else cfg$prevalence_target,
             seed = cfg$seed, verbose = TRUE))
  # Nested CV runs on the training split only, so split$test stays a pure
  # hold-out that never influences selection, fitting, or CV thresholds. It is
  # also the ONLY fold-parallel step (each fold is single-threaded via
  # num_threads=1L in cast_fit). Spin up a PSOCK cluster only when CV actually
  # needs recomputing, then tear it down immediately so no idle worker holds RAM
  # during the memory-heavy raster projection below.
  cv_path <- file.path(cache_dir, "cv_nested.rds")
  cv_workers <- if (is.null(cfg$cv_workers)) 1L else as.integer(cfg$cv_workers)
  n_workers <- max(1L, min(cv_workers, as.integer(cfg$cv_k),
                           as.integer(future::availableCores(omit = 1L))))
  use_parallel <- n_workers > 1L && !cache_is_fresh(cv_path, fp)
  if (use_parallel) {
    # A fresh worker would attach the stale *installed* castSDM and fail to find
    # cast_select(); load the identical dev source on each worker instead.
    cl <- parallelly::makeClusterPSOCK(n_workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterCall(cl, function(root) {
      suppressMessages(pkgload::load_all(root, quiet = TRUE))
      invisible(NULL)
    }, root = cfg$cast_root)
    future::plan(future::cluster, workers = cl)
  }
  cv <- cache_step(cv_path, fp,
    cast_cv(split$train, select_method = cfg$select_method,
            select_args = select_args,
            k = cfg$cv_k, models = cfg$models,
            block_method = cfg$cv_block_method, block_size = cfg$cv_block_size,
            buffer = cfg$cv_buffer,
            buffer_unit = cfg$cv_buffer_unit,
            n_repeat = cfg$cv_n_repeat, tune = cfg$cv_tune,
            threshold_method = cfg$cv_threshold_method,
            rf_ntree = cfg$rf_ntree, rf_mtry = cfg$rf_mtry,
            brt_n_trees = cfg$brt_n_trees, brt_depth = cfg$brt_depth,
            brt_shrinkage = cfg$brt_shrinkage,
            maxent_classes = cfg$maxent_classes, maxent_regmult = cfg$maxent_regmult,
            num_threads = cfg$model_num_threads, tune_folds = cfg$tune_folds,
            aoa = cfg$cv_aoa,
            prevalence_target = if (is.na(cfg$prevalence_target)) NULL else cfg$prevalence_target,
            parallel = use_parallel, seed = cfg$seed, verbose = TRUE))
  if (use_parallel) {
    future::plan(future::sequential)
    parallel::stopCluster(cl)
  }
  # One-shot hold-out evaluation. Thresholds are the pooled out-of-fold
  # operating points from cast_cv() -- selected without touching the hold-out
  # rows -- and are applied unchanged here (protocol rule 2 in castSDM 0.20.0).
  evaluation <- cache_step(file.path(cache_dir, "evaluation.rds"), fp,
    cast_evaluate(fit, split$test,
                  threshold_method = cfg$cv_threshold_method,
                  thresholds = cv$thresholds))
  table_dir <- file.path(out_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(screen$scores, file.path(table_dir, "selection_scores.csv"))
  data.table::fwrite(split$train, file.path(table_dir, "training_data.csv"))
  data.table::fwrite(split$test, file.path(table_dir, "holdout_data.csv"))
  data.table::fwrite(evaluation$metrics, file.path(table_dir, "holdout_metrics.csv"))
  data.table::fwrite(cv$metrics, file.path(table_dir, "cv_metrics.csv"))
  data.table::fwrite(cv$fold_metrics, file.path(table_dir, "cv_fold_metrics.csv"))
  data.table::fwrite(cv$selection_freq, file.path(table_dir, "cv_selection_frequency.csv"))
  data.table::fwrite(data.frame(row = seq_len(nrow(split$train)), fold = cv$folds,
    cv$oof), file.path(table_dir, "cv_oof_predictions.csv"))
  data.table::fwrite(data.frame(fold = seq_len(cv$k), status = cv$fold_status,
    selected = vapply(cv$selections, paste, character(1), collapse = ";")),
    file.path(table_dir, "cv_fold_status.csv"))
  engine_status <- data.frame(model = cfg$models,
    fitted = cfg$models %in% names(fit$models))
  data.table::fwrite(engine_status, file.path(table_dir, "engine_status.csv"))
  # Fold-coverage audit (castSDM >= 0.18.0). cast_cv() records how many folds
  # actually produced predictions per engine; an engine below the coverage
  # floor yields its threshold from a partial out-of-fold surface. Recording it
  # here means a reader can see the coverage instead of trusting a silent
  # fallback to 0.5.
  if (!is.null(cv$pred_folds)) {
    n_folds_total <- length(cv$fold_status)
    coverage <- data.frame(
      model = names(cv$pred_folds),
      folds_with_predictions = as.integer(cv$pred_folds),
      folds_total = n_folds_total,
      coverage = round(as.integer(cv$pred_folds) / max(1L, n_folds_total), 3),
      stringsAsFactors = FALSE)
    coverage$meets_floor <- coverage$folds_with_predictions >=
      max(2L, ceiling(n_folds_total * 0.8))
    data.table::fwrite(coverage, file.path(table_dir, "cv_fold_coverage.csv"))
    if (!any(coverage$meets_floor)) {
      warning(sprintf(
        paste0("No engine met the fold-coverage floor (%d of %d folds); the ",
               "ensemble threshold rests on a partial out-of-fold surface."),
        max(2L, ceiling(n_folds_total * 0.8)), n_folds_total), call. = FALSE)
    }
  }
  # Area of Applicability (castSDM >= 0.18.0). The DI distribution over the
  # held-out rows is the extrapolation evidence for the whole projection: a
  # wide tail means many reported cells sit outside the calibrated region.
  if (isTRUE(cfg$cv_aoa) && !is.null(cv$aoa)) {
    di <- cv$aoa$di_oof %||% cv$di
    if (!is.null(di) && length(di)) {
      probs <- c(0.5, 0.75, 0.9, 0.95, 0.99)
      thr <- cv$aoa$threshold
      aoa_tab <- data.frame(
        statistic = c("n_oof", "di_p50", "di_p75", "di_p90", "di_p95",
                      "di_p99", "di_max", "aoa_threshold",
                      "share_outside_aoa_oof"),
        value = c(length(di), stats::quantile(di, probs, na.rm = TRUE),
                  max(di, na.rm = TRUE), thr,
                  if (is.finite(thr)) mean(di > thr, na.rm = TRUE) else NA_real_),
        stringsAsFactors = FALSE)
      aoa_tab <- aoa_tab[stats::complete.cases(aoa_tab$value) |
                           aoa_tab$statistic %in% c("n_oof", "aoa_threshold"), ]
      data.table::fwrite(aoa_tab, file.path(table_dir, "aoa_summary.csv"))
    }
  }
  data.table::fwrite(data.frame(method = "tramicp", status = "not_applicable",
    reason = "No scientifically prespecified environment groups; spatial folds are not substitutes."),
    file.path(table_dir, "invariant_selection_applicability.csv"))
  if (!all(engine_status$fitted)) {
    stop("Requested engines failed; inspect engine_status.csv before accepting a complete case.")
  }
  holdout_prediction <- cast_predict(fit, split$test, cv = cv,
                                     aoa = cfg$cv_aoa, clamp = cfg$clamp)
  data.table::fwrite(data.frame(presence = split$test$presence,
    holdout_prediction$predictions), file.path(table_dir, "holdout_predictions.csv"))
  # The configured method, not a literal: a hardcoded "weighted" here made
  # holdout_ensemble.csv disagree with rasters/ whenever the configured
  # ensemble method differed, with nothing to flag it.
  holdout_ensemble <- cast_ensemble(fit, cv, split$test,
                                    method = cfg$ensemble_method,
                                    min_score = cfg$ensemble_min_score,
                                    decay = cfg$ensemble_decay,
                                    fallback = cfg$ensemble_fallback,
                                    min_metric_folds = cfg$ensemble_min_metric_folds)
  data.table::fwrite(holdout_ensemble$predictions,
                     file.path(table_dir, "holdout_ensemble.csv"))
  # Ensemble-layer permutation importance (castSDM >= 0.18.0). Distinct from
  # cast_importance() above, which replays the SELECTION path: this permutes
  # each predictor and recomputes the whole ensemble, so it ranks drivers by
  # their effect on the published surface. Reporting-only, so it gets its own
  # fingerprint and never invalidates the cached fit when toggled.
  ensemble_imp <- NULL
  if (isTRUE(cfg$ensemble_permutation)) {
    perm_fp <- digest::digest(list(model = fp, step = "ensemble-permutation",
                                   n_perm = cfg$ensemble_perm_n,
                                   metric = cfg$ensemble_perm_metric,
                                   method = cfg$ensemble_method,
                                   decay = cfg$ensemble_decay,
                                   min_score = cfg$ensemble_min_score,
                                   fallback = cfg$ensemble_fallback))
    ensemble_imp <- cache_step(file.path(cache_dir, "ensemble_importance.rds"),
                               perm_fp, {
      ei <- tryCatch(
        cast_ensemble_importance(fit, cv, split$test,
                                 models = cfg$models,
                                 method = cfg$ensemble_method,
                                 min_score = cfg$ensemble_min_score,
                                 decay = cfg$ensemble_decay,
                                 fallback = cfg$ensemble_fallback,
                                 min_metric_folds = cfg$ensemble_min_metric_folds,
                                 n_perm = cfg$ensemble_perm_n,
                                 metric = cfg$ensemble_perm_metric,
                                 seed = cfg$seed, verbose = TRUE),
        error = function(e) {
          cli::cli_warn("cast_ensemble_importance() failed: {conditionMessage(e)}")
          NULL
        })
      ei
    })
    if (!is.null(ensemble_imp)) {
      data.table::fwrite(ensemble_imp$importance,
                         file.path(table_dir, "ensemble_permutation_importance.csv"))
    }
  }

  # Final models: refit on the FULL data set after validation (ODMAP/Zurell
  # et al. 2020 standard practice). Selection stays decided on the training
  # split; the hold-out evaluation above stays on the training fit; only the
  # final projection and effect products below use every record. With ~25
  # presences, discarding the 30% hold-out from the final model would waste
  # scarce records.
  fit_full <- cache_step(file.path(cache_dir, "fit_full.rds"), fp,
    cast_fit(d, screen = screen, models = cfg$models,
             rf_ntree = cfg$rf_ntree, rf_mtry = .cast_tuned_args(fit)$rf_mtry,
             brt_n_trees = cfg$brt_n_trees,
             brt_depth = .cast_tuned_args(fit)$brt_depth %||% cfg$brt_depth,
             brt_shrinkage = .cast_tuned_args(fit)$brt_shrinkage %||% cfg$brt_shrinkage,
             maxent_classes = .cast_tuned_args(fit)$maxent_classes %||% cfg$maxent_classes,
             maxent_regmult = .cast_tuned_args(fit)$maxent_regmult %||% cfg$maxent_regmult,
             num_threads = cfg$model_num_threads,
             prevalence_target = if (is.na(cfg$prevalence_target)) NULL else cfg$prevalence_target,
             seed = cfg$seed, verbose = TRUE))
  engine_status_full <- data.frame(model = cfg$models,
    fitted = cfg$models %in% names(fit_full$models))
  if (!all(engine_status_full$fitted)) {
    stop("Full-data refit failed for an engine; inspect the refit before accepting the case.")
  }

  current <- wf_load_current_stack(cfg, verbose = TRUE)
  future <- wf_build_future_stacks(
    cfg, current, selected_vars = fit_full$env_vars, verbose = TRUE
  )
  # The fingerprint MUST carry the same method the projection call uses; a
  # literal here left fp and projection_fp unchanged when the configured
  # ensemble method changed, so the run silently reused a stale projection.
  # method must be the CONFIGURED one and must equal what the call below
  # passes. As a literal it left fp and projection_fp unchanged when
  # cfg$ensemble_method changed, so the run silently reused a projection that a
  # different ensemble method had produced.
  projection_inputs <- list(model = fp, current = raster_input_manifest(current),
                            future = raster_input_manifest(future),
                            method = cfg$ensemble_method)
  projection_fp <- digest::digest(projection_inputs)
  projection_manifest <- file.path(out_dir,"projection_manifest.rds")
  old_projection <- if (file.exists(projection_manifest)) readRDS(projection_manifest)$fingerprint else NULL
  projection_overwrite <- !identical(old_projection,projection_fp)
  # The reporting domain is the accessible area, not the whole national grid:
  # a national map extrapolates a model that was fitted inside a buffer.
  # The mask is rebuilt from the boundary polygon instead of read from
  # study_area.rds$mask, because a SpatRaster restored from RDS carries a
  # dangling pointer in a new session and faults on first use.
  sa <- readRDS(file.path(out_dir, "study_area.rds"))
  study_v <- terra::vect(sa$boundary[[1L]][[1L]], type = "polygons",
                         crs = attr(sa$boundary, "crs")$wkt)
  area_mask <- terra::rasterize(study_v, current[[1L]], field = 1, background = NA)
  # Restrict raster processing to the reporting extent before model prediction.
  # This preserves the original 1 km cells and polygon mask while avoiding
  # repeated nationwide evaluations of cells excluded from every output.
  region <- terra::ext(study_v)
  input_dir <- file.path(cache_dir, "scenario_inputs")
  dir.create(input_dir, showWarnings = FALSE)
  crop_to_file <- function(r, stem) terra::crop(r, region, snap = "out",
    filename = file.path(input_dir, paste0(stem, ".tif")), overwrite = TRUE,
    wopt = list(datatype = "FLT8S", gdal = "COMPRESS=LZW"))
  current <- crop_to_file(current[[fit_full$env_vars]], "current")
  future <- lapply(names(future), function(s) crop_to_file(future[[s]], s)) |>
    stats::setNames(names(future))
  area_mask <- crop_to_file(area_mask, "mask")
  projection_workers <- max(1L, as.integer(cfg$projection_workers %||% 1L))
  if (projection_workers > 1L) {
    cl_project <- parallelly::makeClusterPSOCK(projection_workers)
    on.exit(try(parallel::stopCluster(cl_project), silent = TRUE), add = TRUE)
    parallel::clusterCall(cl_project, function(root) {
      options(ranger.num.threads = 1L)
      suppressMessages(pkgload::load_all(root, quiet = TRUE))
      invisible(NULL)
    }, root = cfg$cast_root)
    future::plan(future::cluster, workers = cl_project)
  }
  projection <- cast_project_raster(
    fit_full, cv, current_raster = current[[fit_full$env_vars]],
    future_rasters = future, output_dir = out_dir, mask = area_mask,
    method = cfg$ensemble_method,
    min_score = cfg$ensemble_min_score,
    decay = cfg$ensemble_decay,
    # Same ensemble contract as the hold-out ensemble above: fallback travels
    # with the projection, otherwise the published rasters and the holdout
    # table can disagree on how a missing engine is handled (decay is not a
    # projection-side knob; assert_cast_api guards the combination).
    fallback = cfg$ensemble_fallback,
    min_metric_folds = cfg$ensemble_min_metric_folds,
    clamp = cfg$clamp,
    overwrite = projection_overwrite, verbose = TRUE,
    parallel = projection_workers > 1L, max_memory_mb = 25
  )
  if (length(projection$future) != length(future) ||
      nrow(projection$stats) != length(future) ||
      ("error" %in% names(projection$stats) && any(nzchar(projection$stats$error) & !is.na(projection$stats$error)))) {
    stop("Projection is incomplete; inspect failed scenarios before accepting the case bundle.")
  }
  saveRDS(list(fingerprint=projection_fp, inputs=projection_inputs),projection_manifest)
  packed_current <- terra::wrap(current[[fit_full$env_vars]], proxy = TRUE)
  packed_future <- lapply(future, terra::wrap, proxy = TRUE)
  packed_mask <- terra::wrap(area_mask, proxy = TRUE)
  attribution <- future.apply::future_lapply(names(future), function(scen) {
    path <- file.path(out_dir, "rasters", paste0(scen, "_attribution.tif"))
    obj <- cast_attribution(fit_full, terra::unwrap(packed_current),
      terra::unwrap(packed_future[[scen]]),
      weights = projection$current$weights, drivers = fit_full$env_vars,
      mask = terra::unwrap(packed_mask), filename = path, overwrite = TRUE,
      block_size = cfg$attribution_block_size, support_probs = cfg$effect_support_probs,
      verbose = FALSE)
    obj$summary$scenario <- scen
    obj
  }, future.seed = TRUE)
  if (projection_workers > 1L) {
    future::plan(future::sequential)
    parallel::stopCluster(cl_project)
  }
  names(attribution) <- names(future)
  data.table::fwrite(data.table::rbindlist(lapply(attribution, `[[`, "summary")),
    file.path(table_dir, "scenario_attribution.csv"))
  # Single-shift frozen-weight ensemble response with hard range masking,
  # on the full-data fit and its own training rows (all records). Masked
  # drivers (support < min_support) carry NA estimates.
  audit <- list(); effect_map <- NULL
  if (isTRUE(cfg$single_shift_diagnostics)) {
  audit <- cache_step(file.path(cache_dir, "attribution_audit.rds"), fp, {
    # Out of sample by construction: the audit model is fitted on the 70%
    # training split and the do-shift contrast is measured on the 30%
    # holdout rows. Reading the same contrast off fit_full's own training
    # rows inflates the driver ranking, and that ranking is what selects
    # the mapped driver.
    et <- cast_effect_table(fit, newdata = split$test, drivers = fit$env_vars,
                            shift = cfg$effect_shift,
                            shift_type = cfg$effect_shift_type,
                            support_probs = cfg$effect_support_probs,
                            min_support = cfg$effect_min_support,
                            weights = projection$current$weights, verbose = TRUE)
    et_in <- cast_effect_table(fit_full, newdata = d, drivers = fit_full$env_vars,
                               shift = cfg$effect_shift,
                               shift_type = cfg$effect_shift_type,
                               support_probs = cfg$effect_support_probs,
                               min_support = cfg$effect_min_support,
                               weights = projection$current$weights, verbose = TRUE)
    tab <- as.data.frame(et)[, c("driver", "shift_raw", "mean_abs_dHSS", "mean_signed_dHSS",
                                 "n", "n_supported", "support", "masked")]
    tab <- tab[order(-tab$mean_abs_dHSS, na.last = TRUE), ]
    rownames(tab) <- NULL
    tab_in <- as.data.frame(et_in)
    tab_in <- tab_in[order(-tab_in$mean_abs_dHSS, na.last = TRUE), ]
    rownames(tab_in) <- NULL
    list(effect_table = et, audit = tab, audit_in_sample = tab_in)
  })
  data.table::fwrite(audit$audit,
                     file.path(out_dir, "attribution_audit.csv"))
  data.table::fwrite(audit$audit_in_sample,
                     file.path(out_dir, "tables", "attribution_audit_in_sample.csv"))

  # Map the driver with the LARGEST single-shift mean absolute response
  # (supported first; masked drivers carry NA and sort last). The map carries
  # the same hard masking as the table.
  top_var <- audit$audit$driver[1]
  top_row <- audit$audit[1, , drop = FALSE]
  cli::cli_inform("Mapping {.val {top_var}}: largest single-shift model response (mean_abs_dHSS = {round(top_row$mean_abs_dHSS, 4)}, support = {round(top_row$support, 3)}). Descriptive, not causal.")

  eff_fp <- digest::digest(list(fp, top_var, fit$env_vars, cfg$effect_shift,
                                "step" = "effect-map-v6-holdout",
                                cfg$effect_support_probs, cfg$effect_min_support,
                                cfg$effect_shift_type, projection_inputs$current,
                                step = "effect-map-v5-single-shift"))
  effect_map <- cache_step(file.path(cache_dir, "effect_map.rds"), eff_fp, {
    em <- cast_effect_map(fit_full, terra::mask(current[[fit_full$env_vars]], area_mask), drivers = fit_full$env_vars,
                          shift = cfg$effect_shift,
                          shift_type = cfg$effect_shift_type,
                          support_probs = cfg$effect_support_probs,
                          min_support = cfg$effect_min_support,
                          weights = projection$current$weights,
                          block_rows = max(1L, as.integer(cfg$attribution_block_size /
                            terra::ncol(current))),
                          filename = file.path(out_dir,
                                               paste0("effect_map_", top_var, ".tif")),
                          overwrite = TRUE, verbose = TRUE)
    list(driver = top_var, effect_table = attr(em, "effect_table"),
         raster_file = file.path(out_dir, paste0("effect_map_", top_var, ".tif")))
  })

  data.table::fwrite(as.data.frame(audit$effect_table), file.path(table_dir, "effect_table.csv"))
  data.table::fwrite(effect_map$effect_table, file.path(table_dir, "effect_map_summary.csv"))
  }

  if (!identical(source_info, workflow_source_signature(cfg))) {
    stop("Workflow sources changed during execution; no current-source bundle can be issued.")
  }
  result <- new_cast_result(screen, fit, evaluation, cv = cv, fit_full = fit_full)
  # A weighted ensemble that collapses to a single engine is a scientific
  # fact about the run, not a detail: castSDM records it on the object, so
  # carry it into the bundle instead of leaving it only in the console log.
  ensemble_state <- list(
    method = holdout_ensemble$method,
    weights = holdout_ensemble$weights,
    weights_fallback = holdout_ensemble$weights_fallback,
    fallback_reason = holdout_ensemble$fallback_reason,
    threshold = holdout_ensemble$threshold,
    decay = holdout_ensemble$decay)
  if (isTRUE(holdout_ensemble$weights_fallback)) {
    cli::cli_warn("The hold-out ensemble collapsed to a fallback weighting: ",
                  "{holdout_ensemble$fallback_reason}. See tables/holdout_ensemble.csv ",
                  "and ensemble_state in castSDM_result.rds.")
  }
  saveRDS(list(result = result, projection = projection, config = cfg,
               audit = audit, effect_map = effect_map, attribution = attribution,
               ensemble_state = ensemble_state,
               ensemble_importance = ensemble_imp,
               mode = cfg$run_mode, provenance = source_info),
          file.path(out_dir, "castSDM_result.rds"))
  utils::capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
  write_species_reports(cfg, result, audit)
  figure_dir <- make_species_figures(cfg)
  if (!identical(source_info, workflow_source_signature(cfg))) {
    stop("Workflow sources changed during reporting; completion is not certified.")
  }
  if (!is.null(figure_dir)) {
    manifest <- readRDS(file.path(out_dir, "figure_manifest.rds"))
    status <- manifest$figure_status
    if (!nrow(status) || any(status$status != "written_not_visually_checked") ||
        any(manifest$source_data$status != "present")) {
      stop("Case exports are incomplete; inspect figure_status.csv and source_data.csv.")
    }
    # Validation-chain gate. Two fold outcomes must be distinguished:
    #   selection_error / failed  -> the selector or fitter threw inside the
    #     fold. That is a DEFECT signal (a bug, e.g. the fold-internal
    #     covariate-degeneracy class), and the case is not certified.
    #   single_class -> the fold's train or test split carried one response
    #     class. For clustered occurrences on a national background this is a
    #     DATA property: any blocking, fold count or seed leaves at least one
    #     fold without presences. It is recorded as a caveat and reported in
    #     the tables, but it does not by itself block certification -- the
    #     fold-coverage floor (>= 80% of folds) still has to be met.
    n_folds_total <- as.integer(cfg$cv_k)
    defect_folds <- cv$fold_status %in%
      c("selection_error", "model_error", "no_predictions", "failed")
    if (any(defect_folds)) {
      stop(sprintf(
        paste0("CV folds failed with a defect signal (%s in fold(s) %s).",
               " This is a code/selector failure, not a data property;",
               " completion is not certified. Inspect cv$screens diagnostics",
               " and tables/cv_fold_status.csv."),
        paste(unique(cv$fold_status[defect_folds]), collapse = "/"),
        paste(which(defect_folds), collapse = ",")), call. = FALSE)
    }
    n_eval <- sum(cv$fold_status == "evaluated")
    floor_folds <- max(2L, ceiling(n_folds_total * 0.8))
    if (n_eval < floor_folds) {
      stop(sprintf(
        paste0("Only %d of %d CV folds evaluated (floor %d); ensemble threshold",
               " and CV metrics rest on a partial out-of-fold surface.",
               " Inspect tables/cv_fold_coverage.csv."),
        n_eval, n_folds_total, floor_folds), call. = FALSE)
    }
    # single_class folds (if any) ride along as an explicit caveat.
    skipped <- which(cv$fold_status == "single_class")
    caveats <- if (length(skipped)) {
      sprintf("%d of %d spatial folds carried a single response class and were skipped (fold %s); CV metrics rest on %d folds.",
              length(skipped), n_folds_total, paste(skipped, collapse = ","), n_eval)
    } else character(0)
    artifacts <- c(file.path(out_dir, "castSDM_result.rds"),
      file.path(out_dir, "figure_manifest.rds"), manifest$source_data$path,
      status$path)
    saveRDS(list(status = "computed_and_exported", mode = cfg$run_mode,
      source = source_info, artifacts = tools::md5sum(unique(artifacts)),
      engines = names(fit$models), scenarios = names(projection$future),
      visual_QA = "not performed",
      cv_caveats = caveats), file.path(out_dir, "completion.rds"))
  }
  invisible(result)
}
