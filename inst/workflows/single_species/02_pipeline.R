# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# One package pipeline, followed by one paired scenario comparison.
case_training_year <- function(preparation, split) {
  qc <- data.table::fread(file.path(preparation, "occurrence_qc.csv"), data.table = FALSE)
  if (!all(c("x_albers", "y_albers", "year", "used_in_model") %in% names(qc))) return(NULL)
  qc <- qc[qc$used_in_model %in% TRUE, , drop = FALSE]
  key <- function(x, y) paste(round(x, 3), round(y, 3), sep = "|")
  pres <- split$train[split$train$presence == 1L, ]
  idx <- match(key(pres$lon, pres$lat), key(qc$x_albers, qc$y_albers))
  column <- if ("sampling_year_verified" %in% names(qc)) "sampling_year_verified" else "year"
  year <- suppressWarnings(as.numeric(qc[[column]][idx]))
  year <- year[is.finite(year) & year > 1500 & year < 2100]
  if (!length(year)) return(NULL)
  list(year = stats::median(year), source = if (column == "sampling_year_verified")
    "reported_unambiguous_sampling_years_in_training_presences" else "unverified_record_years_in_training_presences")
}

case_match_training_layers <- function(data, split, meta, reference_year, mode) {
  if (!is.null(reference_year)) {
    if (!is.numeric(reference_year) || length(reference_year) != 1L || !is.finite(reference_year))
      stop("reference_year must be NULL or one finite numeric year.")
    reference_year <- as.numeric(reference_year)
  }
  keep <- rep(TRUE, nrow(meta)); rule <- rep("no_single_year_in_inventory", nrow(meta))
  if (!is.null(reference_year) && mode %in% c("npts", "nts")) {
    for (id in unique(meta$construct_id)) {
      ix <- which(meta$construct_id == id & is.finite(meta$year))
      if (!length(ix)) next
      years <- meta$year[ix]
      chosen <- if (mode == "npts") {
        if (any(years <= reference_year)) max(years[years <= reference_year]) else min(years)
      } else years[which.min(abs(years-reference_year))]
      keep[ix[years != chosen]] <- FALSE
      rule[ix] <- ifelse(years != chosen, "other_time_slice",
        if (mode == "npts" && chosen > reference_year) "no_past_slice_available" else mode)
    }
  } else rule[is.finite(meta$year)] <- "no_reference_year_or_matching_disabled"
  matching <- data.frame(layer = meta$layer, construct = meta$construct_id,
    source_year = meta$year, reference_year = rep(reference_year %||% NA_real_, nrow(meta)),
    retained = keep, rule = rule, stringsAsFactors = FALSE)
  before <- nrow(meta); meta <- meta[keep, , drop = FALSE]
  message("Temporal layers: ", before, " -> ", nrow(meta), "; method=", mode,
    "; reference year=", reference_year %||% "unavailable", ".")
  future_only <- sum(keep & rule == "no_past_slice_available")
  if (future_only) message("Temporal matching: ", future_only,
    " retained layers have no earlier time slice; retained source years are later than the reference. This is not contemporaneous matching; see temporal_layer_matching.csv.")
  vars <- meta$layer; cols <- c("lon", "lat", "presence", vars)
  split$train <- split$train[, cols, drop = FALSE]
  split$test <- split$test[, cols, drop = FALSE]; split$env_vars <- vars
  list(data = data[, cols, drop = FALSE], split = split, meta = meta, matching = matching)
}

run_species <- function(cfg) {
  run_id <- Sys.getenv("CAST_RUN_ID", castSDM:::.cast_workflow_run_id())
  if (!grepl("^[A-Za-z0-9_-]+$", run_id)) stop("CAST_RUN_ID must contain only letters, digits, _ and -.")
  out <- file.path(cfg$output_root, run_id, cfg$species_name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(out, ".case.lock"), timeout = 0)
  if (is.null(lock)) stop("Another process is writing this case: ", out)
  on.exit(filelock::unlock(lock), add = TRUE)
  claim_file <- file.path(out, "run_claim.rds")
  scientific_cfg <- cfg[setdiff(names(cfg), c("runtime", "resources", "example_root",
    "output_root", "cache_root", "preparation_root", "diagnostic_thresholds",
    "figure_response_method", "figure_response_grid", "figure_display_resolution"))]
  scientific_cfg$protocol$num_threads <- NULL
  # Numerical stages check their own code and content. Pin the declared inputs
  # and data adapter here so a repair to execution can resume the same case.
  claim <- list(schema = "cast_case_claim/2", config = scientific_cfg,
    sources = unname(tools::md5sum(file.path(cfg$example_root,
      c("prepare_data.R", "data_adapter.R", "scenario_case.R")))))
  claim$key <- digest::digest(claim, algo = "sha256")
  previous <- if (file.exists(claim_file)) readRDS(claim_file) else NULL
  compatible_previous <- identical(previous$schema, "cast_case_claim/1") &&
    identical(previous$config, claim$config) && identical(previous$sources[-1L], claim$sources)
  if ((!is.null(previous) && !identical(previous$key, claim$key) && !compatible_previous) ||
      (is.null(previous) && file.exists(file.path(out, "case.rds"))))
    stop("This output belongs to a different/older analysis. Set a new CAST_RUN_ID: ", out)
  if (is.null(previous) || compatible_previous) castSDM:::.cast_atomic_rds(claim, claim_file)
  # A failed rerun must never leave a previous completion marker.
  completion <- file.path(out, "completion.rds")
  if (file.exists(completion)) unlink(completion)
  case_start <- Sys.time(); case_clock <- proc.time()[["elapsed"]]
  timings <- list()
  record_time <- function(stage, start, seconds, status) {
    timings[[length(timings)+1L]] <<- data.frame(stage = stage,
      started = format(start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      elapsed_s = seconds, status = status)
    utils::write.csv(do.call(rbind, timings), file.path(out, "workflow_timings.csv"), row.names = FALSE)
  }
  timed <- function(stage, expr) {
    start <- Sys.time(); clock <- proc.time()[["elapsed"]]; status <- "failed_or_interrupted"
    on.exit({
      seconds <- proc.time()[["elapsed"]]-clock
      try(record_time(stage, start, seconds, status), silent = TRUE)
      message("Workflow stage ", stage, ": ", round(seconds, 1), " s; ", status, ".")
    }, add = TRUE)
    value <- force(expr); status <- "completed"; value
  }
  on.exit(try(record_time("workflow_wall", case_start,
    proc.time()[["elapsed"]]-case_clock, if (file.exists(completion)) "completed" else "incomplete"),
    silent = TRUE), add = TRUE, after = FALSE)
  message("Record-year interpretation: ", cfg$occurrence_year_semantics,
    "; spatial validation remains conditional on the declared preprocessing domain.")
  preparation <- timed("preparation", prepare_species(cfg))
  qc <- data.table::fread(file.path(preparation, "occurrence_qc.csv"), data.table = FALSE)
  used <- qc$used_in_model %in% TRUE
  taxon_names <- if ("species" %in% names(qc)) qc$species else rep(cfg$species_name, nrow(qc))
  if ("sampling_date_status" %in% names(qc)) utils::write.csv(as.data.frame(table(
    sampling_date_status = qc$sampling_date_status[used], useNA = "always")),
    file.path(out, "sampling_date_quality.csv"), row.names = FALSE)
  excluded_taxa <- as.data.frame(table(taxon = taxon_names[!used], reason = qc$qc_reason[!used]))
  utils::write.csv(excluded_taxa[excluded_taxa$Freq > 0, ], file.path(out, "excluded_records_by_taxon.csv"), row.names = FALSE)
  years <- suppressWarnings(as.numeric(qc$year[qc$used_in_model %in% TRUE]))
  periods <- cut(years, c(-Inf, 1949, 1979, 1999, 2009, 2019, Inf),
    labels = c("before_1950", "1950_1979", "1980_1999", "2000_2009", "2010_2019", "2020_onward"))
  year_counts <- as.data.frame(table(periods, useNA = "always"))
  names(year_counts) <- c("record_period", "n_model_presences")
  utils::write.csv(year_counts, file.path(out, "occurrence_years.csv"), row.names = FALSE)
  utils::write.csv(data.frame(taxon = cfg$species_name, taxon_level = cfg$taxon_level,
    year_semantics = cfg$occurrence_year_semantics,
    declared_first_year = (cfg$occurrence_year_range %||% c(NA_real_, NA_real_))[1L],
    declared_last_year = (cfg$occurrence_year_range %||% c(NA_real_, NA_real_))[2L],
    preprocessing_scope = cfg$preprocessing_scope,
    domain = cfg$domain_mode, n_taxon_names = length(unique(taxon_names[used])),
    n_possible_coarse_coordinates = if ("qc_tenth_degree_grid" %in% names(qc)) sum(qc$qc_tenth_degree_grid[used], na.rm = TRUE) else NA_integer_,
    source_metadata_available = "sampling_year_verified" %in% names(qc),
    interpretation = "Pooled historical records; spatial hold-out is not external or temporal validation."),
    file.path(out, "case_design.csv"), row.names = FALSE)
  data <- data.table::fread(file.path(preparation,
    paste0(cfg$species_name, "_model_data.csv")), data.table = FALSE)
  p <- cfg$protocol
  split <- timed("spatial_split", castSDM::cast_prepare(data, train_fraction = p$train_fraction,
    n_blocks = p$holdout_blocks, block_method = p$block_method, seed = p$seed,
    coords_crs = p$coords_crs, block_size = p$block_size))
  if (is.data.frame(split$fold_design$class_balance))
    utils::write.csv(split$fold_design$class_balance, file.path(out,"split_class_balance.csv"), row.names=FALSE)
  meta <- case_pool_metadata(cfg, split$env_vars)
  temporal_reference <- case_training_year(preparation, split)
  year <- temporal_reference$year
  if (is.null(year)) message("Usable training sampling years unavailable: temporal matching skipped.") else
    message("Temporal reference source: ", temporal_reference$source, ". A reference year is not record-wise environmental matching.")
  matched <- timed("temporal_matching", case_match_training_layers(data, split, meta, year, cfg$select_temporal_match))
  utils::write.csv(matched$matching, file.path(out, "temporal_layer_matching.csv"), row.names = FALSE)
  prepared_meta <- meta
  meta <- matched$meta
  utils::write.csv(data.frame(reference_year = year %||% NA_real_,
    source = if (is.null(year)) "unavailable" else temporal_reference$source,
    method = cfg$select_temporal_match), file.path(out, "temporal_matching.csv"), row.names = FALSE)
  utils::write.csv(meta, file.path(out, "training_pool.csv"), row.names = FALSE)
  response_data <- timed("response_validation_inputs",
    case_response_validation(cfg, matched$split, meta))
  result <- timed("training_validation", castSDM::cast_run(matched$data, prepared = matched$split, protocol = p,
    select_var_id = stats::setNames(meta$construct_id, meta$layer),
    select_category = stats::setNames(meta$group, meta$layer),
    select_roles = stats::setNames(meta$role, meta$layer),
    response_data = response_data,
    cache_dir = cfg$cache_root, output_dir = file.path(out, "training"),
    report_meta = list(taxon = cfg$species_name,
      purpose = "CHELSA ensemble-mean climate projections; nonclimate predictors held current")))
  selected <- result$screen$selected
  # Preserve the meaning of repeated library layers without inventing an
  # ecological ranking or silently choosing a depth/scale for the user.
  semantics <- meta[, intersect(c("layer", "group", "product", "construct_id", "scale", "stat", "year", "role"), names(meta)), drop = FALSE]
  semantics$selected <- semantics$layer %in% selected
  semantics$construct_variants <- as.integer(table(meta$construct_id)[meta$construct_id])
  semantics$spatial_context <- if (cfg$domain_mode != "freshwater_river")
    rep("Gridded environmental predictor; inspect source, scale and date for ecological relevance", nrow(semantics)) else ifelse(semantics$group == "soil",
    "Terrestrial grid property; not a measured in-stream condition",
    ifelse(semantics$group == "bioclim", "Gridded atmospheric climate; not measured water temperature",
      "See source definition for grid, neighbourhood or upstream aggregation"))
  utils::write.csv(semantics, file.path(out, "predictor_semantics.csv"), row.names = FALSE)
  groups <- sort(unique(prepared_meta$group))
  pool_summary <- data.frame(group = groups,
    prepared = vapply(groups, function(g) sum(prepared_meta$group == g), integer(1)),
    temporally_retained = vapply(groups, function(g) sum(meta$group == g), integer(1)),
    selected = vapply(groups, function(g) sum(meta$group == g & meta$layer %in% selected), integer(1)))
  utils::write.csv(pool_summary, file.path(out, "predictor_pool_summary.csv"), row.names = FALSE)
  message("Retained predictor categories: ", paste(paste0(pool_summary$group, "=", pool_summary$selected), collapse = ", "),
    ". VIF is a collinearity check, not ecological or causal relevance selection.")
  message("Predictor context: ", sum(semantics$selected & semantics$construct_variants > 1L),
    " retained layers share a construct with another candidate; see predictor_semantics.csv for scale and source context. No fixed predictor quota is imposed.")
  excluded <- qc[qc$qc_reason == "missing_environment", , drop = FALSE]
  if (nrow(excluded) && "missing_predictors" %in% names(excluded)) {
    excluded$missing_final_predictors <- vapply(strsplit(excluded$missing_predictors, ";", fixed = TRUE),
      function(z) paste(intersect(z, selected), collapse = ";"), character(1))
    excluded$complete_on_final_predictors <- !nzchar(excluded$missing_final_predictors)
    excluded$cohort_policy <- "Excluded before splitting on the declared candidate pool; not reintroduced after hold-out evaluation"
    utils::write.csv(excluded, file.path(out, "excluded_cohort_final_predictors.csv"), row.names = FALSE)
    message("Frozen validation cohort: ", sum(excluded$complete_on_final_predictors),
      " excluded records are complete on final predictors. They are not added after validation; candidate-pool or missingness policy changes require a new declared run.")
  }
  current <- timed("selected_environment", wf_load_current_stack(cfg, varnames = selected))
  scenarios <- timed("scenario_inputs", case_scenarios(cfg, selected, current, out))
  mask <- terra::rast(file.path(preparation, "analysis_domain_mask.tif"))
  comparison <- timed("scenario_comparison", castSDM::cast_compare(result, current, scenarios$future,
    groups = scenarios$groups, mask = mask, support_probs = cfg$support_probs,
    domain = "common", resources = cfg$resources,
    output_dir = file.path(out, "comparison"), thresholds = cfg$diagnostic_thresholds))
  footprint <- comparison$summary[, c("scenario", "n_valid", "n_supported", "support_share")]
  domain_summary <- utils::read.csv(file.path(preparation, "analysis_domain_summary.csv"))
  footprint$domain_cells <- domain_summary$domain_cells[1L]
  footprint$national_mask_cells <- domain_summary$national_mask_cells[1L]
  footprint$prediction_share_national <- footprint$n_valid/footprint$national_mask_cells
  footprint$n_inside_domain_without_paired_prediction <- footprint$domain_cells-footprint$n_valid
  footprint$scope <- "Declared analysis domain; continuous maps do not discard values below the binary threshold."
  utils::write.csv(footprint, file.path(out, "spatial_output_qc.csv"), row.names = FALSE)
  diagnostic <- timed("diagnostics", castSDM::cast_diagnose(result, comparison = comparison,
    output_dir = file.path(out, "diagnostics")))
  sensitivity <- diagnostic$response_spatial_sensitivity
  if (is.data.frame(sensitivity) && "n_common" %in% names(sensitivity))
    message("Spatial-refit explanation check: ", sensitivity$n_available[1L], " fits on ",
      sensitivity$n_locations[1L], " unchanged hold-out locations; ", sensitivity$n_common[1L],
      " common central-quantile locations (", sensitivity$n_common_presence[1L],
      " presences; ", sensitivity$n_common_background[1L],
      " backgrounds). Observed-range counts are reported separately in diagnostics/response_spatial_sensitivity.csv; neither criterion establishes joint support or confidence.")
  terms <- diagnostic$gam_concurvity
  if (is.data.frame(terms) && nrow(terms)) {
    terms <- terms[terms$term != "para", ]
    message("GAM response diagnostics: ", nrow(terms), " smooths; ",
      sum(terms$estimate > .9, na.rm = TRUE), " with concurvity estimate > 0.9; ",
      sum(terms$edf < .1, na.rm = TRUE), " with EDF < 0.1. These summaries do not remove predictors or identify causal effects.")
  }
  report_scope <- footprint[footprint$scenario == scenarios$report, ]
  message("Reporting map scope: ", report_scope$n_valid, " paired valid cells; ",
    report_scope$n_supported, " inside the declared central marginal-quantile envelope. Outside this envelope does not necessarily mean outside the observed training range. Blank areas are outside the declared analysis domain or lack paired predictions. Grid area must be interpreted according to this domain.")
  case <- list(schema = "cast_single_species/1", species = cfg$species_name,
    cfg = cfg, report_scenario = scenarios$report, preparation = preparation,
    training = "training/result.rds", comparison = "comparison/comparison.rds",
    script_hashes = cfg$runtime$script_hashes %||%
      tools::md5sum(list.files(cfg$example_root, pattern = "\\.R$", full.names = TRUE)),
    session = utils::capture.output(utils::sessionInfo()))
  file.copy(file.path(cfg$example_root, "config", "predictor_metadata.csv"),
    file.path(out, "predictor_metadata.csv"), overwrite = TRUE)
  saveRDS(case, file.path(out, "case.rds"))
  # Main figures use scenario responses and frozen fits. Complete atlases and
  # extra diagnostics remain available through 04_REPLOT.R.
  timed("environmental_response_products", case_figures(comparison, out, scenarios$report,
    result = result, cfg = cfg, current = current,
    future = scenarios$future[[scenarios$report]], prepare_only = TRUE))
  timed("figure_exports", case_figures(comparison, out, scenarios$report,
    result = result, cfg = cfg, render_only = TRUE))
  sensitivity <- summary(comparison, thresholds = cfg$diagnostic_thresholds)
  utils::write.csv(sensitivity, file.path(out, "threshold_sensitivity.csv"), row.names = FALSE)
  saveRDS(list(finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    package_version = as.character(utils::packageVersion("castSDM")), case = "case.rds"), completion)
  message("Outputs: ", out)
  invisible(out)
}
