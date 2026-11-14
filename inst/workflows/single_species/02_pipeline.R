# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# One package pipeline, followed by one paired scenario comparison.
case_training_year <- function(preparation, split) {
  qc <- data.table::fread(file.path(preparation, "occurrence_qc.csv"), data.table = FALSE)
  if (!all(c("x_albers", "y_albers", "year", "used_in_model") %in% names(qc))) return(NULL)
  qc <- qc[qc$used_in_model %in% TRUE, , drop = FALSE]
  key <- function(x, y) paste(round(x, 3), round(y, 3), sep = "|")
  pres <- split$train[split$train$presence == 1L, ]
  idx <- match(key(pres$lon, pres$lat), key(qc$x_albers, qc$y_albers))
  year <- suppressWarnings(as.numeric(qc$year[idx]))
  year <- year[is.finite(year) & year > 1500 & year < 2100]
  if (length(year)) stats::median(year) else NULL
}

case_match_training_layers <- function(data, split, meta, reference_year, mode) {
  keep <- rep(TRUE, nrow(meta))
  if (!is.null(reference_year) && mode %in% c("npts", "nts")) {
    for (id in unique(meta$construct_id)) {
      ix <- which(meta$construct_id == id & is.finite(meta$year))
      if (!length(ix)) next
      years <- meta$year[ix]
      chosen <- if (mode == "npts") {
        if (any(years <= reference_year)) max(years[years <= reference_year]) else min(years)
      } else years[which.min(abs(years-reference_year))]
      keep[ix[years != chosen]] <- FALSE
    }
  }
  meta <- meta[keep, , drop = FALSE]
  vars <- meta$layer
  cols <- c("lon", "lat", "presence", vars)
  split$train <- split$train[, cols, drop = FALSE]
  split$test <- split$test[, cols, drop = FALSE]
  split$env_vars <- vars
  list(data = data[, cols, drop = FALSE], split = split, meta = meta)
}

run_species <- function(cfg) {
  run_id <- Sys.getenv("CAST_RUN_ID", "cast032_chelsa_ensmean")
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
  claim <- list(schema = "cast_case_claim/1", config = scientific_cfg,
    package = as.character(utils::packageVersion("castSDM")),
    sources = unname(tools::md5sum(file.path(cfg$example_root,
      c("02_pipeline.R", "prepare_data.R", "data_adapter.R", "scenario_case.R")))))
  claim$key <- digest::digest(claim, algo = "sha256")
  previous <- if (file.exists(claim_file)) readRDS(claim_file) else NULL
  if ((!is.null(previous) && !identical(previous$key, claim$key)) ||
      (is.null(previous) && file.exists(file.path(out, "case.rds"))))
    stop("This output belongs to a different/older analysis. Set a new CAST_RUN_ID: ", out)
  if (is.null(previous)) saveRDS(claim, claim_file)
  # A failed rerun must never leave a previous completion marker.
  completion <- file.path(out, "completion.rds")
  if (file.exists(completion)) unlink(completion)
  preparation <- prepare_species(cfg)
  qc <- data.table::fread(file.path(preparation, "occurrence_qc.csv"), data.table = FALSE)
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
    interpretation = "Pooled historical records; spatial hold-out is not external or temporal validation."),
    file.path(out, "case_design.csv"), row.names = FALSE)
  data <- data.table::fread(file.path(preparation,
    paste0(cfg$species_name, "_model_data.csv")), data.table = FALSE)
  p <- cfg$protocol
  split <- castSDM::cast_prepare(data, train_fraction = p$train_fraction,
    n_blocks = p$holdout_blocks, block_method = p$block_method, seed = p$seed,
    coords_crs = p$coords_crs, block_size = p$block_size)
  meta <- case_pool_metadata(cfg, split$env_vars)
  year <- case_training_year(preparation, split)
  if (is.null(year)) message("Training occurrence years unavailable: temporal matching skipped.")
  matched <- case_match_training_layers(data, split, meta, year, cfg$select_temporal_match)
  meta <- matched$meta
  utils::write.csv(data.frame(reference_year = year %||% NA_real_,
    source = if (is.null(year)) "unavailable" else "training_presences_only",
    method = cfg$select_temporal_match), file.path(out, "temporal_matching.csv"), row.names = FALSE)
  utils::write.csv(meta, file.path(out, "training_pool.csv"), row.names = FALSE)
  result <- castSDM::cast_run(matched$data, prepared = matched$split, protocol = p,
    select_var_id = stats::setNames(meta$construct_id, meta$layer),
    select_category = stats::setNames(meta$group, meta$layer),
    select_roles = stats::setNames(meta$role, meta$layer),
    cache_dir = cfg$cache_root, output_dir = file.path(out, "training"),
    report_meta = list(taxon = cfg$species_name,
      purpose = "CHELSA ensemble-mean climate projections; nonclimate predictors held current"))
  selected <- result$screen$selected
  current <- wf_load_current_stack(cfg, varnames = selected)
  scenarios <- case_scenarios(cfg, selected, current, out)
  mask <- terra::rast(file.path(preparation, "analysis_domain_mask.tif"))
  comparison <- castSDM::cast_compare(result, current, scenarios$future,
    groups = scenarios$groups, mask = mask, support_probs = cfg$support_probs,
    domain = "common", resources = cfg$resources,
    output_dir = file.path(out, "comparison"), thresholds = cfg$diagnostic_thresholds)
  castSDM::cast_diagnose(result, comparison = comparison,
    output_dir = file.path(out, "diagnostics"))
  case <- list(schema = "cast_single_species/1", species = cfg$species_name,
    cfg = cfg, report_scenario = scenarios$report, preparation = preparation,
    training = "training/result.rds", comparison = "comparison/comparison.rds",
    script_hashes = cfg$runtime$script_hashes %||%
      tools::md5sum(list.files(cfg$example_root, pattern = "\\.R$", full.names = TRUE)),
    session = utils::capture.output(utils::sessionInfo()))
  file.copy(file.path(cfg$example_root, "config", "predictor_metadata.csv"),
    file.path(out, "predictor_metadata.csv"), overwrite = TRUE)
  saveRDS(case, file.path(out, "case.rds"))
  # Every available figure is exported. Driver/PDP/MESS products are computed
  # once from frozen fits, then cached independently of style and rendering.
  case_figures(comparison, out, scenarios$report, result = result, cfg = cfg,
    current = current, future = scenarios$future[[scenarios$report]])
  sensitivity <- summary(comparison, thresholds = cfg$diagnostic_thresholds)
  utils::write.csv(sensitivity, file.path(out, "threshold_sensitivity.csv"), row.names = FALSE)
  saveRDS(list(finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    package_version = as.character(utils::packageVersion("castSDM")), case = "case.rds"), completion)
  message("Outputs: ", out)
  invisible(out)
}
