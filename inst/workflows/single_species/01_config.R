# Resolve one explicitly configured project independently of the taxon name.
castsdm_species_configs <- function(root) {
  settings <- castSDM:::.cast_workflow_selection(root, allow_unconfigured = TRUE)$settings
  if (is.null(settings)) return(list())
  data_root <- settings$data_dir
  output <- Sys.getenv("CAST_OUTPUT_ROOT", file.path(root, "output"))
  product <- settings$lulc_product %||% "glu"
  resources <- settings$resources %||% list(cores = 4L, workers = 1L, memory_mb = 4096)
  seed <- 42L
  protocol <- castSDM::cast_protocol(
    models = settings$models %||% c("rf", "maxent", "gam"),
    select_vif_threshold = settings$vif_threshold %||% 10,
    block_method = settings$block_method %||% "auto",
    block_size = settings$block_size, holdout_block_size = settings$holdout_block_size,
    cv_k = settings$cv_k %||% 5L, tune = FALSE, aoa = TRUE,
    ensemble_min_score = 0, prevalence_target = NULL, gam_backend = "bam",
    maxent_classes = "lq", seed = seed, num_threads = resources$cores,
    coords_crs = terra::crs(terra::rast(file.path(data_root, "mask_china.tif"))))
  cfg <- list(
    species_name = settings$species, occurrence_file = settings$occurrence_file,
    taxon_level = "Taxon as identified in the configured record source",
    example_root = root, nsdm_data_root = data_root, output_root = output,
    preparation_root = file.path(output, ".prepared"),
    cache_root = file.path(output, ".training-cache", settings$species),
    nsdm_pool_csv = file.path(root, "config", paste0("pool_current_eligible_", product, ".csv")),
    future_layers_csv = file.path(root, "config", "future_layers.csv"),
    input_type = settings$input_type %||% "gbif", current_lulc_product = product,
    domain_mode = settings$domain %||% "national", analysis_domain_mask_file = settings$domain_mask,
    validation_target = settings$validation_target %||% "spatial_transfer",
    predictor_file = settings$predictor_file,
    occurrence_uncertainty_m = NULL, occurrence_keep_unknown_uncertainty = TRUE,
    occurrence_year_range = NULL, occurrence_keep_unknown_year = TRUE,
    occurrence_year_semantics = "Source dates and record years retained; matching requires a reported sampling year",
    occurrence_thin_km = settings$thin_km %||% 0, occurrence_thin_order = "random",
    preprocessing_scope = "configured records inside the analysis mask",
    bg_n = settings$background_points, bg_strategy = "random",
    bg_target_group = !is.null(settings$target_group_file), bg_target_group_file = settings$target_group_file,
    predictor_min_coverage = settings$predictor_min_coverage %||% .95,
    predictor_missing = settings$predictor_missing %||% "complete_cases",
    select_temporal_match = "npts", seed = seed, protocol = protocol, resources = resources,
    projection_cores = settings$projection_cores,
    support_probs = c(.01, .99), diagnostic_thresholds = c(.3, .5, .7),
    figure_response_method = "ale", figure_response_grid = 50L, figure_display_resolution = .02,
    scenario_climate_model = "ensmean", climate_baseline = "1981_2010",
    scenario_ssps = c("ssp126", "ssp370", "ssp585"),
    scenario_periods = c("2011_2040", "2041_2070", "2071_2100"),
    scenario_report_ssp = "ssp585", scenario_report_period = "2071_2100")
  stats::setNames(list(cfg), settings$species)
}
