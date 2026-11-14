# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# Scientific settings live in one cast_protocol(); paths and data QC stay here.
castsdm_species_configs <- function(root) {
  settings_file <- file.path(root, "workflow_settings.rds")
  settings <- if (file.exists(settings_file)) readRDS(settings_file) else NULL
  data_root <- Sys.getenv("CAST_DATA_ROOT", settings$data_dir %||% "E:/CHECO26/N-SDM/data")
  scope <- case_covariate_scope(data_root)
  occurrence_root <- Sys.getenv("CAST_OCCURRENCE_ROOT",
    "E:/Package/castSDM_workflows/species_occ")
  seed <- 42L
  resources <- list(cores = 12L, workers = 3L, memory_mb = 8000)
  common <- list(
    example_root = root,
    nsdm_data_root = data_root,
    output_root = Sys.getenv("CAST_OUTPUT_ROOT", file.path(root, "output")),
    preparation_root = file.path(root, "output", ".prepared"),
    cache_root = file.path(root, "output", ".training-cache"),
    future_layers_csv = file.path(root, "config", "future_layers.csv"),
    input_type = "gbif", domain_mode = "national",
    analysis_domain_mask_file = NULL,
    freshwater_river_distance_raster = NULL,
    freshwater_river_buffer_m = 100, freshwater_occurrence_buffer_km = 10,
    domain_min_presence_coverage = .95,
    occurrence_uncertainty_m = 1000, occurrence_keep_unknown_uncertainty = TRUE,
    # Keep the completed study's cohort unless a new analysis is declared.
    occurrence_year_range = NULL, occurrence_keep_unknown_year = TRUE,
    occurrence_year_semantics = "record year; sampling/publication meaning unverified",
    occurrence_thin_km = 1, occurrence_thin_order = "random",
    occurrence_basis_of_record = c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION",
      "OBSERVATION", "PRESERVED_SPECIMEN", "MATERIAL_SAMPLE", "MATERIAL_CITATION",
      "LITERATURE", "OCCURRENCE"),
    preprocessing_scope = "all_records",
    bg_target_group = FALSE,
    bg_target_group_file = NULL,
    bg_n = 10000L, predictor_min_coverage = .95,
    bg_strategy = "random", bg_bin_k = 20L, bg_sre_quantile = .025,
    select_temporal_match = "npts", seed = seed,
    protocol = castSDM::cast_protocol(
      models = c("rf", "maxent", "gam"), select_method = "screen",
      select_inner_k = 5L,
      train_fraction = .7, holdout_blocks = 20L, block_method = "blockcv",
      cv_k = 5L, rf_ntree = 500L, tune = FALSE, aoa = TRUE,
      ensemble_method = "weighted", ensemble_min_score = 0,
      prevalence_target = NULL, gam_backend = "bam",
      coords_crs = terra::crs(terra::rast(file.path(data_root, "mask_china.tif"))), clamp = TRUE, num_threads = resources$cores, seed = seed),
    # One total budget: 3 scenario processes, at most 4 RF threads each.
    resources = resources,
    support_probs = c(.01, .99), diagnostic_thresholds = c(.3, .5, .7),
    figure_response_method = "ale", figure_response_grid = 50L, figure_display_resolution = .02,
    scenario_climate_model = "ensmean", climate_baseline = "1981_2010",
    scenario_ssps = c("ssp126", "ssp370", "ssp585"),
    scenario_periods = c("2011_2040", "2041_2070", "2071_2100"),
    scenario_climate_years = c(`2011_2040` = 2025, `2041_2070` = 2055, `2071_2100` = 2085),
    scenario_report_ssp = "ssp585", scenario_report_period = "2071_2100")
  presets <- list(
    Schizothoracinae = utils::modifyList(common, list(
      species_name = "Schizothoracinae",
      taxon_level = "subfamily-level literature records",
      occurrence_file = file.path(occurrence_root, "Schizothoracinae_schisoft.csv"),
      nsdm_pool_csv = file.path(root, "config", "pool_current_eligible_glu.csv"),
      current_lulc_product = "glu", domain_mode = "freshwater_river",
      freshwater_river_distance_raster = file.path(data_root,
        "covariates", scope, "dist/h90m/stream-dist-proximity",
        paste0(scope, "_dist_h90m_stream-dist-proximity_min.tif")))),
    Callosciurus_erythraeus = utils::modifyList(common, list(
      species_name = "Callosciurus_erythraeus", taxon_level = "species-level GBIF records",
      input_type = "nsdm_occ",
      occurrence_file = file.path(occurrence_root, "Callosciurus_erythraeus.csv"),
      nsdm_pool_csv = file.path(root, "config", "pool_current_eligible_hlu.csv"),
      current_lulc_product = "hlu")))
  if (is.null(settings)) return(presets)
  taxon <- settings$species
  base <- presets[[taxon]] %||% utils::modifyList(common, list(
    species_name = taxon, taxon_level = "user-declared taxon", current_lulc_product = "glu"))
  base$occurrence_file <- settings$occurrence_file
  base$nsdm_data_root <- data_root
  base$domain_mode <- settings$domain %||% base$domain_mode
  base$input_type <- settings$input_type %||% base$input_type
  base$current_lulc_product <- settings$lulc_product %||% base$current_lulc_product
  base$nsdm_pool_csv <- file.path(root, "config", paste0("pool_current_eligible_", base$current_lulc_product, ".csv"))
  base$freshwater_river_distance_raster <- file.path(base$nsdm_data_root, "covariates/reg/dist/h90m/stream-dist-proximity/reg_dist_h90m_stream-dist-proximity_min.tif")
  base$bg_n <- settings$background_points %||% base$bg_n
  base$resources <- settings$resources %||% base$resources
  base$protocol$models <- settings$models %||% base$protocol$models
  base$protocol$cv_k <- settings$cv_k %||% base$protocol$cv_k
  base$protocol$num_threads <- base$resources$cores
  presets[[taxon]] <- base
  presets
}
