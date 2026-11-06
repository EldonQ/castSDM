# Scientific settings live in one cast_protocol(); paths and data QC stay here.
castsdm_species_configs <- function(root) {
  data_root <- Sys.getenv("CAST_DATA_ROOT", "E:/CHECO26/N-SDM/data")
  occurrence_root <- Sys.getenv("CAST_OCCURRENCE_ROOT",
    "E:/Package/castSDM_workflows/species_occ")
  seed <- 42L
  resources <- list(cores = 12L, workers = 3L, memory_mb = 8000)
  common <- list(
    example_root = root,
    nsdm_data_root = data_root,
    output_root = file.path(root, "output"),
    preparation_root = file.path(root, "output", ".prepared"),
    cache_root = file.path(root, "output", ".training-cache"),
    future_layers_csv = file.path(root, "config", "future_layers.csv"),
    input_type = "gbif", domain_mode = "national",
    analysis_domain_mask_file = NULL,
    freshwater_river_distance_raster = NULL,
    freshwater_river_buffer_m = 100, freshwater_occurrence_buffer_km = 10,
    domain_min_presence_coverage = .95,
    occurrence_uncertainty_m = 1000, occurrence_keep_unknown_uncertainty = TRUE,
    occurrence_thin_km = 1, occurrence_thin_order = "random",
    occurrence_basis_of_record = c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION",
      "OBSERVATION", "PRESERVED_SPECIMEN", "MATERIAL_SAMPLE", "MATERIAL_CITATION",
      "LITERATURE", "OCCURRENCE"),
    preprocessing_scope = "all_records",
    bg_target_group = FALSE,
    bg_target_group_file = file.path(data_root, "species/glo/gbif_china_alltaxa_geom.psv"),
    bg_ratio = 10, bg_min = 10000L, bg_max = 20000L,
    bg_strategy = "random", bg_bin_k = 20L, bg_sre_quantile = .025,
    select_temporal_match = "npts", seed = seed,
    protocol = castSDM::cast_protocol(
      models = c("rf", "maxent", "gam"), select_method = "screen",
      select_corcut = .7, select_maxncov = 12L,
      select_ensure_n = c(climate = 2L, lulc = 1L),
      train_fraction = .7, holdout_blocks = 20L, block_method = "blockcv",
      cv_k = 5L, rf_ntree = 500L, tune = FALSE, aoa = TRUE,
      ensemble_method = "weighted", ensemble_min_score = 0,
      prevalence_target = .5, clamp = TRUE, num_threads = resources$cores, seed = seed),
    # One total budget: 3 scenario processes, at most 4 RF threads each.
    resources = resources,
    support_probs = c(.01, .99), driver_attribution = FALSE,
    scenario_climate_model = "bcc_csm2_mr",
    scenario_ssps = c("ssp126", "ssp245", "ssp585"),
    scenario_periods = c("2021_2040", "2041_2070", "2071_2100"),
    scenario_climate_years = c(`2021_2040` = 2030, `2041_2070` = 2055, `2071_2100` = 2085),
    scenario_lulc_years = c(`2021_2040` = 2030, `2041_2070` = 2050, `2071_2100` = 2100),
    scenario_report_ssp = "ssp585", scenario_report_period = "2071_2100")
  list(
    Schizothoracinae = utils::modifyList(common, list(
      species_name = "Schizothoracinae",
      taxon_level = "subfamily-level literature records",
      occurrence_file = file.path(occurrence_root, "Schizothoracinae_schisoft.csv"),
      nsdm_pool_csv = file.path(root, "config", "pool_current_eligible_glu.csv"),
      scenario_lulc_product = "glu", domain_mode = "freshwater_river",
      freshwater_river_distance_raster = file.path(data_root,
        "covariates/glo/dist/h90m/stream-dist-proximity/glo_dist_h90m_stream-dist-proximity.tif"))),
    Callosciurus_erythraeus = utils::modifyList(common, list(
      species_name = "Callosciurus_erythraeus", taxon_level = "species-level GBIF records",
      input_type = "nsdm_occ",
      occurrence_file = file.path(occurrence_root, "Callosciurus_erythraeus.csv"),
      nsdm_pool_csv = file.path(root, "config", "pool_current_eligible_hlu.csv"),
      scenario_lulc_product = "hlu",
      scenario_lulc_years = c(`2021_2040` = 2030, `2041_2070` = 2055, `2071_2100` = 2085))))
}
