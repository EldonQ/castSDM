# Schizopygopsis workflow: 596-layer audited pool, explicit sampling and model settings.
source("E:/Package/castSDM_workflows/R/workflow_provenance.R",encoding="UTF-8")

castsdm_species_configs <- function() {
  common <- list(
    cast_root = "E:/Package/cast",
    run_mode = "production",
    output_root = "E:/Package/castSDM_workflows/single_species/output/case_runs",
    adapter = "E:/Package/castSDM_workflows/R/nsdm_castsdm_adapter.R",
    nsdm_pool_csv = "E:/Package/castSDM_workflows/single_species/pool_layers_full.csv",
    nsdm_data_root = "E:/CHECO26/N-SDM/data",
    n_stack_layers = 596L,
    input_type = "nsdm_occ",
    basemap = "none",
    study_area_method = "buffer",
    study_area_buffer_km = 300,
    occurrence_thin_km = 1,
    occurrence_thin_order = "random",
    bg_target_group = FALSE,
    bg_target_group_file = "E:/CHECO26/N-SDM/data/species/glo/gbif_china_alltaxa_geom.psv",
    background_domain = "national",
    bg_ratio = 10, bg_min = 10000L, bg_max = 20000L,
    bg_strategy = "random",
    bg_bin_k = 20L, bg_sre_quantile = 0.025,
    models = c("rf", "brt", "maxent", "gam"),
    rf_ntree = 500L, rf_mtry = NULL,
    brt_n_trees = 2000L, brt_depth = 5L, brt_shrinkage = 0.005,
    maxent_classes = NULL, maxent_regmult = NULL,
    model_num_threads = 1L, tune_folds = 3L,
    select_method = "covsel",
    select_corcut = 0.7,
    select_algorithms = c("glm", "gam", "rf"),
    select_ncov = NULL,
    select_maxncov = 12L,
    select_keep = character(0),
    select_ensure_climate = TRUE,
    select_varid_strip = c("_f1000", "_f5000", "_f10000", "_f20000",
                           "_max", "_min"),
    select_temporal_match = "npts",
    effect_shift = 1,
    effect_shift_type = "sd",
    require_core_validation = FALSE,
    cv_k = 5L, seed = 42L,
    cv_threshold_method = "max_tss",
    cv_workers = 5L,
    train_fraction = 0.7,
    # Spatial blocking follows castSDM 0.20.0: blockCV square blocks sized
    # from the estimated range of spatial autocorrelation (cv_block_size
    # overrides, in metres); falls back to the dependency-free grid blocking
    # with a warning when the design cannot be built.
    cv_block_method = "blockcv",
    cv_block_size = NULL,
    cv_buffer = 0,
    cv_buffer_unit = "deg",
    cv_aoa = TRUE,
    cv_n_repeat = 1L,
    cv_tune = FALSE,
    ensemble_method = "weighted",
    ensemble_min_score = 0.5,
    ensemble_min_metric_folds = 2L,
    ensemble_decay = 1,
    ensemble_fallback = "best",
    ensemble_permutation = FALSE,
    ensemble_perm_n = 25L,
    ensemble_perm_metric = "pearson",
    min_cast_version = "0.20.0",
    # 0.20.0 default: every engine (including MaxEnt, via a stored
    # logit-scale offset) is calibrated so the weighted training prevalence
    # equals this target; ensemble members then share one scale.
    prevalence_target = 0.5,
    # 0.20.0 default: predictors are clamped to the training range at
    # prediction/projection time for every engine (MESS/AOA flags are still
    # computed on the unclamped input).
    clamp = TRUE,
    effect_support_probs = c(0.01, 0.99),
    effect_min_support = 0.5,
    attribution_block_size = 50000L,
    single_shift_diagnostics = FALSE,
    report_scenario = "ensmean_ssp585_2071_2100",
    projection_workers = 6L,
    future_models = c("bcc_csm2_mr", "ensmean"),
    ssp_scenarios = c("ssp126", "ssp245", "ssp585"),
    periods = c("2021_2040", "2041_2070", "2071_2100")
  )

  occ_dir <- "E:/Package/castSDM_workflows/species_occ"

  list(
    Moschus_berezovskii = utils::modifyList(common, list(
      species_name = "Moschus_berezovskii",
      taxon_level = "species",
      occurrence_file = file.path(occ_dir, "Moschus_berezovskii.csv")
    )),
    Schizothoracinae = utils::modifyList(common, list(
      species_name = "Schizothoracinae",
      taxon_level = "subfamily-level record collection (SchiSOFT literature records)",
      input_type = "gbif",
      occurrence_file = "E:/Package/castSDM_workflows/species_occ/Schizothoracinae_schisoft.csv",
      occurrence_uncertainty_m = 1000,
      occurrence_keep_unknown_uncertainty = TRUE
    )),
    Schizopygopsis = utils::modifyList(common, list(
      species_name = "Schizopygopsis",
      taxon_level = "genus-level record collection (SchiSOFT literature records)",
      input_type = "gbif",
      occurrence_file = "E:/Package/castSDM_workflows/species_occ/Schizopygopsis_schisoft.csv",
      occurrence_uncertainty_m = 1000,
      occurrence_keep_unknown_uncertainty = TRUE
    ))
  )
}
