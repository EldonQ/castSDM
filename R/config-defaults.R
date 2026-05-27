#' Build castSDM Default Configuration
#'
#' Returns a named list of workflow defaults that scripts can override with
#' dataset-specific choices. This keeps common DAG, selection, fitting, CV,
#' CATE, and SHAP parameters aligned across examples without forcing every
#' dataset into the same run plan.
#'
#' @param profile Character. One of `"single"`, `"batch"`, `"disdat"`,
#'   `"fish"`, or `"debug"`.
#' @param overrides Optional named list merged on top of the defaults.
#'
#' @return A named list.
#' @export
cast_default_config <- function(profile = c("single", "batch", "disdat", "fish", "debug"),
                                overrides = list()) {
  profile <- match.arg(profile)

  cfg <- list(
    response = "presence",
    seed = 42L,
    models = c("rf", "brt", "maxent", "gam"),
    train_fraction = 0.7,
    output_dir = "castSDM_output",
    fig_dpi = 300L,
    plot_font_family = "Arial",
    parallel = TRUE,
    resume = TRUE,

    dag_env_vars = NULL,
    dag_R = 100L,
    dag_structure_method = "mb_first",
    dag_include_response = TRUE,
    dag_pc_alpha = 0.05,
    dag_pc_test = NULL,
    dag_mb_method = "fast.iamb",
    dag_mb_alpha = 0.05,
    dag_algorithm = "hc",
    dag_score = NULL,
    dag_strength_threshold = 0.7,
    dag_direction_threshold = 0.6,
    dag_max_rows = 8000L,
    dag_verbose = FALSE,
    learn_shared_dag = FALSE,
    do_backdoor = TRUE,

    select_min_vars = 5L,
    select_min_fraction = 0.3,
    select_num_trees = 300L,
    select_verbose = FALSE,

    fit_rf_ntree = 300L,
    fit_brt_n_trees = 500L,
    fit_brt_depth = 5L,
    fit_verbose = FALSE,

    do_cv = TRUE,
    cv_k = 5L,
    cv_block_method = "grid",
    cv_models = NULL,
    cv_parallel = FALSE,
    cv_verbose = FALSE,

    do_cate = TRUE,
    cate_top_n = 3L,
    cate_n_trees = 1000L,
    cate_variables = NULL,
    cate_verbose = FALSE,
    cate_hss_model = "rf",
    cate_hss_threshold = 0.1,

    do_shap = FALSE,
    shap_nrounds = 200L,
    shap_max_depth = 6L,
    shap_eta = 0.05,
    shap_subsample = 0.8,
    shap_colsample_bytree = 0.8,
    shap_test_fraction = 0.2,
    shap_plot_top_n = 15L,
    shap_fastshap_nsim = 40L,
    shap_max_explain_rows = 50L,

    do_predict = TRUE,
    do_ensemble = TRUE,
    do_spatial_heatmap = FALSE,
    run_future_projection = FALSE,
    future_download_cmip6 = FALSE,
    future_gcms = NULL,
    future_ssps = NULL,
    future_periods = NULL,
    future_var = NULL,
    future_res = NULL,
    future_ensemble = TRUE,
    future_cache_dir = NULL,
    future_save_dir = NULL
  )

  if (identical(profile, "batch")) {
    cfg$output_dir <- "castSDM_multi_species"
    cfg$parallel <- TRUE
    cfg$fig_dpi <- 600L
  } else if (identical(profile, "disdat")) {
    cfg$output_dir <- "castSDM_disdat"
    cfg$models <- c("rf", "gam", "maxent", "brt")
    cfg$fig_dpi <- 300L
    cfg$do_shap <- TRUE
    cfg$cate_n_trees <- 500L
    cfg$learn_shared_dag <- FALSE
  } else if (identical(profile, "fish")) {
    cfg$output_dir <- "castSDM_fish"
    cfg$fig_dpi <- 1200L
    cfg$select_min_fraction <- 0.5
    cfg$fit_rf_ntree <- 300L
  } else if (identical(profile, "debug")) {
    cfg$output_dir <- "castSDM_debug"
    cfg$dag_max_rows <- 1500L
    cfg$dag_R <- 20L
    cfg$select_num_trees <- 100L
    cfg$fit_rf_ntree <- 100L
    cfg$fit_brt_n_trees <- 150L
    cfg$do_cv <- FALSE
    cfg$do_cate <- FALSE
    cfg$do_shap <- FALSE
    cfg$parallel <- FALSE
  }

  cast_merge_config(cfg, overrides)
}

#' Merge castSDM Configuration Lists
#'
#' Thin wrapper around [utils::modifyList()] with a clearer name for scripts.
#'
#' @param defaults Named list of default values.
#' @param overrides Named list of dataset-specific values.
#'
#' @return A named list.
#' @export
cast_merge_config <- function(defaults, overrides = list()) {
  if (is.null(overrides)) return(defaults)
  if (!is.list(defaults) || !is.list(overrides)) {
    cli::cli_abort("{.fn cast_merge_config} expects two lists.")
  }
  utils::modifyList(defaults, overrides, keep.null = TRUE)
}
