# ==============================================================================
# castSDM 多物种并行示例（与 run_ovis_ammon.R 同一套 CONFIG 形参风格）
#
# 与单物种 ovis 脚本对齐之处：
#   - 单一 CONFIG 列表集中配置（cast_prepare / cast_dag / cast_ate /
#     cast_screen / cast_fit / cast_evaluate / cast_predict / cast_cv /
#     cast_cate / SHAP 等与 ovis 字段对应）
#   - 出图 DPI、ggsave(bg=white)、空间图 basemap、一致性图字体等与包内
#     cast_batch() 升级版一致
#
# 多物种增量：
#   - future 按物种并行；每个物种独立子目录
#     <output_dir>/<Species>/figures/ 与 cast_result.rds
#   - 顶层 <output_dir>/figures/multi_species_comparison.png
#
# 运行：在 castSDM 项目下 load_all() 后 source 本文件；或安装包后
#   source(system.file("examples/run_multi_species.R", package = "castSDM"))
#
# 依赖建议：future / future.apply / ggplot2 / patchwork；与 ovis 相同
#   的 DAG / SHAP 可选包见 run_ovis_ammon.R 文首说明。
#
# 说明：本脚本不包含「五种 DAG structure_method」自检；自检仅在
#   run_ovis_ammon.R 中提供。多物种流程在每条物种上跑与 ovis 同构的
#   prepare → DAG → ATE → screen → fit → evaluate →（可选）CV → predict
#   → CATE → 一致性图 →（可选）SHAP，结果写入各自子目录。
# ==============================================================================

if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1))) {
  devtools::load_all()
} else {
  library(castSDM)
}

for (pkg in c("ggplot2", "future", "future.apply", "patchwork", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# ==============================================================================
# 集中配置（字段命名与 inst/examples/run_ovis_ammon.R 的 CONFIG 对齐）
# 多物种专用项：data_dir, species_files, env_grid_path, output_dir,
#   n_workers, parallel_species, do_shap（批量默认 FALSE 以免极慢）
# ==============================================================================
CONFIG <- list(
  # --- multi-species only ---
  data_dir          = NULL, # 非 NULL 时扫描目录内 *.csv，忽略 species_files
  species_files     = list(
    Ovis_ammon           = "CAST_Ovis_ammon_Res9_screened.csv",
    Gazella_subgutturosa = "CAST_Gazella_subgutturosa_Res9_screened.csv",
    Pseudois_nayaur      = "CAST_Pseudois_nayaur_Res9_screened.csv"
  ),
  env_grid_path     = NULL, # 非 NULL 时用该路径；否则尝试包内中国环境栅格
  output_dir        = "castSDM_multi_species",
  fig_dpi           = 600L,
  n_workers         = min(3L, max(1L, parallel::detectCores() - 1L)),
  parallel_species  = TRUE,
  batch_verbose     = TRUE,

  response          = "presence",
  seed              = 42L,
  prepare_train_fraction = 0.7,
  prepare_env_vars       = NULL,
  prepare_verbose        = FALSE,
  dag_env_vars            = NULL,
  dag_R                   = 100L,
  dag_algorithm           = "hc",
  dag_score               = "bic-g",
  dag_strength_threshold  = 0.7,
  dag_direction_threshold = 0.6,
  dag_max_rows            = 8000L,
  dag_verbose             = FALSE,
  dag_structure_method    = "bootstrap_hc",
  dag_pc_alpha            = 0.05,
  dag_pc_test             = "zf",
  dag_fci_alpha           = 0.05,
  dag_fci_test            = "zf",
  dag_bidag_algorithm     = "order",
  dag_bidag_iterations    = NULL,
  dag_notears_lambda      = 0.03,
  dag_notears_max_iter    = 2000L,
  dag_notears_lr          = 0.02,
  dag_notears_tol         = 1e-3,
  dag_notears_rho_init    = 0.1,
  dag_notears_alpha_mult  = 1.01,
  ate_variables     = NULL,
  ate_K             = 5L,
  ate_num_trees     = 300L,
  ate_alpha         = 0.05,
  ate_quantile_cuts = c(0.1, 0.25, 0.5, 0.75, 0.9),
  ate_bonferroni    = TRUE,
  ate_parallel      = TRUE,
  ate_verbose       = FALSE,
  screen_min_vars     = 5L,
  screen_min_fraction = 0.5,
  screen_num_trees    = 500L,
  screen_verbose      = FALSE,
  fit_models            = c("cast", "rf", "maxent", "brt"),
  fit_n_epochs          = 50L,
  fit_n_runs            = 1L,
  fit_patience          = 40L,
  fit_val_fraction      = 0.2,
  fit_focal_gamma       = 2.0,
  fit_focal_alpha_mode  = "sqrt",
  fit_focal_alpha       = 0.75,
  fit_rf_ntree          = 300L,
  fit_brt_n_trees       = 500L,
  fit_brt_depth         = 5L,
  fit_hidden_size       = NULL,
  fit_dropout           = 0.2,
  fit_lr                = 1e-3,
  fit_batch_size        = NULL,
  fit_max_interactions  = 15L,
  fit_tune_grid         = FALSE,
  fit_verbose           = FALSE,
  eval_response   = NULL,
  predict_models  = NULL,
  plot_basemap    = "china",
  run_spatial_cv  = TRUE,
  cv_k             = 5L,
  cv_models        = NULL,
  cv_block_method  = "grid",
  cv_n_epochs      = 120L,
  cv_n_runs        = 2L,
  cv_rf_ntree      = 300L,
  cv_brt_n_trees   = 500L,
  cv_parallel      = TRUE,
  cv_verbose       = FALSE,
  do_cate          = TRUE,
  cate_variables   = NULL,
  cate_top_n       = 3L,
  cate_n_trees     = 300L,
  cate_verbose     = FALSE,
  cate_point_size  = 0.45,
  cate_hss_model      = "cast",
  cate_hss_threshold  = 0.1,
  var_labels          = NULL,
  do_shap             = FALSE,
  shap_nrounds          = 200L,
  shap_max_depth        = 6L,
  shap_eta              = 0.05,
  shap_subsample        = 0.8,
  shap_colsample_bytree = 0.8,
  shap_test_fraction    = 0.2,
  shap_verbose          = FALSE,
  shap_plot_top_n       = 15L,
  shap_fastshap_nsim    = 40L,
  shap_max_explain_rows = 50L
)

# 可选：与 ovis 相同的中文变量标签（仅影响图例）
# CONFIG$var_labels <- c(bio02 = "昼夜温差", ...)

cat("============================================================\n")
cat("castSDM multi-species (CONFIG aligned with run_ovis_ammon.R)\n")
cat("============================================================\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Build species_list + env_grid
# ══════════════════════════════════════════════════════════════════════════════

if (!is.null(CONFIG$data_dir) && dir.exists(CONFIG$data_dir)) {
  csv_paths <- list.files(
    CONFIG$data_dir, pattern = "\\.csv$",
    full.names = TRUE, recursive = FALSE
  )
  if (length(csv_paths) == 0L) {
    stop("No CSV files found in CONFIG$data_dir: ", CONFIG$data_dir)
  }
  sp_names <- tools::file_path_sans_ext(basename(csv_paths))
  species_list <- stats::setNames(
    lapply(csv_paths, function(f) data.table::fread(f, data.table = FALSE)),
    sp_names
  )
  cat(sprintf("[Auto-scan] %d species in %s\n",
              length(sp_names), CONFIG$data_dir))
} else {
  species_list <- lapply(CONFIG$species_files, function(f) {
    path <- system.file("extdata", f, package = "castSDM")
    if (path == "") stop("File not found in package extdata: ", f)
    data.table::fread(path, data.table = FALSE)
  })
}

for (sp in names(species_list)) {
  d <- species_list[[sp]]
  cat(sprintf(
    "  %s: %d rows | presence=%d | absence=%d\n",
    sp, nrow(d), sum(d[[CONFIG$response]] == 1),
    sum(d[[CONFIG$response]] == 0)
  ))
}

env_grid <- NULL
if (!is.null(CONFIG$env_grid_path) && file.exists(CONFIG$env_grid_path)) {
  env_grid <- data.table::fread(CONFIG$env_grid_path, data.table = FALSE)
  cat(sprintf("\nPrediction grid (user file): %d x %d\n",
              nrow(env_grid), ncol(env_grid)))
} else {
  pkg_grid <- system.file(
    "extdata", "China_EnvData_Res9_Screened.csv",
    package = "castSDM"
  )
  if (pkg_grid != "") {
    env_grid <- data.table::fread(pkg_grid, data.table = FALSE)
    cat(sprintf("\nPrediction grid (bundled): %d x %d\n",
                nrow(env_grid), ncol(env_grid)))
  } else {
    cat("\nNo prediction grid; spatial prediction skipped.\n")
  }
}

dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Parallel backend
# ══════════════════════════════════════════════════════════════════════════════

use_par <- isTRUE(CONFIG$parallel_species) && CONFIG$n_workers > 1L
cat(sprintf(
  "\nParallel: %s | workers: %d (cores: %d)\n",
  if (use_par) "ON" else "OFF",
  if (use_par) CONFIG$n_workers else 1L,
  parallel::detectCores()
))
if (use_par) {
  future::plan(future::multisession, workers = CONFIG$n_workers)
} else {
  future::plan(future::sequential)
}

# ══════════════════════════════════════════════════════════════════════════════
# cast_batch (all CONFIG-driven)
# ══════════════════════════════════════════════════════════════════════════════

cat("\nRunning cast_batch()...\n\n")

eval_resp <- CONFIG$eval_response
if (is.null(eval_resp)) eval_resp <- CONFIG$response

batch_result <- cast_batch(
  species_list   = species_list,
  env_data       = env_grid,
  models         = CONFIG$fit_models,
  output_dir     = CONFIG$output_dir,
  fig_dpi        = CONFIG$fig_dpi,
  parallel       = use_par,
  seed           = CONFIG$seed,
  verbose        = CONFIG$batch_verbose,
  fit_verbose    = CONFIG$fit_verbose,
  response       = CONFIG$response,
  train_fraction = CONFIG$prepare_train_fraction,
  prepare_env_vars   = CONFIG$prepare_env_vars,
  prepare_verbose    = CONFIG$prepare_verbose,
  dag_env_vars       = CONFIG$dag_env_vars,
  dag_verbose        = CONFIG$dag_verbose,
  dag_R                   = CONFIG$dag_R,
  dag_structure_method    = CONFIG$dag_structure_method,
  dag_pc_alpha            = CONFIG$dag_pc_alpha,
  dag_pc_test             = CONFIG$dag_pc_test,
  dag_fci_alpha           = CONFIG$dag_fci_alpha,
  dag_fci_test            = CONFIG$dag_fci_test,
  dag_bidag_algorithm     = CONFIG$dag_bidag_algorithm,
  dag_bidag_iterations    = CONFIG$dag_bidag_iterations,
  dag_notears_lambda      = CONFIG$dag_notears_lambda,
  dag_notears_max_iter    = CONFIG$dag_notears_max_iter,
  dag_notears_lr          = CONFIG$dag_notears_lr,
  dag_notears_tol         = CONFIG$dag_notears_tol,
  dag_notears_rho_init    = CONFIG$dag_notears_rho_init,
  dag_notears_alpha_mult  = CONFIG$dag_notears_alpha_mult,
  dag_algorithm           = CONFIG$dag_algorithm,
  dag_score               = CONFIG$dag_score,
  dag_strength_threshold  = CONFIG$dag_strength_threshold,
  dag_direction_threshold = CONFIG$dag_direction_threshold,
  dag_max_rows            = CONFIG$dag_max_rows,
  ate_variables     = CONFIG$ate_variables,
  ate_K             = CONFIG$ate_K,
  ate_alpha         = CONFIG$ate_alpha,
  ate_num_trees     = CONFIG$ate_num_trees,
  ate_quantile_cuts = CONFIG$ate_quantile_cuts,
  ate_bonferroni    = CONFIG$ate_bonferroni,
  ate_parallel      = CONFIG$ate_parallel,
  ate_verbose       = CONFIG$ate_verbose,
  screen_min_vars     = CONFIG$screen_min_vars,
  screen_min_fraction = CONFIG$screen_min_fraction,
  screen_num_trees    = CONFIG$screen_num_trees,
  screen_verbose      = CONFIG$screen_verbose,
  do_cv               = CONFIG$run_spatial_cv,
  cv_k                = CONFIG$cv_k,
  cv_models           = CONFIG$cv_models,
  cv_block_method     = CONFIG$cv_block_method,
  cv_n_epochs         = CONFIG$cv_n_epochs,
  cv_n_runs           = CONFIG$cv_n_runs,
  cv_rf_ntree         = CONFIG$cv_rf_ntree,
  cv_brt_n_trees      = CONFIG$cv_brt_n_trees,
  cv_parallel         = CONFIG$cv_parallel,
  cv_verbose          = CONFIG$cv_verbose,
  do_cate             = CONFIG$do_cate,
  cate_top_n          = CONFIG$cate_top_n,
  cate_n_trees        = CONFIG$cate_n_trees,
  cate_variables      = CONFIG$cate_variables,
  cate_verbose        = CONFIG$cate_verbose,
  cate_point_size     = CONFIG$cate_point_size,
  cate_hss_model      = CONFIG$cate_hss_model,
  cate_hss_threshold  = CONFIG$cate_hss_threshold,
  predict_models      = CONFIG$predict_models,
  plot_basemap        = CONFIG$plot_basemap,
  eval_response       = eval_resp,
  var_labels          = CONFIG$var_labels,
  do_shap             = CONFIG$do_shap,
  shap_nrounds          = CONFIG$shap_nrounds,
  shap_max_depth        = CONFIG$shap_max_depth,
  shap_eta              = CONFIG$shap_eta,
  shap_subsample        = CONFIG$shap_subsample,
  shap_colsample_bytree = CONFIG$shap_colsample_bytree,
  shap_test_fraction    = CONFIG$shap_test_fraction,
  shap_verbose          = CONFIG$shap_verbose,
  shap_plot_top_n       = CONFIG$shap_plot_top_n,
  shap_fastshap_nsim    = CONFIG$shap_fastshap_nsim,
  shap_max_explain_rows = CONFIG$shap_max_explain_rows,
  n_epochs         = CONFIG$fit_n_epochs,
  n_runs           = CONFIG$fit_n_runs,
  patience         = CONFIG$fit_patience,
  val_fraction     = CONFIG$fit_val_fraction,
  focal_gamma      = CONFIG$fit_focal_gamma,
  focal_alpha_mode = CONFIG$fit_focal_alpha_mode,
  focal_alpha      = CONFIG$fit_focal_alpha,
  rf_ntree         = CONFIG$fit_rf_ntree,
  brt_n_trees      = CONFIG$fit_brt_n_trees,
  brt_depth         = CONFIG$fit_brt_depth,
  hidden_size      = CONFIG$fit_hidden_size,
  dropout          = CONFIG$fit_dropout,
  lr               = CONFIG$fit_lr,
  batch_size       = CONFIG$fit_batch_size,
  max_interactions = CONFIG$fit_max_interactions,
  tune_grid        = CONFIG$fit_tune_grid
)

print(batch_result)

# ══════════════════════════════════════════════════════════════════════════════
# Cross-species comparison
# ══════════════════════════════════════════════════════════════════════════════

fig_dir <- file.path(CONFIG$output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n[Plot] Multi-species model performance comparison...\n")

p_compare <- plot(batch_result, metrics = c("auc", "tss", "cbi", "sedi"))
ggplot2::ggsave(
  file.path(fig_dir, "multi_species_comparison.png"),
  p_compare,
  width = 14, height = 6, dpi = CONFIG$fig_dpi,
  bg = "white", limitsize = FALSE
)
cat("  Saved: multi_species_comparison.png\n")

future::plan(future::sequential)

cat("\n============================================================\n")
cat("Multi-species batch complete.\n")
cat("Output: ", normalizePath(CONFIG$output_dir, winslash = "/",
                              mustWork = FALSE), "\n", sep = "")
cat("============================================================\n\n")

for (sp in names(batch_result$results)) {
  sp_figs <- file.path(CONFIG$output_dir, sp, "figures")
  if (dir.exists(sp_figs)) {
    figs <- list.files(sp_figs, pattern = "\\.png$")
    cat(sprintf("  %s: %d PNG(s)\n", sp, length(figs)))
  }
}
