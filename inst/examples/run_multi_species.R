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
#   - 物种并行（有 dev_package_root 时为 PSOCK）；每个物种
#     <output_dir>/<Species>/figures/ 与 cast_result.rds
#   - SHAP：do_shap=TRUE 时写 shap_xgb/rf/cast 各两张图 + shap_panel_2x3.png（2 行 3 列）
#   - 已跑完 batch 仅补 SHAP：CONFIG$only_shap_posthoc <- TRUE 后 source（不重训）
#   - 空间 HSS/CATE：包原生为 geom_point 格点图；CONFIG$only_replot_spatial_heatmap
#     TRUE 时用 terra 热图覆盖各物种 figures/HSS_*.png、CATE_*.png（需 cast_result.rds）
#   - 顶层 <output_dir>/figures/multi_species_comparison.png
#
# 在 RStudio 中运行方式（与 run_ovis_ammon.R 一致）：
#   1. 用 RStudio 打开 castSDM.Rproj（工作目录在包根目录）
#   2. Ctrl+Shift+L (devtools::load_all()) 或运行下面「包加载」代码块
#   3. source() 本脚本；若从子目录 source，代码会从当前目录向上查找
#      含 castSDM 的 DESCRIPTION 的包根并 load_all(路径)，避免用到旧版
#      已安装包里的 cast_dag（否则会出现 structure_method 等 unused arguments）。
#
# 或已安装**当前源码对应版本**的包后：
#   library(castSDM)
#   source(system.file("examples/run_multi_species.R", package = "castSDM"))
#
# # 方式一：开发模式 (推荐) — 若 getwd() 不在包根，可显式指定：
# # devtools::load_all("E:/Package/cast")
#
# # 方式二：安装后再加载
# # devtools::install("E:/Package/cast")
# # library(castSDM)
#
# 并行开发模式：cast_batch(dev_package_root=...) 用 PSOCK，各 worker 先 pkgload::load_all。
# 依赖建议：future / future.apply / ggplot2 / patchwork / pkgload；与 ovis 相同
#   的 DAG / SHAP 可选包见 run_ovis_ammon.R 文首说明。
#
# 说明：本脚本不包含「五种 DAG structure_method」自检；自检仅在
#   run_ovis_ammon.R 中提供。多物种流程在每条物种上跑与 ovis 同构的
#   prepare → DAG → ATE → screen → fit → evaluate →（可选）CV → predict
#   → CATE → 一致性图 →（可选）SHAP，结果写入各自子目录。
# ==============================================================================

.cast_find_package_root <- function(max_up = 8L) {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (.i in seq_len(max_up)) {
    desc <- file.path(d, "DESCRIPTION")
    if (file.exists(desc)) {
      l1 <- tryCatch(readLines(desc, 1L, warn = FALSE), error = function(e) "")
      if (length(l1) && nzchar(l1[1L]) && grepl("castSDM", l1[1L], fixed = TRUE)) {
        return(d)
      }
    }
    desc_cast <- file.path(d, "cast", "DESCRIPTION")
    if (file.exists(desc_cast)) {
      l1 <- tryCatch(readLines(desc_cast, 1L, warn = FALSE), error = function(e) "")
      if (length(l1) && nzchar(l1[1L]) && grepl("castSDM", l1[1L], fixed = TRUE)) {
        return(normalizePath(file.path(d, "cast"), winslash = "/", mustWork = FALSE))
      }
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  env <- Sys.getenv("CASTSDM_ROOT", "")
  if (nzchar(env)) {
    for (env_try in c(env, file.path(env, "cast"))) {
      desc <- file.path(env_try, "DESCRIPTION")
      if (file.exists(desc)) {
        l1 <- tryCatch(readLines(desc, 1L, warn = FALSE), error = function(e) "")
        if (length(l1) && nzchar(l1[1L]) && grepl("castSDM", l1[1L], fixed = TRUE)) {
          return(normalizePath(env_try, winslash = "/", mustWork = FALSE))
        }
      }
    }
  }
  NA_character_
}

# 热图 helper：工作区在 e:/Package、脚本在 cast/inst/examples 时 pkg_root 可能为 NA，
# 故多路径查找（含 getwd()/cast/inst/examples、向上爬目录）。
.cast_find_spatial_heatmap_helpers <- function(pkg_root) {
  bn <- "cast_spatial_heatmap_helpers.R"
  cand <- character(0)
  if (!is.na(pkg_root) && nzchar(pkg_root)) {
    cand <- c(cand, file.path(pkg_root, "inst", "examples", bn))
  }
  cr <- Sys.getenv("CASTSDM_ROOT", "")
  if (nzchar(cr)) {
    cand <- c(
      cand,
      file.path(cr, "inst", "examples", bn),
      file.path(cr, "cast", "inst", "examples", bn)
    )
  }
  cand <- c(cand, system.file(file.path("examples", bn), package = "castSDM"))
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (.i in seq_len(12L)) {
    cand <- c(
      cand,
      file.path(d, "inst", "examples", bn),
      file.path(d, "cast", "inst", "examples", bn)
    )
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  cand <- unique(cand[nzchar(cand)])
  for (p in cand) {
    pp <- tryCatch(
      normalizePath(p, winslash = "/", mustWork = FALSE),
      error = function(e) ""
    )
    if (nzchar(pp) && file.exists(pp)) {
      return(pp)
    }
  }
  NA_character_
}

# 开发模式：与 run_ovis_ammon.R 相同逻辑，并支持从 inst/examples 等子目录 source
pkg_root <- .cast_find_package_root()
if (!is.na(pkg_root)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::load_all(pkg_root)
  # future::multisession 子进程只认 library() 里的包；cast_batch 会读此变量并在各 worker 内 pkgload::load_all
  Sys.setenv(CASTSDM_ROOT = normalizePath(pkg_root, winslash = "/", mustWork = FALSE))
} else {
  Sys.unsetenv("CASTSDM_ROOT")
  library(castSDM)
}

for (pkg in c("ggplot2", "future", "future.apply", "patchwork", "data.table", "pkgload")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# ==============================================================================
# 集中配置（字段命名与 inst/examples/run_ovis_ammon.R 的 CONFIG 对齐）
# 多物种专用项：data_dir, species_files, output_dir, n_workers,
#   only_shap_posthoc, only_replot_spatial_heatmap, spatial_heatmap_res_deg, …
# ==============================================================================
CONFIG <- list(
  # --- multi-species only ---
  only_shap_posthoc  = FALSE, # TRUE：跳过 cast_batch，只读各物种 cast_result.rds 补 SHAP（不重训）
  only_replot_spatial_heatmap = TRUE, # TRUE：仅热图重绘；需 terra 与同仓库 cast_spatial_heatmap_helpers.R
  spatial_heatmap_res_deg      = 0.06,
  spatial_heatmap_interp_method = "nearest",
  spatial_heatmap_display_res   = 0.02,
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
  do_shap             = TRUE,
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
# 仅重绘空间 HSS / CATE 热图（nearest 插值栅格 + 中国边界 shap 遮罩；覆盖各物种 figures/HSS_*.png、CATE_*.png）
# ══════════════════════════════════════════════════════════════════════════════

if (isTRUE(CONFIG$only_replot_spatial_heatmap)) {
  future::plan(future::sequential)
  hlp <- .cast_find_spatial_heatmap_helpers(pkg_root)
  if (is.na(hlp) || !nzchar(hlp)) {
    stop(
      "Could not find cast_spatial_heatmap_helpers.R.\n",
      "  pkg_root = ", deparse(pkg_root), "\n",
      "  getwd()  = ", getwd(), "\n",
      "Set CASTSDM_ROOT to the castSDM source root, or run from e:/Package/cast (or parent with cast/ subfolder).",
      call. = FALSE
    )
  }
  sys.source(hlp, envir = .GlobalEnv)
  spn <- names(species_list)
  cat("\n[only_replot_spatial_heatmap] Re-saving HSS/CATE heatmaps...\n")
  for (sp in spn) {
    rds <- file.path(CONFIG$output_dir, sp, "cast_result.rds")
    if (!file.exists(rds)) {
      warning("Skip ", sp, ": missing ", rds)
      next
    }
    res <- readRDS(rds)
    fig_dir <- file.path(CONFIG$output_dir, sp, "figures")
    cast_spatial_replot_hss_cate_heatmaps(
      pred            = res$predict,
      cate            = res$cate,
      fig_dir         = fig_dir,
      fig_dpi         = CONFIG$fig_dpi,
      var_labels      = CONFIG$var_labels,
      basemap         = CONFIG$plot_basemap,
      res_deg         = CONFIG$spatial_heatmap_res_deg,
      interp_method   = CONFIG$spatial_heatmap_interp_method,
      display_res_deg = CONFIG$spatial_heatmap_display_res,
      hss_model       = CONFIG$cate_hss_model,
      hss_threshold   = CONFIG$cate_hss_threshold,
      species_label   = sp,
      ovis_style      = FALSE
    )
    cat("  ", sp, "\n", sep = "")
  }
  cat("\nSpatial heatmap replot complete.\n\n")
} else if (isTRUE(CONFIG$only_shap_posthoc)) {
  future::plan(future::sequential)
  cfg_shap <- list(
    do_shap               = TRUE,
    response              = CONFIG$response,
    shap_nrounds          = CONFIG$shap_nrounds,
    shap_max_depth        = CONFIG$shap_max_depth,
    shap_eta              = CONFIG$shap_eta,
    shap_subsample        = CONFIG$shap_subsample,
    shap_colsample_bytree = CONFIG$shap_colsample_bytree,
    shap_test_fraction    = CONFIG$shap_test_fraction,
    shap_verbose          = CONFIG$shap_verbose,
    shap_plot_top_n       = CONFIG$shap_plot_top_n,
    shap_fastshap_nsim    = CONFIG$shap_fastshap_nsim,
    shap_max_explain_rows = CONFIG$shap_max_explain_rows
  )
  spn <- names(species_list)
  cat("\n[only_shap_posthoc] Writing SHAP figures (splits replayed via cast_prepare)...\n")
  for (i in seq_along(spn)) {
    sp <- spn[i]
    rds <- file.path(CONFIG$output_dir, sp, "cast_result.rds")
    if (!file.exists(rds)) {
      warning("Skip ", sp, ": missing ", rds)
      next
    }
    res <- readRDS(rds)
    sp_data <- species_list[[sp]]
    seed_i <- if (!is.null(CONFIG$seed)) CONFIG$seed + i else NULL
    split <- cast_prepare(
      sp_data,
      train_fraction = CONFIG$prepare_train_fraction,
      seed = seed_i,
      env_vars = CONFIG$prepare_env_vars,
      verbose = FALSE
    )
    fig_dir <- file.path(CONFIG$output_dir, sp, "figures")
    castSDM::save_cast_batch_shap_outputs(
      train_df = split$train,
      fit = res$fit,
      dag = res$dag,
      screen = res$screen,
      cfg = cfg_shap,
      fig_dir = fig_dir,
      fig_dpi = CONFIG$fig_dpi,
      seed_i = seed_i
    )
    cat("  ", sp, ": SHAP written -> ", fig_dir, "\n", sep = "")
  }
  cat("\nPost-hoc SHAP complete (also see shap_panel_2x3.png per species).\n\n")
} else {

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
  dev_package_root = if (!is.na(pkg_root)) pkg_root else NULL,
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
  bg = "transparent", limitsize = FALSE
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

} # end else (full cast_batch path; skipped when shap-only or heatmap-only)
