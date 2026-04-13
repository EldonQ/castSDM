# ==============================================================================
# castSDM 端到端测试：盘羊 (Ovis ammon)
#
# 在 RStudio 中运行方式：
#   1. 用 RStudio 打开 castSDM.Rproj
#   2. Ctrl+Shift+L (devtools::load_all()) 加载包
#   3. 逐段选中本脚本运行 (Ctrl+Enter)
#
# 或者已安装包后直接：
#   library(castSDM)
#   source(system.file("examples/run_ovis_ammon.R", package = "castSDM"))
# ==============================================================================
#
# # 方式一：开发模式 (推荐)
# devtools::load_all("E:/Package/cast")
#
# # 方式二：重新安装后加载
# devtools::install("E:/Package/cast")
# library(castSDM)

# 开发模式 (在 castSDM 项目目录下) 用 load_all(), 否则用 library()
if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1))) {
  devtools::load_all()
} else {
  library(castSDM)
}

# ── Step 0: 加载示例数据 ─────────────────────────────────────────────────────
data(ovis_ammon)
data(china_env_grid)

cat("盘羊数据:", nrow(ovis_ammon), "行,", ncol(ovis_ammon), "列\n")
cat("存在点:", sum(ovis_ammon$presence == 1),
    "| 缺失点:", sum(ovis_ammon$presence == 0), "\n")

# 变量中文名映射 (可选, 用于美化绘图标签)
var_labels <- c(
  bio02 = "昼夜温差", bio15 = "降水季节性", bio19 = "最冷季降水",
  elevation = "海拔", aridityindexthornthwaite = "干旱指数",
  maxtempcoldest = "最冷月最高温", tri = "地形崎岖度",
  topowet = "地形湿度指数", nontree = "非乔木植被",
  etccdi_cwd = "连续湿日", landcover_igbp = "土地覆盖"
)


# ══════════════════════════════════════════════════════════════════════════════
# 方式 A: 模块化逐步运行 (推荐初次使用, 便于理解每一步)
# ══════════════════════════════════════════════════════════════════════════════

# ── Step 1: 数据准备 ─────────────────────────────────────────────────────────
split <- cast_prepare(ovis_ammon, train_fraction = 0.7, seed = 42)
cat("训练集:", nrow(split$train), "| 测试集:", nrow(split$test), "\n")
cat("环境变量:", paste(split$env_vars, collapse = ", "), "\n")

# ── Step 2: DAG 因果结构学习 ─────────────────────────────────────────────────
# structure_method 可选:
#   "bootstrap_hc"(默认), "pc", "fci", "bidag_bge"(需 BiDAG), "notears_linear"(需 torch)
# 下面展示 cast_dag() 主要参数; R 仅对 bootstrap_hc 有效。
dag <- cast_dag(
  split$train,
  response              = "presence",
  env_vars              = NULL,
  R                     = 50L,
  algorithm             = "hc",
  score                 = "bic-g",
  strength_threshold      = 0.7,
  direction_threshold     = 0.6,
  max_rows              = 8000L,
  seed                  = 42L,
  verbose               = TRUE,
  structure_method      = "bootstrap_hc",
  pc_alpha              = 0.05,
  fci_alpha             = 0.05,
  bidag_algorithm       = "order",
  bidag_iterations      = NULL,
  notears_lambda        = 0.03,
  notears_max_iter      = 2000L,
  notears_lr            = 0.02,
  notears_tol           = 1e-3,
  notears_rho_init      = 0.1,
  notears_alpha_mult    = 1.01
)
print(dag)

# ── Step 3: ATE 因果效应估计 ─────────────────────────────────────────────────
ate <- cast_ate(
  split$train,
  variables      = NULL,
  K              = 5L,
  alpha          = 0.05,
  num_trees      = 300L,
  quantile_cuts  = c(0.25, 0.5, 0.75),
  bonferroni     = TRUE,
  parallel       = TRUE,
  seed           = 42L,
  verbose        = TRUE
)
print(ate)

# ── Step 4: 自适应变量筛选 ───────────────────────────────────────────────────
screen <- cast_screen(
  dag, ate, split$train,
  response     = "presence",
  min_vars      = 5L,
  min_fraction  = 0.5,
  num_trees     = 500L,
  seed          = 42L,
  verbose       = TRUE
)
print(screen)

# ── Step 5: 因果角色分配 ─────────────────────────────────────────────────────
roles <- cast_roles(screen, dag)
print(roles)

# ── Step 6: 模型训练 ─────────────────────────────────────────────────────────
# 如果已安装 torch, 取消下面注释跑完整 CAST:
fit_full <- cast_fit(
  split$train,
  screen = screen, dag = dag, ate = ate,
  models = c("cast", "rf", "maxent", "brt"),
  n_runs = 1, n_epochs = 50, seed = 42,  # 快速测试参数
  tune_grid = FALSE
)

# ── Step 7: 模型评估 ─────────────────────────────────────────────────────────
eval_result <- cast_evaluate(fit_full, split$test)
print(eval_result)

# ── Step 8: 空间预测 ─────────────────────────────────────────────────────────
pred <- cast_predict(fit_full, china_env_grid)
print(pred)

# ── Step 9: 绘图 ─────────────────────────────────────────────────────────────
# 以下绘图需要 ggplot2, ggraph, igraph, sf 等
# install.packages(c("ggplot2", "ggraph", "igraph", "sf", "patchwork"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  # DAG 网络图
  plot(dag, roles = roles, screen = screen, var_labels = var_labels)

  # ATE 森林图
  plot(ate, var_labels = var_labels)

  # 筛选分数分解图
  plot(screen, var_labels = var_labels)

  # 评估对比图
  plot(eval_result)

  # 空间栖息地适宜性地图 (需要 sf)
  if (requireNamespace("sf", quietly = TRUE)) {
    plot(pred, model = "rf", basemap = "china")
  }
    # 模型间空间一致性热力图 (需要 >= 2 个模型)
  if (length(pred$models) >= 2 &&
      requireNamespace("patchwork", quietly = TRUE)) {
    consistency <- cast_consistency(pred)
    print(consistency)
    p_cons <- plot(
      consistency,
      species = "Ovis_ammon",
      font_family = "sans",
      use_bold = TRUE,
      font_base = 11L,
      font_main_title = 14L,
      font_panel_title = 11L,
      font_axis = 8L,
      font_cell_value = 3.5,
      value_decimals = 3L,
      tile_linewidth = 0.5,
      text_white_above = 0.7
    )
    print(p_cons)
    # ggsave("ovis_consistency.png", p_cons, width = 16, height = 5.5, dpi = 300)
  }

  # 空间 CATE + 按 HSS 截断出图 (需 grf; HSS 模型默认 cast, 阈值 0.1)
  if (requireNamespace("grf", quietly = TRUE)) {
    cate <- cast_cate(
      split$train,
      variables    = NULL,
      ate            = ate,
      screen         = screen,
      response       = "presence",
      top_n          = 2L,
      n_trees        = 300L,
      predict_data   = china_env_grid,
      seed           = 42L,
      verbose        = FALSE
    )
    if (!is.null(cate)) {
      for (v in cate$variables) {
        pc <- plot(
          cate,
          variable         = v,
          species          = "Ovis_ammon",
          basemap          = "china",
          var_labels       = var_labels,
          point_size       = 0.45,
          legend_position  = "bottom",
          hss_predict      = pred,
          hss_model        = "cast",
          hss_threshold    = 0.1
        )
        print(pc)
      }
    }
  }

  # XGBoost + SHAP 可解释性 (需 xgboost; 与 SDM 共用 presence + 环境变量)
  if (requireNamespace("xgboost", quietly = TRUE)) {
    sh <- cast_shap_xgb(
      split$train,
      response           = "presence",
      env_vars           = NULL,
      screen             = screen,
      nrounds            = 200L,
      max_depth          = 6L,
      eta                = 0.05,
      subsample          = 0.8,
      colsample_bytree   = 0.8,
      test_fraction      = 0.2,
      seed               = 42L,
      verbose            = FALSE
    )
    print(plot(sh, type = "interaction_network", top_n = 15L))
    print(plot(sh, type = "waterfall", top_n = 15L))
    # cast_shap_write_csv(sh, file.path(getwd(), "shap_export_ovis"))
  }
}




# ══════════════════════════════════════════════════════════════════════════════
# 方式 B: 一键 Pipeline (等价于上面全部步骤)
# ══════════════════════════════════════════════════════════════════════════════
# 注意: 这会跑完整流程, 包含 DAG + ATE + Screening + 模型训练, 需要几分钟
#
# result <- cast(
#   species_data = ovis_ammon,
#   env_data = china_env_grid,
#   models = c("rf", "maxent", "brt"),
#   train_fraction = 0.7,
#   n_bootstrap = 50L,
#   dag_structure_method = "bootstrap_hc",
#   dag_pc_alpha = 0.05, dag_fci_alpha = 0.05,
#   dag_bidag_algorithm = "order", dag_bidag_iterations = NULL,
#   dag_notears_lambda = 0.03, dag_notears_max_iter = 2000L,
#   dag_notears_lr = 0.02, dag_notears_tol = 1e-3,
#   dag_notears_rho_init = 0.1, dag_notears_alpha_mult = 1.01,
#   strength_threshold = 0.7, direction_threshold = 0.6,
#   ate_folds = 2L, ate_alpha = 0.05, screen_min_vars = 5L,
#   do_cv = TRUE, cv_k = 5L, cv_block_method = "grid",
#   do_cate = TRUE, cate_top_n = 3L,
#   seed = 42L, verbose = TRUE
# )
# summary(result)
# plot(result, var_labels = var_labels)


# ══════════════════════════════════════════════════════════════════════════════
# 方式 C: 不同 DAG structure_method 各跑一遍并保存图 (与当次 ATE/screen/roles 衔接)
# ══════════════════════════════════════════════════════════════════════════════
# 将 RUN_DAG_METHOD_SHOWCASE 改为 TRUE 可执行 (除 bootstrap 外部分依赖 BiDAG/torch)
RUN_DAG_METHOD_SHOWCASE <- FALSE

if (isTRUE(RUN_DAG_METHOD_SHOWCASE) && requireNamespace("ggplot2", quietly = TRUE)) {
  out_d <- file.path(getwd(), "figures_ovis_dag_by_method")
  dir.create(out_d, FALSE, TRUE)
  methods <- c("bootstrap_hc", "pc", "fci", "bidag_bge", "notears_linear")
  for (sm in methods) {
    dm <- tryCatch(
      cast_dag(
        split$train,
        R = 40L,
        seed = 42L,
        verbose = FALSE,
        structure_method = sm
      ),
      error = function(e) {
        message("cast_dag skipped (", sm, "): ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(dm)) next
    am <- cast_ate(split$train, variables = dm$nodes, K = 5L,
                   num_trees = 200L, parallel = FALSE, seed = 42L,
                   verbose = FALSE)
    smc <- cast_screen(dm, am, split$train, num_trees = 300L,
                       seed = 42L, verbose = FALSE)
    rm <- cast_roles(smc, dm)
    p <- plot(dm, roles = rm, screen = smc, var_labels = var_labels,
              species = "Ovis_ammon")
    ggplot2::ggsave(file.path(out_d, paste0("dag_", sm, ".png")),
                     p, width = 12, height = 9, dpi = 200)
  }
  message("Saved under: ", normalizePath(out_d, winslash = "/", mustWork = FALSE))
}
