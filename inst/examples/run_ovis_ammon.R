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
# R=50 加速测试; 正式分析用 R=100+
dag <- cast_dag(split$train, R = 50, seed = 42)
print(dag)

# ── Step 3: ATE 因果效应估计 ─────────────────────────────────────────────────
ate <- cast_ate(split$train, K = 2, num_trees = 200, seed = 42)
print(ate)

# ── Step 4: 自适应变量筛选 ───────────────────────────────────────────────────
screen <- cast_screen(dag, ate, split$train, seed = 42)
print(screen)

# ── Step 5: 因果角色分配 ─────────────────────────────────────────────────────
roles <- cast_roles(screen, dag)
print(roles)

# ── Step 6: 模型训练 ─────────────────────────────────────────────────────────
# 仅跑 RF 快速验证 (不需要 torch)
fit_rf <- cast_fit(
  split$train,
  screen = screen, dag = dag, ate = ate,
  models = c("rf"),
  seed = 42
)
print(fit_rf)

# 如果已安装 torch, 取消下面注释跑完整 CAST:
fit_full <- cast_fit(
  split$train,
  screen = screen, dag = dag, ate = ate,
  models = c("cast", "rf", "maxent", "brt"),
  n_runs = 1, n_epochs = 50, seed = 42  # 快速测试参数
)

# ── Step 7: 模型评估 ─────────────────────────────────────────────────────────
eval_result <- cast_evaluate(fit_rf, split$test)
print(eval_result)

# ── Step 8: 空间预测 ─────────────────────────────────────────────────────────
pred <- cast_predict(fit_rf, china_env_grid)
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
}


# ══════════════════════════════════════════════════════════════════════════════
# 方式 B: 一键 Pipeline (等价于上面全部步骤)
# ══════════════════════════════════════════════════════════════════════════════
# 注意: 这会跑完整流程, 包含 DAG + ATE + Screening + 模型训练, 需要几分钟
#
# result <- cast(
#   species_data = ovis_ammon,
#   env_data = china_env_grid,
#   models = c("rf", "maxent", "brt"),  # 不含 cast 以避免 torch 依赖
#   n_bootstrap = 50,
#   seed = 42
# )
# summary(result)
# plot(result, var_labels = var_labels)
