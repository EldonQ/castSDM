# ==============================================================================
# castSDM × disdat 瑞士植物 (SWI) 完整教程
#
# 目标读者: SDM 建模初学者
# 数据来源: R 包 disdat (Elith et al. 2020) — 标准化 SDM 基准数据集
# 区域: 瑞士 (Switzerland, SWI), 30 种高山植物
#
# 本脚本将演示：
#   1. 安装/加载所需 R 包
#   2. 读取 SWI 数据并探索格式
#   3. 选择物种, 构造 presence/absence 训练数据
#   4. 瑞士坐标系 → WGS84 经纬度转换
#   5. cast_prepare / cast_dag / cast_ate / cast_screen / cast_roles
#   6. cast_fit 训练模型 (RF 快速版 + 完整版含 CAST 神经网络)
#   7. cast_evaluate 评估
#   8. cast_predict 空间预测 → 世界地图
#   9. cast_cate 空间异质性因果效应 → CATE 地图
#
# 运行方式 (RStudio):
#   1. 打开 castSDM.Rproj
#   2. Ctrl+Shift+L  → devtools::load_all()  (开发模式)
#   或  library(castSDM)                      (已安装包模式)
#   3. 逐段选中, Ctrl+Enter 运行
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# 第 0 步: 安装所需包
# ══════════════════════════════════════════════════════════════════════════════

# 如果以下包尚未安装, 取消注释后运行一次：
# install.packages(c(
#   "disdat",       # SWI 标准数据集
#   "sf",           # 坐标系转换 + 地图
#   "ggplot2",      # 绘图基础
#   "ggraph",       # DAG 网络图
#   "igraph",       # 图结构
#   "patchwork",    # 多图拼接
#   "bnlearn",      # 贝叶斯网络 (DAG 学习)
#   "ranger",       # 随机森林
#   "gbm",          # 梯度提升树 BRT
#   "maxnet",       # MaxEnt
#   "grf",          # 因果随机森林 (CATE)
#   "pROC",         # AUC 计算
#   "scales",       # 颜色辅助
#   "car"           # VIF 共线性检验
# ))

# torch 安装 (CI-MLP 神经网络需要, 较大约 500MB):
# install.packages("torch")
# torch::install_torch()


# ══════════════════════════════════════════════════════════════════════════════
# 第 1 步: 加载包
# ══════════════════════════════════════════════════════════════════════════════

# 加载 castSDM (开发模式用 load_all, 已安装用 library)
if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1), fixed = TRUE)) {
  devtools::load_all()  # RStudio 项目开发模式
} else {
  library(castSDM)
}

library(disdat)   # SWI 标准数据集
library(sf)       # 坐标系转换


# ══════════════════════════════════════════════════════════════════════════════
# 第 2 步: 读取并探索 SWI 数据
# ══════════════════════════════════════════════════════════════════════════════

cat("=== SWI 数据集概览 ===\n")

# 存在点 (Presence-Only, PO): 物种出现记录 + 环境变量
po <- disdat::disPo("SWI")
cat("PO 数据:", nrow(po), "行,", ncol(po), "列\n")
cat("列名:", paste(names(po)[1:8], collapse = ", "), "...\n")
cat("物种数:", length(unique(po$spid)), "\n\n")

# 背景点 (Background, BG): 伪缺失点 + 环境变量
bg <- disdat::disBg("SWI")
cat("BG 数据:", nrow(bg), "行,", ncol(bg), "列\n\n")

# 存在-缺失测试集 (Presence-Absence, PA): 用于模型评估
pa <- disdat::disPa("SWI")
cat("PA 测试集:", nrow(pa), "行,", ncol(pa), "列\n")
cat("列名 (前10):", paste(names(pa)[1:10], collapse = ", "), "\n\n")

# 环境预测网格 (用于生成空间预测图)
env_grid <- disdat::disEnv("SWI")
cat("预测网格:", nrow(env_grid), "格点\n\n")

# 查看坐标范围 — SWI 用 EPSG:21781 + 常量平移, 不是标准 WGS84 经纬度
# 用 disCRS("SWI") 获取含平移信息的正确 CRS, 不要自己猜 EPSG
cat("坐标范围: x =", round(range(po$x, na.rm=TRUE), 0),
    "| y =", round(range(po$y, na.rm=TRUE), 0), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# 第 3 步: 选择目标物种
# ══════════════════════════════════════════════════════════════════════════════

# 查看各物种存在点数量, 选择记录最多的物种 (保证数据充足)
sp_counts <- sort(table(po$spid), decreasing = TRUE)
cat("存在点数量前5物种:\n")
print(sp_counts[1:5])

# 选择存在点最多的物种 (通常 > 50 条, 适合建模)
target_species <- names(sp_counts)[1]
cat("\n目标物种:", target_species, "\n")
cat("存在点数:", sp_counts[target_species], "\n")

# 如需手动指定物种, 取消下面注释并修改名称:
# target_species <- "Adenostyles_alliariae"


# ══════════════════════════════════════════════════════════════════════════════
# 第 4 步: 构造训练数据 (存在点 + 背景点)
# ══════════════════════════════════════════════════════════════════════════════

# 提取目标物种的存在点
po_sp <- po[po$spid == target_species, ]

# ── 获取官方环境变量列名 ──────────────────────────────────────────────────────
# disPredictors("SWI") 返回 disdat 为该数据集定义的 13 个官方预测变量名
# 直接用它: 比手动排除元数据列更可靠, 不会把 occ/siteid 等误算进去
env_cols <- disdat::disPredictors("SWI")
cat("官方环境变量 (", length(env_cols), "):",
    paste(env_cols, collapse = ", "), "\n\n")

# 存在点: presence = 1
pres_df <- data.frame(
  x        = po_sp$x,
  y        = po_sp$y,
  presence = 1L,
  po_sp[, env_cols, drop = FALSE],
  stringsAsFactors = FALSE
)

# 背景点: presence = 0
# 使用全部背景点 (或随机抽样与存在点1:10比例)
n_bg_sample <- min(nrow(bg), 10 * nrow(pres_df))
set.seed(42)
bg_idx <- sample(nrow(bg), n_bg_sample)
bg_sub <- bg[bg_idx, ]
abs_df <- data.frame(
  x        = bg_sub$x,
  y        = bg_sub$y,
  presence = 0L,
  bg_sub[, env_cols, drop = FALSE],
  stringsAsFactors = FALSE
)

# 合并存在点和背景点
train_raw <- rbind(pres_df, abs_df)
cat("合并后训练数据:", nrow(train_raw), "行\n")
cat("存在点:", sum(train_raw$presence == 1),
    "| 背景点:", sum(train_raw$presence == 0), "\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# 第 5 步: 坐标系转换 → WGS84 (lon/lat)
#
# castSDM 需要 lon (经度) 和 lat (纬度) 列, WGS84 坐标系 (EPSG:4326)
#
# disdat 文档说明: SWI 坐标为 EPSG:21781 (Swiss LV03, Bessel 椭球体),
# 且有常量平移 (constant shift applied). 因此不能直接用 EPSG:21781 转换,
# 必须用 disCRS("SWI") 获取含平移信息的正确 CRS proj4 字符串.
# ══════════════════════════════════════════════════════════════════════════════

# 获取 SWI 数据集的正确 CRS (含平移)
swi_crs <- disdat::disCRS("SWI")
cat("SWI 坐标参考系:\n", swi_crs, "\n\n")

# 通用坐标转换函数: 用 proj4 字符串转到 WGS84
convert_disdat_to_wgs84 <- function(df, proj4_crs) {
  pts_sf    <- sf::st_as_sf(
    df, coords = c("x", "y"), crs = sf::st_crs(proj4_crs), remove = FALSE
  )
  pts_wgs84 <- sf::st_transform(pts_sf, crs = 4326)
  coords    <- sf::st_coordinates(pts_wgs84)
  df$lon    <- coords[, 1]
  df$lat    <- coords[, 2]
  df
}

# 转换训练数据坐标
train_raw <- convert_disdat_to_wgs84(train_raw, swi_crs)
cat("训练数据经度范围:", round(range(train_raw$lon), 2), "\n")
cat("训练数据纬度范围:", round(range(train_raw$lat), 2), "\n")
# 预期: lon ≈ 5.9~10.5°E, lat ≈ 45.8~47.8°N (瑞士范围)

# 转换预测网格坐标
env_grid_df <- convert_disdat_to_wgs84(as.data.frame(env_grid), swi_crs)
cat("预测网格格点数:", nrow(env_grid_df), "\n\n")

# 构造 castSDM 标准格式: lon, lat, presence, [env vars]
species_data <- data.frame(
  lon      = train_raw$lon,
  lat      = train_raw$lat,
  presence = train_raw$presence,
  train_raw[, env_cols, drop = FALSE],
  stringsAsFactors = FALSE
)

# 构造预测网格: lon, lat, [env vars]
pred_grid <- data.frame(
  lon = env_grid_df$lon,
  lat = env_grid_df$lat,
  env_grid_df[, env_cols, drop = FALSE],
  stringsAsFactors = FALSE
)

cat("最终训练数据:", nrow(species_data), "行,",
    ncol(species_data) - 3, "个环境变量\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# 第 6 步: 数据分割 (70% 训练 / 30% 测试)
# ══════════════════════════════════════════════════════════════════════════════

# env_vars 显式传入, 确保与 disPredictors("SWI") 完全一致
# (即使 castSDM 自动检测结果相同, 显式指定更安全可靠)
split <- cast_prepare(
  species_data,
  train_fraction = 0.7,
  seed           = 42,
  env_vars       = env_cols   # 来自 disPredictors("SWI")
)
cat("训练集:", nrow(split$train),
    "(存在:", sum(split$train$presence), ")\n")
cat("测试集:", nrow(split$test),
    "(存在:", sum(split$test$presence), ")\n")
cat("环境变量:", paste(split$env_vars, collapse = ", "), "\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# 第 7 步: 因果 DAG 学习
# cast_dag 使用 bnlearn (Hill-climbing + bootstrap) 学习变量间因果结构
# R = 引导重采样次数: 正式分析用 100+, 本教程用 50 加速
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 7: 学习因果 DAG ===\n")
dag <- cast_dag(
  split$train,
  R           = 50,    # 引导重采样次数 (正式分析建议 100~200)
  seed        = 42,
  verbose     = TRUE
)
print(dag)
# 输出: 边的数量, 各变量在 DAG 中的因果角色 (Root/Mediator/Terminal)


# ══════════════════════════════════════════════════════════════════════════════
# 第 8 步: ATE 因果效应估计 (Double Machine Learning)
# DML 使用交叉拟合 RF 分离混杂, 估计每个环境变量对物种出现的平均因果效应
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 8: 估计 ATE (双机器学习) ===\n")
ate <- cast_ate(
  split$train,
  K         = 2,    # K 折交叉拟合 (正式分析用 5)
  num_trees = 200,  # RF 树数
  seed      = 42,
  verbose   = TRUE
)
print(ate)
# 输出: 每个变量的 ATE 系数、标准误、p 值, 是否显著


# ══════════════════════════════════════════════════════════════════════════════
# 第 9 步: 自适应变量筛选 (CAST Screening)
# 综合 DAG 拓扑分、ATE 因果分、RF 重要性三维度, 选出核心变量
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 9: 自适应变量筛选 ===\n")
screen <- cast_screen(
  dag, ate, split$train,
  seed = 42,
  verbose = TRUE
)
print(screen)
cat("选出变量:", paste(screen$selected, collapse = ", "), "\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# 第 10 步: 因果角色分配
# 根据 DAG 结构为每个变量分配角色: Root (根节点), Mediator (中介), Terminal (终端)
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 10: 因果角色分配 ===\n")
roles <- cast_roles(screen, dag)
print(roles)


# ══════════════════════════════════════════════════════════════════════════════
# 第 11 步: 模型训练
#
# 方案 A: 仅 RF (快速, 无需 torch, 适合初学者验证)
# 方案 B: 完整版 (含 CAST 神经网络, 需要 torch)
# tune_grid = TRUE 会在内部验证集上搜索 hidden_size × dropout × lr,
#   找到最佳超参后再做 n_runs 集成 (大约增加 3× 训练时间)
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 11: 训练模型 ===\n")

## ── 方案 A: 仅随机森林 (推荐初次运行) ───────────────────────────────────────
fit <- cast_fit(
  split$train,
  screen   = screen,
  dag      = dag,
  ate      = ate,
  models   = c("rf"),          # 仅 RF, 不需要 torch
  rf_ntree = 300,
  seed     = 42,
  verbose  = TRUE
)
print(fit)

## ── 方案 B: RF + BRT + MaxEnt + CAST 神经网络 (需要 torch) ──────────────────
## 取消注释后运行; n_epochs=50/n_runs=1 是快速测试参数, 正式分析用默认值
## tune_grid=TRUE 开启超参网格搜索 (正式分析推荐)
# fit <- cast_fit(
#   split$train,
#   screen      = screen,
#   dag         = dag,
#   ate         = ate,
#   models      = c("cast", "rf", "brt", "maxent"),
#   n_epochs    = 50,        # 正式: 200
#   n_runs      = 1,         # 正式: 3
#   rf_ntree    = 300,
#   brt_n_trees = 500,
#   tune_grid   = FALSE,     # TRUE = 开启超参搜索 (需更多时间)
#   seed        = 42,
#   verbose     = TRUE
# )


# ══════════════════════════════════════════════════════════════════════════════
# 第 12 步: 空间交叉验证 (Spatial k-fold CV)
#
# cast_cv() 将数据按空间格网分成 k 个地理块，每次用 k-1 个块训练、
# 第 k 块测试，输出：
#   - 空间独立的 AUC / TSS / CBI / SEDI / Kappa / PRAUC
#   - 每折最优阈值 → 用于后续 HSS 地图二值化
#
# k 建议值：
#   k=5  (默认) ── 适合大多数数据集
#   k=3          ── 存在点 < 100 时
#   k=10         ── 数据量充足且需要更稳定的估计
#
# block_method:
#   "grid"    (默认) ── 均匀空间格网分块，无需额外包
#   "cluster"        ── k-means 聚类分块，更自适应分布
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 12: 空间交叉验证 ===\n")

## ── 方案 A: RF 空间CV (快速, 无需 torch) ─────────────────────────────────────
cv_result <- cast_cv(
  species_data,        # 用全量数据做CV (比只用 split$train 更充分)
  screen       = screen,
  dag          = dag,
  ate          = ate,
  k            = 5L,   # 折数: 用户可修改为 3 / 5 / 10
  models       = c("rf"),
  block_method = "grid",
  rf_ntree     = 300L,
  seed         = 42,
  verbose      = TRUE
)
print(cv_result)

## ── 方案 B: 多模型空间CV (含神经网络, 需要 torch) ────────────────────────────
# cv_result <- cast_cv(
#   species_data,
#   screen       = screen,
#   dag          = dag,
#   ate          = ate,
#   k            = 5L,
#   models       = c("cast", "rf", "brt", "maxent"),
#   block_method = "grid",
#   n_epochs     = 50L,   # 正式: 100
#   n_runs       = 1L,    # 正式: 2
#   rf_ntree     = 300L,
#   seed         = 42,
#   verbose      = TRUE
# )

## ── 可视化: 折叠地图 + 折间指标变化 ─────────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  # 仅指标折线图 (无需提供坐标)
  p_cv_metric <- plot(cv_result, metric = "auc")
  print(p_cv_metric)

  # 含地理折叠分区图 (需传入坐标)
  p_cv_map <- plot(
    cv_result,
    lon     = species_data$lon,
    lat     = species_data$lat,
    metric  = "auc",
    basemap = "world"
  )
  print(p_cv_map)
}


# ══════════════════════════════════════════════════════════════════════════════
# 第 13 步: 保留集评估 (Hold-out AUC/TSS/CBI/SEDI/Kappa/PRAUC)
#
# cast_evaluate() 在随机保留测试集上计算全套指标。
# 配合 cast_cv() 使用：
#   - cast_cv()        → 空间泛化能力 (更可靠, 用于论文主指标)
#   - cast_evaluate()  → 快速参考指标 (补充对比)
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 13: 保留集评估 ===\n")
eval_result <- cast_evaluate(fit, split$test)
print(eval_result)
# 输出: AUC, TSS, CBI, SEDI, Kappa, PRAUC (6 个指标)


# ══════════════════════════════════════════════════════════════════════════════
# 第 14 步: 空间预测 → 栖息地适宜性地图 (HSS)
# 用 pred_grid (瑞士全域环境网格) 生成每个格点的栖息地适宜性评分
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 14: 空间预测 ===\n")
pred <- cast_predict(fit, pred_grid)
print(pred)
cat("HSS 列:", names(pred$predictions), "\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# 第 15 步: CATE 空间异质性因果效应
# 用因果随机森林估计每个格点处每个关键变量对物种出现的条件平均因果效应
# 揭示"哪里的气温变化对物种影响最大"
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 15: CATE 空间估计 ===\n")
cate_result <- cast_cate(
  data         = split$train,     # 用训练数据拟合因果森林
  ate          = ate,             # ATE 结果 (用于选择显著变量)
  screen       = screen,          # 筛选结果 (用于选择核心变量)
  predict_data = pred_grid,       # 在整个预测网格上估计 CATE
  top_n        = 3L,              # 选 top 3 显著变量
  n_trees      = 500L,            # 因果森林树数
  seed         = 42,
  verbose      = TRUE
)
print(cate_result)
cat("CATE 估计变量:", paste(cate_result$variables, collapse = ", "), "\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# 第 16 步: 可视化
# 需要: ggplot2, ggraph, igraph, sf, patchwork, scales
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Step 16: 绘图 ===\n")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("请先安装 ggplot2: install.packages('ggplot2')")
} else {

  ## ── 16.1 DAG 因果网络图 ────────────────────────────────────────────────────
  # 节点颜色: 因果角色 (Root/Mediator/Terminal)
  # 节点形状: 是否被 CAST 选中
  # 边粗细/透明度: bootstrap 强度
  p_dag <- plot(dag, roles = roles, screen = screen)
  print(p_dag)
  # ggsave("swi_dag.png", p_dag, width = 8, height = 6, dpi = 150)


  ## ── 16.2 ATE 森林图 ─────────────────────────────────────────────────────────
  # 红色: 显著 (p < 0.05), 灰色: 不显著
  # 横轴: ATE 系数 (正 = 促进物种出现, 负 = 抑制)
  p_ate <- plot(ate)
  print(p_ate)


  ## ── 16.3 CAST 变量筛选评分分解图 ──────────────────────────────────────────
  # 三色堆叠条形图: DAG 拓扑分 + ATE 因果分 + RF 重要性分
  # 深色: 被选中变量, 浅色: 未被选中
  p_screen <- plot(screen)
  print(p_screen)


  ## ── 16.4 模型性能对比图 (6 指标 facet) ────────────────────────────────────
  # 展示 AUC / TSS / CBI / SEDI / Kappa / PRAUC
  p_eval <- plot(eval_result)
  print(p_eval)

  ## ── 16.4b 空间CV折叠地图 + 折间AUC变化 ────────────────────────────────────
  p_cv <- plot(
    cv_result,
    lon     = species_data$lon,
    lat     = species_data$lat,
    metric  = "auc",
    basemap = "world"
  )
  print(p_cv)


  ## ── 16.5 栖息地适宜性空间地图 (HSS) ──────────────────────────────────────
  # basemap = "world" 自动加载世界底图 (需要包内 basemap/countries.shp)
  # 颜色: turbo 色板, 0 (不适宜) → 1 (高度适宜)
  if (requireNamespace("sf", quietly = TRUE)) {
    p_hss <- plot(
      pred,
      model   = "rf",         # 若训练了 "cast", 改为 model = "cast"
      basemap = "world",
      title   = paste("Habitat Suitability:", gsub("_", " ", target_species))
    )
    print(p_hss)
    # ggsave("swi_hss_rf.png", p_hss, width = 8, height = 5, dpi = 150)
  }
    ## ── 16.5b 模型间空间一致性热力图 ──────────────────────────────────────────
  if (length(pred$models) >= 2 &&
      requireNamespace("patchwork", quietly = TRUE)) {
    consistency <- cast_consistency(pred)
    print(consistency)
    p_cons <- plot(consistency, species = target_species)
    print(p_cons)
    # ggsave("swi_consistency.png", p_cons, width = 16, height = 5.5, dpi = 300)
  }


  ## ── 16.6 CATE 空间异质性因果效应地图 ─────────────────────────────────────
  # 蓝色: 负效应 (该变量增加 → 物种出现概率下降)
  # 红色: 正效应 (该变量增加 → 物种出现概率上升)
  # 揭示"环境驱动力的空间异质性"
  if (requireNamespace("sf", quietly = TRUE) &&
      length(cate_result$variables) > 0) {

    # 绘制第一个变量的 CATE 地图
    p_cate1 <- plot(
      cate_result,
      variable        = cate_result$variables[1],
      species         = target_species,
      basemap         = "world",
      legend_position = "bottom"
    )
    print(p_cate1)
    # ggsave("swi_cate_1.png", p_cate1, width = 8, height = 5, dpi = 150)

    # 若有 3 个变量, 用 patchwork 拼成三联图
    if (requireNamespace("patchwork", quietly = TRUE) &&
        length(cate_result$variables) >= 3) {
      p_cate2 <- plot(cate_result, variable = cate_result$variables[2],
                      species = target_species, basemap = "world",
                      legend_position = "bottom")
      p_cate3 <- plot(cate_result, variable = cate_result$variables[3],
                      species = target_species, basemap = "world",
                      legend_position = "bottom")
      p_cate_all <- p_cate1 / p_cate2 / p_cate3 +
        patchwork::plot_annotation(
          title    = paste("Spatial CATE:", gsub("_", " ", target_species)),
          subtitle = "Conditional Average Treatment Effects (Causal Forest)",
          theme    = ggplot2::theme(
            plot.title    = ggplot2::element_text(face="bold", hjust=0.5),
            plot.subtitle = ggplot2::element_text(hjust=0.5, color="grey40")
          )
        )
      print(p_cate_all)
      # ggsave("swi_cate_all.png", p_cate_all, width=8, height=14, dpi=150)
    }
  }

}  # end ggplot2 block


# ══════════════════════════════════════════════════════════════════════════════
# 第 17 步: 快速汇总
# ══════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("=", 60), "\n")
cat("castSDM SWI 完整流程运行完毕!\n")
cat(strrep("=", 60), "\n\n")
cat("物种:", gsub("_", " ", target_species), "\n")
cat("训练数据:", nrow(split$train), "条 (存在:", sum(split$train$presence), ")\n")
cat("因果 DAG:", nrow(dag$edges), "条边\n")
sig_n <- sum(ate$estimates$significant, na.rm = TRUE)
cat("显著 ATE 变量:", sig_n, "/", nrow(ate$estimates), "\n")
cat("CAST 选出变量:", paste(screen$selected, collapse=", "), "\n")
cat("CATE 变量:", paste(cate_result$variables, collapse=", "), "\n")
cat("预测格点:", nrow(pred_grid), "\n\n")

# 空间CV结果摘要
if (exists("cv_result") && nrow(cv_result$metrics) > 0) {
  cat("── 空间CV结果 (", cv_result$k, "折, ", cv_result$block_method, ") ──\n", sep="")
  for (i in seq_len(nrow(cv_result$metrics))) {
    r <- cv_result$metrics[i, ]
    cat(sprintf(
      "  %s: AUC=%.3f TSS=%.3f CBI=%.3f SEDI=%.3f Kappa=%.3f PRAUC=%.3f\n",
      r$model, r$auc_mean, r$tss_mean,
      r$cbi_mean, r$sedi_mean, r$kappa_mean, r$prauc_mean
    ))
  }
  cat("  最优阈值:", paste(
    names(cv_result$thresholds),
    round(cv_result$thresholds, 3), sep="=", collapse=", "
  ), "\n\n")
}

cat("主要输出对象:\n")
cat("  dag        → 因果 DAG (cast_dag)\n")
cat("  ate        → 因果效应 (cast_ate)\n")
cat("  screen     → 变量筛选 (cast_screen)\n")
cat("  fit        → 拟合模型 (cast_fit)\n")
cat("  cv_result  → 空间交叉验证 (cast_cv)  ← 新增\n")
cat("  eval_result → 保留集评估 (cast_eval)\n")
cat("  pred       → 空间预测 (cast_predict)\n")
cat("  cate_result → CATE 估计 (cast_cate)\n")

