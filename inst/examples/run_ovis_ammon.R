# ==============================================================================
# castSDM 端到端测试：盘羊 (Ovis ammon)
#
# 空间 HSS / CATE：包内 plot.cast_predict / plot.cast_cate 对任意 lon/lat 格点
# 均用 geom_point 绘制（Eco-ISEA3H 只是你的数据形态，非特例）。文末「空间热图
# 重绘」用 terra 栅格化 + focal 得到类 fig7_spatial_prediction_maps.py 的连续热图。
#
# 在 RStudio 中运行方式：
#   1. 用 RStudio 打开 castSDM.Rproj
#   2. Ctrl+Shift+L (devtools::load_all()) 加载包
#   3. 逐段选中本脚本运行 (Ctrl+Enter)
#
# 或者已安装包后直接：
#   library(castSDM)
#   source(system.file("examples/run_ovis_ammon.R", package = "castSDM"))
#
# DAG 四种 structure_method 全量自检默认开启 (见 RUN_DAG_SELFTEST_ALL_METHODS)。
# 依赖 (自检 + 后续 ATE 必用 ranger): 建议一次安装
#   install.packages(c("bnlearn", "BiDAG", "torch", "ranger"))
# 说明: bidag_bge 需 BiDAG;
# notears_linear 需 torch。仅想跑主流程、不做 DAG 五法自检可设
#   RUN_DAG_SELFTEST_ALL_METHODS <- FALSE
# ==============================================================================
#
# # 方式一：开发模式 (推荐)
# devtools::load_all("E:/Package/cast")
#
# # 方式二：重新安装后加载
# devtools::install("E:/Package/cast")
# library(castSDM)

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

# 开发模式：包根 load_all；与 run_multi_species.R 一致。并行 future 子进程需 CASTSDM_ROOT（见该脚本说明）。
pkg_root <- .cast_find_package_root()
if (!is.na(pkg_root)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::load_all(pkg_root)
  Sys.setenv(CASTSDM_ROOT = normalizePath(pkg_root, winslash = "/", mustWork = FALSE))
} else {
  Sys.unsetenv("CASTSDM_ROOT")
  library(castSDM)
}

# 出图统一落盘：当前工作目录下子文件夹，600 dpi（需 ggplot2::ggsave）
OVIS_FIG_DIR <- file.path(getwd(), "cast_ovis_ammon_figures")
OVIS_FIG_DPI <- 1200L
dir.create(OVIS_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ovis_save_plot <- function(plot_obj, filename, width, height) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }
  fp <- file.path(OVIS_FIG_DIR, filename)
  ggplot2::ggsave(
    filename = fp,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = OVIS_FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )
  message("Saved figure: ", normalizePath(fp, winslash = "/", mustWork = FALSE))
  invisible(fp)
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

# ==============================================================================
# cast_dag()：四种 structure_method 的参数与依赖（正式跑前必读）
# ==============================================================================
# 入口：cast_dag(data, response = "presence", env_vars = NULL, ...,
#               structure_method = "<下列之一>", ...)
#
# (1) bootstrap_hc — 自助法 + 打分搜索（默认 CAST 流程）
#     依赖：bnlearn
#     生效参数：R, algorithm, score, strength_threshold, direction_threshold,
#               max_rows, seed, verbose, env_vars, response
#     algorithm：传给 bnlearn::boot.strength() 的学习器，如 "hc","tabu",
#                 "mmhc","pc.stable"（注意：选约束型 PC 本体应设
#                 structure_method="pc"，不要把 algorithm 写成 "pc"）
#     忽略：pc_alpha/pc_test、bidag_*、notears_*（仍会传入，包内忽略）
#
# (2) pc — 约束型 PC（bnlearn::pc.stable）
#     依赖：bnlearn
#     生效参数：pc_alpha, pc_test（如连续高斯常用 "zf"）, max_rows, seed,
#               verbose, env_vars, response
#     忽略：R, algorithm, score, strength_threshold, direction_threshold,
#           bidag_*, notears_*

# (3) bidag_bge — BiDAG + BGe 分数
#     依赖：BiDAG（及数据矩阵）
#     生效参数：bidag_algorithm ("order"|"orderIter"), bidag_iterations,
#               max_rows, seed, verbose, env_vars, response
#     忽略：R, algorithm, score, pc_*, notears_*
#
# (4) notears_linear — 线性 NOTEARS（torch 实现）
#     依赖：torch；变量数 p 不宜过大（包内约 p<=60）
#     生效参数：notears_lambda, notears_max_iter, notears_lr, notears_tol,
#               notears_rho_init, notears_alpha_mult, max_rows, seed, verbose
#     忽略：R, algorithm, score, pc_*, bidag_*
#
# 一键安装常用建议包（按需删减后运行）：
#   install.packages(c("bnlearn", "BiDAG", "torch"))
# ==============================================================================

# ── 集中配置: 与 cast_prepare / cast_dag / cast_ate / cast_screen / cast_fit /
#    cast_evaluate / cast_predict / cast_cv / cast_cate / cast_shap_xgb /
#    cast_shap_fit（需 fastshap）等
#    的形参一一对应, 便于对照帮助文档调参 ───────────────────────────────────
CONFIG <- list(
  response          = "presence",
  seed              = 42L,
  # cast_prepare
  prepare_train_fraction = 0.7,
  prepare_env_vars       = NULL,
  prepare_verbose        = TRUE,
  # cast_dag (注意: 选 PC/NOTEARS/BiDAG 时改 dag_structure_method;
  #   dag_algorithm / dag_score 仅当 structure_method = "bootstrap_hc" 时生效)
  dag_env_vars            = NULL,
  dag_R                   = 100L,  #speed
  dag_algorithm           = "hc",
  dag_score               = "bic-g",
  dag_strength_threshold  = 0.7,
  dag_direction_threshold = 0.6,
  dag_max_rows            = 8000L,
  dag_verbose             = TRUE,
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
  dag_blacklist           = NULL,   # data.frame(from, to) 禁止边
  dag_whitelist           = NULL,   # data.frame(from, to) 强制边
  # cast_ate
  ate_variables     = NULL,
  ate_K             = 5L,
  ate_num_trees     = 300L,
  ate_alpha         = 0.05,
  ate_quantile_cuts = c(0.25, 0.5, 0.75),
  ate_p_adjust     = "fdr",
  ate_parallel      = TRUE,
  ate_verbose       = TRUE,
  # cast_screen
  screen_min_vars     = 5L,
  screen_min_fraction = 0.5,
  screen_num_trees    = 500L,
  screen_verbose      = TRUE,
  # cast_fit (快速演示: n_epochs / n_runs 较小; 正式分析请增大并打开 tune_grid)
  fit_models            = c("cast", "rf", "maxent", "brt"),
  fit_n_epochs          = 50L,  #speed
  fit_n_runs            = 1L, #speed
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
  fit_tune_grid         = FALSE,  #speed
  fit_verbose           = TRUE,
  # cast_evaluate
  eval_response = "presence",
  # cast_predict
  predict_models = NULL,
  # cast_cv (可选 Step 6b)
  cv_k             = 5L,
  cv_models        = c("rf", "maxent", "brt"),
  cv_block_method  = "grid",
  cv_n_epochs      = 120L,
  cv_n_runs        = 2L,
  cv_rf_ntree      = 300L,
  cv_brt_n_trees   = 500L,
  cv_parallel      = TRUE,
  cv_verbose       = TRUE,
  run_spatial_cv   = TRUE,
  # cast_cate
  cate_variables   = NULL,
  cate_top_n       = 2L,
  cate_n_trees     = 300L,
  cate_verbose     = FALSE,
  # cast_shap_xgb
  shap_nrounds          = 200L,
  shap_max_depth        = 6L,
  shap_eta              = 0.05,
  shap_subsample        = 0.8,
  shap_colsample_bytree = 0.8,
  shap_test_fraction    = 0.2,
  shap_verbose          = FALSE,
  shap_plot_top_n       = 15L,
  # cast_shap_fit（RF / CAST，依赖 fastshap；MC 次数越大越稳但更慢）
  shap_fastshap_nsim      = 40L,
  shap_max_explain_rows   = 50L,
  # cast_evalue
  do_evalue               = TRUE,
  evalue_transform        = "RR",
  evalue_p0               = 0.5,
  # cast_backdoor
  do_backdoor             = TRUE,
  # cast_uncertainty (MC Dropout; 仅 ci_mlp / cast 模型可用)
  do_uncertainty          = TRUE,
  uncertainty_n_forward   = 50L,
  # cast_report (需 rmarkdown)
  do_report               = FALSE,
  # 仅重绘 HSS/CATE 热图（需 ovis_spatial_replot_cache.rds；见文末「空间热图重绘」）
  only_replot_spatial_heatmap = FALSE,
  spatial_heatmap_res_deg      = 0.06,
  spatial_heatmap_interp_method = "nearest",
  spatial_heatmap_display_res   = 0.02
)

# 是否对四种 cast_dag structure_method 逐一做「DAG → ATE → screen → roles」闭环自检
RUN_DAG_SELFTEST_ALL_METHODS <- TRUE

#' 盘羊示例：四种 DAG 结构学习 + 下游最小闭环 (内部函数，仅本脚本使用)
#' @noRd
ovis_dag_selftest_all <- function(train_df, cfg) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop(
      "cast_ate 需要 ranger。install.packages(\"ranger\")",
      call. = FALSE
    )
  }

  ncap <- min(1600L, nrow(train_df))
  d <- train_df[seq_len(ncap), , drop = FALSE]
  methods <- c("bootstrap_hc", "pc", "bidag_bge", "notears_linear")

  for (sm in methods) {
    miss <- character(0)
    if (!requireNamespace("bnlearn", quietly = TRUE)) {
      miss <- c(miss, "bnlearn")
    }
    if (identical(sm, "bidag_bge") && !requireNamespace("BiDAG", quietly = TRUE)) {
      miss <- c(miss, "BiDAG")
    }
    if (identical(sm, "notears_linear") && !requireNamespace("torch", quietly = TRUE)) {
      miss <- c(miss, "torch")
    }
    if (length(miss)) {
      stop(
        "DAG 全算法自检未满足依赖 (方法 ", encodeString(sm, quote = "\""), "): ",
        paste(miss, collapse = ", "),
        "\n请安装: install.packages(c(",
        paste(paste0("\"", miss, "\""), collapse = ", "), "))",
        call. = FALSE
      )
    }

    dag <- tryCatch(
      cast_dag(
        data = d,
        response = cfg$response,
        env_vars = cfg$dag_env_vars,
        R = if (identical(sm, "bootstrap_hc")) 15L else cfg$dag_R,
        algorithm = cfg$dag_algorithm,
        score = cfg$dag_score,
        strength_threshold = if (identical(sm, "bootstrap_hc")) {
          0.45
        } else {
          cfg$dag_strength_threshold
        },
        direction_threshold = if (identical(sm, "bootstrap_hc")) {
          0.35
        } else {
          cfg$dag_direction_threshold
        },
        max_rows = min(if (identical(sm, "bidag_bge")) 1800L else 2500L, ncap),
        seed = cfg$seed,
        verbose = FALSE,
        structure_method = sm,
        pc_alpha = cfg$dag_pc_alpha,
        pc_test = cfg$dag_pc_test,
        bidag_algorithm = cfg$dag_bidag_algorithm,
        bidag_iterations = if (identical(sm, "bidag_bge")) {
          if (!is.null(cfg$dag_bidag_iterations)) {
            cfg$dag_bidag_iterations
          } else {
            80L
          }
        } else {
          cfg$dag_bidag_iterations
        },
        notears_lambda = cfg$dag_notears_lambda,
        notears_max_iter = if (identical(sm, "notears_linear")) {
          min(700L, cfg$dag_notears_max_iter)
        } else {
          cfg$dag_notears_max_iter
        },
        notears_lr = cfg$dag_notears_lr,
        notears_tol = if (identical(sm, "notears_linear")) {
          max(0.01, cfg$dag_notears_tol)
        } else {
          cfg$dag_notears_tol
        },
        notears_rho_init = cfg$dag_notears_rho_init,
        notears_alpha_mult = cfg$dag_notears_alpha_mult,
        blacklist = cfg$dag_blacklist,
        whitelist = cfg$dag_whitelist
      ),
      error = function(e) {
        stop(
          "cast_dag(structure_method=",
          encodeString(sm, quote = "\""),
          "): ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (!inherits(dag, "cast_dag") || length(dag$nodes) < 3L) {
      stop("DAG 自检返回异常: ", sm, call. = FALSE)
    }

    ate <- tryCatch(
      cast_ate(
        data = d,
        response = cfg$response,
        variables = dag$nodes,
        K = 3L,
        num_trees = 120L,
        alpha = cfg$ate_alpha,
        quantile_cuts = cfg$ate_quantile_cuts,
        p_adjust = cfg$ate_p_adjust,
        parallel = FALSE,
        seed = cfg$seed,
        verbose = FALSE
      ),
      error = function(e) {
        stop("cast_ate (", sm, "): ", conditionMessage(e), call. = FALSE)
      }
    )

    screen <- tryCatch(
      cast_screen(
        dag = dag,
        ate = ate,
        data = d,
        response = cfg$response,
        min_vars = cfg$screen_min_vars,
        min_fraction = cfg$screen_min_fraction,
        num_trees = 180L,
        seed = cfg$seed,
        verbose = FALSE
      ),
      error = function(e) {
        stop("cast_screen (", sm, "): ", conditionMessage(e), call. = FALSE)
      }
    )

    roles <- cast_roles(screen = screen, dag = dag)
    if (length(screen$selected) < 1L) {
      stop("cast_screen 未保留任何变量 (", sm, ")", call. = FALSE)
    }

    message(
      "✔ DAG 闭环自检: ", sm,
      " | edges=", nrow(dag$edges),
      " | selected=", length(screen$selected),
      " | roles=", nrow(roles$roles)
    )
  }

  message("cast_dag 四种 structure_method 自检均通过 (bootstrap_hc, pc, bidag_bge, notears_linear)。")
  invisible(TRUE)
}


# ══════════════════════════════════════════════════════════════════════════════
# 方式 A: 模块化逐步运行 (推荐初次使用, 便于理解每一步)
# ══════════════════════════════════════════════════════════════════════════════

# ── Step 1: 数据准备 ─────────────────────────────────────────────────────────
split <- cast_prepare(
  data           = ovis_ammon,
  train_fraction = CONFIG$prepare_train_fraction,
  seed           = CONFIG$seed,
  env_vars       = CONFIG$prepare_env_vars,
  verbose        = CONFIG$prepare_verbose
)
cat("训练集:", nrow(split$train), "| 测试集:", nrow(split$test), "\n")
cat("环境变量:", paste(split$env_vars, collapse = ", "), "\n")

# ── DAG 四种算法 + ATE/screen/roles 最小闭环自检 (默认必须全部通过) ─────────
if (isTRUE(RUN_DAG_SELFTEST_ALL_METHODS)) {
  ovis_dag_selftest_all(split$train, CONFIG)
}

# 可选: 多重共线性筛查 (需建议包 car)
if (requireNamespace("car", quietly = TRUE)) {
  vif_tbl <- cast_vif(
    data          = split$train,
    threshold     = 10,
    exclude       = c("lon", "lat", CONFIG$response),
    expert_filter = NULL,
    verbose       = TRUE
  )
  print(vif_tbl)
}

# ── Step 2: DAG 因果结构学习 ─────────────────────────────────────────────────
# structure_method 可选:
#   "bootstrap_hc"(默认), "pc", "bidag_bge"(需 BiDAG), "notears_linear"(需 torch)
# 下面列出 cast_dag() 全部形参; R / algorithm / score 仅用于 bootstrap_hc。
dag <- cast_dag(
  data                 = split$train,
  response             = CONFIG$response,
  env_vars             = CONFIG$dag_env_vars,
  R                    = CONFIG$dag_R,
  algorithm            = CONFIG$dag_algorithm,
  score                = CONFIG$dag_score,
  strength_threshold   = CONFIG$dag_strength_threshold,
  direction_threshold  = CONFIG$dag_direction_threshold,
  max_rows             = CONFIG$dag_max_rows,
  seed                 = CONFIG$seed,
  verbose              = CONFIG$dag_verbose,
  structure_method     = CONFIG$dag_structure_method,
  pc_alpha             = CONFIG$dag_pc_alpha,
  pc_test              = CONFIG$dag_pc_test,
  bidag_algorithm      = CONFIG$dag_bidag_algorithm,
  bidag_iterations     = CONFIG$dag_bidag_iterations,
  notears_lambda       = CONFIG$dag_notears_lambda,
  notears_max_iter     = CONFIG$dag_notears_max_iter,
  notears_lr           = CONFIG$dag_notears_lr,
  notears_tol          = CONFIG$dag_notears_tol,
  notears_rho_init     = CONFIG$dag_notears_rho_init,
  notears_alpha_mult   = CONFIG$dag_notears_alpha_mult,
  blacklist            = CONFIG$dag_blacklist,
  whitelist            = CONFIG$dag_whitelist
)
print(dag)

# ── Step 3: ATE 因果效应估计 ─────────────────────────────────────────────────
ate <- cast_ate(
  data          = split$train,
  response      = CONFIG$response,
  variables     = CONFIG$ate_variables,
  K             = CONFIG$ate_K,
  num_trees     = CONFIG$ate_num_trees,
  alpha         = CONFIG$ate_alpha,
  quantile_cuts = CONFIG$ate_quantile_cuts,
  p_adjust        = CONFIG$ate_p_adjust,
  parallel      = CONFIG$ate_parallel,
  seed          = CONFIG$seed,
  verbose       = CONFIG$ate_verbose
)
print(ate)

# ── Step 3b: E-value 敏感性分析 ──────────────────────────────────────────────
if (isTRUE(CONFIG$do_evalue)) {
  evalue <- cast_evalue(
    ate       = ate,
    transform = CONFIG$evalue_transform,
    p0        = CONFIG$evalue_p0,
    verbose   = TRUE
  )
  print(evalue)
}

# ── Step 3c: 后门准则检查 ────────────────────────────────────────────────────
if (isTRUE(CONFIG$do_backdoor) &&
    requireNamespace("dagitty", quietly = TRUE)) {
  backdoor <- cast_backdoor(
    dag     = dag,
    outcome = CONFIG$response,
    verbose = TRUE
  )
  print(backdoor)
}

# ── Step 4: 自适应变量筛选 ───────────────────────────────────────────────────
screen <- cast_screen(
  dag            = dag,
  ate            = ate,
  data           = split$train,
  response       = CONFIG$response,
  min_vars       = CONFIG$screen_min_vars,
  min_fraction   = CONFIG$screen_min_fraction,
  num_trees      = CONFIG$screen_num_trees,
  seed           = CONFIG$seed,
  verbose        = CONFIG$screen_verbose
)
print(screen)

# ── Step 5: 因果角色分配 ─────────────────────────────────────────────────────
roles <- cast_roles(screen = screen, dag = dag)
print(roles)

# ── Step 6: 模型训练 ─────────────────────────────────────────────────────────
# 若已安装 torch, 可将 fit_models 设为含 "cast" / "mlp_ate" 等并增大 n_epochs。
fit_full <- cast_fit(
  data               = split$train,
  screen             = screen,
  dag                = dag,
  ate                = ate,
  models             = CONFIG$fit_models,
  response           = CONFIG$response,
  n_epochs           = CONFIG$fit_n_epochs,
  n_runs             = CONFIG$fit_n_runs,
  patience           = CONFIG$fit_patience,
  val_fraction       = CONFIG$fit_val_fraction,
  focal_gamma        = CONFIG$fit_focal_gamma,
  focal_alpha_mode   = CONFIG$fit_focal_alpha_mode,
  focal_alpha        = CONFIG$fit_focal_alpha,
  rf_ntree           = CONFIG$fit_rf_ntree,
  brt_n_trees        = CONFIG$fit_brt_n_trees,
  brt_depth          = CONFIG$fit_brt_depth,
  hidden_size        = CONFIG$fit_hidden_size,
  dropout            = CONFIG$fit_dropout,
  lr                 = CONFIG$fit_lr,
  batch_size         = CONFIG$fit_batch_size,
  max_interactions   = CONFIG$fit_max_interactions,
  tune_grid          = CONFIG$fit_tune_grid,
  seed               = CONFIG$seed,
  verbose            = CONFIG$fit_verbose
)

# ── Step 6b (可选): 空间交叉验证 ────────────────────────────────────────────
if (isTRUE(CONFIG$run_spatial_cv)) {
  cv_ovis <- cast_cv(
    data         = ovis_ammon,
    screen       = screen,
    dag          = dag,
    ate          = ate,
    k            = CONFIG$cv_k,
    models       = CONFIG$cv_models,
    block_method = CONFIG$cv_block_method,
    response     = CONFIG$response,
    n_epochs     = CONFIG$cv_n_epochs,
    n_runs       = CONFIG$cv_n_runs,
    rf_ntree     = CONFIG$cv_rf_ntree,
    brt_n_trees  = CONFIG$cv_brt_n_trees,
    parallel     = CONFIG$cv_parallel,
    seed         = CONFIG$seed,
    verbose      = CONFIG$cv_verbose
  )
  print(cv_ovis)
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      all(c("lon", "lat") %in% names(ovis_ammon))) {
    p_cv <- plot(
      cv_ovis,
      lon = ovis_ammon$lon,
      lat = ovis_ammon$lat,
      basemap = "china"
    )
    print(p_cv)
    ovis_save_plot(p_cv, "ovis_cv_spatial.png", width = 12, height = 5.5)
  }
}

# ── Step 7: 模型评估 ─────────────────────────────────────────────────────────
eval_result <- cast_evaluate(
  fit       = fit_full,
  test_data = split$test,
  response  = CONFIG$eval_response
)
print(eval_result)

# ── Step 8: 空间预测 ─────────────────────────────────────────────────────────
pred <- cast_predict(
  fit      = fit_full,
  new_data = china_env_grid,
  models   = CONFIG$predict_models
)
print(pred)

# ── Step 9: 绘图 ─────────────────────────────────────────────────────────────
# 以下绘图需要 ggplot2, ggraph, igraph, sf 等
# install.packages(c("ggplot2", "ggraph", "igraph", "sf", "patchwork", "fastshap"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  # SHAP：网络图 + 瀑布图各保存一张（stem 为文件名前缀）
  ovis_shap_save_pair <- function(sh_obj, stem) {
    p_net <- plot(
      sh_obj,
      type = "interaction_network",
      top_n = CONFIG$shap_plot_top_n
    )
    print(p_net)
    ovis_save_plot(
      p_net,
      paste0(stem, "_interaction_network.png"),
      width = 9,
      height = 9
    )
    bv <- sh_obj$bias_shap
    if (is.null(bv) || length(bv) != nrow(sh_obj$shap)) {
      bv <- rep(as.numeric(sh_obj$base_score)[1], nrow(sh_obj$shap))
    }
    marg <- bv + rowSums(as.matrix(sh_obj$shap))
    wi <- which.min(abs(marg - stats::median(marg)))[1L]
    p_wf <- plot(
      sh_obj,
      type = "waterfall",
      top_n = CONFIG$shap_plot_top_n,
      waterfall_row = wi
    )
    print(p_wf)
    ovis_save_plot(
      p_wf,
      paste0(stem, "_waterfall.png"),
      width = 10,
      height = 6
    )
  }

  # DAG 网络图
  p_dag <- plot(
    dag,
    roles = roles,
    screen = screen,
    var_labels = var_labels,
    species = "Ovis_ammon"
  )
  print(p_dag)
  ovis_save_plot(p_dag, "ovis_dag.png", width = 12, height = 9)

  # ATE 森林图
  p_ate <- plot(ate, var_labels = var_labels)
  print(p_ate)
  ovis_save_plot(p_ate, "ovis_ate.png", width = 8, height = 10)

  # E-value 敏感性分析图
  if (isTRUE(CONFIG$do_evalue) && exists("evalue")) {
    p_evalue <- plot(evalue, var_labels = var_labels, type = "evalue")
    print(p_evalue)
    ovis_save_plot(p_evalue, "ovis_evalue.png", width = 9, height = 8)
  }

  # Backdoor 识别性检验表格图
  if (isTRUE(CONFIG$do_backdoor) && exists("backdoor") &&
      is.data.frame(backdoor) && nrow(backdoor) > 0) {
    bd <- backdoor
    bd$variable <- if (!is.null(var_labels)) {
      ifelse(bd$variable %in% names(var_labels), var_labels[bd$variable], bd$variable)
    } else {
      bd$variable
    }
    bd$status <- ifelse(bd$identifiable, "\u2713 Identifiable", "\u2717 Not identifiable")
    p_bd <- ggplot2::ggplot(bd, ggplot2::aes(
      x = stats::reorder(.data$variable, .data$n_paths),
      y = .data$n_paths,
      fill = .data$identifiable
    )) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$status),
        hjust = -0.1, size = 3
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#4DBBD5", "FALSE" = "#E64B35"),
        guide = "none"
      ) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::labs(
        title = "Backdoor Criterion Check",
        subtitle = sprintf(
          "%d/%d variables identifiable (DAG-implied)",
          sum(bd$identifiable), nrow(bd)
        ),
        x = NULL, y = "Number of backdoor paths"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.margin = ggplot2::margin(5, 40, 5, 5))
    print(p_bd)
    ovis_save_plot(p_bd, "ovis_backdoor.png", width = 9, height = 7)
  }

  # 筛选分数分解图
  p_screen <- plot(screen, var_labels = var_labels)
  print(p_screen)
  ovis_save_plot(p_screen, "ovis_screen.png", width = 10, height = 8)

  # 评估对比图
  p_eval <- plot(eval_result)
  print(p_eval)
  ovis_save_plot(p_eval, "ovis_evaluate.png", width = 10, height = 7)

  # 空间栖息地适宜性地图 (需要 sf)：已训练模型各保存一张 HSS 图
  if (requireNamespace("sf", quietly = TRUE)) {
    for (md in pred$models) {
      p_pred <- plot(pred, model = md, basemap = "china")
      print(p_pred)
      ovis_save_plot(
        p_pred,
        sprintf("ovis_predict_%s_china.png", md),
        width = 10,
        height = 7
      )
    }
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
    ovis_save_plot(p_cons, "ovis_consistency.png", width = 16, height = 5.5)
  }

  # 空间 CATE + 按 HSS 截断出图 (需 grf; HSS 模型 fallback: cast → rf → 第一个可用)
  if (requireNamespace("grf", quietly = TRUE)) {
    # 选择 HSS 掩膜模型：优先 cast，若 cast 预测全 NA 则回退到 rf
    cate_hss_model <- "cast"
    if ("cast" %in% pred$models) {
      hss_vals <- pred$predictions[["HSS_cast"]]
      if (all(is.na(hss_vals))) {
        cate_hss_model <- if ("rf" %in% pred$models) "rf" else pred$models[1]
        message("CATE HSS mask: cast 预测全 NA，回退到 ", cate_hss_model)
      }
    } else {
      cate_hss_model <- if ("rf" %in% pred$models) "rf" else pred$models[1]
    }

    cate <- cast_cate(
      data           = split$train,
      variables      = CONFIG$cate_variables,
      ate            = ate,
      screen         = screen,
      response       = CONFIG$response,
      top_n          = CONFIG$cate_top_n,
      n_trees        = CONFIG$cate_n_trees,
      predict_data   = china_env_grid,
      seed           = CONFIG$seed,
      verbose        = CONFIG$cate_verbose
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
          hss_model        = cate_hss_model,
          hss_threshold    = 0.1
        )
        print(pc)
        ovis_save_plot(
          pc,
          sprintf("ovis_cate_%s.png", v),
          width = 10,
          height = 7
        )
      }
    }
  }

  # ── SHAP 三种路径（变量个数与含义不同，勿混读一张图）────────────────────────
  # 1) XGBoost：单独训练的 surrogate，TreeSHAP；自 dag 传入后与 RF 共用 **dag$nodes**
  #    原始环境列（非 cast_features 的 int_*）。瀑布图纵轴为 **logit**。
  # 2) RF：cast_fit 里的 ranger，fastshap，**原始环境列**，概率尺度。
  # 3) CAST：cast_fit 里的神经网络，fastshap，输入为 **cast_features()**（含 int_A_B
  #    = DAG 边 A-B 上的标准化乘积×边强度）；与 XGB/RF 的列空间不同，属预期行为。
  # XGBoost（TreeSHAP）
  if (requireNamespace("xgboost", quietly = TRUE)) {
    sh_xgb <- cast_shap_xgb(
      data               = split$train,
      response           = CONFIG$response,
      env_vars           = NULL,
      screen             = screen,
      dag                = dag,
      nrounds            = CONFIG$shap_nrounds,
      max_depth          = CONFIG$shap_max_depth,
      eta                = CONFIG$shap_eta,
      subsample          = CONFIG$shap_subsample,
      colsample_bytree   = CONFIG$shap_colsample_bytree,
      test_fraction      = CONFIG$shap_test_fraction,
      seed               = CONFIG$seed,
      verbose            = CONFIG$shap_verbose
    )
    ovis_shap_save_pair(sh_xgb, "ovis_shap_xgb")
    # cast_shap_write_csv(sh_xgb, file.path(OVIS_FIG_DIR, "shap_export_xgb"))
  }

  # RF（fastshap + 已拟合 ranger，概率尺度）
  if (requireNamespace("fastshap", quietly = TRUE) &&
      requireNamespace("ranger", quietly = TRUE) &&
      "rf" %in% names(fit_full$models) &&
      !is.null(fit_full$models$rf$model)) {
    sh_rf <- tryCatch(
      cast_shap_fit(
        fit                  = fit_full,
        which                = "rf",
        data                 = split$train,
        response             = CONFIG$response,
        test_fraction        = CONFIG$shap_test_fraction,
        seed                 = CONFIG$seed,
        fastshap_nsim        = CONFIG$shap_fastshap_nsim,
        max_explain_rows     = CONFIG$shap_max_explain_rows,
        verbose              = FALSE
      ),
      error = function(e) {
        message("cast_shap_fit(rf) 跳过: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(sh_rf)) {
      ovis_shap_save_pair(sh_rf, "ovis_shap_rf")
    }
  }

  # CAST 神经网络（fastshap + torch；列为 cast_features 工程特征）
  torch_ok <- requireNamespace("torch", quietly = TRUE) &&
    tryCatch(torch::torch_is_installed(), error = function(e) FALSE)
  if (requireNamespace("fastshap", quietly = TRUE) && isTRUE(torch_ok) &&
      "cast" %in% names(fit_full$models) &&
      !is.null(fit_full$models$cast$model)) {
    sh_cast <- tryCatch(
      cast_shap_fit(
        fit                  = fit_full,
        which                = "cast",
        data                 = split$train,
        response             = CONFIG$response,
        test_fraction        = CONFIG$shap_test_fraction,
        seed                 = CONFIG$seed,
        fastshap_nsim        = CONFIG$shap_fastshap_nsim,
        max_explain_rows     = CONFIG$shap_max_explain_rows,
        verbose              = FALSE
      ),
      error = function(e) {
        message("cast_shap_fit(cast) 跳过: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(sh_cast)) {
      ovis_shap_save_pair(sh_cast, "ovis_shap_cast")
    }
  }

  # ── MC Dropout 不确定性地图 (仅 cast / ci_mlp 模型) ─────────────────────────
  if (isTRUE(CONFIG$do_uncertainty)) {
    torch_ok2 <- requireNamespace("torch", quietly = TRUE) &&
      tryCatch(torch::torch_is_installed(), error = function(e) FALSE)
    has_ci_mlp <- isTRUE(torch_ok2) &&
      ("cast" %in% names(fit_full$models) &&
       !is.null(fit_full$models$cast$model))
    if (has_ci_mlp) {
      unc <- tryCatch(
        cast_uncertainty(
          fit       = fit_full,
          new_data  = china_env_grid,
          n_forward = CONFIG$uncertainty_n_forward,
          verbose   = TRUE
        ),
        error = function(e) {
          message("cast_uncertainty 跳过: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(unc) && requireNamespace("sf", quietly = TRUE)) {
        # 绘制 CV (变异系数) 不确定性地图
        p_unc <- ggplot2::ggplot(
          unc,
          ggplot2::aes(x = lon, y = lat, colour = cv)
        ) +
          ggplot2::geom_point(size = 0.3) +
          ggplot2::scale_colour_viridis_c(
            name = "CV", option = "inferno", direction = -1
          ) +
          ggplot2::coord_sf() +
          ggplot2::labs(
            title = "Ovis ammon — MC Dropout Uncertainty (CV)",
            x = "Longitude", y = "Latitude"
          ) +
          ggplot2::theme_minimal(base_size = 11)
        print(p_unc)
        ovis_save_plot(p_unc, "ovis_uncertainty_cv.png", width = 10, height = 7)
      }
    } else {
      message("cast_uncertainty: 无可用 cast/ci_mlp 模型，跳过 MC Dropout。")
    }
  }

  # 供「仅重绘空间热图」：castSDM 原生 plot.cast_predict / plot.cast_cate 用 geom_point
  # 画 Eco-ISEA3H 格点；此处缓存 pred/cate 供 terra 栅格热图后处理（见文末）。
  cate_hss_model_final <- if (exists("cate_hss_model")) cate_hss_model else "rf"
  tryCatch(
    saveRDS(
      list(
        pred         = pred,
        cate         = if (exists("cate")) cate else NULL,
        var_labels   = var_labels,
        config_flags = list(
          cate_hss_model     = cate_hss_model_final,
          cate_hss_threshold = 0.1,
          species            = "Ovis_ammon"
        )
      ),
      file.path(OVIS_FIG_DIR, "ovis_spatial_replot_cache.rds")
    ),
    error = function(e) {
      warning("ovis_spatial_replot_cache.rds: ", conditionMessage(e))
    }
  )

  # ── 插值渲染：用 terra 栅格热图覆盖 geom_point 的 HSS / CATE 图 ──────────
  hlp_path <- .cast_find_spatial_heatmap_helpers(pkg_root)
  if (!is.na(hlp_path) && nzchar(hlp_path) &&
      requireNamespace("terra", quietly = TRUE) &&
      requireNamespace("sf", quietly = TRUE)) {
    message("── 正在用 terra 插值生成出版级 HSS / CATE 热图 ──")
    sys.source(hlp_path, envir = environment())
    tryCatch(
      cast_spatial_replot_hss_cate_heatmaps(
        pred            = pred,
        cate            = if (exists("cate")) cate else NULL,
        fig_dir         = OVIS_FIG_DIR,
        fig_dpi         = OVIS_FIG_DPI,
        var_labels      = var_labels,
        basemap         = "china",
        res_deg         = CONFIG$spatial_heatmap_res_deg,
        interp_method   = CONFIG$spatial_heatmap_interp_method,
        display_res_deg = CONFIG$spatial_heatmap_display_res,
        hss_model       = cate_hss_model_final,
        hss_threshold   = 0.1,
        species_label   = "Ovis_ammon",
        ovis_style      = TRUE
      ),
      error = function(e) {
        warning("Spatial heatmap replot failed: ", conditionMessage(e))
      }
    )
  } else {
    message("Skip terra heatmap replot: missing helpers or packages (terra/sf).")
  }
}




# ══════════════════════════════════════════════════════════════════════════════
# 方式 B: 一键 Pipeline cast() (与上面模块化等价; 列出 cast() 全部形参)
# ══════════════════════════════════════════════════════════════════════════════
# 注意: 完整流程含空间 CV 时耗时较长。将 RUN_CAST_PIPELINE 改为 TRUE 再运行本段。
RUN_CAST_PIPELINE <- FALSE

if (isTRUE(RUN_CAST_PIPELINE)) {
  result <- cast(
    species_data           = ovis_ammon,
    env_data               = china_env_grid,
    models                 = c("rf", "maxent", "brt"),
    train_fraction         = CONFIG$prepare_train_fraction,
    n_bootstrap            = CONFIG$dag_R,
    dag_structure_method   = CONFIG$dag_structure_method,
    dag_pc_alpha           = CONFIG$dag_pc_alpha,
    dag_pc_test            = CONFIG$dag_pc_test,
    dag_bidag_algorithm    = CONFIG$dag_bidag_algorithm,
    dag_bidag_iterations   = CONFIG$dag_bidag_iterations,
    dag_notears_lambda     = CONFIG$dag_notears_lambda,
    dag_notears_max_iter   = CONFIG$dag_notears_max_iter,
    dag_notears_lr         = CONFIG$dag_notears_lr,
    dag_notears_tol        = CONFIG$dag_notears_tol,
    dag_notears_rho_init   = CONFIG$dag_notears_rho_init,
    dag_notears_alpha_mult = CONFIG$dag_notears_alpha_mult,
    strength_threshold     = CONFIG$dag_strength_threshold,
    direction_threshold    = CONFIG$dag_direction_threshold,
    ate_folds                = 2L,
    ate_alpha                = CONFIG$ate_alpha,
    screen_min_vars          = CONFIG$screen_min_vars,
    do_cv                    = TRUE,
    cv_k                     = CONFIG$cv_k,
    cv_block_method          = CONFIG$cv_block_method,
    do_predict               = NULL,
    do_cate                  = TRUE,
    cate_top_n               = 3L,
    blacklist                = CONFIG$dag_blacklist,
    whitelist                = CONFIG$dag_whitelist,
    do_evalue                = CONFIG$do_evalue,
    do_backdoor              = CONFIG$do_backdoor,
    seed                     = CONFIG$seed,
    verbose                  = TRUE
  )
  print(summary(result))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_cast <- plot(result, var_labels = var_labels)
    print(p_cast)
    ovis_save_plot(p_cast, "ovis_cast_pipeline.png", 12, 10)
  }
  # Pipeline 模式下生成 HTML 报告
  if (isTRUE(CONFIG$do_report) &&
      requireNamespace("rmarkdown", quietly = TRUE)) {
    cast_report(
      result      = result,
      output_file = file.path(OVIS_FIG_DIR, "ovis_ammon_report.html"),
      species     = "Ovis ammon",
      var_labels  = var_labels,
      open        = FALSE,
      verbose     = TRUE
    )
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# 方式 C (可选): 将五种 DAG 网络图落盘 (依赖已通过上面自检)
# ══════════════════════════════════════════════════════════════════════════════
# 与「全算法自检」重复时可关；需要对比图时改为 TRUE。
RUN_DAG_METHOD_SHOWCASE <- TRUE

if (isTRUE(RUN_DAG_METHOD_SHOWCASE) && requireNamespace("ggplot2", quietly = TRUE)) {
  methods <- c("bootstrap_hc", "pc", "bidag_bge", "notears_linear")
  for (sm in methods) {
    dm <- tryCatch(
      cast_dag(
        data = split$train,
        response = CONFIG$response,
        env_vars = CONFIG$dag_env_vars,
        R = CONFIG$dag_R,
        algorithm = CONFIG$dag_algorithm,
        score = CONFIG$dag_score,
        strength_threshold = CONFIG$dag_strength_threshold,
        direction_threshold = CONFIG$dag_direction_threshold,
        max_rows = CONFIG$dag_max_rows,
        seed = CONFIG$seed,
        verbose = FALSE,
        structure_method = sm,
        pc_alpha = CONFIG$dag_pc_alpha,
        pc_test = CONFIG$dag_pc_test,
        bidag_algorithm = CONFIG$dag_bidag_algorithm,
        bidag_iterations = CONFIG$dag_bidag_iterations,
        notears_lambda = CONFIG$dag_notears_lambda,
        notears_max_iter = CONFIG$dag_notears_max_iter,
        notears_lr = CONFIG$dag_notears_lr,
        notears_tol = CONFIG$dag_notears_tol,
        notears_rho_init = CONFIG$dag_notears_rho_init,
        notears_alpha_mult = CONFIG$dag_notears_alpha_mult,
        blacklist = CONFIG$dag_blacklist,
        whitelist = CONFIG$dag_whitelist
      ),
      error = function(e) {
        message("cast_dag (", sm, ") 出图跳过: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(dm)) next
    am <- cast_ate(
      data = split$train,
      response = CONFIG$response,
      variables = dm$nodes,
      K = 5L,
      num_trees = 200L,
      parallel = FALSE,
      seed = CONFIG$seed,
      verbose = FALSE
    )
    smc <- cast_screen(
      dag = dm,
      ate = am,
      data = split$train,
      response = CONFIG$response,
      num_trees = 300L,
      seed = CONFIG$seed,
      verbose = FALSE
    )
    rm <- cast_roles(screen = smc, dag = dm)
    p <- plot(
      dm,
      roles = rm,
      screen = smc,
      var_labels = var_labels,
      species = "Ovis_ammon"
    )
    ovis_save_plot(p, paste0("ovis_dag_showcase_", sm, ".png"), 12, 9)
  }
  message(
    "DAG showcase figures saved under: ",
    normalizePath(OVIS_FIG_DIR, winslash = "/", mustWork = FALSE)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 空间 HSS / CATE 仅热图重绘（lon/lat 规则网格 + nearest 插值，对齐 fig7 Python 的
# griddata + pcolormesh 思路；并按中国边界 shap 做遮罩）
# ══════════════════════════════════════════════════════════════════════════════
# 将 CONFIG$only_replot_spatial_heatmap 设为 TRUE 后，选中本段到结尾运行（或 source
# 全脚本）。需要 ovis_spatial_replot_cache.rds（全量跑 Step 9 后自动写入 OVIS_FIG_DIR）。
# 会覆盖 ovis_predict_*_china.png、ovis_cate_*.png。

if (isTRUE(CONFIG[["only_replot_spatial_heatmap"]])) {
  hlp <- .cast_find_spatial_heatmap_helpers(pkg_root)
  if (is.na(hlp) || !nzchar(hlp)) {
    stop(
      "Could not find cast_spatial_heatmap_helpers.R.\n",
      "  pkg_root = ", deparse(pkg_root), "\n",
      "  getwd()  = ", getwd(), "\n",
      "Set CASTSDM_ROOT to the castSDM source root, or use getwd() under e:/Package/cast.",
      call. = FALSE
    )
  }
  sys.source(hlp, envir = .GlobalEnv)
  cache <- file.path(OVIS_FIG_DIR, "ovis_spatial_replot_cache.rds")
  if (!file.exists(cache)) {
    stop(
      "Missing cache file:\n  ", cache,
      "\nRun the full script once (Step 9 writes this file), then set ",
      "CONFIG$only_replot_spatial_heatmap <- TRUE and re-run this block.",
      call. = FALSE
    )
  }
  z <- readRDS(cache)
  cf <- z$config_flags
  cast_spatial_replot_hss_cate_heatmaps(
    pred            = z$pred,
    cate            = z$cate,
    fig_dir         = OVIS_FIG_DIR,
    fig_dpi         = OVIS_FIG_DPI,
    var_labels      = z$var_labels,
    basemap         = "china",
    res_deg         = CONFIG$spatial_heatmap_res_deg,
    interp_method   = CONFIG$spatial_heatmap_interp_method,
    display_res_deg = CONFIG$spatial_heatmap_display_res,
    hss_model       = if (!is.null(cf) && !is.null(cf$cate_hss_model)) cf$cate_hss_model else "cast",
    hss_threshold   = if (!is.null(cf) && !is.null(cf$cate_hss_threshold)) cf$cate_hss_threshold else 0.1,
    species_label   = if (!is.null(cf) && !is.null(cf$species)) cf$species else "Ovis_ammon",
    ovis_style      = TRUE
  )
  message("Spatial heatmap replot finished (ovis_predict_* / ovis_cate_*).")
}
