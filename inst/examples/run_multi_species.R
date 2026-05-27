# ==============================================================================
# castSDM v0.3.0 Multi-Species Batch Example
#
# This script auto-loads all 32 bundled species datasets (data/*.rda) plus
# china_env_grid and runs the full castSDM pipeline for each species.
# Results are saved to <output_dir>/<Species>/ sub-directories.
#
# v0.3.0 pipeline per species:
#   cast_prepare -> cast_dag (PC/MB-First, presence as node) -> cast_select (MB + RF)
#   -> cast_roles -> cast_fit (RF/BRT/MaxEnt/GAM) -> cast_evaluate
#   -> cast_cv -> cast_predict -> cast_ensemble -> (optional) cast_cate
#   -> (optional) SHAP
#
# Features:
#   - cast_worker_budget(): auto-allocate species x intra parallelism
#   - cast_batch_resume(): crash-safe checkpoint/resume
#   - GAM algorithm in default model set
#   - Per-step checkpoints + resource_log.csv
#   - DAG-guided variable selection via Markov Blanket (first in SDM)
#   - Performance-weighted ensemble prediction
#
# Multi-species additions:
#   - Species-level parallelism (PSOCK via future::multisession)
#   - SHAP: do_shap=TRUE writes shap_xgb/rf figures + shap_panel_2x2.png
#   - Post-hoc SHAP only: CONFIG$only_shap_posthoc <- TRUE
#   - Spatial HSS/CATE heatmap replot: CONFIG$only_replot_spatial_heatmap <- TRUE
#   - Top-level <output_dir>/figures/multi_species_comparison.png
#
# Running in RStudio:
#   1. Open castSDM.Rproj (sets working directory to package root)
#   2. Ctrl+Shift+L (devtools::load_all())
#   3. source() this script
#
# Or after installing:
#   library(castSDM)
#   source(system.file("examples/run_multi_species.R", package = "castSDM"))
#
# Dependencies: future / future.apply / ggplot2 / patchwork / pkgload
#   DAG: bnlearn (required), BiDAG / dagitty / pcalg (optional)
#   Models: ranger, gbm, maxnet, mgcv
#   SHAP: xgboost, fastshap
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

# Dev mode: find package root and load_all
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
if (exists("cast_set_plot_defaults", mode = "function")) {
  cast_set_plot_defaults("Arial")
}

for (pkg in c("ggplot2", "future", "future.apply", "patchwork", "data.table", "pkgload")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# ==============================================================================
# CONFIG — centralised settings (aligned with v0.3.0 API)
# ==============================================================================
CONFIG <- cast_merge_config(cast_default_config("batch"), list(
  # --- multi-species only ---
  plot_font_family  = "Arial",
  only_shap_posthoc  = FALSE,
  only_replot_spatial_heatmap = FALSE,
  run_future_projection = TRUE,
  future_download_cmip6 = TRUE,
  future_gcms = c("ACCESS-CM2", "CMCC-ESM2", "MIROC6",
                  "MRI-ESM2-0", "IPSL-CM6A-LR"),
  future_ssps = c("245", "585"),
  future_periods = c("2041-2060", "2061-2080"),
  future_var = "bioc",
  future_res = 2.5,
  future_ensemble = TRUE,
  future_cache_dir = file.path("castSDM_multi_species", "cmip6_cache"),
  future_envs_path = NULL,
  future_save_dir = NULL,
  spatial_heatmap_res_deg      = 0.06,
  spatial_heatmap_interp_method = "nearest",
  spatial_heatmap_display_res   = 0.02,

  # --- species source ---
  data_dir          = NULL,
  env_grid_path     = NULL,

  output_dir        = "castSDM_multi_species",
  fig_dpi           = 600L,
  n_workers         = min(6L, max(1L, parallel::detectCores() - 1L)),
  parallel_species  = TRUE,
  resume            = TRUE,
  batch_verbose     = TRUE,

  response          = "presence",
  seed              = 42L,

  # -- cast_prepare --
  prepare_train_fraction = 0.7,
  prepare_env_vars       = NULL,
  prepare_verbose        = FALSE,

  # -- cast_dag (v0.3.0: response-focused MB, presence as node) --
  dag_env_vars            = NULL,
  dag_R                   = 100L,
  learn_shared_dag        = FALSE,
  dag_algorithm           = "hc",
  dag_score               = NULL,   # auto: bic-cg for mixed data
  dag_strength_threshold  = 0.7,
  dag_direction_threshold = 0.6,
  dag_max_rows            = 8000L,
  dag_verbose             = FALSE,
  dag_structure_method    = "mb_first", # pc | bootstrap_hc | mb_first | bidag_bge
  dag_include_response    = TRUE,    # v0.3.0: presence as DAG node
  dag_pc_alpha            = 0.05,
  dag_pc_test             = NULL,    # auto: mi-cg for mixed data
  dag_mb_method           = "fast.iamb", # MB discovery (mb_first only): fast.iamb | iamb | inter.iamb | gs
  dag_mb_alpha            = 0.05,        # significance for MB discovery phase
  dag_bidag_algorithm     = "order",
  dag_bidag_iterations    = NULL,


  # -- cast_select (v0.3.0: replaces ATE + screen) --
  select_min_vars     = 5L,
  select_min_fraction = 0.3,
  select_num_trees    = 300L,
  select_verbose      = FALSE,

  # -- cast_fit (RF / BRT / MaxEnt / GAM) --
  fit_models            = c("rf", "brt", "maxent", "gam"),
  fit_rf_ntree          = 300L,
  fit_brt_n_trees       = 500L,
  fit_brt_depth         = 5L,
  fit_verbose           = FALSE,

  # -- cast_evaluate --
  eval_response   = NULL,

  # -- cast_predict --
  predict_models  = NULL,
  plot_basemap    = "china",

  # -- cast_cv --
  run_spatial_cv  = TRUE,
  cv_k             = 5L,
  cv_models        = NULL,
  cv_block_method  = "grid",
  cv_rf_ntree      = 300L,
  cv_brt_n_trees   = 500L,
  cv_parallel      = TRUE,
  cv_verbose       = FALSE,

  # -- cast_cate (DAG-guided confounders) --
  do_cate          = TRUE,
  cate_variables   = NULL,
  cate_top_n       = 3L,
  cate_n_trees     = 300L,
  cate_verbose     = FALSE,
  cate_point_size  = 0.45,
  cate_hss_model      = "rf",       # v0.3.0: "cast" removed
  cate_hss_threshold  = 0.1,

  var_labels          = NULL,

  # -- SHAP --
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
))

if (exists("cast_set_plot_defaults", mode = "function")) {
  cast_set_plot_defaults(CONFIG$plot_font_family)
}

cat("============================================================\n")
cat("castSDM v0.3.0 multi-species batch\n")
cat("============================================================\n\n")

# ==============================================================================
# Build species_list + env_grid
# ==============================================================================

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
  data_dir <- if (!is.na(pkg_root)) {
    file.path(pkg_root, "data")
  } else {
    system.file("data", package = "castSDM")
  }
  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  rda_names <- tools::file_path_sans_ext(basename(rda_files))
  exclude <- c("china_env_grid")
  keep <- !rda_names %in% exclude
  rda_files <- rda_files[keep]
  rda_names <- rda_names[keep]

  species_list <- list()
  for (i in seq_along(rda_files)) {
    e <- new.env(parent = emptyenv())
    load(rda_files[i], envir = e)
    obj_name <- ls(e)[1]
    species_list[[rda_names[i]]] <- get(obj_name, envir = e)
  }
  cat(sprintf("[Package data] %d species loaded from %s\n",
              length(species_list), data_dir))
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
  grid_rda <- if (!is.na(pkg_root)) {
    file.path(pkg_root, "data", "china_env_grid.rda")
  } else {
    ""
  }
  if (nzchar(grid_rda) && file.exists(grid_rda)) {
    e <- new.env(parent = emptyenv())
    load(grid_rda, envir = e)
    env_grid <- get(ls(e)[1], envir = e)
  } else {
    tryCatch({
      env_grid <- get("china_env_grid", envir = asNamespace("castSDM"))
    }, error = function(err) NULL)
  }
  if (!is.null(env_grid)) {
    cat(sprintf("\nPrediction grid (china_env_grid): %d x %d\n",
                nrow(env_grid), ncol(env_grid)))
  } else {
    pkg_grid <- system.file(
      "extdata", "China_EnvData_Res9_Screened.csv",
      package = "castSDM"
    )
    if (pkg_grid != "") {
      env_grid <- data.table::fread(pkg_grid, data.table = FALSE)
      cat(sprintf("\nPrediction grid (extdata CSV): %d x %d\n",
                  nrow(env_grid), ncol(env_grid)))
    } else {
      cat("\nNo prediction grid; spatial prediction skipped.\n")
    }
  }
}

dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Spatial heatmap replot only
# ==============================================================================

if (isTRUE(CONFIG$only_replot_spatial_heatmap)) {
  future::plan(future::sequential)
  hlp <- .cast_find_spatial_heatmap_helpers(pkg_root)
  if (is.na(hlp) || !nzchar(hlp)) {
    stop(
      "Could not find cast_spatial_heatmap_helpers.R.\n",
      "  pkg_root = ", deparse(pkg_root), "\n",
      "  getwd()  = ", getwd(), "\n",
      "Set CASTSDM_ROOT to the castSDM source root.",
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
  cat("\n[only_shap_posthoc] Writing SHAP figures...\n")
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
  cat("\nPost-hoc SHAP complete (also see shap_panel_2x2.png per species).\n\n")
} else {

# ==============================================================================
# Parallel backend + worker budget
# ==============================================================================

use_par <- isTRUE(CONFIG$parallel_species) && CONFIG$n_workers > 1L

worker_budget <- cast_worker_budget(
  total_workers = CONFIG$n_workers,
  n_species     = length(species_list)
)
cat(sprintf(
  "\nParallel: %s | budget: total=%d species=%d intra=%d (cores: %d)\n",
  if (use_par) "ON" else "OFF",
  worker_budget$total, worker_budget$species, worker_budget$intra,
  parallel::detectCores()
))
if (use_par) {
  future::plan(future::multisession, workers = worker_budget$species)
} else {
  future::plan(future::sequential)
}

# ==============================================================================
# cast_batch (v0.3.0 API: no ATE, no NN params)
# ==============================================================================

cat("\nRunning cast_batch()...\n\n")

eval_resp <- CONFIG$eval_response
if (is.null(eval_resp)) eval_resp <- CONFIG$response

do_resume <- isTRUE(CONFIG$resume) && dir.exists(CONFIG$output_dir) &&
  any(vapply(names(species_list), function(sp)
    file.exists(file.path(CONFIG$output_dir, sp, "cast_result.rds")),
    logical(1)))

batch_args <- list(
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
  learn_shared_dag        = CONFIG$learn_shared_dag,
  shared_dag_data         = env_grid,
  dag_structure_method    = CONFIG$dag_structure_method,
  dag_include_response    = CONFIG$dag_include_response,
  dag_pc_alpha            = CONFIG$dag_pc_alpha,
  dag_pc_test             = CONFIG$dag_pc_test,
  dag_mb_method           = CONFIG$dag_mb_method,
  dag_mb_alpha            = CONFIG$dag_mb_alpha,
  dag_bidag_algorithm     = CONFIG$dag_bidag_algorithm,
  dag_bidag_iterations    = CONFIG$dag_bidag_iterations,
  dag_algorithm           = CONFIG$dag_algorithm,
  dag_score               = CONFIG$dag_score,
  dag_strength_threshold  = CONFIG$dag_strength_threshold,
  dag_direction_threshold = CONFIG$dag_direction_threshold,
  dag_max_rows            = CONFIG$dag_max_rows,
  select_min_vars     = CONFIG$select_min_vars,
  select_min_fraction = CONFIG$select_min_fraction,
  select_num_trees    = CONFIG$select_num_trees,
  select_verbose      = CONFIG$select_verbose,
  do_cv               = CONFIG$run_spatial_cv,
  cv_k                = CONFIG$cv_k,
  cv_models           = CONFIG$cv_models,
  cv_block_method     = CONFIG$cv_block_method,
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
  rf_ntree         = CONFIG$fit_rf_ntree,
  brt_n_trees      = CONFIG$fit_brt_n_trees,
  brt_depth         = CONFIG$fit_brt_depth
)

if (do_resume) {
  cat("[Resume mode] Skipping species with existing cast_result.rds...\n")
  resume_args <- batch_args
  resume_args$output_dir <- NULL
  batch_result <- do.call(cast_batch_resume, c(
    list(output_dir = CONFIG$output_dir),
    resume_args
  ))
} else {
  batch_result <- do.call(cast_batch, batch_args)
}

if (!is.null(batch_result)) {
  print(batch_result)
} else {
  cat("\nAll species already complete (resume mode, nothing to re-run).\n")
}

# ==============================================================================
# Cross-species comparison
# ==============================================================================

fig_dir <- file.path(CONFIG$output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!is.null(batch_result) && !is.null(batch_result$species_metrics) &&
    nrow(batch_result$species_metrics) > 0) {
  cat("\n[Plot] Multi-species model performance comparison...\n")
  p_compare <- plot(batch_result, metrics = c("auc", "tss", "cbi"))
  if (exists("cast_safe_ggsave", mode = "function")) {
    cast_safe_ggsave(
      file.path(fig_dir, "multi_species_comparison.png"),
      p_compare,
      width = 14, height = 6, dpi = CONFIG$fig_dpi,
      bg = "white", limitsize = FALSE
    )
  } else {
    ggplot2::ggsave(
      file.path(fig_dir, "multi_species_comparison.png"),
      p_compare,
      width = 14, height = 6, dpi = CONFIG$fig_dpi,
      bg = "white", limitsize = FALSE
    )
  }
  cat("  Saved: multi_species_comparison.png\n")
}

if (isTRUE(CONFIG$run_future_projection)) {
  cat("\n[Future] Multi-species projection...\n")
  if (!is.null(CONFIG$future_envs_path) && nzchar(CONFIG$future_envs_path)) {
    future_envs <- cast_load_future_envs(CONFIG$future_envs_path)
  } else if (isTRUE(CONFIG$future_download_cmip6) &&
             !is.null(env_grid) &&
             all(c("lon", "lat") %in% names(env_grid)) &&
             any(grepl("^bio\\d{2}$", names(env_grid)))) {
    cmip6_data <- cast_download_cmip6(
      gcms     = CONFIG$future_gcms,
      ssps     = CONFIG$future_ssps,
      periods  = CONFIG$future_periods,
      var      = CONFIG$future_var,
      res      = CONFIG$future_res,
      ensemble = CONFIG$future_ensemble,
      path     = CONFIG$future_cache_dir,
      verbose  = TRUE
    )
    future_envs <- cast_prepare_future_env(
      rasters     = cmip6_data,
      coords      = env_grid[, c("lon", "lat")],
      static_vars = env_grid,
      save_dir    = CONFIG$future_save_dir %||% file.path(CONFIG$output_dir, "future_envs"),
      verbose     = TRUE
    )
  } else {
    warning("Future projection skipped: set CONFIG$future_envs_path or enable CMIP6 with bio01-bio19 current grid columns.")
    future_envs <- NULL
  }

  if (is.null(future_envs)) {
    invisible(NULL)
  } else {
    future_results <- if (!is.null(batch_result) && !is.null(batch_result$results) &&
                          length(batch_result$results)) {
      batch_result$results
    } else {
      rds_paths <- file.path(CONFIG$output_dir, names(species_list), "cast_result.rds")
      names(rds_paths) <- names(species_list)
      loaded <- lapply(rds_paths[file.exists(rds_paths)], readRDS)
      names(loaded) <- names(rds_paths)[file.exists(rds_paths)]
      loaded
    }
    if (!length(future_results)) {
      warning("Future projection skipped: no in-memory or saved cast_result.rds species results found.")
    } else {
    future_root <- CONFIG$future_save_dir %||% file.path(CONFIG$output_dir, "future_projection")
    dir.create(future_root, recursive = TRUE, showWarnings = FALSE)
    for (sp in names(future_results)) {
      r <- future_results[[sp]]
      if (is.null(r) || is.null(r$fit) || is.null(r$cv) || is.null(r$split)) next
      env_vars <- r$split$env_vars %||% get_env_vars(species_list[[sp]], response = CONFIG$response)
      current_env <- env_grid[, unique(c("lon", "lat", env_vars)), drop = FALSE]
      future_envs_use <- tryCatch(
        lapply(future_envs, function(x) {
          need <- unique(c("lon", "lat", env_vars))
          miss <- setdiff(need, names(x))
          if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
          x[, need, drop = FALSE]
        }),
        error = function(e) {
          warning(sprintf("%s future projection skipped: %s", sp, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(future_envs_use)) next
      tryCatch(
        cast_save_future_projection(
          fit = r$fit,
          cv = r$cv,
          current_env = current_env,
          future_envs = future_envs_use,
          save_dir = file.path(future_root, sp),
          basemap = CONFIG$plot_basemap,
          fig_dpi = CONFIG$fig_dpi,
          method = "weighted",
          prefix = paste0(sp, "_future")
        ),
        error = function(e) warning(sprintf("%s future projection failed: %s", sp, conditionMessage(e)))
      )
    }
    }
  }
}

res_log <- file.path(CONFIG$output_dir, "resource_log.csv")
if (file.exists(res_log)) {
  rl <- utils::read.csv(res_log, stringsAsFactors = FALSE)
  cat(sprintf("\n[Resource log] %d entries across %d species\n",
              nrow(rl), length(unique(rl$species))))
  cat(sprintf("  Total wall time: %.1f min\n",
              sum(rl$elapsed_sec, na.rm = TRUE) / 60))
}

future::plan(future::sequential)

cat("\n============================================================\n")
cat("Multi-species batch complete.\n")
cat("Output: ", normalizePath(CONFIG$output_dir, winslash = "/",
                              mustWork = FALSE), "\n", sep = "")
cat("============================================================\n\n")

for (sp in names(species_list)) {
  sp_figs <- file.path(CONFIG$output_dir, sp, "figures")
  sp_rds  <- file.path(CONFIG$output_dir, sp, "cast_result.rds")
  if (dir.exists(sp_figs)) {
    figs <- list.files(sp_figs, pattern = "\\.png$")
    cat(sprintf("  %s: %d PNG(s) %s\n", sp, length(figs),
                if (file.exists(sp_rds)) "[done]" else "[partial]"))
  } else if (file.exists(sp_rds)) {
    cat(sprintf("  %s: [done]\n", sp))
  }
}

} # end else (full cast_batch path)
