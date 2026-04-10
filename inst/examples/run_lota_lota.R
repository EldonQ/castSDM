# ==============================================================================
# castSDM Freshwater Fish Example: Lota lota (Burbot / 江鳕)
#
# This script demonstrates how to use castSDM for freshwater fish species
# distribution modeling using external CSV data (e.g., from RiverATLAS).
#
# ── INPUT CSV FORMAT ──────────────────────────────────────────────────────────
# Your fish occurrence CSV must contain these columns:
#
#   Column 1: lon        — Longitude (decimal degrees, WGS84)
#   Column 2: lat        — Latitude  (decimal degrees, WGS84)
#   Column 3: presence   — 1 = presence, 0 = background/pseudo-absence
#   Columns 4+: <env_vars> — Numeric environmental predictor variables
#                            (e.g., discharge, temperature, elevation, ...)
#
# Optional columns (will be auto-excluded from modeling):
#   hyriv_id, species, region, occ_lon, occ_lat, snap_dist, etc.
#
# If you have a column_metadata.csv that flags which columns are environmental
# variables, set META_CSV below. Otherwise, castSDM will auto-detect numeric
# columns (excluding lon, lat, presence, and known metadata names).
#
# ── WORKFLOW ─────────────────────────────────────────────────────────────────
#   VIF → cast_prepare → cast_dag → cast_ate → cast_screen
#   → cast_roles → cast_fit → cast_evaluate → cast_cv → cast_predict
#   → cast_cate → Plotting & Export
#
# ── CHECKPOINT SYSTEM ────────────────────────────────────────────────────────
#   Each step saves results to MODEL_DIR/*.rds.
#   Re-running auto-skips completed steps; delete an .rds to re-compute.
#
# Date: 2026-04-10
# ==============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 0. CONFIGURATION — edit this section for your species
# ─────────────────────────────────────────────────────────────────────────────

# ── Paths ────────────────────────────────────────────────────────────────────
# ROOT_DIR: working directory containing your CSV and output/
# TRAIN_CSV: path to your fish occurrence CSV (see format above)
# META_CSV: path to column metadata CSV (optional; set NULL to auto-detect)
ROOT_DIR     <- "."
OUT_DIR      <- file.path(ROOT_DIR, "output")
MODEL_DIR    <- file.path(OUT_DIR, "model")
PRED_DIR     <- file.path(OUT_DIR, "predictions")
FIG_DIR      <- file.path(OUT_DIR, "figures")

TRAIN_CSV    <- file.path(OUT_DIR, "CAST_Lota_lota.csv")
META_CSV     <- file.path(OUT_DIR, "column_metadata.csv")

# ── Model selection ──────────────────────────────────────────────────────────
# Options: "cast"(needs torch), "mlp_ate"(needs torch), "mlp"(needs torch),
#          "rf", "maxent", "brt"
MODELS <- c("cast", "rf", "maxent", "brt")

# ── Spatial cross-validation ─────────────────────────────────────────────────
DO_SPATIAL_CV <- TRUE
CV_FOLDS      <- 5L
CV_METHOD     <- "grid"      # "grid" = uniform grid | "cluster" = k-means

# ── Spatial prediction regions ───────────────────────────────────────────────
# RiverATLAS region codes: eu=Europe, na=North America, si=Siberia, ar=Arctic
# Set NULL or empty to skip spatial prediction
PREDICT_REGIONS <- c("eu", "na", "si", "ar")

# ── Manual environment variable pre-selection (before VIF) ───────────────────
# Set NULL to use ALL environmental variables from the CSV.
# Otherwise, list the variable names you want to consider.
MANUAL_ENV_VARS <- c(
  # Hydrology
  "dis_av_cms",  # Mean annual discharge (m³/s)
  "lka_pc_use",  # Lake area percentage (%)
  "rev_mc_usu",  # Reservoir volume (million m³)
  "ria_ha_usu",  # River area (hectares)
  "riv_tc_usu",  # River volume (thousand m³)
  # Topography
  "ele_mt_uav",  # Elevation (m)
  "slp_dg_uav",  # Terrain slope (degrees)
  "sgr_dk_rav",  # Stream gradient (dm/km)
  # Climate
  "tmp_dc_uyr",  # Annual mean temperature (°C × 10)
  "pre_mm_uyr",  # Annual precipitation (mm)
  "aet_mm_uyr",  # Actual evapotranspiration (mm)
  "ari_ix_uav",  # Global aridity index
  "cmi_ix_uyr",  # Climate moisture index
  "snw_pc_uyr",  # Snow cover extent (%)
  # Land cover
  "for_pc_use",  # Forest cover (%)
  "crp_pc_use",  # Cropland cover (%)
  "pst_pc_use",  # Pasture cover (%)
  # Soil
  "cly_pc_uav",  # Clay content (%)
  "slt_pc_uav",  # Silt content (%)
  "snd_pc_uav",  # Sand content (%)
  "soc_th_uav",  # Soil organic carbon (t/ha)
  "swc_pc_uyr",  # Soil water content (%)
  "ero_kh_uav",  # Soil erosion intensity (t/ha/yr)
  # Human activity
  "pop_ct_usu",  # Population count
  "ppd_pk_uav",  # Population density (ppl/km²)
  "urb_pc_use",  # Urban cover (%)
  "nli_ix_uav",  # Nighttime light index
  "rdd_mk_uav",  # Road density (km/km²)
  "gdp_ud_usu"   # GDP PPP (million USD)
)

# ── DAG parameters ───────────────────────────────────────────────────────────
DAG_R                   <- 50L
DAG_ALGORITHM           <- "hc"
DAG_SCORE               <- "bic-g"
DAG_STRENGTH_THRESHOLD  <- 0.7
DAG_DIRECTION_THRESHOLD <- 0.6
DAG_MAX_ROWS            <- 8000L

SEED <- 2024L

# ── Figure output ────────────────────────────────────────────────────────────
FIG_DPI        <- 600L
FIG_DPI_GLOBAL <- 600L


# ─────────────────────────────────────────────────────────────────────────────
# Checkpoint helper: skip completed steps on re-run
# ─────────────────────────────────────────────────────────────────────────────
checkpoint <- function(path, expr) {
  if (file.exists(path)) {
    message(sprintf("  [skip] Checkpoint exists: %s", basename(path)))
    return(readRDS(path))
  }
  result <- force(expr)
  saveRDS(result, path)
  message(sprintf("  [save] Checkpoint: %s", basename(path)))
  result
}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Package loading
# ─────────────────────────────────────────────────────────────────────────────

for (d in c(MODEL_DIR, PRED_DIR, FIG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Load castSDM (installed package or dev mode)
if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(castSDM)
}

for (pkg in c("data.table", "dplyr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

set.seed(SEED)
cat("=============================================================\n")
cat("castSDM Freshwater Fish Modeling: Lota lota (Burbot)\n")
cat("=============================================================\n")


# ─────────────────────────────────────────────────────────────────────────────
# 2. Load training data
# ─────────────────────────────────────────────────────────────────────────────

message("\n[1/8] Loading training data...")
if (!file.exists(TRAIN_CSV))
  stop("Training CSV not found: ", TRAIN_CSV,
       "\nPlease prepare your fish occurrence data first.")

cast_data <- data.table::fread(TRAIN_CSV, data.table = FALSE)

# Determine environment variables from metadata or auto-detection
if (!is.null(META_CSV) && file.exists(META_CSV)) {
  col_meta <- read.csv(META_CSV, stringsAsFactors = FALSE)
  env_vars_all <- col_meta$column_name[col_meta$is_env_var]
  message("  Variable detection: using column_metadata.csv")
} else {
  env_vars_all <- get_env_vars(cast_data)
  message("  Variable detection: auto-detected from numeric columns")
}

if (is.null(MANUAL_ENV_VARS)) {
  env_vars_input <- env_vars_all
  message("  Variable selection: all environment variables")
} else {
  bad <- setdiff(MANUAL_ENV_VARS, env_vars_all)
  if (length(bad) > 0)
    stop("MANUAL_ENV_VARS not found in data: ", paste(bad, collapse = ", "))
  env_vars_input <- MANUAL_ENV_VARS
  message(sprintf("  Variable selection: manual pre-selection (%d vars)",
                  length(env_vars_input)))
}

message(sprintf("  %d rows | presence %d | background %d",
                nrow(cast_data), sum(cast_data$presence == 1L),
                sum(cast_data$presence == 0L)))


# ─────────────────────────────────────────────────────────────────────────────
# 3. VIF collinearity screening
# ─────────────────────────────────────────────────────────────────────────────

message("\n[2/8] VIF screening (threshold = 10)...")

vif_result <- checkpoint(file.path(MODEL_DIR, "vif_result.rds"), {
  if (exists("col_meta")) {
    exclude_for_vif <- unique(c(
      col_meta$column_name[col_meta$role != "environment"],
      setdiff(env_vars_all, env_vars_input)
    ))
  } else {
    exclude_for_vif <- setdiff(names(cast_data), c(env_vars_input, "lon", "lat", "presence"))
  }
  cast_vif(cast_data, exclude = exclude_for_vif, threshold = 10)
})

env_vars_screened <- vif_result$selected
env_vars_screened <- intersect(env_vars_screened, names(cast_data))
if (length(env_vars_screened) < 3)
  stop("VIF retained < 3 variables. Check MANUAL_ENV_VARS or data quality.")

message(sprintf("  VIF: %d -> %d variables retained",
                length(env_vars_input), length(env_vars_screened)))


# ─────────────────────────────────────────────────────────────────────────────
# 4. Data split (cast_prepare)
# ─────────────────────────────────────────────────────────────────────────────

message("\n[3/8] cast_prepare(): data split...")

split <- checkpoint(file.path(MODEL_DIR, "split.rds"), {
  cast_prepare(cast_data, train_fraction = 0.7, seed = SEED,
               env_vars = env_vars_screened)
})

train_data <- split$train
test_data  <- split$test
message(sprintf("  Train %d rows | Test %d rows | Env vars %d",
                nrow(train_data), nrow(test_data), length(split$env_vars)))


# ─────────────────────────────────────────────────────────────────────────────
# 5. Causal DAG learning
# ─────────────────────────────────────────────────────────────────────────────

message("\n[4/8] cast_dag(): learning causal DAG...")

dag_result <- checkpoint(file.path(MODEL_DIR, "dag_result.rds"), {
  cast_dag(
    train_data,
    env_vars            = split$env_vars,
    R                   = DAG_R,
    algorithm           = DAG_ALGORITHM,
    score               = DAG_SCORE,
    strength_threshold  = DAG_STRENGTH_THRESHOLD,
    direction_threshold = DAG_DIRECTION_THRESHOLD,
    max_rows            = DAG_MAX_ROWS,
    seed                = SEED,
    verbose             = TRUE
  )
})

message(sprintf("  DAG: %d nodes, %d edges",
                length(dag_result$nodes), nrow(dag_result$edges)))


# ─────────────────────────────────────────────────────────────────────────────
# 6. ATE causal effect estimation
# ─────────────────────────────────────────────────────────────────────────────

message("\n[5/8] cast_ate(): DML ATE estimation...")

ate_result <- checkpoint(file.path(MODEL_DIR, "ate_result.rds"), {
  cast_ate(train_data,
           variables     = dag_result$nodes,
           K             = 5L,
           quantile_cuts = c(0.25, 0.50, 0.75),
           bonferroni    = TRUE,
           num_trees     = 300L,
           seed          = SEED,
           verbose       = TRUE)
})

n_sig <- sum(ate_result$estimates$significant, na.rm = TRUE)
message(sprintf("  ATE: %d / %d variables significant (Bonferroni)",
                n_sig, nrow(ate_result$estimates)))


# ─────────────────────────────────────────────────────────────────────────────
# 7. Adaptive variable screening + causal role assignment
# ─────────────────────────────────────────────────────────────────────────────

message("\n[6/8] cast_screen() + cast_roles()...")

screened <- checkpoint(file.path(MODEL_DIR, "screened_result.rds"), {
  cast_screen(dag_result, ate_result, train_data, seed = SEED, verbose = TRUE)
})

roles_data <- checkpoint(file.path(MODEL_DIR, "roles_result.rds"), {
  cast_roles(screened, dag_result)
})

message(sprintf("  Screened variables: %d", length(screened$selected)))
message(sprintf("  Variables: %s", paste(screened$selected, collapse = ", ")))


# ─────────────────────────────────────────────────────────────────────────────
# 8. Model fitting
# ─────────────────────────────────────────────────────────────────────────────

message("\n[7/8] cast_fit(): fitting [", paste(MODELS, collapse=", "), "]...")

fit_result <- checkpoint(file.path(MODEL_DIR, "fit_result.rds"), {
  cast_fit(
    train_data,
    screen  = screened,
    dag     = dag_result,
    ate     = ate_result,
    models  = MODELS,
    seed    = SEED,
    verbose = TRUE,
    n_epochs = 400L, n_runs = 5L, patience = 40L,
    tune_grid = TRUE
  )
})

message("  Model fitting complete.")

# Hold-out evaluation
eval_result <- checkpoint(file.path(MODEL_DIR, "eval_result.rds"), {
  cast_evaluate(fit_result, test_data)
})

cat("\n--- Hold-out evaluation metrics ---\n")
{
  m <- eval_result$metrics
  all_metric_cols <- c("model","auc_mean","tss_mean","cbi_mean",
                       "sedi_mean","kappa_mean","prauc_mean")
  present_cols <- intersect(all_metric_cols, names(m))
  print(format(m[, present_cols, drop=FALSE], digits=4), row.names=FALSE)

  if ("auc_mean" %in% names(m)) {
    best_idx <- which.max(m$auc_mean)
    cat(sprintf("\nBest AUC: %s = %.4f\n", m$model[best_idx], m$auc_mean[best_idx]))
  }
}

if (!is.null(eval_result$metrics))
  write.csv(eval_result$metrics,
            file.path(OUT_DIR, "results_summary.csv"), row.names = FALSE)


# ─────────────────────────────────────────────────────────────────────────────
# 9. Spatial cross-validation
# ─────────────────────────────────────────────────────────────────────────────

if (DO_SPATIAL_CV) {
  message("\n[8/8] cast_cv(): spatial CV (", CV_METHOD, ", k=", CV_FOLDS, ")...")

  cv_result <- checkpoint(file.path(MODEL_DIR, "cv_result.rds"), {
    cast_cv(
      cast_data,
      screen       = screened,
      dag          = dag_result,
      ate          = ate_result,
      k            = CV_FOLDS,
      models       = MODELS,
      block_method = CV_METHOD,
      seed         = SEED,
      verbose      = TRUE
    )
  })

  cat("\n--- Spatial CV metrics ---\n")
  {
    m_cv <- cv_result$metrics
    all_cv_cols <- c("model","auc_mean","tss_mean","cbi_mean",
                     "sedi_mean","kappa_mean","prauc_mean")
    present_cv <- intersect(all_cv_cols, names(m_cv))
    print(format(m_cv[, present_cv, drop=FALSE], digits=4), row.names=FALSE)
  }

  if (!is.null(cv_result$metrics))
    write.csv(cv_result$metrics,
              file.path(OUT_DIR, "cv_summary.csv"), row.names = FALSE)
}


# ─────────────────────────────────────────────────────────────────────────────
# 10. Plots (castSDM native plot functions)
# ─────────────────────────────────────────────────────────────────────────────

message("\n[Plots] Saving castSDM figures...")

save_plot <- function(plot_obj, filename, width = 12, height = 8, dpi = FIG_DPI) {
  path <- file.path(FIG_DIR, filename)
  tryCatch({
    ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
    message(sprintf("  [OK] %s", filename))
  }, error = function(e) {
    message(sprintf("  [skip] %s: %s", filename, conditionMessage(e)))
  })
}

# DAG network
if (requireNamespace("ggraph", quietly=TRUE) && requireNamespace("igraph", quietly=TRUE)) {
  p_dag <- plot(dag_result, roles = roles_data, screen = screened)
  save_plot(p_dag, "causal_dag.png", width = 16, height = 12)
}

# ATE forest plot
p_ate <- plot(ate_result)
save_plot(p_ate, "ate_forest_plot.png", width = 10, height = 7)

# Variable screening scores
p_screen <- plot(screened)
save_plot(p_screen, "variable_screening.png", width = 10, height = 7)

# Model evaluation comparison
p_eval <- plot(eval_result)
save_plot(p_eval, "model_evaluation.png", width = 10, height = 6)

# Spatial CV map
if (DO_SPATIAL_CV && exists("cv_result") && !is.null(cv_result)) {
  p_cv <- tryCatch(
    plot(cv_result,
         lon     = cast_data$lon,
         lat     = cast_data$lat,
         metric  = "auc",
         basemap = "world"),
    error = function(e) { message("  [skip] CV map: ", conditionMessage(e)); NULL }
  )
  if (!is.null(p_cv)) save_plot(p_cv, "spatial_cv_map.png", width = 14, height = 8)
}


# ─────────────────────────────────────────────────────────────────────────────
# 11. Spatial prediction (per region) + HSS maps
# ─────────────────────────────────────────────────────────────────────────────

fitted_models <- names(fit_result$models)

if (!is.null(PREDICT_REGIONS) && length(PREDICT_REGIONS) > 0) {
  message("\n[Spatial prediction] Regions: ", paste(PREDICT_REGIONS, collapse=", "))

  for (reg in PREDICT_REGIONS) {
    out_pred     <- file.path(PRED_DIR, paste0("Lota_lota_HSS_", reg, ".csv"))
    out_pred_rds <- file.path(PRED_DIR, paste0("Lota_lota_pred_", reg, ".rds"))

    if (!file.exists(out_pred)) {
      env_grid_file <- file.path(OUT_DIR, paste0("RiverATLAS_env_", reg, ".csv"))
      if (!file.exists(env_grid_file)) {
        message(sprintf("  [skip] %s: prediction grid not found", reg))
        next
      }

      message(sprintf("  Predicting region '%s'...", reg))
      t0 <- proc.time()
      env_grid <- data.table::fread(env_grid_file, data.table = FALSE)

      pred <- tryCatch(
        cast_predict(fit_result, new_data = env_grid),
        error = function(e) {
          message(sprintf("  [fail] %s: %s", reg, conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(pred)) {
        pred_df <- if (inherits(pred, "cast_predict")) pred$predictions else pred
        data.table::fwrite(pred_df, out_pred)
        saveRDS(pred, out_pred_rds)
        elapsed <- (proc.time() - t0)[["elapsed"]]
        message(sprintf("  [OK] %s -> %s (%s rows, %.1f sec)",
                        reg, basename(out_pred),
                        format(nrow(pred_df), big.mark=","), elapsed))
      }
    } else {
      message(sprintf("  [skip] %s: prediction file exists", reg))
    }

    # HSS maps per model
    pred_obj <- if (file.exists(out_pred_rds)) {
      readRDS(out_pred_rds)
    } else if (file.exists(out_pred)) {
      pred_df <- data.table::fread(out_pred, data.table = FALSE)
      hss_models <- sub("^HSS_", "", names(pred_df)[startsWith(names(pred_df), "HSS_")])
      new_cast_predict(predictions = pred_df, models = hss_models)
    } else {
      NULL
    }

    if (!is.null(pred_obj) && requireNamespace("sf", quietly = TRUE)) {
      for (mdl in fitted_models) {
        fig_path <- file.path(FIG_DIR, sprintf("HSS_%s_%s.png", reg, mdl))
        if (file.exists(fig_path)) next

        p_hss <- tryCatch(
          plot(pred_obj,
               model   = mdl,
               basemap = "world",
               title   = sprintf("Lota lota HSS - %s (%s)", toupper(reg), mdl)),
          error = function(e) {
            message(sprintf("  [skip] HSS_%s_%s: %s", reg, mdl, conditionMessage(e)))
            NULL
          }
        )
        if (!is.null(p_hss))
          save_plot(p_hss, sprintf("HSS_%s_%s.png", reg, mdl), width = 14, height = 8)
      }
    }
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# 12. CATE spatial heterogeneous causal effects
# ─────────────────────────────────────────────────────────────────────────────

message("\n[Optional] cast_cate(): spatial CATE (requires grf)...")

cate_rds <- file.path(MODEL_DIR, "cate_result.rds")

if (file.exists(cate_rds)) {
  message("  [skip] CATE checkpoint exists")
  cate_result <- readRDS(cate_rds)
} else {
  cate_grid_list <- list()
  for (reg in PREDICT_REGIONS) {
    f <- file.path(OUT_DIR, paste0("RiverATLAS_env_", reg, ".csv"))
    if (file.exists(f)) {
      g <- data.table::fread(f, data.table = FALSE)
      g$region <- reg
      cate_grid_list[[reg]] <- g
    }
  }
  first_region_grid <- if (length(cate_grid_list) > 0) {
    message(sprintf("  Merging %d region grids for CATE prediction (%s)",
                    length(cate_grid_list), paste(names(cate_grid_list), collapse=", ")))
    do.call(rbind, cate_grid_list)
  } else NULL

  cate_result <- tryCatch({
    if (!requireNamespace("grf", quietly = TRUE))
      stop("grf package required: install.packages('grf')")
    if (is.null(first_region_grid))
      stop("No prediction grid available, skipping CATE")

    cate_env <- screened$selected
    cate_train_slim <- train_data[,
      intersect(c("lon", "lat", "presence", cate_env), names(train_data)),
      drop = FALSE
    ]
    cate_pred_slim <- first_region_grid[,
      intersect(c("lon", "lat", cate_env), names(first_region_grid)),
      drop = FALSE
    ]

    cast_cate(
      data         = cate_train_slim,
      ate          = ate_result,
      screen       = screened,
      predict_data = cate_pred_slim,
      top_n        = 3L,
      n_trees      = 500L,
      seed         = SEED,
      verbose      = TRUE
    )
  }, error = function(e) {
    message("  CATE skipped: ", conditionMessage(e))
    NULL
  })

  if (!is.null(cate_result)) {
    saveRDS(cate_result, cate_rds)
    message("  CATE done, variables: ", paste(cate_result$variables, collapse=", "))
  }
}

if (!is.null(cate_result) && requireNamespace("sf", quietly = TRUE)) {
  message("  Plotting CATE maps...")
  for (cv in cate_result$variables) {
    fig_name <- paste0("cate_", cv, ".png")
    if (file.exists(file.path(FIG_DIR, fig_name))) next

    p_cate <- tryCatch(
      plot(cate_result,
           variable          = cv,
           species           = "Lota_lota",
           basemap           = "world",
           legend_position   = "bottom"),
      error = function(e) {
        message(sprintf("  [skip] CATE %s: %s", cv, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(p_cate))
      save_plot(p_cate, fig_name, width = 14, height = 8)
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# 13. Global merged maps (HSS + CATE)
# ─────────────────────────────────────────────────────────────────────────────

if (requireNamespace("torch", quietly = TRUE) &&
    exists("fit_result") && !is.null(fit_result)) {
  fit_result <- NULL
  invisible(gc(verbose = FALSE))
}

message("\n[Global maps] Merging regional predictions...")

hss_parts <- lapply(PREDICT_REGIONS, function(reg) {
  f <- file.path(PRED_DIR, paste0("Lota_lota_HSS_", reg, ".csv"))
  if (!file.exists(f)) return(NULL)
  data.table::fread(f, data.table = FALSE)
})
hss_parts <- Filter(Negate(is.null), hss_parts)

if (length(hss_parts) > 0 && requireNamespace("sf", quietly = TRUE)) {
  global_pred_df  <- do.call(rbind, hss_parts)
  hss_cols        <- names(global_pred_df)[startsWith(names(global_pred_df), "HSS_")]
  hss_models_glob <- sub("^HSS_", "", hss_cols)
  global_pred_obj <- new_cast_predict(predictions = global_pred_df,
                                      models      = hss_models_glob)

  for (mdl in hss_models_glob) {
    fig_name <- sprintf("HSS_GLOBAL_%s.png", mdl)
    p_glob <- tryCatch(
      plot(global_pred_obj,
           model   = mdl,
           basemap = "world",
           title   = sprintf("Lota lota - Global HSS (%s)", toupper(mdl))),
      error = function(e) {
        message(sprintf("  [skip] %s: %s", fig_name, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(p_glob))
      save_plot(p_glob, fig_name, width = 20, height = 10, dpi = FIG_DPI_GLOBAL)
  }
}

if (!is.null(cate_result) && requireNamespace("sf", quietly = TRUE)) {
  for (cv in cate_result$variables) {
    fig_name <- sprintf("CATE_GLOBAL_%s.png", cv)
    p_cate_glob <- tryCatch(
      plot(cate_result,
           variable        = cv,
           species         = "Lota_lota",
           basemap         = "world",
           legend_position = "bottom"),
      error = function(e) {
        message(sprintf("  [skip] %s: %s", fig_name, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(p_cate_glob))
      save_plot(p_cate_glob, fig_name, width = 20, height = 10, dpi = FIG_DPI_GLOBAL)
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# 14. Summary
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n")
cat("castSDM modeling complete: Lota lota (Burbot)\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("  Models          : %s\n", paste(MODELS, collapse=", ")))
cat(sprintf("  Active env vars : %d (after VIF + cast_screen)\n",
            length(screened$selected)))
cat(sprintf("  DAG edges       : %d\n", nrow(dag_result$edges)))
cat(sprintf("  Significant ATE : %d / %d\n",
            sum(ate_result$estimates$significant, na.rm=TRUE),
            nrow(ate_result$estimates)))
cat(sprintf("  Output dir      : %s\n\n", OUT_DIR))

if (DO_SPATIAL_CV && exists("cv_result") && !is.null(cv_result) &&
    !is.null(cv_result$metrics)) {
  cat(sprintf("-- Spatial CV (%d-fold, %s) --\n", cv_result$k, cv_result$block_method))
  for (i in seq_len(nrow(cv_result$metrics))) {
    r <- cv_result$metrics[i, ]
    cat(sprintf("  %s: AUC=%.3f  TSS=%.3f  CBI=%.3f  SEDI=%.3f\n",
                r$model, r$auc_mean, r$tss_mean,
                r$cbi_mean, r$sedi_mean))
  }
}

cat("\nCheckpoint files (delete to recompute):\n")
for (f in c("vif_result.rds","split.rds","dag_result.rds",
            "ate_result.rds","screened_result.rds","roles_result.rds",
            "fit_result.rds","eval_result.rds","cv_result.rds")) {
  mark <- if (file.exists(file.path(MODEL_DIR, f))) "+" else "-"
  cat(sprintf("  [%s] %s\n", mark, f))
}
