# ==============================================================================
# castSDM Freshwater Fish Example: Lota lota (Burbot / 江鳕)
#
# This script demonstrates the COMPLETE workflow for freshwater fish species
# distribution modeling with castSDM, starting from raw GBIF occurrence data.
#
# ── INPUT ────────────────────────────────────────────────────────────────────
#   1. A GBIF-style occurrence CSV with 4 columns:
#        taxon, longitude, latitude, occurrenceStatus
#   2. RiverATLAS v1.0 regional shapefiles (one per study region)
#      Download from: https://www.hydrosheds.org/products/riveratlas
#
# ── WORKFLOW ─────────────────────────────────────────────────────────────────
#   Phase 1 — Data Preparation (from raw CSV):
#     Load occurrences → Load RiverATLAS shapefiles → Snap to reaches
#     → Sample background → Assemble training data → Save checkpoint
#
#   Phase 2 — castSDM Modeling:
#     VIF → cast_prepare → cast_dag → cast_ate → cast_screen
#     → cast_roles → cast_fit → cast_evaluate → cast_cv → cast_predict
#     → cast_cate → cast_consistency → Plotting & Export
#
# ── CHECKPOINT SYSTEM ────────────────────────────────────────────────────────
#   Each step saves results to MODEL_DIR/*.rds (or .csv).
#   Re-running auto-skips completed steps; delete a file to re-compute.
#
# Date: 2026-04-10
# ==============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 0. CONFIGURATION — edit this section for your species / data
# ─────────────────────────────────────────────────────────────────────────────

# ── Raw occurrence data ──────────────────────────────────────────────────────
# Your GBIF CSV must have: taxon, longitude, latitude, occurrenceStatus
OCC_CSV  <- "E:/Package/testFish/Lota_lota.csv"

# ── RiverATLAS shapefiles ───────────────────────────────────────────────────
# Each region file: RiverATLAS_v10_<region>.shp
ATLAS_DIR <- "E:/Package/testFish/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp"

# ── Study regions (RiverATLAS codes) ─────────────────────────────────────────
# eu=Europe, na=North America, si=Siberia, ar=Arctic, as=Asia
TARGET_REGIONS <- c("eu", "na", "si", "ar")

# ── Output directories ──────────────────────────────────────────────────────
OUT_DIR      <- "output"
MODEL_DIR    <- file.path(OUT_DIR, "model")
PRED_DIR     <- file.path(OUT_DIR, "predictions")
FIG_DIR      <- file.path(OUT_DIR, "figures")

# ── Background sampling ─────────────────────────────────────────────────────
BG_PER_REGION   <- 20000L
BG_TOTAL_CAP    <- 100000L
MAX_SNAP_DIST   <- 0.5     # degrees (~55 km)

# ── Model selection ──────────────────────────────────────────────────────────
MODELS <- c("cast", "rf", "maxent", "brt")

# ── Spatial cross-validation ─────────────────────────────────────────────────
DO_SPATIAL_CV <- TRUE
CV_FOLDS      <- 5L
CV_METHOD     <- "grid"

# ── Spatial prediction regions (set NULL to skip) ────────────────────────────
PREDICT_REGIONS <- c("eu", "na", "si", "ar")

# ── Manual environment variable pre-selection (before VIF) ───────────────────
# These are RiverATLAS attribute names. Set NULL to use ALL env vars.
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

# ── RiverATLAS non-environmental columns (topology / ID) ─────────────────────
NON_ENV_COLS <- c(
  "hyriv_id", "next_down", "main_riv", "length_km", "dist_dn_km",
  "dist_up_km", "order_", "class", "endorheic", "ord_stra",
  "ord_clas", "ord_flow", "hybas_l12", "upland_skm",
  "region", "species", "occ_lon", "occ_lat", "snap_dist"
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

# ── RiverATLAS nodata sentinels ──────────────────────────────────────────────
NODATA_VALS <- c(-9999L, -9998L, -999L, -1L)


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

if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(castSDM)
}

for (pkg in c("data.table", "dplyr", "ggplot2", "sf", "FNN")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

set.seed(SEED)
cat("=============================================================\n")
cat("castSDM Freshwater Fish Modeling: Lota lota (Burbot)\n")
cat("=============================================================\n")


# ==========================================================================
#  PHASE 1: DATA PREPARATION  (raw GBIF CSV → castSDM training data)
# ==========================================================================

TRAIN_CSV <- file.path(OUT_DIR, "CAST_Lota_lota.csv")

if (file.exists(TRAIN_CSV)) {
  # ── Fast path: training data already prepared ─────────────────────────────
  message("\n[Phase 1] Training data found, loading checkpoint...")
  cast_data <- data.table::fread(TRAIN_CSV, data.table = FALSE)

  # Detect env vars
  meta_csv <- file.path(OUT_DIR, "column_metadata.csv")
  if (file.exists(meta_csv)) {
    col_meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
    env_vars_all <- col_meta$column_name[col_meta$is_env_var]
  } else {
    env_vars_all <- get_env_vars(cast_data)
  }

} else {
  # ── Full data preparation from raw occurrence CSV ─────────────────────────
  message("\n[Phase 1] Preparing training data from raw occurrences...")

  # ── 1a. Load and clean occurrence records ────────────────────────────────
  message("\n  [1/5] Loading occurrence data: ", OCC_CSV)
  if (!file.exists(OCC_CSV))
    stop("Occurrence CSV not found: ", OCC_CSV)

  occ_raw <- data.table::fread(OCC_CSV, data.table = FALSE)

  occ <- occ_raw |>
    dplyr::filter(occurrenceStatus == "PRESENT") |>
    dplyr::rename(lon = longitude, lat = latitude) |>
    dplyr::filter(!is.na(lon), !is.na(lat),
                  lon >= -180, lon <= 180,
                  lat >= -90,  lat <= 90) |>
    dplyr::distinct(lon, lat, .keep_all = TRUE) |>
    dplyr::mutate(species = sub(" ", "_", taxon)) |>
    dplyr::select(species, lon, lat)

  message(sprintf("    Raw: %d records | Clean: %d unique points",
                  nrow(occ_raw), nrow(occ)))

  # ── 1b. Load RiverATLAS shapefiles ──────────────────────────────────────
  message(sprintf("\n  [2/5] Loading RiverATLAS regions: %s",
                  paste(TARGET_REGIONS, collapse = ", ")))

  load_river_region <- function(region_code) {
    shp_path <- file.path(ATLAS_DIR,
                          paste0("RiverATLAS_v10_", region_code, ".shp"))
    if (!file.exists(shp_path)) {
      warning("Shapefile not found, skipping: ", shp_path)
      return(NULL)
    }
    message(sprintf("    Loading '%s' ...", region_code))
    t0 <- proc.time()
    riv_sf <- sf::st_read(shp_path, quiet = TRUE)
    suppressWarnings(
      centroids <- sf::st_centroid(sf::st_geometry(riv_sf))
    )
    coords <- sf::st_coordinates(centroids)
    riv_tbl <- sf::st_drop_geometry(riv_sf)
    names(riv_tbl) <- tolower(names(riv_tbl))
    riv_tbl$lon    <- coords[, 1]
    riv_tbl$lat    <- coords[, 2]
    riv_tbl$region <- region_code
    elapsed <- (proc.time() - t0)[["elapsed"]]
    message(sprintf("      %s reaches | %.1f sec",
                    format(nrow(riv_tbl), big.mark = ","), elapsed))
    riv_tbl
  }

  river_list <- lapply(TARGET_REGIONS, load_river_region)
  names(river_list) <- TARGET_REGIONS
  river_list <- Filter(Negate(is.null), river_list)
  if (length(river_list) == 0)
    stop("No RiverATLAS regions loaded. Check ATLAS_DIR.")

  river_all <- dplyr::bind_rows(river_list)
  message(sprintf("    Combined: %s reaches",
                  format(nrow(river_all), big.mark = ",")))

  # ── 1c. Snap occurrences to nearest river reach ─────────────────────────
  message("\n  [3/5] Snapping occurrences to nearest reach...")

  nn_result <- FNN::get.knnx(
    data  = as.matrix(river_all[, c("lon", "lat")]),
    query = as.matrix(occ[, c("lon", "lat")]),
    k     = 1L
  )
  occ$nn_row  <- nn_result$nn.index[, 1]
  occ$nn_dist <- nn_result$nn.dist[, 1]

  occ_ok <- dplyr::filter(occ, nn_dist <= MAX_SNAP_DIST)
  message(sprintf("    Kept: %d / %d (dropped %d > %.2f deg)",
                  nrow(occ_ok), nrow(occ),
                  nrow(occ) - nrow(occ_ok), MAX_SNAP_DIST))

  matched <- river_all[occ_ok$nn_row, ]
  matched$occ_lon   <- occ_ok$lon
  matched$occ_lat   <- occ_ok$lat
  matched$snap_dist <- occ_ok$nn_dist
  matched$presence  <- 1L
  matched$species   <- occ_ok$species

  occ_unique <- matched |>
    dplyr::arrange(snap_dist) |>
    dplyr::distinct(hyriv_id, .keep_all = TRUE)

  message(sprintf("    Unique presence reaches: %d", nrow(occ_unique)))

  # ── 1d. Sample background reaches ───────────────────────────────────────
  message("\n  [4/5] Sampling background reaches...")
  set.seed(SEED)
  pres_ids <- occ_unique$hyriv_id

  bg_list <- lapply(names(river_list), function(reg) {
    reg_data <- river_all[river_all$region == reg, ]
    pool     <- reg_data[!(reg_data$hyriv_id %in% pres_ids), ]
    n_sample <- min(BG_PER_REGION, nrow(pool))
    if (n_sample == 0L) return(NULL)
    pool[sample(nrow(pool), n_sample, replace = FALSE), ]
  })
  bg_data <- dplyr::bind_rows(bg_list)

  if (nrow(bg_data) > BG_TOTAL_CAP)
    bg_data <- bg_data[sample(nrow(bg_data), BG_TOTAL_CAP), ]

  bg_data$presence  <- 0L
  bg_data$occ_lon   <- NA_real_
  bg_data$occ_lat   <- NA_real_
  bg_data$snap_dist <- NA_real_
  bg_data$species   <- unique(occ$species)[1]

  message(sprintf("    Background: %s reaches",
                  format(nrow(bg_data), big.mark = ",")))

  # ── 1e. Assemble training dataset ───────────────────────────────────────
  message("\n  [5/5] Assembling training dataset...")

  meta_cols <- c("hyriv_id", "lon", "lat", "presence", "species", "region",
                 "occ_lon", "occ_lat", "snap_dist")
  all_cols  <- union(names(occ_unique), names(bg_data))
  env_cols  <- setdiff(all_cols, c(meta_cols, tolower(NON_ENV_COLS)))
  topo_cols <- intersect(tolower(NON_ENV_COLS), all_cols)
  col_order <- unique(c(meta_cols, topo_cols, env_cols))

  align_cols <- function(df, target_cols) {
    for (m in setdiff(target_cols, names(df))) df[[m]] <- NA
    df[, intersect(target_cols, names(df))]
  }

  cast_data <- dplyr::bind_rows(
    align_cols(occ_unique, col_order),
    align_cols(bg_data,    col_order)
  )
  cast_data$presence <- as.integer(cast_data$presence)

  # Replace RiverATLAS nodata sentinels
  num_cols <- names(cast_data)[
    vapply(cast_data, is.numeric, logical(1)) &
    !(names(cast_data) %in% c("lon", "lat", "presence", meta_cols, topo_cols))
  ]
  for (col in num_cols)
    cast_data[[col]][cast_data[[col]] %in% NODATA_VALS] <- NA

  # Save training data
  data.table::fwrite(cast_data, TRAIN_CSV)
  message(sprintf("    Saved: %s (%d presence + %d background = %d rows)",
                  TRAIN_CSV, sum(cast_data$presence == 1L),
                  sum(cast_data$presence == 0L), nrow(cast_data)))

  # Save column metadata
  determine_role <- function(cn) {
    if (cn %in% c("lon", "lat"))                                 return("spatial")
    if (cn == "presence")                                        return("response")
    if (cn %in% c("hyriv_id", "region", "species"))              return("identifier")
    if (cn %in% c("occ_lon", "occ_lat", "snap_dist"))            return("qc_metadata")
    if (cn %in% intersect(tolower(NON_ENV_COLS), names(cast_data))) return("network_topo")
    "environment"
  }
  col_meta <- data.frame(
    column_name = names(cast_data),
    role        = vapply(names(cast_data), determine_role, character(1)),
    stringsAsFactors = FALSE
  )
  col_meta$is_env_var <- col_meta$role == "environment"
  write.csv(col_meta, file.path(OUT_DIR, "column_metadata.csv"), row.names = FALSE)

  env_vars_all <- col_meta$column_name[col_meta$is_env_var]

  # Save per-region env grids for spatial prediction
  message("    Saving per-region prediction grids...")
  for (reg in names(river_list)) {
    out_grid <- file.path(OUT_DIR, paste0("RiverATLAS_env_", reg, ".csv"))
    if (!file.exists(out_grid)) {
      reg_tbl <- river_list[[reg]]
      names(reg_tbl) <- tolower(names(reg_tbl))
      for (col in names(reg_tbl)) {
        if (is.numeric(reg_tbl[[col]]))
          reg_tbl[[col]][reg_tbl[[col]] %in% NODATA_VALS] <- NA
      }
      data.table::fwrite(reg_tbl, out_grid)
      message(sprintf("      [OK] %s (%s rows)",
                      basename(out_grid),
                      format(nrow(reg_tbl), big.mark = ",")))
    }
  }

  # Release large objects
  rm(river_list, river_all, bg_data, occ_unique, matched, occ_ok, occ)
  invisible(gc(verbose = FALSE))
}

message(sprintf("\n  Training data: %d rows | presence %d | background %d",
                nrow(cast_data), sum(cast_data$presence == 1L),
                sum(cast_data$presence == 0L)))


# ==========================================================================
#  PHASE 2: castSDM MODELING
# ==========================================================================

# ── Select environment variables ─────────────────────────────────────────────
if (is.null(MANUAL_ENV_VARS)) {
  env_vars_input <- env_vars_all
  message("  Variable selection: all environment variables")
} else {
  env_vars_input <- intersect(MANUAL_ENV_VARS, env_vars_all)
  if (length(env_vars_input) < 3)
    stop("Fewer than 3 MANUAL_ENV_VARS found in data. Check variable names.")
  message(sprintf("  Variable selection: manual pre-selection (%d vars)",
                  length(env_vars_input)))
}


# ─────────────────────────────────────────────────────────────────────────────
# 2. VIF collinearity screening
# ─────────────────────────────────────────────────────────────────────────────

message("\n[2/9] VIF screening (threshold = 10)...")

vif_result <- checkpoint(file.path(MODEL_DIR, "vif_result.rds"), {
  exclude_for_vif <- setdiff(
    names(cast_data),
    c(env_vars_input, "lon", "lat", "presence")
  )
  cast_vif(cast_data, exclude = exclude_for_vif, threshold = 10)
})

env_vars_screened <- intersect(vif_result$selected, names(cast_data))
if (length(env_vars_screened) < 3)
  stop("VIF retained < 3 variables. Relax threshold or check data.")

message(sprintf("  VIF: %d -> %d variables retained",
                length(env_vars_input), length(env_vars_screened)))


# ─────────────────────────────────────────────────────────────────────────────
# 3. Data split
# ─────────────────────────────────────────────────────────────────────────────

message("\n[3/9] cast_prepare(): data split...")

split <- checkpoint(file.path(MODEL_DIR, "split.rds"), {
  cast_prepare(cast_data, train_fraction = 0.7, seed = SEED,
               env_vars = env_vars_screened)
})

train_data <- split$train
test_data  <- split$test
message(sprintf("  Train %d | Test %d | Env vars %d",
                nrow(train_data), nrow(test_data), length(split$env_vars)))


# ─────────────────────────────────────────────────────────────────────────────
# 4. Causal DAG learning
# ─────────────────────────────────────────────────────────────────────────────

message("\n[4/9] cast_dag(): learning causal DAG...")

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
# 5. ATE causal effect estimation
# ─────────────────────────────────────────────────────────────────────────────

message("\n[5/9] cast_ate(): DML ATE estimation...")

ate_result <- checkpoint(file.path(MODEL_DIR, "ate_result.rds"), {
  cast_ate(train_data,
           variables     = dag_result$nodes,
           K             = 5L,
           quantile_cuts = c(0.25, 0.50, 0.75),
           bonferroni    = TRUE,
           num_trees     = 300L,
           parallel      = TRUE,
           seed          = SEED,
           verbose       = TRUE)
})

n_sig <- sum(ate_result$estimates$significant, na.rm = TRUE)
message(sprintf("  ATE: %d / %d variables significant (Bonferroni)",
                n_sig, nrow(ate_result$estimates)))


# ─────────────────────────────────────────────────────────────────────────────
# 6. Adaptive variable screening + causal role assignment
# ─────────────────────────────────────────────────────────────────────────────

message("\n[6/9] cast_screen() + cast_roles()...")

screened <- checkpoint(file.path(MODEL_DIR, "screened_result.rds"), {
  cast_screen(dag_result, ate_result, train_data, seed = SEED, verbose = TRUE)
})

roles_data <- checkpoint(file.path(MODEL_DIR, "roles_result.rds"), {
  cast_roles(screened, dag_result)
})

message(sprintf("  Screened variables: %d (%s)",
                length(screened$selected),
                paste(screened$selected, collapse = ", ")))


# ─────────────────────────────────────────────────────────────────────────────
# 7. Model fitting
# ─────────────────────────────────────────────────────────────────────────────

message("\n[7/9] cast_fit(): fitting [", paste(MODELS, collapse=", "), "]...")

fit_result <- checkpoint(file.path(MODEL_DIR, "fit_result.rds"), {
  cast_fit(
    train_data,
    screen    = screened,
    dag       = dag_result,
    ate       = ate_result,
    models    = MODELS,
    seed      = SEED,
    verbose   = TRUE,
    n_epochs  = 400L,
    n_runs    = 5L,
    patience  = 40L,
    tune_grid = TRUE
  )
})

message("  Model fitting complete.")

eval_result <- checkpoint(file.path(MODEL_DIR, "eval_result.rds"), {
  cast_evaluate(fit_result, test_data)
})

cat("\n--- Hold-out evaluation ---\n")
print(format(eval_result$metrics, digits = 4), row.names = FALSE)

if (!is.null(eval_result$metrics))
  write.csv(eval_result$metrics,
            file.path(OUT_DIR, "results_summary.csv"), row.names = FALSE)


# ─────────────────────────────────────────────────────────────────────────────
# 8. Spatial cross-validation
# ─────────────────────────────────────────────────────────────────────────────

if (DO_SPATIAL_CV) {
  message("\n[8/9] cast_cv(): spatial CV (", CV_METHOD, ", k=", CV_FOLDS, ")...")

  cv_result <- checkpoint(file.path(MODEL_DIR, "cv_result.rds"), {
    cast_cv(
      cast_data,
      screen       = screened,
      dag          = dag_result,
      ate          = ate_result,
      k            = CV_FOLDS,
      models       = MODELS,
      block_method = CV_METHOD,
      parallel     = TRUE,
      seed         = SEED,
      verbose      = TRUE
    )
  })

  cat("\n--- Spatial CV ---\n")
  print(format(cv_result$metrics, digits = 4), row.names = FALSE)

  if (!is.null(cv_result$metrics))
    write.csv(cv_result$metrics,
              file.path(OUT_DIR, "cv_summary.csv"), row.names = FALSE)
}


# ─────────────────────────────────────────────────────────────────────────────
# 9. Plots (all castSDM native plot functions)
# ─────────────────────────────────────────────────────────────────────────────

message("\n[9/9] Saving castSDM figures...")

save_plot <- function(plot_obj, filename, width = 12, height = 8,
                      dpi = FIG_DPI) {
  path <- file.path(FIG_DIR, filename)
  tryCatch({
    ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
    message(sprintf("  [OK] %s", filename))
  }, error = function(e) {
    message(sprintf("  [skip] %s: %s", filename, conditionMessage(e)))
  })
}

# 9a. DAG network (nodes = roles, edges = bootstrap strength)
if (requireNamespace("ggraph", quietly = TRUE) &&
    requireNamespace("igraph", quietly = TRUE)) {
  p_dag <- plot(dag_result, roles = roles_data, screen = screened,
                species = "Lota_lota")
  save_plot(p_dag, "causal_dag.png", width = 16, height = 12)
}

# 9b. ATE forest plot
p_ate <- plot(ate_result)
save_plot(p_ate, "ate_forest_plot.png", width = 10, height = 7)

# 9c. Variable screening scores
p_screen <- plot(screened)
save_plot(p_screen, "variable_screening.png", width = 10, height = 7)

# 9d. Model evaluation comparison (6 metrics)
p_eval <- plot(eval_result)
save_plot(p_eval, "model_evaluation.png", width = 10, height = 6)

# 9e. Spatial CV map
if (DO_SPATIAL_CV && exists("cv_result") && !is.null(cv_result)) {
  p_cv <- tryCatch(
    plot(cv_result,
         lon     = cast_data$lon,
         lat     = cast_data$lat,
         metric  = "auc",
         basemap = "world"),
    error = function(e) {
      message("  [skip] CV map: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p_cv))
    save_plot(p_cv, "spatial_cv_map.png", width = 14, height = 8)
}


# ─────────────────────────────────────────────────────────────────────────────
# 10. Spatial prediction (per region) + HSS maps
# ─────────────────────────────────────────────────────────────────────────────

fitted_models <- names(fit_result$models)

if (!is.null(PREDICT_REGIONS) && length(PREDICT_REGIONS) > 0) {
  message("\n[Spatial prediction] Regions: ",
          paste(PREDICT_REGIONS, collapse = ", "))

  for (reg in PREDICT_REGIONS) {
    out_pred     <- file.path(PRED_DIR,
                              paste0("Lota_lota_HSS_", reg, ".csv"))
    out_pred_rds <- file.path(PRED_DIR,
                              paste0("Lota_lota_pred_", reg, ".rds"))

    if (!file.exists(out_pred)) {
      env_grid_file <- file.path(OUT_DIR,
                                 paste0("RiverATLAS_env_", reg, ".csv"))
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
                        format(nrow(pred_df), big.mark = ","), elapsed))
      }
    } else {
      message(sprintf("  [skip] %s: prediction file exists", reg))
    }

    # HSS maps per model
    pred_obj <- if (file.exists(out_pred_rds)) {
      readRDS(out_pred_rds)
    } else if (file.exists(out_pred)) {
      pred_df <- data.table::fread(out_pred, data.table = FALSE)
      hss_models <- sub("^HSS_", "",
                        names(pred_df)[startsWith(names(pred_df), "HSS_")])
      new_cast_predict(predictions = pred_df, models = hss_models)
    } else {
      NULL
    }

    if (!is.null(pred_obj) && requireNamespace("sf", quietly = TRUE)) {
      for (mdl in fitted_models) {
        fig_path <- file.path(FIG_DIR,
                              sprintf("HSS_%s_%s.png", reg, mdl))
        if (file.exists(fig_path)) next
        p_hss <- tryCatch(
          plot(pred_obj, model = mdl, basemap = "world",
               title = sprintf("Lota lota HSS - %s (%s)",
                               toupper(reg), mdl)),
          error = function(e) {
            message(sprintf("  [skip] HSS_%s_%s: %s",
                            reg, mdl, conditionMessage(e)))
            NULL
          }
        )
        if (!is.null(p_hss))
          save_plot(p_hss, sprintf("HSS_%s_%s.png", reg, mdl),
                    width = 14, height = 8)
      }
    }
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# 11. CATE spatial heterogeneous causal effects
# ─────────────────────────────────────────────────────────────────────────────

message("\n[CATE] cast_cate(): spatial CATE (requires grf)...")

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
    do.call(rbind, cate_grid_list)
  } else NULL

  cate_result <- tryCatch({
    if (!requireNamespace("grf", quietly = TRUE))
      stop("grf package required: install.packages('grf')")
    if (is.null(first_region_grid))
      stop("No prediction grid available")

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
    message("  CATE done: ", paste(cate_result$variables, collapse = ", "))
  }
}

if (!is.null(cate_result) && requireNamespace("sf", quietly = TRUE)) {
  for (cv in cate_result$variables) {
    fig_name <- paste0("cate_", cv, ".png")
    if (file.exists(file.path(FIG_DIR, fig_name))) next
    p_cate <- tryCatch(
      plot(cate_result, variable = cv, species = "Lota_lota",
           basemap = "world", legend_position = "bottom"),
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
# 12. Global merged maps + Inter-model Consistency
# ─────────────────────────────────────────────────────────────────────────────

# Release torch objects before heavy plotting
if (requireNamespace("torch", quietly = TRUE) &&
    exists("fit_result") && !is.null(fit_result)) {
  fit_result <- NULL
  invisible(gc(verbose = FALSE))
}

message("\n[Global] Merging regional predictions...")

# ── 12a. Global HSS maps ─────────────────────────────────────────────────────
hss_parts <- lapply(PREDICT_REGIONS, function(reg) {
  f <- file.path(PRED_DIR, paste0("Lota_lota_HSS_", reg, ".csv"))
  if (!file.exists(f)) return(NULL)
  data.table::fread(f, data.table = FALSE)
})
hss_parts <- Filter(Negate(is.null), hss_parts)

global_pred_obj <- NULL

if (length(hss_parts) > 0 && requireNamespace("sf", quietly = TRUE)) {
  global_pred_df  <- do.call(rbind, hss_parts)
  hss_cols        <- names(global_pred_df)[startsWith(names(global_pred_df),
                                                       "HSS_")]
  hss_models_glob <- sub("^HSS_", "", hss_cols)
  global_pred_obj <- new_cast_predict(predictions = global_pred_df,
                                      models      = hss_models_glob)

  for (mdl in hss_models_glob) {
    fig_name <- sprintf("HSS_GLOBAL_%s.png", mdl)
    p_glob <- tryCatch(
      plot(global_pred_obj, model = mdl, basemap = "world",
           title = sprintf("Lota lota - Global HSS (%s)", toupper(mdl))),
      error = function(e) {
        message(sprintf("  [skip] %s: %s", fig_name, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(p_glob))
      save_plot(p_glob, fig_name, width = 20, height = 10,
                dpi = FIG_DPI_GLOBAL)
  }
}

# ── 12b. Global CATE maps ────────────────────────────────────────────────────
if (!is.null(cate_result) && requireNamespace("sf", quietly = TRUE)) {
  for (cv in cate_result$variables) {
    fig_name <- sprintf("CATE_GLOBAL_%s.png", cv)
    p_cate_glob <- tryCatch(
      plot(cate_result, variable = cv, species = "Lota_lota",
           basemap = "world", legend_position = "bottom"),
      error = function(e) {
        message(sprintf("  [skip] %s: %s", fig_name, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(p_cate_glob))
      save_plot(p_cate_glob, fig_name, width = 20, height = 10,
                dpi = FIG_DPI_GLOBAL)
  }
}

# ── 12c. Inter-model Spatial Consistency Heatmaps ─────────────────────────────
if (!is.null(global_pred_obj) && length(global_pred_obj$models) >= 2 &&
    requireNamespace("patchwork", quietly = TRUE)) {
  message("\n[Consistency] Computing inter-model spatial consistency...")

  global_cons <- tryCatch(
    cast_consistency(global_pred_obj),
    error = function(e) {
      message("  [skip] Consistency: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(global_cons)) {
    print(global_cons)
    p_cons <- plot(global_cons, species = "Lota_lota")
    save_plot(p_cons, "model_consistency_GLOBAL.png",
              width = 16, height = 5.5, dpi = FIG_DPI_GLOBAL)
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# 13. Summary
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n")
cat("castSDM modeling complete: Lota lota (Burbot)\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("  Models          : %s\n", paste(MODELS, collapse = ", ")))
cat(sprintf("  Active env vars : %d (after VIF + cast_screen)\n",
            length(screened$selected)))
cat(sprintf("  DAG edges       : %d\n", nrow(dag_result$edges)))
cat(sprintf("  Significant ATE : %d / %d\n",
            sum(ate_result$estimates$significant, na.rm = TRUE),
            nrow(ate_result$estimates)))
cat(sprintf("  Output dir      : %s\n\n", OUT_DIR))

if (DO_SPATIAL_CV && exists("cv_result") && !is.null(cv_result) &&
    !is.null(cv_result$metrics)) {
  cat(sprintf("-- Spatial CV (%d-fold, %s) --\n",
              cv_result$k, cv_result$block_method))
  for (i in seq_len(nrow(cv_result$metrics))) {
    r <- cv_result$metrics[i, ]
    cat(sprintf("  %s: AUC=%.3f  TSS=%.3f  CBI=%.3f  SEDI=%.3f\n",
                r$model, r$auc_mean, r$tss_mean,
                r$cbi_mean, r$sedi_mean))
  }
}

cat("\nOutput figures:\n")
for (f in sort(list.files(FIG_DIR, pattern = "\\.png$")))
  cat(sprintf("  %s\n", f))

cat("\nCheckpoint files (delete to recompute):\n")
for (f in c("vif_result.rds", "split.rds", "dag_result.rds",
            "ate_result.rds", "screened_result.rds", "roles_result.rds",
            "fit_result.rds", "eval_result.rds", "cv_result.rds",
            "cate_result.rds")) {
  mark <- if (file.exists(file.path(MODEL_DIR, f))) "+" else "-"
  cat(sprintf("  [%s] %s\n", mark, f))
}
