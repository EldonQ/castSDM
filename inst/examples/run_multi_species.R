# ==============================================================================
# castSDM Multi-Species Parallel Modeling Example
#
# Demonstrates parallel batch processing of multiple species using
# cast_batch().  Supports two modes:
#   A) Auto-scan: point DATA_DIR at a folder of CSV files; the script
#      auto-detects every *.csv and treats each as one species.
#   B) Explicit list: enumerate species_files manually (see below).
#
# Output structure (one subfolder per species + shared comparison):
#   <OUTPUT_DIR>/
#   ├── figures/                           # cross-species comparison
#   │   └── multi_species_comparison.png
#   ├── <Species_1>/
#   │   ├── cast_result.rds
#   │   └── figures/
#   │       ├── causal_dag.png
#   │       ├── ate_forest_plot.png
#   │       ├── variable_screening.png
#   │       ├── model_evaluation.png
#   │       ├── spatial_cv_map.png
#   │       ├── HSS_rf.png, HSS_maxent.png, ...
#   │       ├── CATE_<var>.png
#   │       └── model_consistency.png
#   └── <Species_2>/  ...
#
# Requirements:
#   install.packages(c("future", "future.apply", "ggplot2", "patchwork"))
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — USER CONFIGURATION (edit this section)
# ══════════════════════════════════════════════════════════════════════════════

# ── 1a. Data source ──────────────────────────────────────────────────────────
# Option A: Auto-scan — set DATA_DIR to a folder containing CSV files.
#   Each CSV becomes one species. Species name is derived from filename.
#   Set to NULL to disable auto-scan and use the explicit list below.
DATA_DIR <- NULL     # e.g.  "D:/MyProject/species_csv/"

# Option B: Explicit list — files bundled with castSDM
species_files <- list(
  Ovis_ammon            = "CAST_Ovis_ammon_Res9_screened.csv",
  Gazella_subgutturosa  = "CAST_Gazella_subgutturosa_Res9_screened.csv",
  Pseudois_nayaur       = "CAST_Pseudois_nayaur_Res9_screened.csv"
)

# ── 1b. Prediction grid (optional) ──────────────────────────────────────────
# Shared environmental grid for spatial prediction.
# Set to NULL to skip spatial prediction entirely.
ENV_GRID_PATH <- NULL   # e.g.  "D:/MyProject/EnvGrid.csv"

# ── 1c. Models to fit ───────────────────────────────────────────────────────
MODELS <- c("rf", "maxent", "brt", "cast")   # add "cast" if torch is available

# ── 1d. Output ───────────────────────────────────────────────────────────────
OUTPUT_DIR <- "castSDM_multi_species"

# ── 1e. Parallel workers ────────────────────────────────────────────────────
# How many species run simultaneously. Guidelines:
#   - Each worker loads the full env_data into memory, so RAM is the
#     primary constraint.  Rule of thumb:
#       available_RAM_GB / (env_grid_rows * n_vars * 8 / 1e9 + 2)
#   - For 3-10 species   → N_WORKERS = min(n_species, parallel::detectCores() - 1)
#   - For 100+ species   → N_WORKERS = physical_cores / 2   (leave headroom)
#   - For 1000+ species  → N_WORKERS = physical_cores / 2
#     and consider reducing DAG bootstrap replicates (dag_R = 50)
#     and ATE trees (ate_num_trees = 200) to reduce per-species time.
#   - If env_data is very large (>1M rows), use fewer workers to avoid OOM.
#   - Set to 1 to disable parallelism (sequential processing).
N_WORKERS <- min(3L, parallel::detectCores() - 1L)

# ── 1f. Reproducibility & figure quality ────────────────────────────────────
SEED    <- 42L
FIG_DPI <- 1200L       # unified DPI for ALL saved figures


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — LOAD PACKAGES
# ══════════════════════════════════════════════════════════════════════════════

if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(castSDM)
}

for (pkg in c("ggplot2", "future", "future.apply", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — BUILD SPECIES LIST
# ══════════════════════════════════════════════════════════════════════════════

cat("============================================================\n")
cat("castSDM Multi-Species Parallel Modeling\n")
cat("============================================================\n\n")

if (!is.null(DATA_DIR) && dir.exists(DATA_DIR)) {
  # ── Auto-scan mode ──
  csv_paths <- list.files(DATA_DIR, pattern = "\\.csv$",
                          full.names = TRUE, recursive = FALSE)
  if (length(csv_paths) == 0) stop("No CSV files found in DATA_DIR: ", DATA_DIR)

  sp_names <- tools::file_path_sans_ext(basename(csv_paths))
  species_list <- setNames(
    lapply(csv_paths, function(f) data.table::fread(f, data.table = FALSE)),
    sp_names
  )
  cat(sprintf("[Auto-scan] Found %d species in %s\n", length(sp_names), DATA_DIR))
} else {
  # ── Explicit list mode (bundled data) ──
  species_list <- lapply(species_files, function(f) {
    path <- system.file("extdata", f, package = "castSDM")
    if (path == "") stop("File not found in package: ", f)
    data.table::fread(path, data.table = FALSE)
  })
}

for (sp in names(species_list)) {
  d <- species_list[[sp]]
  cat(sprintf("  %s: %d rows | presence=%d | absence=%d\n",
              sp, nrow(d), sum(d$presence == 1), sum(d$presence == 0)))
}

# ── Load shared prediction grid ──
env_grid <- NULL
if (!is.null(ENV_GRID_PATH) && file.exists(ENV_GRID_PATH)) {
  env_grid <- data.table::fread(ENV_GRID_PATH, data.table = FALSE)
  cat(sprintf("\nPrediction grid: %d rows x %d cols\n",
              nrow(env_grid), ncol(env_grid)))
} else {
  pkg_grid <- system.file("extdata", "China_EnvData_Res9_Screened.csv",
                           package = "castSDM")
  if (pkg_grid != "") {
    env_grid <- data.table::fread(pkg_grid, data.table = FALSE)
    cat(sprintf("\nPrediction grid (bundled): %d rows x %d cols\n",
                nrow(env_grid), ncol(env_grid)))
  } else {
    cat("\nNo prediction grid; spatial prediction will be skipped.\n")
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — SET UP PARALLEL BACKEND
# ══════════════════════════════════════════════════════════════════════════════

cat(sprintf("\nParallel workers: %d  (available cores: %d)\n",
            N_WORKERS, parallel::detectCores()))
future::plan(future::multisession, workers = N_WORKERS)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — RUN BATCH PIPELINE
# ══════════════════════════════════════════════════════════════════════════════
# cast_batch() is the one-stop interface:
#   - ALL pipeline parameters are configurable here
#   - Per-species results + ALL figures saved automatically
#   - Per-species folders: <OUTPUT_DIR>/<species>/figures/

cat("\nRunning cast_batch()...\n\n")

batch_result <- cast_batch(
  species_list = species_list,
  env_data     = env_grid,
  models       = MODELS,
  output_dir   = OUTPUT_DIR,
  fig_dpi      = FIG_DPI,
  parallel     = TRUE,
  seed         = SEED,

  # ── Data splitting ──────────────────────────────────────────────────
  train_fraction = 0.7,

  # ── DAG structure learning ──────────────────────────────────────────
  dag_R                   = 100L,     # 100 for publication; 50 for speed
  dag_strength_threshold  = 0.7,
  dag_direction_threshold = 0.6,

  # ── ATE (Double Machine Learning) ──────────────────────────────────
  ate_K          = 5L,
  ate_alpha      = 0.05,
  ate_num_trees  = 300L,
  ate_bonferroni = TRUE,
  ate_parallel   = TRUE,

  # ── Variable screening ─────────────────────────────────────────────
  screen_min_vars  = 5L,

  # ── Spatial cross-validation ───────────────────────────────────────
  do_cv           = TRUE,
  cv_k            = 10L,
  cv_block_method = "grid",

  # ── CATE (Causal forests) ──────────────────────────────────────────
  do_cate    = TRUE,
  cate_top_n = 3L,

  # ── Model hyper-parameters (passed through ... to cast_fit) ────────
  n_epochs  = 300L,    # CI-MLP training epochs
  n_runs    = 5L,      # CI-MLP ensemble seeds
  rf_ntree  = 300L,    # RF trees
  brt_n_trees = 500L,  # BRT trees
  tune_grid = TRUE    # grid search for RF/BRT
)

print(batch_result)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — CROSS-SPECIES COMPARISON PLOTS
# ══════════════════════════════════════════════════════════════════════════════

fig_dir <- file.path(OUTPUT_DIR, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n[Plot] Multi-species model performance comparison...\n")

p_compare <- plot(batch_result, metrics = c("auc", "tss", "cbi", "sedi"))
ggplot2::ggsave(
  file.path(fig_dir, "multi_species_comparison.png"),
  p_compare, width = 14, height = 6, dpi = FIG_DPI
)
cat("  Saved: multi_species_comparison.png\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — CLEANUP & SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

future::plan(future::sequential)

cat("\n============================================================\n")
cat("Multi-species batch modeling complete!\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("============================================================\n\n")

for (sp in names(batch_result$results)) {
  sp_figs <- file.path(OUTPUT_DIR, sp, "figures")
  if (dir.exists(sp_figs)) {
    figs <- list.files(sp_figs, pattern = "\\.png$")
    cat(sprintf("  %s: %d figures saved\n", sp, length(figs)))
    for (f in figs) cat(sprintf("    - %s\n", f))
  }
}
