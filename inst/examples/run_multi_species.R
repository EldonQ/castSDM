# ==============================================================================
# castSDM Multi-Species Parallel Modeling Example
#
# Demonstrates parallel processing of multiple species using cast_batch().
# Uses three ungulate species from Central/Western China:
#   - Ovis ammon (Argali / 盘羊)
#   - Gazella subgutturosa (Goitered Gazelle / 鹅喉羚)
#   - Pseudois nayaur (Blue Sheep / 岩羊)
#
# Output structure:
#   castSDM_multi_species/
#   ├── figures/                    # Cross-species comparison plots
#   │   └── multi_species_comparison.png
#   ├── Ovis_ammon/                 # Per-species results
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
#   ├── Gazella_subgutturosa/       # (same structure)
#   └── Pseudois_nayaur/            # (same structure)
#
# Requirements:
#   install.packages(c("future", "future.apply", "ggplot2", "patchwork"))
#
# Run:
#   library(castSDM)
#   source(system.file("examples/run_multi_species.R", package = "castSDM"))
# ==============================================================================

# ── Configuration ────────────────────────────────────────────────────────────

MODELS      <- c("rf", "maxent", "brt")   # add "cast" if torch is installed
OUTPUT_DIR  <- "castSDM_multi_species"
N_WORKERS   <- 3L
SEED        <- 42L
FIG_DPI     <- 300L    # unified DPI for all saved figures

# ── Load castSDM ─────────────────────────────────────────────────────────────

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

# ── Step 1: Load species data ────────────────────────────────────────────────

cat("============================================================\n")
cat("castSDM Multi-Species Parallel Modeling\n")
cat("============================================================\n\n")

species_files <- list(
  Ovis_ammon            = "CAST_Ovis_ammon_Res9_screened.csv",
  Gazella_subgutturosa  = "CAST_Gazella_subgutturosa_Res9_screened.csv",
  Pseudois_nayaur       = "CAST_Pseudois_nayaur_Res9_screened.csv"
)

species_list <- lapply(species_files, function(f) {
  path <- system.file("extdata", f, package = "castSDM")
  if (path == "") stop("File not found in package: ", f)
  data.table::fread(path, data.table = FALSE)
})

for (sp in names(species_list)) {
  d <- species_list[[sp]]
  cat(sprintf("  %s: %d rows | presence=%d | absence=%d\n",
              sp, nrow(d), sum(d$presence == 1), sum(d$presence == 0)))
}

# Load shared prediction grid
env_grid_path <- system.file(
  "extdata", "China_EnvData_Res9_Screened.csv", package = "castSDM"
)
env_grid <- if (env_grid_path != "") {
  data.table::fread(env_grid_path, data.table = FALSE)
} else {
  message("Prediction grid not found; spatial prediction will be skipped.")
  NULL
}


# ── Step 2: Set up parallel backend ──────────────────────────────────────────

cat(sprintf("\nSetting up parallel backend: %d workers\n", N_WORKERS))
future::plan(future::multisession, workers = N_WORKERS)


# ── Step 3: Run batch pipeline ───────────────────────────────────────────────
# cast_batch() runs the full pipeline per species:
#   prepare → DAG → ATE → screen → roles → fit → evaluate → CV
#   → predict → CATE → consistency
# All per-species plots are saved automatically to <output>/<species>/figures/

cat("\nRunning cast_batch()...\n\n")

batch_result <- cast_batch(
  species_list = species_list,
  env_data     = env_grid,
  models       = MODELS,
  output_dir   = OUTPUT_DIR,
  do_cv        = TRUE,
  cv_k         = 5L,
  n_bootstrap  = 50L,       # increase to 100+ for publication
  fig_dpi      = FIG_DPI,
  parallel     = TRUE,
  seed         = SEED
)

print(batch_result)


# ── Step 4: Cross-species comparison plots (shared figures/) ─────────────────
# These go in the top-level figures/ directory, not per-species folders.

fig_dir <- file.path(OUTPUT_DIR, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n[Plot] Multi-species model performance comparison...\n")

p_compare <- plot(batch_result, metrics = c("auc", "tss", "cbi", "sedi"))
ggplot2::ggsave(
  file.path(fig_dir, "multi_species_comparison.png"),
  p_compare, width = 14, height = 6, dpi = FIG_DPI
)
cat("  [OK] multi_species_comparison.png\n")


# ── Cleanup ──────────────────────────────────────────────────────────────────

future::plan(future::sequential)

cat("\n============================================================\n")
cat("Multi-species batch modeling complete!\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("============================================================\n\n")

# List all generated per-species plots
for (sp in names(batch_result$results)) {
  sp_figs <- file.path(OUTPUT_DIR, sp, "figures")
  if (dir.exists(sp_figs)) {
    figs <- list.files(sp_figs, pattern = "\\.png$")
    cat(sprintf("  %s: %d figures\n", sp, length(figs)))
    for (f in figs) cat(sprintf("    %s\n", f))
  }
}
