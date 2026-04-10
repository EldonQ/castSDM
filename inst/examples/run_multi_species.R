# ==============================================================================
# castSDM Multi-Species Parallel Modeling Example
#
# Demonstrates parallel processing of multiple species using cast_batch().
# Uses three ungulate species from Central/Western China:
#   - Ovis ammon (Argali / 盘羊)
#   - Gazella subgutturosa (Goitered Gazelle / 鹅喉羚)
#   - Pseudois nayaur (Blue Sheep / 岩羊)
#
# Workflow:
#   1. Load species data from inst/extdata/ CSVs
#   2. Set up parallel backend (future)
#   3. Run cast_batch() for all species
#   4. Compare model performance across species (boxplot)
#   5. Compute inter-model spatial consistency per species
#   6. Save all outputs
#
# Requirements:
#   install.packages(c("future", "future.apply"))
#   # For full CAST model: install torch
#   # For MaxEnt: install.packages("maxnet")
#
# Run:
#   library(castSDM)
#   source(system.file("examples/run_multi_species.R", package = "castSDM"))
# ==============================================================================

# ── Configuration ────────────────────────────────────────────────────────────

# Models: remove "cast" if torch is not installed
MODELS     <- c("rf", "maxent", "brt")
OUTPUT_DIR <- "castSDM_multi_species"
N_WORKERS  <- 3L
SEED       <- 42L

# ── Load castSDM ─────────────────────────────────────────────────────────────

if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(castSDM)
}

for (pkg in c("ggplot2", "future", "future.apply")) {
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

cat("\nRunning cast_batch()...\n\n")

batch_result <- cast_batch(
  species_list = species_list,
  env_data     = env_grid,
  models       = MODELS,
  output_dir   = OUTPUT_DIR,
  do_cv        = TRUE,
  cv_k         = 5L,
  parallel     = TRUE,
  seed         = SEED
)

print(batch_result)


# ── Step 4: Model performance comparison boxplot ─────────────────────────────

cat("\n[Plot] Multi-species model performance comparison...\n")

fig_dir <- file.path(OUTPUT_DIR, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p_compare <- plot(batch_result, metrics = c("auc", "tss", "cbi", "sedi"))
ggplot2::ggsave(
  file.path(fig_dir, "multi_species_comparison.png"),
  p_compare, width = 14, height = 6, dpi = 300
)
cat("  [OK] multi_species_comparison.png\n")


# ── Step 5: Inter-model consistency per species ──────────────────────────────

cat("\n[Plot] Inter-model spatial consistency heatmaps...\n")

for (sp in names(batch_result$results)) {
  r <- batch_result$results[[sp]]
  if (is.null(r) || is.null(r$predict)) next

  consistency <- tryCatch(
    cast_consistency(r$predict),
    error = function(e) {
      message(sprintf("  [skip] %s: %s", sp, e$message))
      NULL
    }
  )

  if (!is.null(consistency)) {
    print(consistency)

    p_cons <- tryCatch(
      plot(consistency, species = sp),
      error = function(e) NULL
    )
    if (!is.null(p_cons)) {
      ggplot2::ggsave(
        file.path(fig_dir, sprintf("consistency_%s.png", sp)),
        p_cons, width = 16, height = 5.5, dpi = 300
      )
      cat(sprintf("  [OK] consistency_%s.png\n", sp))
    }
  }
}


# ── Step 6: Per-species DAG plots ────────────────────────────────────────────

cat("\n[Plot] Per-species causal DAG plots...\n")

for (sp in names(batch_result$results)) {
  r <- batch_result$results[[sp]]
  if (is.null(r) || is.null(r$dag)) next

  if (requireNamespace("ggraph", quietly = TRUE) &&
      requireNamespace("igraph", quietly = TRUE)) {
    p_dag <- tryCatch(
      plot(r$dag, roles = r$roles, screen = r$screen, species = sp),
      error = function(e) NULL
    )
    if (!is.null(p_dag)) {
      ggplot2::ggsave(
        file.path(fig_dir, sprintf("dag_%s.png", sp)),
        p_dag, width = 14, height = 10, dpi = 300
      )
      cat(sprintf("  [OK] dag_%s.png\n", sp))
    }
  }
}


# ── Cleanup ──────────────────────────────────────────────────────────────────

future::plan(future::sequential)

cat("\n============================================================\n")
cat("Multi-species batch modeling complete!\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("============================================================\n")
