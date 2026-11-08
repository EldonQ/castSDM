# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# source("E:/Package/castSDM_workflows/single_species/00_RUN_ALL.R", encoding = "UTF-8")
# or Rscript 00_RUN_ALL.R [--dry-run]
local({
  sources <- vapply(sys.frames(), function(frame)
    if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else
      if (is.character(frame$file) && length(frame$file) == 1L &&
          basename(frame$file) == "00_RUN_ALL.R") frame$file else "", character(1))
  script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
  override <- Sys.getenv("CAST_EXAMPLE_ROOT", "")
  roots <- if (nzchar(override)) override else c(
    dirname(tail(sources[nzchar(sources)], 1L)), dirname(script), getwd())
  roots <- roots[file.exists(file.path(roots, "runtime.R"))]
  if (!length(roots)) stop("Cannot locate runtime.R; source 00_RUN_ALL.R using its full path.")
  root <- normalizePath(roots[1L], winslash = "/", mustWork = TRUE)
  sys.source(file.path(root, "runtime.R"), envir = environment())
  runtime <- case_runtime_setup(root, packages = c("castSDM", "terra", "sf", "data.table",
    "digest", "filelock", "ranger", "maxnet", "mgcv", "glmnet", "pROC", "blockCV", "FNN",
    "ggplot2", "patchwork", "scales", "ragg"), minimum_cast = "0.28.1")
  for (file in c("01_config.R", "prepare_data.R", "data_adapter.R", "inventory.R",
                 "scenario_case.R", "02_pipeline.R", "03_figures.R"))
    sys.source(file.path(root, file), envir = environment(), keep.source = FALSE)
  taxon <- Sys.getenv("CAST_SPECIES", "Schizothoracinae")
  cfg <- castsdm_species_configs(root)[[taxon]]
  if (is.null(cfg)) stop("Unknown CAST_SPECIES: ", taxon)
  cfg$runtime <- runtime
  cfg$runtime$script_hashes <- tools::md5sum(list.files(root, pattern = "\\.R$", full.names = TRUE))
  for (path in c(cfg$occurrence_file, cfg$nsdm_pool_csv, cfg$future_layers_csv))
    if (!file.exists(path)) stop("Missing input: ", path, ". Check paths or run 00_INVENTORY.R.")
  if (cfg$input_type == "gbif" && !nzchar(Sys.which("cs2cs")))
    stop("Add PROJ cs2cs to PATH for this data adapter before running.")
  if (identical(Sys.getenv("CAST_DRY_RUN"), "1") || "--dry-run" %in% commandArgs(TRUE)) {
    message("Preflight OK: ", taxon, "; no preparation, training or raster projection started.")
    invisible(cfg)
  } else {
    case_runtime_resume_patch(cfg)
    run_species(cfg)
  }
})
