# Rscript 00_RUN_ALL.R [--dry-run]
case_main <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", args[grepl("^--file=", args)])
  root <- Sys.getenv("CAST_EXAMPLE_ROOT", if (length(script)) dirname(script[1L]) else getwd())
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("castSDM", quietly = TRUE) || utils::packageVersion("castSDM") < "0.26.0")
    stop("Install castSDM >= 0.26.0 first; this example uses the installed package.")
  needed <- c("terra", "sf", "data.table", "digest", "filelock", "ranger", "maxnet", "mgcv", "glmnet", "blockCV", "FNN")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  for (file in c("01_config.R", "prepare_data.R", "data_adapter.R", "inventory.R",
                 "scenario_case.R", "02_pipeline.R", "03_figures.R"))
    sys.source(file.path(root, file), envir = environment(), keep.source = FALSE)
  taxon <- Sys.getenv("CAST_SPECIES", "Schizothoracinae")
  cfg <- castsdm_species_configs(root)[[taxon]]
  if (is.null(cfg)) stop("Unknown CAST_SPECIES: ", taxon)
  for (path in c(cfg$occurrence_file, cfg$nsdm_pool_csv, cfg$future_layers_csv))
    if (!file.exists(path)) stop("Missing input: ", path, ". Check paths or run 00_INVENTORY.R.")
  if (cfg$input_type == "gbif" && !nzchar(Sys.which("cs2cs")))
    stop("Add PROJ cs2cs to PATH for this data adapter before running.")
  if (identical(Sys.getenv("CAST_DRY_RUN"), "1") || "--dry-run" %in% commandArgs(TRUE)) {
    message("Preflight OK: ", taxon, "; no preparation, training or raster projection started.")
    return(invisible(cfg))
  }
  run_species(cfg)
}
case_main()
