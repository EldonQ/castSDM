# AI-assisted development (2026-11-11); scope: inst/AI_ASSISTANCE.md.
# source("E:/Package/castSDM_workflows/single_species/05_SHINY.R", encoding = "UTF-8")
# Rscript 05_SHINY.R [case-directory]
local({
  sources <- vapply(sys.frames(), function(frame)
    if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else
      if (is.character(frame$file) && length(frame$file) == 1L &&
          basename(frame$file) == "05_SHINY.R") frame$file else "", character(1))
  script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
  override <- Sys.getenv("CAST_EXAMPLE_ROOT", "")
  roots <- if (nzchar(override)) override else c(
    dirname(tail(sources[nzchar(sources)], 1L)), dirname(script), getwd())
  roots <- roots[file.exists(file.path(roots, "runtime.R"))]
  if (!length(roots)) stop("Cannot locate runtime.R; source 05_SHINY.R using its full path.")
  root <- normalizePath(roots[1L], winslash = "/", mustWork = TRUE)
  sys.source(file.path(root, "runtime.R"), envir = environment())
  case_runtime_setup(root, packages = c("castSDM", "terra", "ranger", "maxnet", "mgcv",
    "shiny", "bslib", "leaflet", "callr", "rpart", "zip"), minimum_cast = "0.29.0")
  sys.source(file.path(root, "01_config.R"), envir = environment())
  sys.source(file.path(root, "inventory.R"), envir = environment())
  configs <- castsdm_species_configs(root)
  taxon <- Sys.getenv("CAST_SPECIES", "Schizothoracinae")
  cfg <- configs[[taxon]]
  if (is.null(cfg)) stop("Unknown CAST_SPECIES: ", taxon)
  args <- commandArgs(TRUE)
  positional <- args[!startsWith(args, "--")]
  directory <- Sys.getenv("CAST_CASE_DIR", if (length(positional)) positional[1L] else
    file.path(cfg$output_root, Sys.getenv("CAST_RUN_ID", paste0("cast028_", cfg$scenario_lulc_product)), taxon))
  dry <- identical(Sys.getenv("CAST_DRY_RUN"), "1") || "--dry-run" %in% args
  app <- castSDM::cast_app(directory, workflow = file.path(root, "00_RUN_ALL.R"), run = !dry)
  if (dry) message("Shiny entry ready; no server, training or scenario computation started.")
  invisible(app)
})
