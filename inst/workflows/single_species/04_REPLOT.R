# Rscript 04_REPLOT.R <species-output-directory>
local({
  sources <- vapply(sys.frames(), function(frame)
    if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else
      if (is.character(frame$file) && length(frame$file) == 1L &&
          basename(frame$file) == "04_REPLOT.R") frame$file else "", character(1))
  script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
  override <- Sys.getenv("CAST_EXAMPLE_ROOT", "")
  roots <- if (nzchar(override)) override else c(
    dirname(tail(sources[nzchar(sources)], 1L)), dirname(script), getwd())
  roots <- roots[file.exists(file.path(roots, "runtime.R"))]
  if (!length(roots)) stop("Cannot locate runtime.R; source 04_REPLOT.R using its full path.")
  root <- normalizePath(roots[1L], winslash = "/", mustWork = TRUE)
  sys.source(file.path(root, "runtime.R"), envir = environment())
  case_runtime_setup(root, c("castSDM", "terra", "sf", "ggplot2", "patchwork", "scales", "ragg"),
    minimum_cast = "0.33.0")
  args <- commandArgs(TRUE)
  positional <- args[!startsWith(args, "--")]
  directory <- Sys.getenv("CAST_CASE_DIR", if (length(positional) == 1L) positional[1L] else "")
  if (!nzchar(directory)) stop("Supply the case directory: Rscript 04_REPLOT.R <directory>, or set CAST_CASE_DIR before source().")
  out <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(out, "case.rds")))
    stop("This case has not finished its scientific stages yet; choose a completed output directory.")
  case <- readRDS(file.path(out, "case.rds"))
  if (!identical(case$schema, "cast_single_species/1")) stop("Unsupported case bundle.")
  sys.source(file.path(root, "prepare_data.R"), envir = environment())
  sys.source(file.path(root, "03_figures.R"), envir = environment())
  comparison <- readRDS(file.path(out, case$comparison))
  # Rebase relative scenario artifacts if the output bundle was moved.
  comparison$results <- lapply(comparison$results, function(x) {
    if (!is.null(x$artifact)) x$raster_file <- file.path(out, "comparison", "rasters", x$artifact)
    x
  })
  run <- readRDS(file.path(out, case$training))
  render_only <- identical(Sys.getenv("CAST_RENDER_ONLY"), "1") || "--render-only" %in% args
  include <- strsplit(Sys.getenv("CAST_FIGURES", "all"), ",", fixed = TRUE)[[1L]]
  case_figures(comparison, out, case$report_scenario, result = run, cfg = case$cfg,
    render_only = render_only, include = include)
})
