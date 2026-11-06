# Rscript 04_REPLOT.R <species-output-directory>
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Supply the directory containing case.rds.")
if (!requireNamespace("castSDM", quietly = TRUE) || utils::packageVersion("castSDM") < "0.26.0")
  stop("Install castSDM >= 0.26.0 before replotting.")
out <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
case <- readRDS(file.path(out, "case.rds"))
if (!identical(case$schema, "cast_single_species/1")) stop("Unsupported case bundle.")
source(file.path(case$cfg$example_root, "03_figures.R"), encoding = "UTF-8")
comparison <- readRDS(file.path(out, case$comparison))
# Rebase relative scenario artifacts if the output bundle was moved.
comparison$results <- lapply(comparison$results, function(x) {
  if (!is.null(x$artifact)) x$raster_file <- file.path(out, "comparison", "rasters", x$artifact)
  x
})
case_figures(comparison, out, case$report_scenario)
