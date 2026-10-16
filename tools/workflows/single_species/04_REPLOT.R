# Replot a completed case into a separate directory, preserving its receipt.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript 04_REPLOT.R <completed-case-directory> [figure-directory]")
case_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
bundle <- readRDS(file.path(case_dir, "castSDM_result.rds"))
cfg <- bundle$config
pkgload::load_all(cfg$cast_root, quiet = TRUE)
source(file.path(dirname(cfg$adapter), "workflow_provenance.R"), encoding = "UTF-8")
source(file.path(dirname(dirname(cfg$adapter)), "single_species", "03_figures.R"), encoding = "UTF-8")
output <- if (length(args) >= 2L) args[2] else file.path(case_dir,
  paste0("replot_", format(Sys.time(), "%Y%m%dT%H%M%S")))
if (dir.exists(output)) stop("Choose a new figure directory to preserve existing exports")
dir.create(output, recursive = TRUE)
receipt <- bundle$computation_provenance %||% bundle$provenance
if (is.null(receipt)) stop("Completed bundle has no computation source identity")
receipt_path <- file.path(output, "computation_source.rds")
saveRDS(receipt, receipt_path)
make_species_figures(cfg, computation_source = receipt_path, output_dir = output)
message("Figures exported to ", output, "; visual inspection is still required")
