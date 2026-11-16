# Rscript tests/real-data/temporal.R <prepare_stage.rds> <occurrence_qc.csv> <eligible_pool.csv> <QA_directory>
# Uses saved real rows/labels and extracts only missing source time slices.
args <- commandArgs(TRUE)
if (length(args) != 4L || any(!file.exists(args[1:3])))
  stop("Supply real prepared rows, occurrence QC, eligible raster inventory and a separate QA directory.")
library(castSDM)
env <- new.env(parent = asNamespace("castSDM"))
sys.source(system.file("workflows/single_species/02_pipeline.R", package = "castSDM"), env)
split <- readRDS(args[1L])
meta <- utils::read.csv(args[3L])
dated <- unique(meta$construct_id[meta$layer %in% split$env_vars & is.finite(meta$year)])
meta <- meta[meta$layer %in% split$env_vars | (is.finite(meta$year) & meta$construct_id %in% dated), ]
data <- rbind(split$train, split$test)
missing <- setdiff(meta$layer, names(data))
if (length(missing)) {
  paths <- meta$path[match(missing, meta$layer)]
  stopifnot(all(file.exists(paths)))
  raster <- terra::rast(paths); names(raster) <- missing
  values <- terra::extract(raster, as.matrix(data[, c("lon", "lat")]))
  data[missing] <- values
}
data <- data[, c("lon", "lat", "presence", meta$layer)]
n_train <- nrow(split$train)
split$train <- data[seq_len(n_train), ]
split$test <- data[n_train+seq_len(nrow(split$test)), ]
split$env_vars <- meta$layer
dir.create(args[4L], recursive = TRUE, showWarnings = FALSE)
qc <- file.path(args[4L], "occurrence_qc.csv")
if (normalizePath(args[2L], winslash = "/") != normalizePath(qc, winslash = "/", mustWork = FALSE))
  file.copy(args[2L], qc, overwrite = TRUE)
reference <- env$case_training_year(args[4L], split)
stopifnot(is.list(reference), is.numeric(reference$year), length(reference$year) == 1L,
  is.null(attributes(reference$year)), is.character(reference$source))
matched <- env$case_match_training_layers(data, split, meta, reference$year, "npts")
stopifnot(nrow(matched$matching) == nrow(meta), nrow(matched$data) == nrow(data),
  nrow(matched$split$train) == nrow(split$train), nrow(matched$split$test) == nrow(split$test),
  identical(matched$meta$layer, matched$split$env_vars),
  all(matched$matching$reference_year == reference$year))
# Replay the exact old scalar shape using the same observed year/provenance.
legacy <- reference$year; attr(legacy, "source") <- reference$source
replayed <- env$case_match_training_layers(data, split, meta, legacy, "npts")
stopifnot(identical(matched, replayed))
utils::write.csv(matched$matching, file.path(args[4L], "temporal_layer_matching.csv"), row.names = FALSE)
utils::write.csv(data.frame(reference_year = reference$year, source = reference$source),
  file.path(args[4L], "temporal_matching.csv"), row.names = FALSE)
roundtrip <- utils::read.csv(file.path(args[4L], "temporal_layer_matching.csv"))
stopifnot(nrow(roundtrip) == nrow(meta), all(roundtrip$reference_year == reference$year))
cat("Real-data reference -> layer matching -> CSV chain passed:", nrow(data), "rows;",
  nrow(meta), "input layers ->", nrow(matched$meta), "retained; reference", reference$year, "\n")
