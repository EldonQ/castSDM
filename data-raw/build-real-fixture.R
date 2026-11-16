# Rscript data-raw/build-real-fixture.R <completed Schizothoracinae case> <package root>
# Reads existing observations, background locations and actual CHELSA rasters.
# Does not generate observations, fit models, resample grids or alter measurements.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
library(terra)
case <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
root <- normalizePath(args[2L], winslash = "/", mustWork = TRUE)
split <- readRDS(file.path(case, "training/stages/01_prepare.rds"))
run <- readRDS(file.path(case, "training/result.rds"))
comparison <- readRDS(file.path(case, "comparison/comparison.rds"))
info <- readRDS(file.path(case, "case.rds"))
drivers <- as.character(unlist(comparison$groups, use.names = FALSE))
stopifnot(length(drivers) > 0L, all(grepl("^bio[0-9]+$", drivers)))
data <- rbind(split$train, split$test)[, c("lon", "lat", "presence", drivers)]
partition <- rep(c("train", "holdout"), c(nrow(split$train), nrow(split$test)))
fold <- c(run$cv$fold_info$fold, rep(NA_integer_, nrow(split$test)))
stopifnot(length(fold) == nrow(data), all(data$presence %in% 0:1), !anyNA(data))
index <- c(split$train_index, split$test_index)
order <- order(index)
data <- data[order, ]; partition <- partition[order]; fold <- fold[order]
rownames(data) <- NULL
pool <- read.csv(file.path(case, "training_pool.csv"))
mapping <- read.csv(file.path(case, "scenario_table.csv"))
mapping <- mapping[mapping$scenario == info$report_scenario, ]
mapping <- mapping[match(drivers, mapping$layer), ]
paths <- pool$path[match(drivers, pool$layer)]
current_grid <- rast(paths); names(current_grid) <- drivers
future_grid <- rast(mapping$path); names(future_grid) <- drivers
future <- extract(future_grid, as.matrix(data[, c("lon", "lat")]))
stopifnot(identical(names(future), drivers))
reference <- data[partition == "train", drivers]
outside <- Reduce(`|`, lapply(drivers, function(v)
  future[[v]] < min(reference[[v]]) | future[[v]] > max(reference[[v]])))
centre <- which(outside & data$presence == 1L)[1L]
stopifnot(is.finite(centre))
cell <- cellFromXY(current_grid, as.matrix(data[centre, c("lon", "lat")]))
rc <- rowColFromCell(current_grid, cell)
# A small window of original grid cells around an actual occurrence. Window
# size bounds CI cost and is not a modelling or predictor-selection setting.
rows <- max(1L, rc[1, 1]-2L):min(nrow(current_grid), rc[1, 1]+2L)
cols <- max(1L, rc[1, 2]-2L):min(ncol(current_grid), rc[1, 2]+2L)
window <- ext(current_grid[rows, cols, drop = FALSE])
directory <- file.path(root, "inst/extdata/schizothoracinae")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
writeRaster(crop(current_grid, window), file.path(directory, "current.tif"),
  datatype = "FLT8S", overwrite = TRUE)
writeRaster(crop(future_grid, window), file.path(directory, "future.tif"),
  datatype = "FLT8S", overwrite = TRUE)
saveRDS(list(data = data, future = future, partition = partition, fold = fold,
  source_row = index[order], variables = drivers,
  scenario = info$report_scenario, crs = crs(current_grid),
  background = "Saved presence/background design; zero denotes background, not observed absence.",
  population = "Pooled Schizothoracinae records in the declared river/locality analysis domain.",
  fit_settings = run$protocol[c("models", "rf_ntree", "maxent_classes", "gam_backend", "seed")],
  centre_source_row = index[order][centre]), file.path(directory, "records.rds"), compress = "xz")
sources <- unique(c(file.path(case, c("training/stages/01_prepare.rds", "training/result.rds",
  "training_pool.csv", "scenario_table.csv")), paths, mapping$path))
write.csv(data.frame(file = basename(sources), md5 = unname(tools::md5sum(sources))),
  file.path(directory, "sources.csv"), row.names = FALSE)
cat(nrow(data), "unaltered saved rows;", sum(data$presence == 1), "presences;",
  sum(data$presence == 0), "backgrounds;", length(drivers), "actual changing climate predictors.\n")
