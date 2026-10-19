.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
pkgload::load_all("E:/Package/cast", quiet = TRUE)
base_root <- "E:/Package/_toolaudit/castSDM_raster_20261019"
root <- file.path(base_root, "benchmark_v3")
dir.create(root, showWarnings = FALSE)
output <- file.path(root, "evidence")
dir.create(output, showWarnings = FALSE)
old <- new.env(parent = asNamespace("castSDM"))
for (f in c("ensemble-raster", "predict-tiled", "effect-products", "scenario-attribution"))
  sys.source(file.path(base_root, "baseline/R", paste0(f, ".R")), envir = old)
reference <- data.frame(x1 = seq(-2, 2, length.out = 100), z1 = seq(-1, 1, length.out = 100))
fit <- cast_as_fit(list(a = 1, b = 2), reference, list(
  a = function(model, newdata) plogis(newdata$x1 - .4 * newdata$z1),
  b = function(model, newdata) plogis(.8 * newdata$x1 - .2 * newdata$z1)))
cv <- new_cast_cv(metrics = data.frame(model = c("a", "b"), auc_mean = .8, tss_mean = .6,
  n_folds = 2), fold_metrics = data.frame(), folds = 1:2, k = 2,
  block_method = "external", thresholds = c(a = .5, b = .5),
  oof = data.frame(obs = c(0L, 1L), HSS_a = c(.2, .8), HSS_b = c(.3, .7)))
r <- terra::rast(nrows = 600, ncols = 400, xmin = 0, xmax = 400, ymin = 0, ymax = 600, nlyrs = 2)
names(r) <- c("x1", "z1")
set.seed(19)
terra::values(r) <- cbind(runif(terra::ncell(r), -1, 1), runif(terra::ncell(r), -.5, .5))
r <- terra::writeRaster(r, file.path(root, "input.tif"), overwrite = TRUE)
profile <- function(name, code) {
  gc()
  file <- file.path(root, paste0(name, "_Rprofmem.txt"))
  Rprofmem(file)
  started <- proc.time()[[3]]
  result <- force(code)
  elapsed <- proc.time()[[3]] - started
  Rprofmem(NULL)
  lines <- readLines(file, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  list(value = result, metrics = data.frame(run = name, elapsed_seconds = elapsed,
    largest_recorded_R_allocation_bytes = max(bytes, na.rm = TRUE),
    total_recorded_R_allocations_bytes = sum(bytes, na.rm = TRUE)))
}
old_ensemble <- profile("old_ensemble", old$cast_ensemble_raster(fit, cv, r,
  file.path(root, "old_ensemble"), method = "equal", max_memory_mb = 1, verbose = FALSE))
new_ensemble <- profile("new_ensemble", cast_ensemble_raster(fit, cv, r,
  file.path(root, "new_ensemble"), method = "equal", max_memory_mb = 1, verbose = FALSE))
old_effect <- profile("old_effect", old$cast_effect_map(fit, r, shift = .1,
  block_rows = 20, min_support = 0, filename = file.path(root, "old_effect.tif"), verbose = FALSE))
new_effect <- profile("new_effect", cast_effect_map(fit, r, shift = .1,
  block_rows = 20, min_support = 0, filename = file.path(root, "new_effect.tif"), verbose = FALSE))
old_tiled <- profile("old_tiled", old$cast_predict_tiled(fit, r,
  file.path(root, "old_tiled"), tile_size = c(20, 100), verbose = FALSE))
new_tiled <- profile("new_tiled", cast_predict_tiled(fit, r,
  file.path(root, "new_tiled"), tile_size = c(20, 100), verbose = FALSE))
compare <- function(a, b, name, tolerance) {
  va <- terra::values(a); vb <- terra::values(b)
  stopifnot(terra::compareGeom(a, b), identical(is.na(va), is.na(vb)))
  difference <- max(abs(va - vb), na.rm = TRUE)
  stopifnot(difference <= tolerance)
  data.frame(comparison = name, cells = nrow(va), layers = ncol(va),
    max_abs_difference = difference, tolerance = tolerance)
}
checks <- list()
for (layer in c("hss_path", "hss_sd_path", "binary_path", "mess_path")) checks[[layer]] <- compare(
  terra::rast(old_ensemble$value[[layer]]), terra::rast(new_ensemble$value[[layer]]), layer, 1e-7)
checks$effect <- compare(old_effect$value, new_effect$value, "single_shift", 1e-7)
stopifnot(isTRUE(all.equal(attr(old_effect$value, "effect_table"),
  attr(new_effect$value, "effect_table"), tolerance = 1e-12)))
for (model in names(fit$models)) checks[[model]] <- compare(
  terra::rast(old_tiled$value$rasters[[model]]), terra::rast(new_tiled$value$rasters[[model]]),
  paste0("tiled_", model), 1e-7)
future <- r + .1; names(future) <- names(r)
before <- old$cast_attribution(fit, r, future, block_size = 8000, verbose = FALSE)
after <- cast_attribution(fit, r, future, block_size = 8000, verbose = FALSE)
checks$attribution <- compare(terra::rast(before$raster_file), terra::rast(after$raster_file), "paired_attribution", 1e-12)
stopifnot(isTRUE(all.equal(before$summary, after$summary, tolerance = 1e-12)))
metrics <- do.call(rbind, lapply(list(old_ensemble, new_ensemble, old_effect, new_effect, old_tiled, new_tiled), `[[`, "metrics"))
write.csv(metrics, file.path(output, "R_allocation_profile.csv"), row.names = FALSE)
write.csv(do.call(rbind, checks), file.path(output, "numerical_equivalence.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
writeLines(c("Baseline is f5cbb7d (0.20.0) raster entry points, using identical current engine callbacks.",
  "600 x 400 equal-geometry grid, two external engines, finite inputs, same row/tile limits.",
  "Rprofmem measures R allocations only; largest single allocation is not process peak memory.",
  "GDAL/terra C++ allocations and resident fitted models are not captured; timings are one run, not a speed benchmark.",
  "Effect map uses FLT8S rather than FLT4S storage to preserve table/raster numerical agreement.",
  "No ecological, universal performance or new source-matched historical-case claim."), file.path(output, "PROFILE_LIMITS.txt"))
print(metrics)
print(do.call(rbind, checks))
