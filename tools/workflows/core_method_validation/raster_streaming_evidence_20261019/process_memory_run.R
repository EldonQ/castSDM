.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
root <- "E:/Package/_toolaudit/castSDM_raster_20261019"
writeLines(as.character(Sys.getpid()), file.path(root, paste0("memory_", mode, ".pid")))
pkgload::load_all("E:/Package/cast", quiet = TRUE)
old <- new.env(parent = asNamespace("castSDM"))
sys.source(file.path(root, "baseline/R/ensemble-raster.R"), envir = old)
ref <- data.frame(x1 = c(-2, 0, 2), z1 = c(-1, 0, 1))
fit <- cast_as_fit(list(a = 1, b = 2), ref, list(
  a = function(model, newdata) plogis(newdata$x1 - .4 * newdata$z1),
  b = function(model, newdata) plogis(.8 * newdata$x1 - .2 * newdata$z1)))
cv <- new_cast_cv(metrics = data.frame(model = c("a", "b"), auc_mean = .8,
  tss_mean = .6, n_folds = 2), fold_metrics = data.frame(), folds = 1:2, k = 2,
  block_method = "external", thresholds = c(a = .5, b = .5),
  oof = data.frame(obs = rep(0:1, 10), HSS_a = rep(c(.2, .8), 10), HSS_b = rep(c(.3, .7), 10)))
path <- file.path(root, "memory_input.tif")
if (mode == "prepare") {
  r <- terra::rast(nrows = 1500, ncols = 1500, nlyrs = 2, xmin = 0, xmax = 1500, ymin = 0, ymax = 1500)
  writer <- .cast_raster_writer(r, path, c("x1", "z1"), overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  set.seed(19)
  for (row in seq.int(1L, 1500L, by = 20L)) .cast_writer_block(writer,
    cbind(runif(30000, -1, 1), runif(30000, -.5, .5)), row, 20L)
  .cast_writer_finish(writer)
} else {
  r <- terra::rast(path)
  fun <- if (mode == "old") old$cast_ensemble_raster else cast_ensemble_raster
  fun(fit, cv, r, file.path(root, paste0("memory_", mode)), method = "equal",
    max_memory_mb = 2, extrapolation = FALSE, verbose = FALSE, overwrite = TRUE)
}
