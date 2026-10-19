# Regression verification of raster contracts; historical audit outputs are retained.
.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
pkgload::load_all("E:/Package/cast", quiet = TRUE)
root <- "E:/Package/castSDM_workflows/core_method_validation"
output <- file.path(root, "raster_contract_validation_20261019")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
reference <- data.frame(x = c(-2, 0, 2))
fit <- cast_as_fit(list(external = 1), reference,
  function(model, newdata) stats::plogis(newdata$x))
cv <- new_cast_cv(metrics = data.frame(model = "external", auc_mean = .8,
  tss_mean = .6, n_folds = 2), fold_metrics = data.frame(),
  folds = 1:2, k = 2L, block_method = "external", thresholds = c(external = .5))
r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2,
  crs = "EPSG:4326")
names(r) <- "x"; terra::values(r) <- c(0, .5, 1, Inf)
mask <- r; terra::values(mask) <- c(1, 0, 1, 1)
legacy <- cast_ensemble_raster(fit, cv, r, file.path(output, "rasters"),
  method = "equal", mask = mask, extrapolation = FALSE, verbose = FALSE, overwrite = TRUE)
future <- r + .1; names(future) <- "x"
paired <- cast_attribution(fit, r, future, mask = mask, verbose = FALSE)
legacy_values <- as.numeric(terra::values(terra::rast(legacy$hss_path)))
paired_values <- as.numeric(terra::values(terra::rast(paired$raster_file)[["current"]]))
stopifnot(is.na(legacy_values[2]), is.na(paired_values[2]),
  is.na(legacy_values[4]), is.na(paired_values[4]),
  isTRUE(all.equal(legacy_values, paired_values, tolerance = 1e-7)))
unknown_first <- castSDM:::.cast_thin_points(c(0, 0), c(0, 0), 1,
  order = "year_desc", year = c(NA_real_, 2020))
stopifnot(identical(unknown_first, c(FALSE, TRUE)))
results <- data.frame(
  check = c("zero_mask_cell", "infinite_predictor_cell", "year_desc_missing_vs_known"),
  observed = c("ensemble and paired attribution exclude zero mask", "ensemble and paired attribution drop Inf",
    "year 2020 retained before unknown year"),
  classification = rep("regression passed", 3))
utils::write.csv(results, file.path(output, "verified_contracts.csv"), row.names = FALSE)
utils::write.csv(data.frame(cell = 1:4, mask = c(1, 0, 1, 1),
  predictor = c(0, .5, 1, Inf), legacy = legacy_values, paired = paired_values),
  file.path(output, "mask_nonfinite_example.csv"), row.names = FALSE)
files <- c("E:/Package/cast/R/utils.R", "E:/Package/cast/R/raster-io.R", "E:/Package/cast/R/ensemble-raster.R",
  "E:/Package/cast/R/scenario-attribution.R", file.path(root, "06_legacy_contract_audit.R"))
utils::write.csv(data.frame(file = files, sha256 = vapply(files, function(f)
  digest::digest(file = f, algo = "sha256"), character(1))),
  file.path(output, "source_manifest.csv"), row.names = FALSE)
utils::capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
writeLines(c("# Raster contract regressions verified", "",
  "Zero/non-finite mask cells and Inf predictors are excluded consistently.",
  "year_desc puts finite known years before unknown years.",
  "Historical legacy_contract_audit artifacts are preserved and describe the previous source snapshot.",
  "This is an engineering verification, not a new ecological experiment."),
  file.path(output, "RESULTS.md"))
print(results)
