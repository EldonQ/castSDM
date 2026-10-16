# Executable audit of remaining legacy contracts; does not alter package sources.
.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
pkgload::load_all("E:/Package/cast", quiet = TRUE)
root <- "E:/Package/castSDM_workflows/core_method_validation"
output <- file.path(root, "legacy_contract_audit")
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
  method = "equal", mask = mask, extrapolation = FALSE, verbose = FALSE)
future <- r + .1; names(future) <- "x"
paired <- cast_attribution(fit, r, future, mask = mask, verbose = FALSE)
legacy_values <- as.numeric(terra::values(terra::rast(legacy$hss_path)))
paired_values <- as.numeric(terra::values(terra::rast(paired$raster_file)[["current"]]))
stopifnot(is.finite(legacy_values[2]), is.na(paired_values[2]),
  is.finite(legacy_values[4]), is.na(paired_values[4]))
unknown_first <- castSDM:::.cast_thin_points(c(0, 0), c(0, 0), 1,
  order = "year_desc", year = c(NA_real_, 2020))
stopifnot(identical(unknown_first, c(TRUE, FALSE)))
results <- data.frame(
  check = c("zero_mask_cell", "infinite_predictor_cell", "year_desc_missing_vs_known"),
  observed = c("legacy predicts; paired attribution masks", "legacy predicts; paired attribution drops",
    "missing year retained before year 2020"),
  classification = c("documented legacy contract differs from paired contract",
    "nonfinite handling inconsistent across prediction interfaces", "confirmed ordering defect"))
utils::write.csv(results, file.path(output, "remaining_contract_findings.csv"), row.names = FALSE)
utils::write.csv(data.frame(cell = 1:4, mask = c(1, 0, 1, 1),
  predictor = c(0, .5, 1, Inf), legacy = legacy_values, paired = paired_values),
  file.path(output, "mask_nonfinite_example.csv"), row.names = FALSE)
files <- c("E:/Package/cast/R/utils.R", "E:/Package/cast/R/ensemble-raster.R",
  "E:/Package/cast/R/scenario-attribution.R", file.path(root, "06_legacy_contract_audit.R"))
utils::write.csv(data.frame(file = files, sha256 = vapply(files, function(f)
  digest::digest(file = f, algo = "sha256"), character(1))),
  file.path(output, "source_manifest.csv"), row.names = FALSE)
utils::capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
writeLines(c("# Remaining legacy contracts reproduced", "",
  "Zero mask values are valid under the documented legacy non-NA mask rule, but not under paired attribution.",
  "Legacy raster complete.cases passes Inf to a callback; paired attribution rejects that cell.",
  "year_desc incorrectly retains an unknown-date record before a collocated 2020 record.",
  "These findings are not fixed in this source snapshot; existing completed evidence is preserved."),
  file.path(output, "RESULTS.md"))
print(results)
