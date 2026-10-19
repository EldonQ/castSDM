raster_contract_fixture <- function() {
  r <- terra::rast(nrows = 3, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 3)
  names(r) <- "x"
  terra::values(r) <- c(0, .5, 1, Inf, NA, -.5, 1.5, 0, .1, .2, .3, .4)
  fit <- cast_as_fit(list(external = 1), data.frame(x = c(-2, 0, 2)),
    function(model, newdata) plogis(newdata$x))
  cv <- new_cast_cv(metrics = data.frame(model = "external", auc_mean = .8,
    tss_mean = .6, n_folds = 2), fold_metrics = data.frame(), folds = 1:2,
    k = 2L, block_method = "external", thresholds = c(external = .5))
  cv$oof <- data.frame(obs = c(0L, 1L), HSS_external = c(.2, .8))
  list(r = r, fit = fit, cv = cv)
}

test_that("raster mask and finite-cell semantics agree with paired attribution", {
  skip_if_not_installed("terra")
  fx <- raster_contract_fixture()
  mask <- fx$r; terra::values(mask) <- c(1, 0, -1, 1, 1, NA, Inf, 1, 1, 1, 1, 1)
  td <- tempfile(); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res <- cast_ensemble_raster(fx$fit, fx$cv, fx$r, td, method = "equal", mask = mask,
    max_memory_mb = .001, extrapolation = FALSE, verbose = FALSE)
  future <- fx$r + .1; names(future) <- "x"
  at <- cast_attribution(fx$fit, fx$r, future, mask = mask, verbose = FALSE)
  hss <- terra::values(terra::rast(res$hss_path), mat = FALSE)
  expect_equal(hss, terra::values(terra::rast(at$raster_file)[["current"]], mat = FALSE), tolerance = 1e-7)
  expect_true(all(is.na(hss[c(2, 4, 5, 6, 7)])))
  expect_equal(res$n_valid_cells, 7)
  shifted <- terra::shift(mask, dx = .5)
  expect_error(cast_ensemble_raster(fx$fit, fx$cv, fx$r, td, mask = shifted), "Mask geometry")
  expect_error(cast_attribution(fx$fit, fx$r, future, mask = shifted), "Mask geometry")
  expect_error(cast_ensemble_raster(fx$fit, fx$cv, fx$r, td, mask = c(mask, mask)), "one layer")
  expect_error(cast_ensemble_raster(fx$fit, fx$cv, fx$r, td, max_memory_mb = 0), "max_memory_mb")
})

test_that("tiled outputs preserve cell order and reject stale input or callback", {
  skip_if_not_installed("terra")
  fx <- raster_contract_fixture()
  td <- tempfile(); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  out <- cast_predict_tiled(fx$fit, fx$r, td, tile_size = c(2, 3), verbose = FALSE)
  vals <- terra::values(terra::rast(out$rasters[[1]]), mat = FALSE)
  raw <- terra::values(fx$r, mat = FALSE)
  expected <- plogis(raw); expected[!is.finite(raw)] <- NA_real_
  expect_equal(vals, expected, tolerance = 1e-7)
  expect_equal(out$n_tiles, 4)
  expect_equal(cast_predict_tiled(fx$fit, fx$r, td, tile_size = c(2, 3), verbose = FALSE)$rasters, out$rasters)
  changed <- fx$r; terra::values(changed) <- raw[c(3, 2, 1, 4:12)]
  expect_error(cast_predict_tiled(fx$fit, changed, td, tile_size = c(2, 3)), "does not match")
  fx$fit$models$external$predict_function <- function(model, newdata) plogis(-newdata$x)
  expect_error(cast_predict_tiled(fx$fit, fx$r, td, tile_size = c(2, 3)), "does not match")
  overwrite <- cast_predict_tiled(fx$fit, changed, td, tile_size = c(2, 3), overwrite = TRUE, verbose = FALSE)
  expect_equal(terra::values(terra::rast(overwrite$rasters[[1]]), mat = FALSE)[1:3], plogis(-raw[3:1]), tolerance = 1e-7)
  expect_error(cast_predict_tiled(fx$fit, fx$r, td, tile_size = 1.5), "whole numbers")
  expect_error(cast_predict_tiled(fx$fit, fx$r, td, models = c("external", "unknown")), "unique fitted")
})

test_that("failed tiles resume without publishing a partial output", {
  skip_if_not_installed("terra")
  fx <- raster_contract_fixture()
  fail <- TRUE; calls <- 0L
  local_mocked_bindings(predict_single_model = function(model, X, ...) {
    calls <<- calls + 1L
    if (fail && calls == 2L) stop("tile engine failure")
    plogis(X$x)
  })
  td <- tempfile(); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  expect_error(cast_predict_tiled(fx$fit, fx$r, td, tile_size = c(1, 4), verbose = FALSE), "tile engine failure")
  expect_false(file.exists(file.path(td, "HSS_external.tif")))
  expect_length(list.files(td, pattern = "^\\.cast-raster-", all.files = TRUE), 0)
  expect_length(list.files(td, pattern = "[.]rds$", recursive = TRUE, all.files = TRUE), 1)
  fail <- FALSE; before <- calls
  out <- cast_predict_tiled(fx$fit, fx$r, td, tile_size = c(1, 4), verbose = FALSE)
  expect_equal(calls - before, 2L)
  expect_true(file.exists(out$rasters[[1]]))
  # A readable but altered checkpoint must be recomputed, not replayed.
  files <- list.files(td, pattern = "[.]rds$", recursive = TRUE, all.files = TRUE, full.names = TRUE)
  tile_file <- files[grepl("1_1[.]rds$", files)]
  checkpoint <- readRDS(tile_file); checkpoint$values[] <- .9
  saveRDS(checkpoint, tile_file)
  unlink(out$rasters[[1]])
  before <- calls
  cast_predict_tiled(fx$fit, fx$r, td, tile_size = c(1, 4), verbose = FALSE)
  expect_equal(calls - before, 1L)
})

test_that("mask read failure aborts before publication and preserves existing output", {
  skip_if_not_installed("terra")
  fx <- raster_contract_fixture()
  td <- tempfile(); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res <- cast_ensemble_raster(fx$fit, fx$cv, fx$r, td, method = "equal", extrapolation = FALSE, verbose = FALSE)
  hash <- tools::md5sum(res$hss_path)
  local_mocked_bindings(.cast_mask_rows = function(...) stop("Mask read failed"))
  expect_error(cast_ensemble_raster(fx$fit, fx$cv, fx$r, td, method = "equal",
    overwrite = TRUE, extrapolation = FALSE, verbose = FALSE), "Mask read failed")
  expect_equal(tools::md5sum(res$hss_path), hash)
  expect_length(list.files(td, pattern = "^\\.cast-raster-", all.files = TRUE), 0)
})

test_that("a failed replacement restores the original file without copying", {
  td <- tempfile(); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  old <- file.path(td, "old.tif"); new <- file.path(td, "new.tif")
  writeLines("old complete output", old); writeLines("new complete output", new)
  rename <- base::file.rename
  local_mocked_bindings(file.rename = function(from, to) {
    if (identical(from, new)) return(FALSE)
    rename(from, to)
  }, .package = "base")
  expect_error(.cast_publish_raster(new, old, overwrite = TRUE), "Cannot publish")
  expect_equal(readLines(old), "old complete output")
  expect_equal(readLines(new), "new complete output")
  expect_length(list.files(td, pattern = "cast-backup", all.files = TRUE), 0)
})

test_that("year ordering puts finite dates first and keeps tie order", {
  expect_equal(.cast_thin_points(rep(0, 4), rep(0, 4), 1, order = "year_desc",
    year = c(NA, 1990, 2020, Inf)), c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(.cast_thin_points(rep(0, 3), rep(0, 3), 1, order = "year_desc",
    year = rep(NA_real_, 3)), c(TRUE, FALSE, FALSE))
  expect_equal(.cast_thin_points(c(0, 0), c(0, 0), 1, order = "year_desc",
    year = c(2020, 2020)), c(TRUE, FALSE))
  expect_error(.cast_thin_points(0, 0, 1, order = "year_desc", year = "2020"), "year vector")
  data <- data.frame(lon = c(NA, NA, 100), lat = c(NA, NA, 30))
  expect_equal(cast_thin(data, seed = 1, verbose = FALSE), data[3, , drop = FALSE], ignore_attr = TRUE)
})

test_that("callback fingerprints track used closure values only", {
  make <- function(slope) function(model, newdata) plogis(slope * newdata$x)
  a <- make(1); b <- make(2)
  expect_false(identical(.cast_digest(.cast_callback_identity(a)), .cast_digest(.cast_callback_identity(b))))
  before <- .cast_digest(.cast_callback_identity(a))
  environment(a)$unrelated <- "saved output"
  expect_equal(.cast_digest(.cast_callback_identity(a)), before)
})
