# AI-assisted maintenance (2026-11-09); scope and verification: inst/AI_ASSISTANCE.md.
test_that("missing physical inputs fail before comparison creates output", {
  skip_if_not_installed("terra")
  a <- data.frame(x = seq(-1, 1, length.out = 12), z = rep(0:1, 6))
  fit <- cast_as_fit(list(one = 1), a, function(model, newdata) stats::plogis(newdata$x))
  memory <- terra::rast(nrows = 3, ncols = 4, nlyrs = 2)
  terra::values(memory) <- as.matrix(a); names(memory) <- names(a)
  root <- tempfile("missing-input-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "relocated.tif")
  disk <- terra::writeRaster(memory, path, datatype = "FLT8S")
  unlink(path)
  groups <- list(first = "x", second = "z")
  output <- file.path(root, "output")
  expect_error(cast_compare(fit, disk, memory, groups, output_dir = output),
    "Missing file-backed raster source for current")
  expect_false(dir.exists(output))
  expect_error(cast_compare(fit, memory, list(moved_future = disk), groups, output_dir = output),
    "future scenario moved_future")
  expect_error(cast_compare(fit, memory, memory, groups, mask = disk[[1]], output_dir = output),
    "Missing file-backed raster source for mask")
  expect_false(dir.exists(output))
})
