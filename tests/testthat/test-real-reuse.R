test_that("sparse native reads preserve actual cell values, order and duplicates", {
  grid <- real_grids()$current
  i <- which(is.finite(terra::values(grid[[1L]], mat = FALSE)))
  indices <- c(rev(i), i)
  expected <- castSDM:::.cast_raster_block(grid, 1L, terra::nrow(grid))[indices, , drop = FALSE]
  actual <- castSDM:::.cast_raster_cells(grid, 1L, terra::nrow(grid), indices)
  expect_equal(actual, expected, ignore_attr = TRUE)
  expect_identical(names(actual), names(grid))
})

test_that("display requests prepare only the selected products", {
  core <- castSDM:::.cast_figure_display_request("core")
  expect_false(core$range)
  expect_identical(core$other_layers, "joint")
  expect_identical(core$driver_mode, "focus")
  expect_null(castSDM:::.cast_figure_display_request("all"))
  expect_true(castSDM:::.cast_figure_display_request("08")$range)
})

test_that("saved real response maps require matching data and numerical policy", {
  skip_if_not_installed("filelock")
  fit <- real_fit(); g <- real_grids()
  root <- tempfile("real-response-reuse-"); dir.create(file.path(root,"figures",".data"),recursive=TRUE)
  directory <- file.path(root,"figures",".data")
  attr <- cast_attribution(fit,g$current,g$future,drivers=fit$env_vars,
    models="rf",weights=c(rf=1),filename=file.path(directory,"driver_attribution.tif"),verbose=FALSE)
  focus <- attr$drivers[1L]
  display <- castSDM:::.cast_figure_driver_display(attr, file.path(directory,"focus.tif"),
    resolution=.02, drivers=focus, include_diagnostics=FALSE)
  expect_setequal(names(terra::rast(display)),paste0(c("effect_","retained_"),focus))
  preview <- castSDM:::.cast_explore_preview(attr$raster_file, file.path(directory,"preview.tif"))
  expect_true(all(c("total_change", "supported_total_change", "support_total") %in% names(terra::rast(preview))))
  expect_true(terra::ncell(terra::rast(preview)) <= 250000L)
  project <- list(directory=root,fit=fit,comparison=list(weights=c(rf=1),clamp=FALSE,
    support_probs=c(.01,.99),results=list(report=list(raster_file=attr$raster_file))))
  memo <- new.env(parent=emptyenv())
  contract <- list(fit=castSDM:::.cast_fit_identity(fit),scenario="report",weights=c(rf=1),
    clamp=FALSE,support_probs=c(.01,.99),
    current=castSDM:::.cast_data_identity(g$current,memo),future=castSDM:::.cast_data_identity(g$future,memo),
    comparison=unname(tools::md5sum(attr$raster_file)))
  receipt <- file.path(directory,"driver_catalog.rds")
  castSDM:::.cast_cache_write(list(contract=contract,attribution=attr),receipt,
    castSDM:::.cast_digest(contract),attr$raster_file)
  expect_equal(castSDM:::.cast_explore_saved(project,"report",g$current,g$future,memo)$summary,attr$summary)
  project$comparison$clamp <- TRUE
  expect_null(castSDM:::.cast_explore_saved(project,"report",g$current,g$future,memo))
  expect_true(all(attr$model_summary_common$n == unique(attr$common_summary$n)))
})
