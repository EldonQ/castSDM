test_that("actual raster and table contrasts agree and resume without recomputation", {
  fit <- real_fit(); grids <- real_grids(); x <- real_case()
  a <- as.data.frame(terra::values(grids$current)); b <- as.data.frame(terra::values(grids$future))
  valid <- complete.cases(a) & complete.cases(b)
  drivers <- x$variables
  w <- setNames(rep(1/length(fit$models), length(fit$models)), names(fit$models))
  folder <- tempfile("real_compare_")
  for (clamp in c(FALSE, TRUE)) {
    output <- paste0(folder, "_", clamp)
    raster <- cast_compare(fit, grids$current, list(actual = grids$future),
      groups = list(climate = drivers), weights = w, clamp = clamp,
      output_dir = output, verbose = FALSE)
    again <- cast_compare(fit, grids$current, list(actual = grids$future),
      groups = list(climate = drivers), weights = w, clamp = clamp,
      output_dir = output, verbose = FALSE)
    expect_true(all(again$events$cache_hit))
    expect_equal(raster$summary, again$summary)
    tabular <- cast_attribution(fit, a[valid, ], b[valid, ], drivers = drivers,
      weights = w, clamp = clamp, verbose = FALSE, diagnostics = TRUE)
    mapped <- cast_attribution(fit, grids$current, grids$future, drivers = drivers,
      weights = w, clamp = clamp, verbose = FALSE, diagnostics = TRUE,
      filename = file.path(output, "responses.tif"))
    values <- as.data.frame(terra::values(terra::rast(mapped$raster_file)))
    for (name in names(tabular$data)) expect_equal(values[[name]][valid],
      tabular$data[[name]], tolerance = 1e-12)
    projected <- as.data.frame(terra::values(terra::rast(raster$results$actual$raster_file)))
    expect_equal(projected$joint_change[valid], tabular$data$total_change, tolerance = 1e-12)
    expect_equal(tabular$data$total_change, tabular$data$future-tabular$data$current,
      tolerance = 1e-12)
  }
})

test_that("tiled predictions use the same explicit raw clamp", {
  fit <- real_fit(); grids <- real_grids()
  data <- as.data.frame(terra::values(grids$future)); valid <- complete.cases(data)
  for (clamp in c(FALSE, TRUE)) {
    paths <- cast_predict_tiled(fit, grids$future, output_dir = tempfile("real_tiles_"),
      tile_size = terra::nrow(grids$future), clamp = clamp, verbose = FALSE)
    expected <- cast_predict(fit, data[valid, ], clamp = clamp, extrapolation = FALSE)$predictions
    for (model in names(fit$models)) expect_equal(
      terra::values(terra::rast(paths$rasters[[model]]), mat = FALSE)[valid],
      # Tiled exports explicitly use IEEE float32; kernels are checked to
      # double precision above. This bound is the float32 machine epsilon.
      expected[[paste0("HSS_", model)]], tolerance = 2^(-23))
  }
})
