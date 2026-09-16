test_that("ensemble cancellation is preserved between effect tables and rasters", {
  skip_if_not_installed("terra")
  # Opposing fitted response surfaces have exactly zero ensemble response.
  # This catches mean(abs(delta_model)) versus abs(mean(delta_model)).
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    cbind(plogis(X$x), plogis(-X$x))
  })
  fit <- new_cast_fit(models = list(a = list(), b = list()),
    cast_vars = "x", env_vars = "x",
    scaling = list(sds = c(x = 1), reference = data.frame(x = c(-2, 2))))
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2)
  names(r) <- "x"
  terra::values(r) <- c(-1, 0, 1, NA, 0.5, -0.5)
  tab <- cast_effect_table(fit, data.frame(x = terra::values(r)[, 1]), verbose = FALSE)
  em <- cast_effect_map(fit, r, block_rows = 1, verbose = FALSE)
  expect_equal(tab$mean_abs_dHSS, 0, tolerance = 1e-14)
  expect_equal(mean(terra::values(em[["absdHSS_x"]]), na.rm = TRUE),
               tab$mean_abs_dHSS, tolerance = 1e-14)
  expect_true(is.na(terra::values(em)[4, 1]))
  expect_equal(attr(tab, "rows_complete"), 5L)
  expect_gt(tab$outside_range_fraction, 0)
})

test_that("effect input failures are explicit rather than all-NA results", {
  fit <- new_cast_fit(models = list(a = list()), cast_vars = "x", env_vars = "x",
    scaling = list(sds = c(x = 1)))
  expect_error(cast_effect_table(fit, data.frame(x = 1), shifts = Inf), "finite")
  expect_error(cast_effect_table(fit, data.frame(x = 1), steps = list(x = Inf)), "finite")
  expect_error(cast_effect_table(fit, data.frame(x = 1), shift_type = "typo"), "arg")
  expect_error(cast_effect_table(fit, data.frame(y = 1)), "Missing fitted")
})

test_that("real multiengine maps agree with tabular ensemble magnitudes", {
  skip_if_not_installed("terra")
  skip_if_not_installed("ranger")
  skip_if_not_installed("gbm")
  set.seed(91)
  x <- rnorm(220); z <- rnorm(220)
  d <- data.frame(lon = runif(220), lat = runif(220), x = x, z = z,
                  presence = rbinom(220, 1, plogis(x - z)))
  fit <- cast_fit(d, models = c("rf", "brt"), rf_ntree = 40,
                  brt_n_trees = 60, seed = 92, verbose = FALSE)
  r <- terra::rast(nrows = 4, ncols = 5, nlyrs = 2)
  names(r) <- c("x", "z")
  terra::values(r) <- as.matrix(d[1:20, c("x", "z")])
  tab <- cast_effect_table(fit, d[1:20, c("x", "z")], verbose = FALSE)
  em <- cast_effect_map(fit, r, block_rows = 2, verbose = FALSE)
  for (v in c("x", "z")) {
    expect_equal(mean(terra::values(em[[paste0("absdHSS_", v)]])),
      tab$mean_abs_dHSS[tab$driver == v], tolerance = 1e-10)
  }
})
