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

test_that("effect tables carry a positivity-linked support column", {
  skip_if_not_installed("ranger")
  set.seed(93)
  x <- rnorm(220); z <- rnorm(220)
  d <- data.frame(lon = runif(220), lat = runif(220), x = x, z = z,
                  presence = rbinom(220, 1, plogis(x - z)))
  fit <- cast_fit(d, models = "rf", rf_ntree = 40, seed = 94, verbose = FALSE)
  tab <- cast_effect_table(fit, d[1:60, c("x", "z")], verbose = FALSE)
  expect_true("support" %in% names(tab))
  expect_true(all(is.na(tab$support) |
                    (tab$support >= 0 & tab$support <= 1)))
})

test_that("effect-table and necessity plots render as ggplots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(95)
  x <- rnorm(200); z <- rnorm(200)
  d <- data.frame(lon = runif(200), lat = runif(200), x = x, z = z,
                  presence = rbinom(200, 1, plogis(x - z)))
  fit <- cast_fit(d, models = "rf", rf_ntree = 30, seed = 96, verbose = FALSE)
  tab <- cast_effect_table(fit, d[1:50, c("x", "z")], verbose = FALSE)
  expect_s3_class(plot(tab), "ggplot")
  nec <- cast_necessity(d, variables = c("x", "z"), k = 2, num_trees = 30,
                        seed = 97, verbose = FALSE)
  expect_s3_class(plot(nec), "ggplot")
})

test_that("prediction registers engine methods in a bare session", {
  skip_if_not_installed("ranger")
  set.seed(98)
  x <- rnorm(120)
  d <- data.frame(lon = runif(120), lat = runif(120), x = x,
                  presence = rbinom(120, 1, plogis(x)))
  fit <- cast_fit(d, models = "rf", rf_ntree = 20, seed = 99, verbose = FALSE)
  # Simulate a fresh session (e.g. a fit reloaded via readRDS): without the
  # engine namespace, S3 dispatch has no predict method and every contrast
  # silently degrades to NA.
  tryCatch(unloadNamespace("ranger"), error = function(e) NULL)
  skip_if_not_installed("ranger")  # still installed, just unloaded
  p <- predict_single_model(fit$models[["rf"]], d[1:5, "x", drop = FALSE])
  expect_true(all(is.finite(p)))
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
