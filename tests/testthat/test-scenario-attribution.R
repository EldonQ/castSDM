attribution_test_fit <- function() {
  new_cast_fit(models = list(a = list(), b = list()), cast_vars = c("x", "z"),
    env_vars = c("x", "z"), scaling = list(reference = data.frame(x = c(-2, 2), z = c(-2, 2))))
}

test_that("paired scenarios use exact frozen weights and preserve joint response", {
  local_mocked_bindings(predict_single_model = function(model, X) {
    if (identical(model$id, "b")) 0.2 + .1 * X$x else .3 + .1 * X$x + .05 * X$x * X$z
  })
  f <- attribution_test_fit(); f$models$b$id <- "b"
  a <- data.frame(x = c(0, .5), z = c(0, 1)); b <- a + .2
  z <- cast_attribution(f, a, b, weights = c(a = .8, b = .2), verbose = FALSE)
  oracle <- function(X) .8 * (.3 + .1 * X$x + .05 * X$x * X$z) + .2 * (.2 + .1 * X$x)
  only <- a; only$x <- b$x
  held <- b; held$x <- a$x
  expect_equal(z$data$total_change, oracle(b) - oracle(a))
  expect_equal(z$data$effect_x, oracle(only) - oracle(a))
  expect_equal(z$data$retained_x, oracle(b) - oracle(held))
  expect_equal(z$data$remainder, rep(.8 * .05 * .2 * .2, 2))
  expect_equal(z$summary$mean_signed[1], mean(z$data$total_change))
})

test_that("raster and table attribution agree across block sizes and missing cells", {
  skip_if_not_installed("terra")
  local_mocked_bindings(predict_single_model = function(model, X) plogis(X$x - X$z))
  f <- attribution_test_fit()
  r <- terra::rast(nrows = 3, ncols = 2, nlyrs = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 3)
  names(r) <- c("x", "z")
  a <- data.frame(x = c(-1, 0, 1, NA, .5, -.5), z = c(0, .2, .4, 0, .8, -1))
  b <- a; b$x <- b$x + .3
  terra::values(r) <- as.matrix(a); s <- r; terra::values(s) <- as.matrix(b)
  table <- cast_attribution(f, a, b, verbose = FALSE)
  table_clamped <- cast_attribution(f, a, b, clamp = TRUE, verbose = FALSE)
  expect_equal(table_clamped$summary, table$summary, tolerance = 1e-12)
  m1 <- cast_attribution(f, r, s, block_size = 2, verbose = FALSE)
  m2 <- cast_attribution(f, r, s, block_size = 20, verbose = FALSE)
  expect_equal(terra::values(terra::rast(m1$raster_file)), as.matrix(table$data), tolerance = 1e-12)
  expect_equal(terra::values(terra::rast(m1$raster_file)), terra::values(terra::rast(m2$raster_file)))
  expect_equal(m1$summary, table$summary, tolerance = 1e-12)
  expect_true(all(is.na(table$data[4, ])))
  expect_equal(table$data$effect_z[is.finite(table$data$effect_z)], rep(0, 5))
  expect_equal(table$data$effect_x, table$data$total_change)
  shifted <- terra::shift(s, dx = 1)
  expect_error(cast_attribution(f, r, shifted, verbose = FALSE))
})

test_that("alignment, range flags and failure checks are explicit", {
  local_mocked_bindings(predict_single_model = function(model, X) plogis(X$x))
  f <- attribution_test_fit()
  a <- data.frame(lon = c(1, 2), lat = c(1, 2), x = c(0, 1), z = 0)
  b <- a; b$x <- c(.2, 3)
  z <- cast_attribution(f, a, b[2:1, ], verbose = FALSE)
  expect_equal(z$data$future, plogis(b$x))
  expect_equal(z$data$support_total, c(1, 0))
  expect_true(all(is.finite(z$data$total_change)))
  expect_error(cast_attribution(f, a, b, weights = c(a = 0, b = 0)), "positive")
  expect_error(cast_attribution(f, a, b, weights = c(a = 1)), "weights")
  expect_error(cast_attribution(f, a, b, drivers = "q"), "drivers")
  expect_error(cast_attribution(f, a, b, block_size = 0), "block_size")
  local_mocked_bindings(predict_single_model = function(model, X) stop("broken engine"))
  expect_error(cast_attribution(f, a, b), "broken engine")
})
