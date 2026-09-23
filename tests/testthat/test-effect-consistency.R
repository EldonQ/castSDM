# Single-shift effects agree between tables and rasters; masking is hard ----

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
  tab <- cast_effect_table(fit, data.frame(x = terra::values(r)[, 1]),
                           shift = 1, verbose = FALSE)
  em <- cast_effect_map(fit, r, shift = 1, block_rows = 1, verbose = FALSE)
  expect_equal(tab$mean_abs_dHSS, 0, tolerance = 1e-14)
  expect_equal(mean(terra::values(em[["absdHSS_x"]]), na.rm = TRUE),
               tab$mean_abs_dHSS, tolerance = 1e-14)
  expect_true(is.na(terra::values(em)[4, 1]))
  expect_equal(attr(tab, "rows_complete"), 5L)
  # x = 1 shifted by +1 leaves the box, so it is masked, not ranked.
  expect_equal(tab$n_supported, 4L)
  expect_equal(tab$support, 0.8)
  expect_equal(unname(terra::values(em)[, "support_x"]), c(1, 1, 0, NA, 1, 1))
})

test_that("tables pair engines before averaging missing predictions", {
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    cbind(ifelse(X$x <= 0, 0.1, NA_real_), 0.8 + 0.01 * X$x)
  })
  fit <- new_cast_fit(models = list(a = list(), b = list()), cast_vars = "x",
    env_vars = "x", scaling = list(sds = c(x = 1),
      reference = data.frame(x = c(-2, 2))))
  X <- data.frame(x = 0)
  tab <- cast_effect_table(fit, X, "x", shift = 1, shift_type = "raw",
                           verbose = FALSE)
  # Engine a is NA after the shift; engine b moves 0.01. Paired averaging
  # over finite engines gives 0.01, not an all-NA contrast.
  expect_equal(tab$mean_signed_dHSS, 0.01)
  expect_equal(tab$mean_abs_dHSS, 0.01)
})

test_that("effect tables carry a quantile-box coverage column", {
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
  expect_true("masked" %in% names(tab))
})

test_that("effect-table plots render as ggplots, masked as crosses", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ranger")
  set.seed(95)
  x <- rnorm(200); z <- rnorm(200)
  d <- data.frame(lon = runif(200), lat = runif(200), x = x, z = z,
                  presence = rbinom(200, 1, plogis(x - z)))
  fit <- cast_fit(d, models = "rf", rf_ntree = 30, seed = 96, verbose = FALSE)
  tab <- cast_effect_table(fit, d[1:50, c("x", "z")], verbose = FALSE)
  expect_s3_class(plot(tab), "ggplot")
  tab$masked[] <- TRUE
  tab$mean_abs_dHSS[] <- NA_real_
  tab$mean_signed_dHSS[] <- NA_real_
  p <- plot(tab)
  expect_true(all(p$data$mask == "masked"))
  expect_true(all(ggplot2::ggplot_build(p)$data[[2]]$shape == 4))
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
  expect_error(cast_effect_table(fit, data.frame(x = 1), shift = Inf), "finite")
  expect_error(cast_effect_table(fit, data.frame(x = 1),
                                 shift = list(x = Inf)), "finite")
  expect_error(cast_effect_table(fit, data.frame(x = 1), shift_type = "typo"), "arg")
  expect_error(cast_effect_table(fit, data.frame(y = 1)), "Missing fitted")
})

test_that("per-driver shifts keep their names through SD conversion", {
  fit <- new_cast_fit(models = list(a = list()), cast_vars = c("x", "z"),
    env_vars = c("x", "z"),
    scaling = list(sds = c(x = 2, z = 0.5),
      reference = data.frame(x = 0:10, z = 0:10)))
  expect_identical(.shift_from_fit(fit, c("x", "z"), 0.5, "sd"),
                   list(x = 1, z = 0.25))
  expect_identical(.shift_from_fit(fit, c("x", "z"), c(x = 1, z = 2), "raw"),
                   list(x = 1, z = 2))
  expect_identical(.shift_from_fit(fit, c("x", "z"), 3, "raw"),
                   list(x = 3, z = 3))
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
  tab <- cast_effect_table(fit, d[1:20, c("x", "z")], shift = 1,
                           min_support = 0, verbose = FALSE)
  em <- cast_effect_map(fit, r, shift = 1, block_rows = 2, min_support = 0,
                        verbose = FALSE)
  for (v in c("x", "z")) {
    expect_equal(mean(terra::values(em[[paste0("absdHSS_", v)]]), na.rm = TRUE),
      tab$mean_abs_dHSS[tab$driver == v], tolerance = 1e-10)
  }
  expect_equal(attr(em, "effect_table")$support, tab$support)
})

test_that("map coverage checks full baselines and masks per cell", {
  skip_if_not_installed("terra")
  bounds_calls <- 0L
  rows_calls <- 0L
  support_bounds <- .cast_support_bounds
  support_rows <- .cast_support_rows
  local_mocked_bindings(
    .cast_predict_matrix = function(fit, X, models, ...) {
      cbind(0.4 + 0.02 * X$x - 0.01 * X$z,
            0.3 + 0.04 * X$x - 0.03 * X$z)
    },
    .cast_support_bounds = function(ref, probs) {
      bounds_calls <<- bounds_calls + 1L
      support_bounds(ref, probs)
    },
    .cast_support_rows = function(X, driver, steps, bounds) {
      rows_calls <<- rows_calls + 1L
      support_rows(X, driver, steps, bounds)
    })
  fit <- new_cast_fit(models = list(a = list(), b = list()),
    cast_vars = c("x", "z"), env_vars = c("x", "z"),
    scaling = list(sds = c(x = 1, z = 1),
      reference = data.frame(x = 0:10, z = 0:10), impute = c(x = 5, z = 5)))
  X <- data.frame(x = c(5, 5, -1, NA, 10, 5, 5, 5),
                  z = c(5, 11, 5, 5, 5, NA, Inf, NaN))
  r <- terra::rast(nrows = 4, ncols = 2, nlyrs = 2)
  names(r) <- c("x", "z")
  terra::values(r) <- as.matrix(X)
  em <- cast_effect_map(fit, r, shift = list(x = 1, z = 2), block_rows = 1,
                        support_probs = c(0, 1), verbose = FALSE)
  expect_equal(bounds_calls, 1L)
  expect_equal(rows_calls, 6L)  # Two drivers, three nonempty blocks.
  expect_equal(terra::nlyr(em), 6L)
  expect_equal(names(em), c("dHSS_x", "dHSS_z", "absdHSS_x", "absdHSS_z",
                            "support_x", "support_z"))
  vals <- terra::values(em)
  # Row 2 fails on a non-driver baseline; row 3 starts outside the box.
  expect_equal(unname(vals[, "support_x"]), c(1, 0, 0, NA, 0, NA, NA, NA))
  expect_equal(unname(vals[, "support_z"]), c(1, 0, 0, NA, 1, NA, NA, NA))
  expect_true(all(is.na(vals[c(4, 6, 7, 8), ])))
  # Masked cells never contribute: only supported cells carry deltas.
  expect_equal(unname(vals[c(1, 5), "dHSS_x"]), c(0.03, NA))
  expect_equal(unname(vals[c(1, 5), "dHSS_z"]), c(-0.04, -0.04))
  expect_equal(attr(em, "support_probs"), c(0, 1))
  expect_equal(attr(em, "steps"), list(x = 1, z = 2))
  expect_equal(attr(em, "min_support"), 0.5)

  tab <- cast_effect_table(fit, X, shift = list(x = 1, z = 2),
                           support_probs = c(0, 1), verbose = FALSE)
  map_tab <- attr(em, "effect_table")
  expect_equal(map_tab$support, tab$support)
  expect_equal(map_tab$support, c(0.25, 0.5))
  # x is masked at the default threshold, z is not.
  expect_true(map_tab$masked[map_tab$driver == "x"])
  expect_false(map_tab$masked[map_tab$driver == "z"])
  expect_true(is.na(map_tab$mean_abs_dHSS[map_tab$driver == "x"]))
  expect_equal(map_tab$mean_signed_dHSS[map_tab$driver == "z"], -0.04)
  whole <- cast_effect_map(fit, r, shift = list(x = 1, z = 2), block_rows = 4,
                           support_probs = c(0, 1), verbose = FALSE)
  expect_equal(terra::values(whole), vals)
  expect_equal(attr(whole, "effect_table"), map_tab)
})

test_that("one-cell maps mask shifts that leave the box", {
  skip_if_not_installed("terra")
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    matrix(0.1 + 0.01 * X$x, ncol = 1)
  })
  fit <- new_cast_fit(models = list(a = list()), cast_vars = "x", env_vars = "x",
    scaling = list(sds = c(x = 2), reference = data.frame(x = 0:10)))
  r <- terra::rast(nrows = 1, ncols = 1)
  names(r) <- "x"
  terra::values(r) <- 10
  # -0.5 stays inside [0, 10]; +0.5 leaves it and is masked.
  em_in <- cast_effect_map(fit, r, shift = -0.5, support_probs = c(0, 1),
                           verbose = FALSE)
  expect_equal(unname(terra::values(em_in)[1, ]), c(-0.005, 0.005, 1))
  expect_false(attr(em_in, "effect_table")$masked)
  em_out <- cast_effect_map(fit, r, shift = 0.5, support_probs = c(0, 1),
                            verbose = FALSE)
  expect_equal(unname(terra::values(em_out)[1, ]), c(NA_real_, NA_real_, 0))
  expect_true(attr(em_out, "effect_table")$masked)
  # SD units convert through the stored scale: 0.5 SD = 1 raw unit.
  em_sd <- cast_effect_map(fit, r, shift = 0.5, shift_type = "sd",
                           support_probs = c(0, 1), verbose = FALSE)
  expect_equal(attr(em_sd, "steps")$x, 1)
  expect_true(attr(em_sd, "effect_table")$masked)
  terra::values(r) <- 0
  em <- cast_effect_map(fit, r, shift = 1, support_probs = c(0, 1),
                        verbose = FALSE)
  expect_equal(unname(terra::values(em)[1, ]), c(0.01, 0.01, 1))
})

test_that("map support counts complete cells even when paired effects fail", {
  skip_if_not_installed("terra")
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    cbind(ifelse(X$x <= 0, 0.1, NA_real_),
          ifelse(X$x < 8, 0.8 + 0.01 * X$x, NA_real_))
  })
  fit <- new_cast_fit(models = list(a = list(), b = list()),
    cast_vars = "x", env_vars = "x",
    scaling = list(reference = data.frame(x = -2:10)))
  r <- terra::rast(nrows = 3, ncols = 1)
  names(r) <- "x"
  terra::values(r) <- c(0, 10, NA)
  em <- cast_effect_map(fit, r, shift = 1, block_rows = 1,
                        support_probs = c(0, 1), verbose = FALSE)
  expect_equal(unname(terra::values(em)[, "support_x"]), c(1, 0, NA))
  expect_equal(unname(terra::values(em)[, "dHSS_x"]), c(0.01, NA, NA))
  expect_equal(attr(em, "effect_table")$n, 2L)
  expect_equal(attr(em, "effect_table")$support, 0.5)
  tab <- cast_effect_table(fit, data.frame(x = c(0, 10, NA)), shift = 1,
                           support_probs = c(0, 1), verbose = FALSE)
  expect_equal(attr(em, "effect_table")$support, tab$support)
  expect_equal(attr(em, "effect_table")$mean_abs_dHSS, tab$mean_abs_dHSS)
})

test_that("missing reference or degenerate bounds mask instead of ranking", {
  skip_if_not_installed("terra")
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    matrix(0.1 + 0.01 * X$x, ncol = 1)
  })
  fit <- new_cast_fit(models = list(a = list()), cast_vars = "x", env_vars = "x",
    scaling = list())
  r <- terra::rast(nrows = 2, ncols = 1)
  names(r) <- "x"
  terra::values(r) <- c(5, NA)
  for (ref in list(NULL, data.frame(z = 0:10), data.frame(x = NA_real_),
                   data.frame(x = numeric()))) {
    fit$scaling$reference <- ref
    em <- cast_effect_map(fit, r, shift = 1, block_rows = 1, verbose = FALSE)
    expect_true(all(is.na(terra::values(em)[, "support_x"])))
    expect_true(is.na(attr(em, "effect_table")$support))
    expect_true(attr(em, "effect_table")$masked)
    expect_true(all(is.na(terra::values(em)[2, ])))
  }
  fit$scaling$reference <- data.frame(x = 0:10)
  terra::values(r) <- NA_real_
  em <- cast_effect_map(fit, r, shift = 1, block_rows = 1, verbose = FALSE)
  expect_true(all(is.na(terra::values(em))))
  expect_identical(attr(em, "effect_table")$support, NA_real_)
  expect_equal(attr(em, "effect_table")$n, 0L)
})

test_that("file-backed maps preserve all TIFF layers and positional arguments", {
  skip_if_not_installed("terra")
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    cbind(0.4 + 0.02 * X$x - 0.01 * X$z,
          0.3 + 0.04 * X$x - 0.03 * X$z)
  })
  fit <- new_cast_fit(models = list(a = list(), b = list()),
    cast_vars = c("x", "z"), env_vars = c("x", "z"),
    scaling = list(reference = data.frame(x = 0:10, z = 0:10)))
  r <- terra::rast(nrows = 3, ncols = 2, nlyrs = 2)
  names(r) <- c("x", "z")
  terra::values(r) <- cbind(c(5, 10, -1, NA, 0, 5), c(5, 5, 5, 5, 0, 11))
  input <- tempfile(fileext = ".tif")
  output <- tempfile(fileext = ".tif")
  on.exit(unlink(c(input, output)), add = TRUE)
  terra::writeRaster(r, input)
  backed <- terra::rast(input)
  expect_false(terra::inMemory(backed))
  em <- cast_effect_map(fit, backed, c("x", "z"), list(x = 1, z = 2), "raw",
                        1L, output, FALSE, FALSE, c(0, 1), 0.5)
  disk <- terra::rast(output)
  expect_false(terra::inMemory(disk))
  expect_equal(names(disk), names(em))
  expect_equal(terra::nlyr(disk), 6L)
  expect_equal(terra::values(disk), terra::values(em), tolerance = 1e-6)
  expect_equal(unname(terra::values(disk)[, "support_x"]),
                c(1, 0, 0, NA, 1, 0), tolerance = 1e-6)
  expect_equal(unname(terra::values(disk)[, "support_z"]),
                c(1, 1, 0, NA, 1, 0), tolerance = 1e-6)
  expect_equal(attr(em, "effect_table")$support, c(0.4, 0.6))
  expect_equal(attr(em, "support_probs"), c(0, 1))
})

test_that("maps reject malformed support probabilities without a reference", {
  skip_if_not_installed("terra")
  fit <- new_cast_fit(models = list(a = list()), cast_vars = "x", env_vars = "x",
    scaling = list())
  r <- terra::rast(nrows = 1, ncols = 1)
  names(r) <- "x"
  terra::values(r) <- NA_real_
  for (ref in list(data.frame(x = 0:10), NULL)) {
    fit$scaling$reference <- ref
    for (bad in list(NULL, numeric(), c(0.9, 0.1), c(0.5, 0.5), c(-1, 1),
                     c(0, 2), c(0, Inf), c(NA, 1), c(NaN, 1), 0.5,
                     c(0, 0.5, 1), c("0", "1"), c(FALSE, TRUE))) {
      expect_error(cast_effect_map(fit, r, shift = 1, support_probs = bad,
                                   verbose = FALSE), "support_probs")
    }
  }
})
