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

test_that("curves and tables pair engines before averaging missing predictions", {
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    cbind(ifelse(X$x <= 0, 0.1, NA_real_), 0.8 + 0.01 * X$x)
  })
  fit <- new_cast_fit(models = list(a = list(), b = list()), cast_vars = "x",
    env_vars = "x", scaling = list(sds = c(x = 1),
      reference = data.frame(x = c(-2, 2))))
  X <- data.frame(x = 0)
  dr <- cast_dose_response(fit, "x", 1, "raw", X)
  tab <- cast_effect_table(fit, X, "x", steps = 1, verbose = FALSE)
  expect_equal(dr$curve$mean_delta, 0.01)
  expect_equal(dr$curve$mean_abs_delta, tab$mean_abs_dHSS)
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
  tab$support[] <- NA_real_
  p <- plot(tab)
  expect_true(all(p$data$box_coverage == "box coverage unknown"))
  expect_true(all(ggplot2::ggplot_build(p)$data[[2]]$shape == 4))
  tab$support <- NULL
  expect_true(all(plot(tab)$data$box_coverage == "box coverage unknown"))
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
  expect_equal(attr(em, "effect_table")$support, tab$support)
})

test_that("map coverage checks full baselines and uses per-shift global totals", {
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
  steps <- list(x = c(-1, 1), z = c(-1, 1, 6))
  em <- cast_effect_map(fit, r, steps = steps, block_rows = 1,
                        support_probs = c(0, 1), verbose = FALSE)
  expect_equal(bounds_calls, 1L)
  expect_equal(rows_calls, 6L)  # Two drivers, three nonempty blocks.
  expect_equal(terra::nlyr(em), 3L * length(steps))
  expect_equal(names(em), c("dHSS_x", "dHSS_z", "absdHSS_x", "absdHSS_z",
                            "support_x", "support_z"))
  vals <- terra::values(em)
  # Row 2 fails on a non-driver; row 3 cannot enter the box from outside it.
  expect_equal(unname(vals[, "support_x"]), c(1, 0, 0, NA, 0.5, NA, NA, NA))
  expect_equal(unname(vals[, "support_z"]), c(2/3, 0, 0, NA, 2/3, NA, NA, NA))
  expect_true(all(is.na(vals[c(4, 6, 7, 8), ])))
  expect_equal(unname(vals[c(1, 2, 3, 5), "dHSS_x"]), rep(0.03, 4))
  expect_equal(unname(vals[c(1, 2, 3, 5), "dHSS_z"]), rep(-0.02 * 8/3, 4))
  expect_equal(attr(em, "support_probs"), c(0, 1))
  expect_equal(attr(em, "steps"), steps)

  tab <- cast_effect_table(fit, X, steps = steps, support_probs = c(0, 1),
                           verbose = FALSE)
  map_tab <- attr(em, "effect_table")
  expect_equal(map_tab, as.data.frame(tab)[, names(map_tab)])
  expect_equal(map_tab$support, c(0.25, 0))
  # Mean cell coverage averages across shifts, not the worst shift.
  expect_equal(mean(vals[, "support_x"], na.rm = TRUE), 0.375)
  expect_false(isTRUE(all.equal(map_tab$support[1],
                                mean(vals[, "support_x"], na.rm = TRUE))))
  for (v in names(steps)) {
    dr <- cast_dose_response(fit, v, steps[[v]], "raw", X, c(0, 1))
    sp <- cast_effect_support(fit, v, steps[[v]], "raw", X, c(0, 1))
    idx <- match(v, map_tab$driver)
    expect_equal(dr$curve$support, sp$support$support)
    expect_equal(map_tab$support[idx], min(dr$curve$support))
    expect_equal(map_tab$mean_abs_dHSS[idx], mean(dr$curve$mean_abs_delta))
    expect_equal(map_tab$mean_signed_dHSS[idx],
                  mean(dr$curve$mean_delta * sign(steps[[v]])))
  }
  whole <- cast_effect_map(fit, r, steps = steps, block_rows = 4,
                           support_probs = c(0, 1), verbose = FALSE)
  expect_equal(terra::values(whole), vals)
  expect_equal(attr(whole, "effect_table"), map_tab)
})

test_that("one-cell maps retain shape for one or multiple shifts and SD units", {
  skip_if_not_installed("terra")
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    matrix(0.1 + 0.01 * X$x, ncol = 1)
  })
  fit <- new_cast_fit(models = list(a = list()), cast_vars = "x", env_vars = "x",
    scaling = list(sds = c(x = 2), reference = data.frame(x = 0:10)))
  r <- terra::rast(nrows = 1, ncols = 1)
  names(r) <- "x"
  terra::values(r) <- 10
  shifts <- list(-0.5, 0.5, c(-0.5, 0.5))
  coverage <- c(1, 0, 0.5)
  worst_shift <- c(1, 0, 0)
  for (i in seq_along(shifts)) {
    em <- cast_effect_map(fit, r, shifts = shifts[[i]], support_probs = c(0, 1),
                          verbose = FALSE)
    expect_equal(dim(terra::values(em)), c(1L, 3L))
    expect_equal(unname(terra::values(em)[1, ]), c(0.01, 0.01, coverage[i]))
    expect_equal(attr(em, "effect_table")$support, worst_shift[i])
    expect_equal(attr(em, "steps")$x, shifts[[i]] * 2)
  }
  terra::values(r) <- 0
  em <- cast_effect_map(fit, r, steps = 1, verbose = FALSE)
  expect_equal(attr(em, "support_probs"), c(0.01, 0.99))
  # The baseline endpoint is outside the default quantile box, not the range.
  expect_equal(unname(terra::values(em)[1, ]), c(0.01, 0.01, 0))
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
  em <- cast_effect_map(fit, r, steps = 1, block_rows = 1,
                        support_probs = c(0, 1), verbose = FALSE)
  expect_equal(unname(terra::values(em)[, "support_x"]), c(1, 0, NA))
  expect_equal(unname(terra::values(em)[, "dHSS_x"]), c(0.01, NaN, NA))
  expect_equal(attr(em, "effect_table")$n, 1L)
  expect_equal(attr(em, "effect_table")$support, 0.5)
  tab <- cast_effect_table(fit, data.frame(x = c(0, 10, NA)), steps = 1,
                           support_probs = c(0, 1), verbose = FALSE)
  expect_equal(attr(em, "effect_table")$support, tab$support)
  expect_equal(attr(em, "effect_table")$mean_abs_dHSS, tab$mean_abs_dHSS)
})

test_that("unknown reference bounds and empty grids never imply coverage", {
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
    em <- cast_effect_map(fit, r, steps = 1, block_rows = 1, verbose = FALSE)
    expect_true(all(is.na(terra::values(em)[, "support_x"])))
    expect_true(is.na(attr(em, "effect_table")$support))
    expect_equal(unname(terra::values(em)[1, 1:2]), c(0.01, 0.01))
    expect_true(all(is.na(terra::values(em)[2, ])))
  }
  fit$scaling$reference <- data.frame(x = 0:10)
  terra::values(r) <- NA_real_
  em <- cast_effect_map(fit, r, steps = 1, block_rows = 1, verbose = FALSE)
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
  em <- cast_effect_map(fit, backed, c("x", "z"), c(-1, 1, 6), "raw", NULL,
                        1L, output, FALSE, FALSE, c(0, 1))
  disk <- terra::rast(output)
  expect_false(terra::inMemory(disk))
  expect_equal(names(disk), names(em))
  expect_equal(terra::nlyr(disk), 6L)
  expect_equal(terra::values(disk), terra::values(em), tolerance = 1e-6)
  expect_equal(unname(terra::values(disk)[, "support_x"]),
                c(2/3, 1/3, 0, NA, 2/3, 0), tolerance = 1e-6)
  expect_equal(unname(terra::values(disk)[, "support_z"]),
                c(2/3, 2/3, 0, NA, 2/3, 0), tolerance = 1e-6)
  expect_equal(attr(em, "effect_table")$support, c(0.2, 0.2))
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
      expect_error(cast_effect_map(fit, r, steps = 1, support_probs = bad,
                                   verbose = FALSE), "support_probs")
    }
  }
})

test_that("sensitivity GAM scenario backends agree for signed shifts and row order", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("marginaleffects")
  set.seed(104)
  d <- data.frame(lon = runif(180), lat = runif(180),
                  x = rnorm(180), z = rnorm(180))
  d$presence <- rbinom(180, 1, plogis(d$x - 0.4 * d$x^2 + d$z))
  fit <- cast_fit(d, models = "gam", verbose = FALSE)
  # An unselected engine must not change the single-GAM estimand.
  fit$models$rf <- list()
  X <- data.frame(easting = c(9, 2, 7, 9, 4, 1),
                  northing = c(3, 8, 6, 3, 5, 2),
                  x = c(-1.2, 0, 1.3, -1.2, NA, 2),
                  z = c(0.2, NA, -0.4, 0.2, 0.8, -0.7),
                  row.names = c("r9", "r2", "r7", "duplicate", "r4", "r1"))
  base <- .cast_impute(X[, c("x", "z")], fit$scaling$impute)
  p0 <- as.numeric(predict(fit$models$gam$model, base, type = "response"))
  for (unit in c("raw", "sd", "percent")) {
    for (shift in c(-0.6, 0, 0.6)) {
      if (unit == "percent") shift <- shift * 100
      step <- switch(unit, raw = shift, sd = shift * fit$scaling$sds[["x"]],
                     percent = shift / 100 * base$x)
      cf <- base
      cf$x <- cf$x + step
      p1 <- as.numeric(predict(fit$models$gam$model, cf, type = "response"))
      native <- cast_sensitivity(fit, X, "x", shift, unit, "gam",
                                  c("easting", "northing"))
      expect_identical(native, cast_sensitivity(fit, X, "x", shift, unit, "gam",
                         c("easting", "northing"), backend = "native"))
      for (direction in c("forward", "backward", "center")) {
        withr::local_options(marginaleffects_contrast_direction = direction)
        result <- cast_sensitivity(fit, X, "x", shift, unit, "gam",
                         c("easting", "northing"), backend = "marginaleffects")
        expect_s3_class(result, "cast_sensitivity")
        expect_equal(result, native, tolerance = 1e-10)
        expect_identical(result$predictions$lon, X$easting)
        expect_identical(result$predictions$lat, X$northing)
        expect_equal(result$predictions$baseline, p0, tolerance = 1e-10)
        expect_equal(result$predictions$counterfactual, p1, tolerance = 1e-10)
        expect_equal(result$predictions$delta_hss, p1 - p0, tolerance = 1e-10)
        expect_true(all(is.na(result$predictions$delta_sd)))
        expect_false(any(c("std.error", "conf.low", "conf.high") %in%
                           names(result$predictions)))
        if (shift == 0) expect_equal(result$predictions$delta_hss, rep(0, nrow(X)))
      }
    }
  }
  # One row exercises a length-one scenario, not a global numeric shortcut.
  one <- cast_sensitivity(fit, X[1, ], "x", -1, "raw", "gam",
                          c("easting", "northing"), backend = "marginaleffects")
  expect_equal(one, cast_sensitivity(fit, X[1, ], "x", -1, "raw", "gam",
                                     c("easting", "northing")), tolerance = 1e-10)
})

test_that("sensitivity marginaleffects requires an explicitly selected supported GAM", {
  fit <- new_cast_fit(models = list(gam = list(name = "gam",
    model = structure(list(), class = "gam")), rf = list(name = "rf", model = list())),
    cast_vars = "x", env_vars = "x", scaling = list(sds = c(x = 1)))
  X <- data.frame(lon = 1, lat = 2, x = 3)
  expect_error(cast_sensitivity(fit, X, "x", backend = "unknown"), "arg")
  for (model in list(NULL, character(), c("gam", "rf"), c("gam", "missing"),
                     c("gam", "gam"), "missing", "rf", NA_character_, 1)) {
    expect_error(cast_sensitivity(fit, X, "x", model = model,
                                  backend = "marginaleffects"), "explicitly selected\\s+fitted GAM")
  }
  fit$models$gam$model <- NULL
  expect_error(cast_sensitivity(fit, X, "x", model = "gam",
                                backend = "marginaleffects"), "fitted GAM")
  fit$models$gam$model <- structure(list(), class = "glm")
  expect_error(cast_sensitivity(fit, X, "x", model = "gam",
                                backend = "marginaleffects"), "fitted GAM")
})

test_that("sensitivity rejects malformed shifts and nonnumeric predictors before computation", {
  fit <- new_cast_fit(models = list(gam = list(name = "gam",
    model = structure(list(), class = "gam"))), cast_vars = c("x", "z"),
    env_vars = c("x", "z"), scaling = list(sds = c(x = 1, z = 1)))
  X <- data.frame(lon = 1, lat = 2, x = 3, z = 4)
  local_mocked_bindings(
    check_suggested = function(...) invisible(NULL),
    .cast_impute = function(...) stop("unsafe computation reached"),
    .package = "castSDM")
  for (backend in c("native", "marginaleffects")) {
    for (bad in list(NA_real_, NaN, Inf, -Inf, NULL, numeric(), c(1, 2),
                     "1", TRUE, list(1), 1i, matrix(1))) {
      expect_error(cast_sensitivity(fit, X, "x", shift = bad, model = "gam",
                                    backend = backend), "shift.*finite numeric")
    }
    for (v in c("x", "z")) {
      for (bad in list("3", factor("3"), TRUE)) {
        invalid <- X
        invalid[[v]] <- bad
        expect_error(cast_sensitivity(fit, invalid, "x", model = "gam",
                                      backend = backend), "Non-numeric predictor")
      }
    }
  }
})

test_that("sensitivity rejects nonfinite scenarios before prediction", {
  fit <- new_cast_fit(models = list(gam = list(name = "gam",
    model = structure(list(), class = "gam"))), cast_vars = "x", env_vars = "x",
    scaling = list(sds = c(x = 1)))
  X <- data.frame(lon = 1, lat = 2, x = 3)
  local_mocked_bindings(check_suggested = function(...) invisible(NULL),
    .cast_predict_matrix = function(...) stop("unsafe prediction reached"),
    .package = "castSDM")
  for (backend in c("native", "marginaleffects")) {
    for (bad in list(Inf, NA_real_, "1", c(1, 2), -1)) {
      fit$scaling$sds <- list(x = bad)
      expect_error(cast_sensitivity(fit, X, "x", model = "gam", backend = backend),
                   "predictor SD")
    }
    fit$scaling$sds <- c(x = 2)
    expect_error(cast_sensitivity(fit, X, "x", shift = .Machine$double.xmax,
                                  model = "gam", backend = backend), "finite scenario")
    huge <- X
    huge$x <- .Machine$double.xmax
    expect_error(cast_sensitivity(fit, huge, "x", shift = .Machine$double.xmax,
                     shift_type = "raw", model = "gam", backend = backend), "finite scenario")
    expect_error(cast_sensitivity(fit, huge, "x", shift = 200,
                     shift_type = "percent", model = "gam", backend = backend), "finite scenario")
    huge$x <- Inf
    expect_error(cast_sensitivity(fit, huge, "x", model = "gam", backend = backend),
                   "finite predictors")
  }
})

test_that("sensitivity marginaleffects dependency failures do not fall back", {
  fit <- new_cast_fit(models = list(gam = list(name = "gam",
    model = structure(list(), class = "gam"))), cast_vars = "x", env_vars = "x",
    scaling = list(sds = c(x = 1)))
  X <- data.frame(lon = 1, lat = 2, x = 3)
  local_mocked_bindings(
    check_suggested = function(pkg, ...) stop(paste("Package", pkg, "is required")),
    .cast_predict_matrix = function(...) stop("native fallback reached"),
    .package = "castSDM")
  expect_error(cast_sensitivity(fit, X, "x", model = "gam",
                                backend = "marginaleffects"), "Package marginaleffects is required")
})

test_that("sensitivity marginaleffects computation errors propagate without fallback", {
  skip_if_not_installed("marginaleffects")
  fit <- new_cast_fit(models = list(gam = list(name = "gam",
    model = structure(list(), class = "gam"))), cast_vars = "x", env_vars = "x",
    scaling = list(sds = c(x = 1)))
  X <- data.frame(lon = 1, lat = 2, x = 3)
  local_mocked_bindings(comparisons = function(...) stop("scenario computation failed"),
                        .package = "marginaleffects")
  local_mocked_bindings(.cast_predict_matrix = function(...) stop("native fallback reached"),
                        .package = "castSDM")
  expect_error(cast_sensitivity(fit, X, "x", model = "gam",
                                backend = "marginaleffects"), "scenario computation failed")
})
