# Quantile-box coverage is not a joint-positivity test.

make_grid_data <- function(n = 900, seed = 31) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  d <- data.frame(
    lon = stats::runif(n, 70, 130), lat = stats::runif(n, 20, 50),
    presence = stats::rbinom(n, 1, stats::plogis(1.4 * x1 - 0.8 * x2)),
    x1 = x1, x2 = x2)
  d
}

grid_fit <- function(d, models = "rf") {
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  cast_fit(d, screen = screen, models = models, rf_ntree = 120,
           seed = 77, verbose = FALSE)
}

test_that("cast_dose_response reports box coverage and drops zero shifts", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  dr <- cast_dose_response(fit, "x1", shift = seq(-2, 2, by = 0.5))
  expect_s3_class(dr, "cast_dose_response")
  expect_true(all(c("shift", "shift_raw", "mean_delta", "mean_abs_delta",
                    "support", "range_supported") %in% names(dr$curve)))
  # A zero shift is not an intervention and is dropped from the scan.
  expect_false(any(dr$curve$shift == 0))
  expect_true(all(dr$curve$mean_abs_delta >= 0))
  # The true driver must move the fit more than a pure noise variable.
  noise <- data.frame(presence = d$presence, x1 = d$x1, x2 = d$x2,
                      noise = stats::rnorm(nrow(d)))
  screen_n <- new_cast_select(c("x1", "x2", "noise"),
                              data.frame(variable = c("x1", "x2", "noise")),
                              method = "manual")
  fit2 <- cast_fit(noise, screen = screen_n, models = "rf", rf_ntree = 120,
                   seed = 77, verbose = FALSE)
  dr2 <- cast_dose_response(fit2, "noise", shift = 1)
  dr1 <- cast_dose_response(fit, "x1", shift = 1)
  expect_gt(dr1$curve$mean_abs_delta[1], dr2$curve$mean_abs_delta[1])
})

test_that("support decays as the shift grows and is bounded in [0, 1]", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  sp <- cast_effect_support(fit, c("x1", "x2"), shift = c(0.5, 2, 5))
  expect_s3_class(sp, "cast_support")
  expect_true(all(sp$support$support >= 0 & sp$support$support <= 1))
  s <- sp$support
  for (v in c("x1", "x2")) {
    z <- s[s$driver == v, ]
    # A larger shift can only remove observations from the training box.
    expect_true(z$support[z$shift == 0.5] >= z$support[z$shift == 5] - 1e-12)
  }
  expect_gt(s$support[s$driver == "x1" & s$shift == 0.5], 0.9)
})

test_that("an absurd shift is refused rather than answered by extrapolation", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  expect_error(cast_dose_response(fit, "x1", shift = 1e6), "training range")
})



test_that("grid inputs are validated against the fitted predictors", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  expect_error(cast_dose_response(fit, "nope"), "one fitted predictor")
  expect_error(cast_effect_support(fit, "nope"), "Unknown")
  expect_error(cast_dose_response(fit, "x1", shift = 0), "non-zero")
})

test_that("a fit without a training reference cannot claim support", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  fit$scaling$reference <- NULL
  expect_error(cast_effect_support(fit, "x1"), "reference")
})

test_that("effect grid objects plot", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("ggplot2")
  d <- make_grid_data()
  fit <- grid_fit(d)
  expect_s3_class(plot(cast_dose_response(fit, "x1", shift = c(-1, 1))), "ggplot")

  expect_s3_class(plot(cast_effect_support(fit, c("x1", "x2"))), "ggplot")
})

test_that("range coverage uses supplied rows and every baseline predictor", {
  ref <- data.frame(x = 0:10, z = 0:10, constant = 3)
  X <- data.frame(x = c(5, 5, -1, 10), z = c(5, 11, 5, 5), constant = 3)
  bounds <- .cast_support_bounds(ref, c(0, 1))
  rows <- .cast_support_rows(X, "x", c(-1, 1), bounds)
  expect_equal(rows, cbind(c(TRUE, FALSE, FALSE, TRUE),
                           c(TRUE, FALSE, FALSE, FALSE)))
  expect_equal(.cast_support_fraction(ref, "x", c(-1, 1),
                                      probs = c(0, 1), newdata = X), c(0.5, 0.25))
  expect_equal(.cast_support_rows(X[1, ], "x", 1, bounds), matrix(TRUE, 1, 1))
  expect_equal(.cast_support_rows(X[1, ], "constant", 1, bounds), matrix(FALSE, 1, 1))
  ref$x[] <- NA_real_
  expect_true(all(is.na(.cast_support_fraction(ref, "x", 1))))
})

test_that("quantile-box coverage does not certify conditional positivity", {
  ref <- data.frame(x = 0:10, z = 0:10)
  X <- data.frame(x = 5, z = 5)
  expect_equal(.cast_support_fraction(ref, "x", 1, newdata = X), 1)
  expect_false(any(ref$x == X$x + 1 & ref$z == X$z))
})

test_that("effect diagnostics use the same complete reference population", {
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    matrix(plogis(X$x - X$z), ncol = 1)
  })
  fit <- new_cast_fit(models = list(a = list()), cast_vars = c("x", "z"),
    env_vars = c("x", "z"), scaling = list(sds = c(x = 1, z = 1),
      reference = data.frame(x = 0:10, z = 0:10), impute = c(x = 5, z = 5)))
  X <- data.frame(x = c(5, 5, -1, 10, NA), z = c(5, 11, 5, 5, 5))
  sp <- cast_effect_support(fit, "x", c(-1, 1), "raw", X, c(0, 1))
  dr <- cast_dose_response(fit, "x", c(-1, 1), "raw", X, c(0, 1))
  tab <- cast_effect_table(fit, X, "x", steps = c(-1, 1),
                           support_probs = c(0, 1), verbose = FALSE)
  expect_equal(sp$support$support, c(0.5, 0.25))
  expect_equal(dr$curve$support, sp$support$support)
  expect_equal(dr$curve$range_supported, c(TRUE, FALSE))
  expect_equal(tab$support, min(sp$support$support))
  expect_equal(tab$mean_abs_dHSS, mean(dr$curve$mean_abs_delta))
  expect_equal(attr(tab, "rows_complete"), 4L)
  train <- cast_effect_support(fit, "x", 1, "raw", support_probs = c(0, 1))
  expect_equal(train$support$support, 10 / 11)
  limited <- cast_effect_support(fit, "x", 1, "raw", X, c(0, 1), max_rows = 1)
  expect_equal(limited$support$support, 1)
  unlimited <- cast_effect_support(fit, "x", 1, "raw", X, c(0, 1), max_rows = Inf)
  expect_equal(unlimited$support$support, 0.25)
  for (bad in list(NULL, NA_real_, 0, -1, 1.5, c(1, 2), "2", -Inf)) {
    expect_error(cast_effect_support(fit, "x", max_rows = bad), "max_rows")
    expect_error(cast_dose_response(fit, "x", max_rows = bad), "max_rows")
  }
  for (bad in list(NULL, c(0.9, 0.1), c(0.5, 0.5), c(-1, 1), c(0, 2),
                   c(0, Inf), c(NA, 1), 0.5, c("0", "1"))) {
    expect_error(cast_effect_support(fit, "x", support_probs = bad), "support_probs")
    expect_error(cast_dose_response(fit, "x", support_probs = bad), "support_probs")
    expect_error(cast_effect_table(fit, X, "x", support_probs = bad,
                                  verbose = FALSE), "support_probs")
  }
  expect_error(cast_effect_support(fit, "x", newdata = X[FALSE, ]), "No complete")
  expect_error(cast_effect_support(fit, "x", newdata = data.frame(x = "5", z = 5)),
               "numeric")
})
