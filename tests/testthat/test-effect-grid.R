# Effect grids and the support (positivity) rule ----------------------------

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

test_that("cast_dose_response reports a curve, a support column and a zero anchor", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  dr <- cast_dose_response(fit, "x1", shift = seq(-2, 2, by = 0.5))
  expect_s3_class(dr, "cast_dose_response")
  expect_true(all(c("shift", "shift_raw", "mean_delta", "mean_abs_delta",
                    "support", "estimable") %in% names(dr$curve)))
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
