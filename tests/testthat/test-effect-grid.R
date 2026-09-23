# Single-shift effect with hard range-masking ----

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

test_that("removed products abort with a forwarding error", {
  expect_error(cast_dose_response(), "removed in 0.12.0")
  expect_error(cast_effect_support(), "removed in 0.12.0")
  expect_error(cast_necessity(), "removed in 0.12.0")
  expect_error(cast_sensitivity(), "removed in 0.12.0")
  expect_error(plot(structure(list(), class = "cast_dose_response")),
               "removed in 0.12.0")
})

test_that("the shift table carries one raw shift per driver and masks", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  tab <- cast_effect_table(fit, d, shift = list(x1 = 1, x2 = 1), verbose = FALSE)
  expect_s3_class(tab, "cast_effect_table")
  expect_true(all(c("driver", "shift_raw", "mean_abs_dHSS", "mean_signed_dHSS",
                    "n", "n_supported", "support", "masked") %in% names(tab)))
  expect_equal(tab$shift_raw, c(1, 1))
  # A positive raw shift moves with the data-generating signs.
  expect_lt(tab$mean_signed_dHSS[tab$driver == "x2"], 0)
  expect_true(all(tab$mean_abs_dHSS >= 0, na.rm = TRUE))
  expect_false(any(tab$masked))
  expect_equal(tab$support, tab$n_supported / tab$n)
  # The true drivers must move the fit more than a pure noise variable.
  noise <- data.frame(presence = d$presence, x1 = d$x1, x2 = d$x2,
                      noise = stats::rnorm(nrow(d)))
  screen_n <- new_cast_select(c("x1", "x2", "noise"),
                              data.frame(variable = c("x1", "x2", "noise")),
                              method = "manual")
  fit2 <- cast_fit(noise, screen = screen_n, models = "rf", rf_ntree = 120,
                   seed = 77, verbose = FALSE)
  tab2 <- cast_effect_table(fit2, noise, verbose = FALSE)
  expect_gt(tab2$mean_abs_dHSS[tab2$driver == "x1"],
            tab2$mean_abs_dHSS[tab2$driver == "noise"])
})

test_that("an absurd shift is refused rather than answered by extrapolation", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  expect_error(cast_effect_table(fit, d, shift = list(x1 = 1e6, x2 = 1),
                                 verbose = FALSE), "training range")
})

test_that("table inputs are validated", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  expect_error(cast_effect_table(fit, d, drivers = "nope", verbose = FALSE), "Unknown")
  expect_error(cast_effect_table(fit, d, shift = 0, verbose = FALSE), "non-zero")
  expect_error(cast_effect_table(fit, d, shift = list(x1 = 1), verbose = FALSE), "lacks entries")
  expect_error(cast_effect_table(fit, d, shift = c(1, 2), verbose = FALSE), "length 1")
  expect_error(cast_effect_table(fit, d, min_support = 2, verbose = FALSE), "min_support")
  expect_error(cast_effect_table(fit, d, support_probs = c(0.9, 0.1),
                                 verbose = FALSE), "support_probs")
})

test_that("a fit without a training reference masks instead of ranking", {
  skip_if_not_installed("ranger")
  d <- make_grid_data()
  fit <- grid_fit(d)
  fit$scaling$reference <- NULL
  tab <- cast_effect_table(fit, d, verbose = FALSE)
  expect_true(all(tab$masked))
  expect_true(all(is.na(tab$mean_abs_dHSS)))
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

test_that("masking uses the same complete rows as the estimates", {
  local_mocked_bindings(.cast_predict_matrix = function(fit, X, models, ...) {
    matrix(plogis(X$x - X$z), ncol = 1)
  })
  fit <- new_cast_fit(models = list(a = list()), cast_vars = c("x", "z"),
    env_vars = c("x", "z"), scaling = list(sds = c(x = 1, z = 1),
      reference = data.frame(x = 0:10, z = 0:10), impute = c(x = 5, z = 5)))
  X <- data.frame(x = c(5, 5, -1, 10, NA), z = c(5, 11, 5, 5, 5))
  tab <- cast_effect_table(fit, X, "x", shift = 1, shift_type = "raw",
                           support_probs = c(0, 1), verbose = FALSE)
  # Rows: (5,5) covered, (5,11) baseline outside, (-1,5) shifted outside,
  # (10,5) shifted to 11 outside; NA row dropped.
  expect_equal(tab$n, 4L)
  expect_equal(tab$n_supported, 1L)
  expect_equal(tab$support, 0.25)
  expect_true(tab$masked)
  expect_true(is.na(tab$mean_abs_dHSS))
  # A permissive threshold unmasks the same supported row.
  tab2 <- cast_effect_table(fit, X, "x", shift = 1, shift_type = "raw",
                            support_probs = c(0, 1), min_support = 0.2,
                            verbose = FALSE)
  expect_false(tab2$masked)
  expect_equal(tab2$n_supported, 1L)
  expect_error(cast_effect_table(fit, X[FALSE, ], "x", verbose = FALSE), "No complete")
  expect_error(cast_effect_table(fit, data.frame(x = "5", z = 5), "x",
                                 verbose = FALSE), "numeric")
})

test_that("the shift table plots", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("ggplot2")
  d <- make_grid_data()
  fit <- grid_fit(d)
  expect_s3_class(plot(cast_effect_table(fit, d, verbose = FALSE)), "ggplot")
})
