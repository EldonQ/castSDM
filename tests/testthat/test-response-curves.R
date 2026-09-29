# P1 review round, W7: cast_response_curves() — fixed-profile and partial-
# dependence univariate curves plus bivariate interaction grids.

.w7_fit <- function() {
  skip_if_not_installed("ranger")
  set.seed(71)
  dat <- data.frame(
    presence = c(rep(1L, 60), rep(0L, 300)),
    x1 = c(runif(60, 0.6, 1), runif(300, 0, 1)),
    x2 = c(runif(60, 0.6, 1), runif(300, 0, 1)),
    x3 = runif(360)
  )
  cast_fit(dat, models = "rf", seed = 7, verbose = FALSE)
}

test_that("fixed response curves return a tidy long table with a signal", {
  fit <- .w7_fit()
  rc <- cast_response_curves(fit, variables = "x1", grid_size = 10)
  expect_s3_class(rc, "cast_response")
  expect_setequal(names(rc), c("variable", "value", "model", "prediction"))
  expect_equal(nrow(rc), 10L)
  expect_true(all(rc$variable == "x1"))
  expect_true(all(is.finite(rc$prediction)))
  # x1 drives presence: the high end of the curve must beat the low end.
  expect_gt(max(rc$prediction), min(rc$prediction))
})

test_that("pdp curves average over the background and stay finite", {
  fit <- .w7_fit()
  rc <- cast_response_curves(fit, variables = c("x1", "x2"),
                             type = "pdp", grid_size = 8)
  expect_equal(nrow(rc), 16L)  # 2 variables x 8 grid points x 1 model
  expect_setequal(unique(rc$variable), c("x1", "x2"))
  expect_true(all(is.finite(rc$prediction)))
})

test_that("bivariate pair grids cover the full cross product", {
  skip_if_not_installed("ranger")
  fit <- .w7_fit()
  rc <- cast_response_curves(fit, pair = c("x1", "x2"), pair_size = 6)
  expect_equal(nrow(rc), 36L)
  expect_setequal(names(rc),
                  c("var1", "value1", "var2", "value2", "model", "prediction"))
  expect_equal(length(unique(rc$value1)), 6L)
  expect_equal(length(unique(rc$value2)), 6L)
})

test_that("argument validation aborts informatively", {
  fit <- .w7_fit()
  expect_error(cast_response_curves(fit, variables = "nope"), "predictors")
  expect_error(cast_response_curves(fit, pair = c("x1", "x1")),
              "two distinct predictors")
  expect_error(cast_response_curves(fit, pair = c("x1", "nope")),
              "two distinct predictors")
  expect_error(cast_response_curves(fit, models = "brt"), "subset")
  expect_error(cast_response_curves(fit, grid_size = 1), "integer >= 2")
  expect_error(cast_response_curves(list()), "cast_fit")
})

test_that("plot.cast_response builds a ggplot for both layouts", {
  skip_if_not_installed("ggplot2")
  fit <- .w7_fit()
  rc <- cast_response_curves(fit, variables = "x1", grid_size = 10)
  p1 <- plot(rc)
  expect_s3_class(p1, "ggplot")
  rb <- cast_response_curves(fit, pair = c("x1", "x2"), pair_size = 6)
  p2 <- plot(rb)
  expect_s3_class(p2, "ggplot")
})
