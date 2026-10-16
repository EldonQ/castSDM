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

test_that("each response row agrees with its named model and grid value", {
  bg <- data.frame(x1 = c(0, 1, 2), x2 = c(0, 1, 5))
  fit <- new_cast_fit(
    models = list(a = list(model = "a"), b = list(model = "b")),
    cast_vars = names(bg), env_vars = names(bg), scaling = list(reference = bg))
  local_mocked_bindings(predict_single_model = function(model, X, ...) {
    if (model$model == "a") 0.1 + 0.05 * X$x1 + 0.02 * X$x2
    else 0.8 - 0.04 * X$x1 - 0.03 * X$x2
  })
  for (type in c("fixed", "pdp")) {
    for (models in list(c("a", "b"), c("b", "a"))) {
      for (size in c(3L, 5L)) {
        rc <- cast_response_curves(fit, type = type, grid_size = size,
                                   models = models)
        x1 <- ifelse(rc$variable == "x1", rc$value,
                     if (type == "fixed") median(bg$x1) else mean(bg$x1))
        x2 <- ifelse(rc$variable == "x2", rc$value,
                     if (type == "fixed") median(bg$x2) else mean(bg$x2))
        expected <- ifelse(rc$model == "a", 0.1 + 0.05 * x1 + 0.02 * x2,
                           0.8 - 0.04 * x1 - 0.03 * x2)
        expect_equal(rc$prediction, expected)
        expect_identical(rc$model, rep(rep(models, each = size), 2L))
        expect_equal(nrow(rc), 2L * size * length(models))
      }
    }
  }
})

test_that("fixed profiles predict one grid per engine and variable", {
  bg <- data.frame(x1 = c(0, 1, 2), x2 = c(0, 1, 5))
  fit <- new_cast_fit(
    models = list(a = list(model = "a"), b = list(model = "b")),
    cast_vars = names(bg), env_vars = names(bg), scaling = list(reference = bg))
  calls <- list()
  local_mocked_bindings(predict_single_model = function(model, X, ...) {
    calls[[length(calls) + 1L]] <<- X
    plogis(X$x1 - X$x2)
  })
  cast_response_curves(fit, grid_size = 5L)
  expect_length(calls, 4L)
  expect_true(all(vapply(calls, nrow, integer(1)) == 5L))
  expect_equal(calls[[1]]$x1, seq(0, 2, length.out = 5L))
  expect_equal(calls[[1]]$x2, rep(median(bg$x2), 5L))
})
