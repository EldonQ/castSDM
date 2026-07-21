test_that("RF benchmark selector returns importance scores", {
  skip_if_not_installed("ranger")
  set.seed(11)
  n <- 300
  x1 <- rnorm(n); x2 <- rnorm(n)
  noise <- replicate(6, rnorm(n))
  y <- rbinom(n, 1, plogis(1.8 * x1 - 1.2 * x2))
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = y, x1 = x1, x2 = x2, noise
  )
  out <- cast_select(
    dat, method = "rf", max_candidates = 6, min_vars = 2,
    num_trees = 60, seed = 12, verbose = FALSE
  )
  expect_s3_class(out, "cast_select")
  expect_identical(out$method, "rf")
  expect_gte(length(out$selected), 2)
  expect_true(all(out$selected %in% names(dat)))
  expect_true("rf_importance" %in% names(out$scores))
})

test_that("DML selector keeps FDR-significant predictors and exposes effects", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")
  set.seed(11)
  n <- 400
  x1 <- rnorm(n); x2 <- rnorm(n)
  noise <- replicate(5, rnorm(n))
  y <- rbinom(n, 1, plogis(2 * x1 - 1.5 * x2))
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = y, x1 = x1, x2 = x2, noise
  )
  out <- cast_select(
    dat, method = "dml", alpha = 0.05, dml_folds = 3L,
    num_trees = 60, min_vars = 1, seed = 12, verbose = FALSE
  )
  expect_s3_class(out, "cast_select")
  expect_identical(out$method, "dml")
  expect_true(all(c("estimate", "p_adjusted", "abs_statistic") %in% names(out$scores)))
  expect_true(any(c("x1", "x2") %in% out$selected))

  eff <- cast_effect(out, conf_level = 0.95)
  expect_s3_class(eff, "cast_effect")
  expect_true(all(c("conf_low", "conf_high", "p_adjusted") %in% names(eff$effects)))
})

test_that("cast_counterfactual maps a what-if shift on the current data", {
  skip_if_not_installed("ranger")
  set.seed(13)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1)), x1 = x1, x2 = x2
  )
  screen <- new_cast_select(c("x1", "x2"), data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 60,
                  seed = 14, verbose = FALSE)
  cf <- cast_counterfactual(fit, dat, variable = "x1", shift = 1,
                            shift_type = "sd")
  expect_s3_class(cf, "cast_counterfactual")
  expect_true(all(c("lon", "lat", "baseline", "counterfactual", "delta_hss")
                  %in% names(cf$predictions)))
  expect_equal(nrow(cf$predictions), n)
  expect_true(is.finite(cf$summary$mean_delta))
})
