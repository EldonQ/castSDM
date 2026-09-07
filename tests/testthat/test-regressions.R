# Regression tests for defects fixed in 0.7.0 -------------------------------

test_that("change classes mark non-binary cells as NA, never empty string", {
  ch <- .cast_change_classes(
    cur_bin = c(1, 0, NA, 1, 0),
    fut_bin = c(1, NA, 1, NA, 0)
  )
  expect_identical(ch, c("stable_present", NA, NA, NA, "stable_absent"))
  expect_false(any(ch == "", na.rm = TRUE))
})

test_that("cast_vif rejects non-finite predictors with a targeted error", {
  skip_if_not_installed("car")
  dat <- data.frame(
    a = rnorm(50), b = rnorm(50), c = rnorm(50),
    d = c(1, rep(Inf, 49))
  )
  expect_error(cast_vif(dat, threshold = 10, verbose = FALSE),
               "Non-finite")
})

test_that("cast_predict reports missing predictors by name", {
  skip_if_not_installed("ranger")
  set.seed(1)
  n <- 100
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 2, verbose = FALSE)
  bad_grid <- data.frame(x1 = rnorm(5))
  expect_error(cast_predict(fit, bad_grid), "missing fitted predictor")
})



test_that("ensemble excludes models with non-finite predictions and warns", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(3)
  n <- 120
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, plogis(rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 40, seed = 4, verbose = FALSE)
  # Sabotage the gam model so its predictions are NA.
  fit$models$gam$model <- NULL

  cv <- new_cast_cv(
    metrics = data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.8, 0.7), auc_sd = c(0.1, 0.1),
      tss_mean = c(0.5, 0.4), tss_sd = c(0.1, 0.1),
      cbi_mean = c(0.4, 0.3), cbi_sd = c(0.1, 0.1),
      n_folds = c(3, 3), n_selected_mean = c(3, 3)
    ),
    fold_metrics = data.frame(), folds = integer(n), k = 3L,
    block_method = "grid",
    thresholds = c(rf = 0.5, gam = 0.5)
  )

  expect_warning(ens <- cast_ensemble(fit, cv, dat, method = "weighted"),
                 "non-finite")
  expect_true(all(is.finite(ens$predictions$hss_ensemble)))
  # rf keeps the full renormalised weight.
  expect_equal(unname(ens$weights["rf"]), 1)
})


