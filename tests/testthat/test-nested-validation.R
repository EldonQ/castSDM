test_that("spatial CV stores fold-specific selections", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ranger")
  set.seed(21)
  n <- 300
  x1 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(x1)), x1 = x1,
    x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n), x5 = rnorm(n)
  )
  cv <- cast_cv(
    dat, select_method = "stable_no_invariance",
    select_args = list(max_vars = 4, min_vars = 2, num_trees = 40),
    k = 3, models = "rf", rf_ntree = 40, seed = 22, verbose = FALSE
  )
  expect_s3_class(cv, "cast_cv")
  expect_length(cv$selections, 3)
  expect_true(all(vapply(cv$selections, length, integer(1)) >= 2))
  expect_true("n_selected" %in% names(cv$fold_metrics))
})
