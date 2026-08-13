test_that("cast_cv stores fold selection frequency and cast_consensus aggregates", {
  skip_if_not_installed("ranger")
  set.seed(80)
  n <- 300
  x1 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.2 * x1)), x1 = x1,
    x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n), x5 = rnorm(n)
  )
  cv <- cast_cv(
    dat, select_method = "rf",
    select_args = list(max_candidates = 4, min_vars = 2, num_trees = 40),
    k = 3, models = "rf", rf_ntree = 40, seed = 81, verbose = FALSE
  )
  expect_true(all(c("variable", "freq") %in% names(cv$selection_freq)))
  expect_true(all(cv$selection_freq$freq >= 0 & cv$selection_freq$freq <= 1))

  cons <- cast_consensus(cv, threshold = 0.5)
  expect_s3_class(cons, "cast_select")
  expect_identical(cons$method, "consensus")
  expect_true(all(cons$selected %in% cv$selection_freq$variable))
  # a manual cast_cv without selection_freq falls back to cv$selections
  cv_manual <- new_cast_cv(
    metrics = cv$metrics, fold_metrics = cv$fold_metrics, folds = cv$folds,
    k = cv$k, block_method = cv$block_method, thresholds = cv$thresholds,
    selections = list(c("x1", "x2"), c("x1", "x3"), c("x1", "x2"))
  )
  cons2 <- cast_consensus(cv_manual, threshold = 2 / 3)
  expect_setequal(cons2$selected, c("x1", "x2"))
})
