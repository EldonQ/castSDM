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
    dat, select_method = "two_stage",
    select_args = list(num_trees = 40, metric = "brier"),
    k = 3, models = "rf", rf_ntree = 40, seed = 81, verbose = FALSE
  )
  expect_true(all(c("variable", "freq") %in% names(cv$selection_freq)))
  expect_true(all(cv$selection_freq$freq >= 0 & cv$selection_freq$freq <= 1))
  expect_true(all(c("obs", "HSS_rf") %in% names(cv$oof)))
  expect_identical(nrow(cv$oof), nrow(dat))

  cons <- cast_consensus(cv, threshold = 0.5)
  expect_s3_class(cons, "cast_select")
  expect_identical(cons$method, "consensus")
  expect_true("x1" %in% cons$selected)
  # a manual cast_cv without selection_freq falls back to cv$selections
  cv_manual <- new_cast_cv(
    metrics = cv$metrics, fold_metrics = cv$fold_metrics, folds = cv$folds,
    k = cv$k, block_method = cv$block_method, thresholds = cv$thresholds,
    selections = list(c("x1", "x2"), c("x1", "x3"), c("x1", "x2"))
  )
  cons2 <- cast_consensus(cv_manual, threshold = 2 / 3)
  expect_setequal(cons2$selected, c("x1", "x2"))
})

test_that("CV reports finite fold counts separately for each metric", {
  skip_if_not_installed("ranger")
  metric_rows <- list(
    c(auc = 0.7, tss = 0.2, cbi = -0.3),
    c(auc = 0.6, tss = NA_real_, cbi = NA_real_),
    c(auc = NA_real_, tss = NA_real_, cbi = NA_real_)
  )
  evaluated <- 0L
  local_mocked_bindings(
    make_spatial_folds = function(lon, lat, k, method, seed) rep(1:4, each = 6),
    evaluate_model_full = function(pred, obs) {
      evaluated <<- evaluated + 1L
      metric_rows[[evaluated]]
    }
  )
  dat <- data.frame(lon = seq_len(24), lat = rep(1:6, 4),
                    presence = c(rep(0:1, 9), rep(0, 6)),
                    x1 = seq_len(24), x2 = rep(1:4, 6), x3 = rep(1:6, 4))
  expect_warning(
    cv <- cast_cv(dat, select_method = "full", k = 4, models = "rf",
                  rf_ntree = 10, seed = 87, verbose = FALSE),
    "single response class"
  )
  expect_equal(cv$metrics$n_folds, 3L)
  expect_equal(cv$metrics$auc_n_folds, 2L)
  expect_equal(cv$metrics$tss_n_folds, 1L)
  expect_equal(cv$metrics$cbi_n_folds, 1L)
  expect_equal(cv$metrics$auc_mean, 0.65)
  expect_equal(cv$metrics$auc_sd, sd(c(0.7, 0.6)))
  expect_equal(cv$metrics$cbi_mean, -0.3)
  expect_true(is.na(cv$metrics$cbi_sd))
  expect_identical(cv$fold_status[4], "single_class")
  expect_true(all(is.na(cv$oof$HSS_rf[19:24])))

  evaluated <- 0L
  metric_rows <- lapply(metric_rows, function(x) {
    x["cbi"] <- NA_real_
    x
  })
  expect_warning(
    cv <- cast_cv(dat, select_method = "full", k = 4, models = "rf",
                  rf_ntree = 10, seed = 87, verbose = FALSE),
    "single response class"
  )
  expect_equal(cv$metrics$cbi_n_folds, 0L)
  expect_identical(cv$metrics$cbi_mean, NA_real_)
  expect_identical(cv$metrics$cbi_sd, NA_real_)
})
