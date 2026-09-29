# P1 review round, W3: prevalence-correcting case weights. Weight = target/n
# per class, so the weighted prevalence is exactly the target. Only the
# fitting objective is re-prevalenced; evaluation metrics stay unweighted.
# MaxEnt has no row-weight interface and is always fitted unweighted.

.w3_test_data <- function(seed = 41, n_pres = 20, n_bg = 300) {
  set.seed(seed)
  data.frame(
    presence = c(rep(1L, n_pres), rep(0L, n_bg)),
    x1 = c(runif(n_pres, 0.6, 1), runif(n_bg, 0, 1)),
    x2 = c(runif(n_pres, 0.6, 1), runif(n_bg, 0, 1)),
    x3 = runif(n_pres + n_bg)
  )
}

test_that("invalid prevalence_target aborts informatively", {
  dat <- .w3_test_data()
  expect_error(cast_fit(dat, models = "rf", prevalence_target = 0),
               "in \\(0, 1\\)")
  expect_error(cast_fit(dat, models = "rf", prevalence_target = 1),
               "in \\(0, 1\\)")
  expect_error(cast_fit(dat, models = "rf", prevalence_target = 1.5),
               "in \\(0, 1\\)")
  expect_error(cast_fit(dat, models = "rf", prevalence_target = "half"),
               "in \\(0, 1\\)")
  # A single-class response cannot define a prevalence.
  one_class <- .w3_test_data()
  one_class$presence <- rep(1L, nrow(one_class))
  expect_error(cast_fit(one_class, models = "rf", prevalence_target = 0.5),
               "both response classes")
})

test_that("prevalence_target shifts RF predictions toward the target prevalence", {
  skip_if_not_installed("ranger")
  # Strong signal, heavily imbalanced sample (20:300). Weighting to 0.5 gives
  # presences a ~15x relative weight, so the weighted model must score the
  # high-suitability corner noticeably higher than the unweighted one.
  dat <- .w3_test_data()
  fit_uw <- cast_fit(dat, models = "rf", seed = 7, verbose = FALSE)
  fit_w <- cast_fit(dat, models = "rf", prevalence_target = 0.5,
                    seed = 7, verbose = FALSE)
  expect_equal(fit_w$scaling$prevalence_target, 0.5)
  expect_null(fit_uw$scaling$prevalence_target)

  hi <- data.frame(x1 = c(0.9, 0.85, 0.95), x2 = c(0.9, 0.85, 0.95),
                   x3 = c(0.5, 0.5, 0.5))
  p_uw <- predict_single_model(fit_uw$models$rf, hi)
  p_w <- predict_single_model(fit_w$models$rf, hi)
  expect_gt(mean(p_w), mean(p_uw))
})

test_that("MaxEnt stays unweighted regardless of prevalence_target", {
  skip_if_not_installed("maxnet")
  dat <- .w3_test_data(42)
  fit_uw <- cast_fit(dat, models = "maxent", seed = 3, verbose = FALSE)
  fit_w <- cast_fit(dat, models = "maxent", prevalence_target = 0.5,
                    seed = 3, verbose = FALSE)
  p_uw <- predict_single_model(fit_uw$models$maxent, dat[, c("x1", "x2", "x3")])
  p_w <- predict_single_model(fit_w$models$maxent, dat[, c("x1", "x2", "x3")])
  expect_identical(p_uw, p_w)
})

test_that("BRT and GAM accept the weights", {
  skip_if_not_installed("gbm")
  skip_if_not_installed("mgcv")
  dat <- .w3_test_data(43)
  fit_brt <- cast_fit(dat, models = "brt", prevalence_target = 0.3,
                      seed = 5, verbose = FALSE)
  expect_s3_class(fit_brt$models$brt$model, "gbm")
  fit_gam <- cast_fit(dat, models = "gam", prevalence_target = 0.3,
                      seed = 5, verbose = FALSE)
  expect_s3_class(fit_gam$models$gam$model, "gam")
})

test_that("tuned RF with prevalence_target runs the grid (weight-aware scoring)", {
  skip_if_not_installed("ranger")
  dat <- .w3_test_data(44)
  fit <- cast_fit(dat, models = "rf", tune = TRUE,
                  prevalence_target = 0.5, seed = 9, verbose = FALSE)
  t <- fit$models$rf$tune
  expect_false(is.null(t))
  expect_true(is.finite(t$best_score))
  expect_named(t$best, "mtry")
})

test_that("cast_cv passes prevalence_target into every fold", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(45)
  n <- 320
  dat <- data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = c(rep(1L, 20), rep(0L, n - 20)),
    x1 = c(runif(20, 0.6, 1), runif(n - 20, 0, 1)),
    x2 = c(runif(20, 0.6, 1), runif(n - 20, 0, 1)),
    x3 = runif(n)
  )
  cv <- cast_cv(dat, models = "rf", k = 3L, select_method = "full",
                prevalence_target = 0.5, seed = 21, verbose = FALSE)
  expect_true("evaluated" %in% cv$fold_status)
  expect_true(nrow(cv$fold_metrics) >= 1L)
  expect_true(all(is.finite(cv$fold_metrics$tss)))
})
