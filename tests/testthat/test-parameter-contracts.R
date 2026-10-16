test_that("imputed predictions do not claim observed environmental support", {
  skip_if_not_installed("ranger")
  set.seed(1601)
  d <- data.frame(presence = rep(0:1, 40), x = rnorm(80), z = rnorm(80))
  fit <- cast_fit(d, models = "rf", rf_ntree = 20, seed = 1601, verbose = FALSE)
  new <- d[1:4, c("x", "z")]
  new$x[2] <- NA; new$z[3] <- Inf
  p <- cast_predict(fit, new)$predictions
  expect_identical(p$n_imputed, c(0L, 1L, 1L, 0L))
  expect_true(all(is.finite(p$HSS_rf)))
  expect_true(all(is.na(p$mess[2:3])))
  expect_true(all(is.na(p$extrapolating[2:3])))
  expect_true(all(is.finite(p$mess[c(1, 4)])))
})

test_that("explicit RF configuration survives tuning and refitting", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(1602)
  d <- data.frame(presence = rep(0:1, 60), x = rnorm(120), z = rnorm(120))
  fit <- cast_fit(d, models = "rf", rf_ntree = 20, rf_mtry = 1,
                  tune = TRUE, seed = 1602, verbose = FALSE)
  expect_equal(fit$models$rf$model$mtry, 1)
  expect_equal(fit$models$rf$tune$grid$mtry, 1)
  expect_equal(.cast_tuned_args(fit)$rf_mtry, 1)
  plain <- cast_fit(d, models = "rf", rf_ntree = 20, rf_mtry = 2,
                    seed = 1602, verbose = FALSE)
  expect_equal(.cast_tuned_args(plain)$rf_mtry, 2)
})

test_that("every spatial fold fits the requested engine configuration", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(1603)
  d <- data.frame(lon = rep(1:3, each = 40), lat = rep(seq_len(40), 3),
                  presence = rep(0:1, 60), x = rnorm(120), z = rnorm(120))
  observed <- list()
  real_fit <- cast_fit
  testthat::local_mocked_bindings(
    make_spatial_folds = function(...) rep(1:3, each = 40),
    cast_fit = function(...) {
      fit <- real_fit(...)
      observed[[length(observed) + 1L]] <<- fit$models$rf$model$mtry
      fit
    })
  cv <- cast_cv(d, select_method = "full", models = "rf", k = 3,
                rf_ntree = 20, rf_mtry = 1, seed = 1603, verbose = FALSE)
  expect_equal(unlist(observed), rep(1, 3))
  expect_equal(cv$pred_folds["rf"], c(rf = 3L))
})

test_that("integer parameters reject fractional or missing values early", {
  d <- data.frame(lon = 1:20, lat = 0, presence = rep(0:1, 10), x = 1:20)
  expect_error(cast_fit(d, models = character()), "at least one unique")
  expect_error(cast_fit(d, models = c("rf", "RF")), "at least one unique")
  expect_error(cast_fit(d, models = "rf", rf_mtry = 1.5), "rf_mtry")
  expect_error(cast_fit(d, models = "rf", rf_ntree = NA_real_), "rf_ntree")
  expect_error(cast_cv(d, k = 2.5), "k")
  expect_error(cast_cv(d, n_repeat = 1.5), "n_repeat")
  d$presence <- 1
  expect_error(cast_fit(d, models = "rf"), "both response classes")
  d$presence <- rep(0:1, 10); d$x <- Inf
  expect_error(cast_fit(d, screen = new_cast_select("x", data.frame()), models = "rf"), "no finite training")
})

test_that("effect coverage captions use evaluated rather than supported records", {
  skip_if_not_installed("ggplot2")
  tab <- data.frame(driver = c("x", "z"), mean_abs_dHSS = c(.1, .2),
                     mean_signed_dHSS = c(.1, -.2), n = 100,
                     n_supported = c(60, 80), support = c(.6, .8), masked = FALSE)
  class(tab) <- c("cast_effect_table", "data.frame")
  p <- plot(tab)
  expect_match(p$labels$subtitle, "60-80%")
  expect_match(p$labels$subtitle, "Evaluated records: 100")
})

test_that("external spatial folds remain aligned and are not regenerated", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(1605)
  d <- data.frame(lon = 1:120, lat = 0, presence = rep(0:1, 60),
                  x = rnorm(120), z = rnorm(120))
  ids <- rep(1:3, each = 40)
  local_mocked_bindings(make_spatial_folds = function(...) stop("must not regenerate"))
  cv <- cast_cv(d, select_method = "full", models = "rf", rf_ntree = 20,
    k = 3, fold_ids = ids, seed = 1605, verbose = FALSE)
  expect_identical(cv$folds, ids)
  expect_identical(cv$block_method, "external")
  expect_error(cast_cv(d, k = 3, fold_ids = ids[-1]), "every data row")
  expect_error(cast_cv(d, k = 3, fold_ids = ids, n_repeat = 2), "n_repeat = 1")
})
