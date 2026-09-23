# AUC direction regression (M1): reversed predictors must report AUC < 0.5 ---

test_that("evaluate_model_full does not mirror reversed predictions", {
  skip_if_not_installed("pROC")
  set.seed(1)
  obs <- rbinom(300, 1, 0.4)
  pred_bad <- ifelse(obs == 1, 0.05, 0.95)  # perfectly reversed ranking
  met <- evaluate_model_full(pred_bad, obs)
  expect_lt(unname(met["auc"]), 0.5)

  pred_good <- 1 - pred_bad
  met_good <- evaluate_model_full(pred_good, obs)
  expect_gt(unname(met_good["auc"]), 0.9)
  expect_gt(unname(met_good["tss"]), 0.9)
})

test_that("compute_auc (Wilcoxon) is directional too", {
  set.seed(2)
  y <- rbinom(200, 1, 0.5)
  pred <- ifelse(y == 1, 0.1, 0.9)
  expect_lt(compute_auc(y, pred), 0.5)
  expect_gt(compute_auc(y, 1 - pred), 0.5)
})

test_that("full evaluation shares one ROC between AUC and TSS", {
  skip_if_not_installed("pROC")
  obs <- rep(0:1, 20)
  pred <- seq(0.01, 0.99, length.out = length(obs))
  roc <- pROC::roc
  reference <- roc(obs, pred, quiet = TRUE, direction = "<")
  threshold <- pROC::coords(reference, "best", ret = c("sensitivity", "specificity"))
  calls <- 0L
  local_mocked_bindings(roc = function(...) {
    calls <<- calls + 1L
    roc(...)
  }, .package = "pROC")
  metrics <- evaluate_model_full(pred, obs)
  expect_identical(calls, 1L)
  expect_equal(unname(metrics["auc"]), as.numeric(pROC::auc(reference)))
  expect_equal(unname(metrics["tss"]),
               threshold$sensitivity[1] + threshold$specificity[1] - 1)
})

test_that("failed ROC evaluation retains missing metrics", {
  skip_if_not_installed("pROC")
  metrics <- evaluate_model_full(c(0.2, 0.4, 0.6), rep(1, 3))
  expect_named(metrics, c("auc", "tss", "cbi"))
  expect_true(all(is.na(metrics)))
})
