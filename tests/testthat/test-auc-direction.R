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
