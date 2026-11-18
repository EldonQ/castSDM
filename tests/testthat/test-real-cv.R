test_that("real spatial folds retain all OOF locations and report every candidate", {
  skip_if_not_installed("ranger"); skip_if_not_installed("pROC")
  x <- real_case(); at <- x$partition == "train"
  data <- x$data[at, c("lon", "lat", "presence", x$variables)]
  held <- x$partition == "holdout"
  response <- list(current = x$data[held, x$variables], future = x$future[held, ],
    drivers = x$variables, scenario = x$scenario, support_probs = c(.01, .99))
  cache <- tempfile("real-cv-folds-")
  result <- cast_cv(data, models = "rf", select_method = "vif", k = 5L,
    fold_ids = x$fold[at], coords_crs = x$crs, rf_ntree = x$fit_settings$rf_ntree,
    seed = x$fit_settings$seed, clamp = TRUE, response_data = response, num_threads = 1L, verbose = FALSE,
    cache_dir = cache)
  reused <- cast_cv(data, models = "rf", select_method = "vif", k = 5L,
    fold_ids = x$fold[at], coords_crs = x$crs, rf_ntree = x$fit_settings$rf_ntree,
    seed = x$fit_settings$seed, clamp = TRUE, response_data = response, num_threads = 1L, verbose = FALSE,
    cache_dir = cache)
  expect_true(all(reused$timings$cache_hit))
  expect_identical(reused$oof, result$oof)
  expect_identical(reused$response_validation, result$response_validation)
  expect_equal(nrow(result$oof), nrow(data))
  expect_equal(result$fold_info$fold, x$fold[at])
  expect_setequal(result$selection_freq$variable, x$variables)
  expect_true(result$clamp)
  expect_length(result$response_validation$runs, 5L)
  sensitivity <- castSDM:::.cast_response_validation_summary(result, c(rf = 1))
  expect_true(all(sensitivity$n_available == 5L))
  expect_true(all(sensitivity$n_common > 0L))
  expect_true(all(sensitivity$n_locations == sum(held)))
  expect_true(all(sensitivity$direction_unanimous_fraction >= 0 & sensitivity$direction_unanimous_fraction <= 1))
  # Refit the first actual spatial fold independently to check row order,
  # selection, training-only range bounds, and full future-variable changes.
  train <- data[x$fold[at] != 1L, ]
  screen <- cast_select(train, verbose = FALSE)
  fit <- cast_fit(train, screen = screen, models = "rf",
    rf_ntree = x$fit_settings$rf_ntree, seed = x$fit_settings$seed+101L,
    num_threads = 1L, verbose = FALSE)
  direct <- cast_attribution(fit, response$current, response$future,
    drivers = intersect(x$variables, screen$selected), weights = c(rf = 1),
    clamp = TRUE, support_probs = response$support_probs, verbose = FALSE)
  first <- result$response_validation$runs[[1L]]
  expect_equal(as.numeric(first$values[, "rf", "total"]), direct$data$total_change)
  expect_equal(first$supported, direct$data$support_total == 1)
  expect_true(all(first$observed_range[first$supported]))
  expect_setequal(sensitivity$criterion, c("marginal_training_quantiles", "observed_training_range"))
  expect_true(all(sensitivity$status == "descriptive_sensitivity"))
  for (v in intersect(x$variables, screen$selected))
    expect_equal(as.numeric(first$values[, "rf", paste0("only:", v)]), direct$data[[paste0("effect_", v)]])
  incomplete <- result; incomplete$response_validation$runs[1L] <- list(NULL)
  expect_equal(castSDM:::.cast_response_validation_summary(incomplete,c(rf=1))$status, "incomplete_spatial_fits")
  for (v in x$variables) expect_equal(result$selection_freq$freq[
    result$selection_freq$variable == v], mean(vapply(result$selections,
      function(selected) v %in% selected, logical(1))))
  expect_true(result$complete)
  expect_true(all(result$fold_balance$can_fit))
})

test_that("single-class actual test records retain proper losses without discrimination scores", {
  x <- real_case(); fit <- real_fit()
  at <- x$partition == "holdout" & x$data$presence == 0L
  data <- x$data[at, x$variables]
  pred <- castSDM:::predict_single_model(fit$models$rf, data)
  obs <- x$data$presence[at]
  value <- castSDM:::evaluate_model_full(pred, obs)
  clipped <- pmin(pmax(pred, 1e-7), 1 - 1e-7)
  expect_equal(unname(value["brier"]), mean((clipped - obs)^2))
  expect_equal(unname(value["logloss"]), -mean(log1p(-clipped)))
  expect_true(all(is.na(value[c("auc", "tss", "sedi", "kappa", "boyce")])))
  expect_true(is.na(value["tss_threshold"]))
})

test_that("unfit real-record designs are identified before model fitting", {
  x <- real_case(); at <- x$partition == "train" & x$data$presence == 0L
  data <- x$data[at, c("lon", "lat", "presence", x$variables)]
  plan <- castSDM:::.cast_cv_plan(data, 5L, "blockcv", NULL, x$crs,
    0, "deg", 1L, x$fit_settings$seed, fold_ids = x$fold[at])
  expect_true(all(!plan$balance$can_fit))
  expect_error(castSDM:::.cast_cv_check_plan(plan, verbose = FALSE), class = "cast_cv_invalid_design")
})


test_that("spatial hold-out reports class coverage and keeps physical blocks intact", {
  skip_if_not_installed("blockCV"); skip_if_not_installed("sf")
  x <- real_case()
  first <- cast_prepare(x$data, coords_crs=x$crs, block_method="blockcv", seed=42L, verbose=FALSE)
  balance <- first$fold_design$class_balance
  expect_equal(balance$total, as.numeric(table(x$data$presence)))
  expect_equal(balance$holdout, as.numeric(table(factor(first$test$presence,levels=c(0,1)))))
  expect_equal(balance$training+balance$holdout,balance$total)
  expect_equal(balance$achieved_holdout_fraction,balance$holdout/balance$total)
  ids <- first$fold_design$block_ids
  expect_length(intersect(unique(ids[first$train_index]),unique(ids[first$test_index])),0L)
  second <- cast_prepare(x$data,coords_crs=x$crs,block_method="blockcv",seed=42L,verbose=FALSE)
  expect_identical(first$test_index,second$test_index)
})
