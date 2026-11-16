test_that("real spatial folds retain all OOF locations and report every candidate", {
  skip_if_not_installed("ranger"); skip_if_not_installed("pROC")
  x <- real_case(); at <- x$partition == "train"
  data <- x$data[at, c("lon", "lat", "presence", x$variables)]
  result <- cast_cv(data, models = "rf", select_method = "vif", k = 5L,
    fold_ids = x$fold[at], coords_crs = x$crs, rf_ntree = x$fit_settings$rf_ntree,
    seed = x$fit_settings$seed, clamp = TRUE, num_threads = 1L, verbose = FALSE)
  expect_equal(nrow(result$oof), nrow(data))
  expect_equal(result$fold_info$fold, x$fold[at])
  expect_setequal(result$selection_freq$variable, x$variables)
  expect_true(result$clamp)
  for (v in x$variables) expect_equal(result$selection_freq$freq[
    result$selection_freq$variable == v], mean(vapply(result$selections,
      function(selected) v %in% selected, logical(1))))
})
