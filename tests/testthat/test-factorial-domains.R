test_that("all factorial counts use the same finite domain and mask", {
  skip_if_not_installed("terra")
  a <- terra::rast(matrix(.2, 2, 2))
  mask <- terra::rast(matrix(c(1, 0, 0, 0), 2, 2))
  out <- cast_factorial_attribution(a, a+.5, a+.4, a+.1,
    threshold = .5, mask = mask, verbose = FALSE)
  expect_true(all(out$summary$n_valid == 1))
  expect_true(all(out$binary$n_valid == 1))
  expect_equal(out$binary$estimate[3], 1)
  c <- a+.4; terra::values(c)[1] <- NA
  out <- cast_factorial_attribution(a, a+.5, c, a+.1, threshold = .5, verbose = FALSE)
  expect_equal(out$binary$estimate[3:4], c(3, 0))
})

test_that("factorial tables reject invalid binary values and exclude Inf", {
  expect_error(cast_factorial_attribution(rep(.2, 4), rep(.8, 4), rep(.4, 4), rep(.5, 4),
    current_binary = rep(2, 4), joint_binary = rep(1, 4), climate_binary = rep(0, 4),
    lulc_binary = rep(0, 4), verbose = FALSE), "only 0, 1 or NA")
  a <- c(.2, Inf, NA, .2)
  out <- cast_factorial_attribution(a, rep(.8, 4), rep(.4, 4), rep(.5, 4), threshold = .5,
    support = c(1, 1, 1, 0), block_size = 1, verbose = FALSE)
  expect_true(all(out$summary$n_valid == 2))
  expect_true(all(out$binary$n_valid == 2))
  expect_true(all(out$binary_supported$n_valid == 1))
  expect_equal(out$binary$estimate[3], 2)
  expect_equal(out$binary_supported$estimate[3], 1)
})

test_that("AOA training neighbours exclude self exactly once", {
  skip_if_not_installed("FNN")
  expect_equal(.cast_aoa_di_train(matrix(c(0, 1, 10), ncol = 1)), c(1, 1, 9))
  expect_equal(.cast_aoa_di_train(matrix(c(0, 1), ncol = 1)), c(1, 1))
  expect_true(is.na(.cast_aoa_di_train(matrix(0, ncol = 1))))
})
