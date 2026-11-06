test_that("training stages resume across output directories and invalidate scientific changes", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(835)
  d <- data.frame(lon = runif(120, 70, 130), lat = runif(120, 20, 50),
    x = rnorm(120), z = rnorm(120))
  d$presence <- rbinom(120, 1, plogis(d$x))
  p <- cast_protocol(models = "rf", select_method = "full", rf_ntree = 15,
    cv_k = 2, block_method = "grid", aoa = FALSE, seed = 8)
  root <- tempfile("training-cache-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  run <- function(where, protocol = p) cast_run(d, protocol = protocol,
    output_dir = file.path(root, where), cache_dir = file.path(root, "shared"), verbose = FALSE)
  first <- run("a"); second <- run("b")
  hits <- function(x) vapply(x$receipts, function(r) isTRUE(r$cache_hit), logical(1))
  expect_true(all(hits(second)))
  expect_equal(second$cv$oof, first$cv$oof)
  p$rf_ntree <- 20
  third <- run("c", p)
  expect_equal(hits(third)[1:2], c(TRUE, TRUE))
  expect_false(any(hits(third)[3:6]))
  split <- cast_prepare(d, block_method = "grid", seed = 8, verbose = FALSE)
  expect_identical(.cast_verify_prepared(d, split), split)
  split$train$x[1] <- split$train$x[1] + 1
  expect_error(.cast_verify_prepared(d, split), "do not match")
})
