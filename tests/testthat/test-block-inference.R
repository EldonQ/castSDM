# Spatial-block CPI inference (option a: cluster-robust t on spatial blocks) --

make_block_data <- function(n = 400, seed = 42) {
  set.seed(seed)
  data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(2.5 * rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n)
  )
}

skip_cpi <- function() {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")
  skip_if_not_installed("knockoff")
}

test_that("block inference runs on coordinate data and reports block df", {
  skip_cpi()
  dat <- make_block_data()
  sel <- cast_select(dat, response = "presence", method = "cpi",
                     n_folds = 5L, num_trees = 50L,
                     inference = "block", n_blocks = 12L, min_vars = 1L,
                     seed = 7, verbose = FALSE, num_threads = 1L)
  expect_s3_class(sel, "cast_select")
  expect_identical(sel$diagnostics$inference, "block")
  expect_true(sel$diagnostics$n_blocks >= 2L)
  expect_true(all(is.finite(sel$scores$p_value)))
  expect_true(all(sel$scores$p_value >= 0 & sel$scores$p_value <= 1))
})

test_that("block inference is the default when lon/lat are present", {
  skip_cpi()
  dat <- make_block_data()
  sel <- cast_select(dat, response = "presence", method = "cpi",
                     n_folds = 5L, num_trees = 50L,
                     min_vars = 1L, seed = 7, verbose = FALSE,
                     num_threads = 1L)
  expect_identical(sel$diagnostics$inference, "block")
})

test_that("block inference falls back to fold when coordinates are absent", {
  skip_cpi()
  set.seed(1)
  n <- 200
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  dat$presence <- rbinom(n, 1, plogis(2 * dat$x1))
  expect_warning(
    sel <- cast_select(dat, response = "presence", method = "cpi",
                       n_folds = 5L, num_trees = 40L,
                       inference = "block", min_vars = 1L, seed = 3,
                       verbose = FALSE, num_threads = 1L),
    "falling back"
  )
  expect_identical(sel$diagnostics$inference, "fold")
})
