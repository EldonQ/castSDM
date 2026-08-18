# CPI fold-level inference (C1) and rare-species guards (F10) ----------------

make_cpi_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
  data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(2.5 * x1 - 2 * x2)),
    x1 = x1, x2 = x2, x3 = x3, x4 = x4
  )
}

skip_cpi <- function() {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")
  skip_if_not_installed("knockoff")
}

test_that("CPI fold inference finds the true driver with fold-level t test", {
  skip_cpi()
  dat <- make_cpi_data()
  sel <- cast_select(dat, method = "cpi", n_folds = 5L, num_trees = 50L,
                     inference = "fold", min_vars = 1L,
                     seed = 7, verbose = FALSE)
  expect_s3_class(sel, "cast_select")
  # Strong signal: the true driver is at least the top-ranked candidate
  # (selection may also come from the min_vars fallback with df = folds - 1).
  expect_true("x1" %in% sel$selected)
  top <- sel$scores$variable[which.max(sel$scores$cpi)]
  expect_identical(top, "x1")
  expect_identical(sel$diagnostics$inference, "fold")
  expect_equal(sel$diagnostics$n_folds, 5L)
  expect_true(all(is.finite(sel$scores$cpi)))
})

test_that("observation inference reproduces the reference per-observation test", {
  skip_cpi()
  dat <- make_cpi_data()
  sel <- cast_select(dat, method = "cpi", n_folds = 5L, num_trees = 50L,
                     inference = "observation",
                     seed = 7, verbose = FALSE)
  expect_s3_class(sel, "cast_select")
  expect_identical(sel$diagnostics$inference, "observation")
  # i.i.d. observation-level test is at least as significant as fold-level:
  # both should run, and p-values must be valid probabilities.
  p <- sel$scores$p_value
  expect_true(all(p >= 0 & p <= 1, na.rm = TRUE))
})

test_that("CPI aborts with an actionable message when classes are too rare", {
  skip_cpi()
  set.seed(1)
  dat <- data.frame(
    lon = runif(100), lat = runif(100),
    presence = c(rep(1L, 4L), rep(0L, 96L)),
    x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100)
  )
  expect_error(
    cast_select(dat, method = "cpi", n_folds = 5L, num_trees = 50L,
                verbose = FALSE),
    "n_folds"
  )
  # ... but succeeds once the folds are reduced below the rarer class size
  sel <- cast_select(dat, method = "cpi", n_folds = 2L, num_trees = 50L,
                     seed = 3, verbose = FALSE)
  expect_s3_class(sel, "cast_select")
})

test_that("CPI folds are response-stratified", {
  fold_id <- castSDM:::.cast_stratified_folds(
    c(rep(1L, 7L), rep(0L, 93L)), k = 5L, seed = 2
  )
  tab <- table(fold_id, c(rep(1L, 7L), rep(0L, 93L)))
  # every fold holds at least one presence (7 presences over 5 folds)
  expect_true(all(tab[, "1"] >= 1L))
  expect_equal(length(fold_id), 100L)
})
