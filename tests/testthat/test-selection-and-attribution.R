# Two-stage selection, effect table, and the necessity audit ----------------

# x3 is a near-copy of x1 with no effect of its own: stage 1 must drop it,
# and the necessity audit must find x1 substitutable when x3 is kept.
make_collinear_data <- function(n = 300, r = 0.98, seed = 21) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- r * x1 + sqrt(1 - r^2) * rnorm(n)
  data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - 1.2 * x2)),
    x1 = x1, x2 = x2, x3 = x3, x4 = rnorm(n), x5 = rnorm(n)
  )
}

test_that("stage 1 drops a near-duplicate predictor", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data()
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 9,
                     seed = 22, verbose = FALSE)
  expect_s3_class(scr, "cast_select")
  expect_identical(scr$method, "two_stage")
  # One of the two collinear partners must be thinned, never both kept.
  expect_true(any(scr$scores$collinear_thinned[scr$scores$variable %in% c("x1", "x3")]))
  expect_false(all(c("x1", "x3") %in% scr$diagnostics$stage1_kept))
})

test_that("stage 2 thresholds against the permutation null, not zero", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data()
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 19,
                     seed = 23, verbose = FALSE)
  expect_true(is.finite(scr$diagnostics$null_threshold))
  expect_equal(scr$diagnostics$null_quantile, 0.95)
  expect_equal(scr$diagnostics$n_perm, 19L)
  # Pure noise predictors have importance above zero but not above the null.
  noise <- scr$scores[scr$scores$variable %in% c("x4", "x5"), ]
  expect_true(all(noise$p_value > 0.05 | !noise$selected))
  # The retained set must be a subset of the stage-1 survivors.
  expect_true(all(scr$selected %in% scr$diagnostics$stage1_kept))
})

test_that("method = 'full' keeps every predictor and skips both stages", {
  dat <- make_collinear_data(n = 120)
  scr <- cast_select(dat, method = "full", verbose = FALSE)
  expect_setequal(scr$selected, c("x1", "x2", "x3", "x4", "x5"))
  expect_identical(scr$method, "full")
})

test_that("retired selection arguments warn instead of silently changing the screen", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 150)
  expect_warning(
    cast_select(dat, method = "two_stage", num_trees = 40, n_perm = 5,
                alpha = 0.01, min_vars = 3L, seed = 24, verbose = FALSE),
    "Deprecated"
  )
})

test_that("cast_effect_table averages a symmetric shift set and splits magnitude from sign", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 200)
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 60,
                  seed = 25, verbose = FALSE)
  eff <- cast_effect_table(fit, newdata = dat, verbose = FALSE)
  expect_s3_class(eff, "cast_effect_table")
  expect_setequal(eff$driver, c("x1", "x2"))
  expect_true(all(c("mean_abs_dHSS", "mean_signed_dHSS") %in% names(eff)))
  # The frozen shift set is +/- 1, 2 SD: four interventions per driver.
  expect_true(all(eff$n_shifts == 4L))
  # Magnitude is non-negative and never smaller than |direction|.
  expect_true(all(eff$mean_abs_dHSS >= 0))
  expect_true(all(eff$mean_abs_dHSS >= abs(eff$mean_signed_dHSS)))
  # x1 raises suitability, x2 lowers it (the data-generating signs).
  expect_gt(eff$mean_signed_dHSS[eff$driver == "x1"], 0)
  expect_lt(eff$mean_signed_dHSS[eff$driver == "x2"], 0)
})

test_that("cast_effect_table rejects a zero shift", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 120)
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 26, verbose = FALSE)
  expect_error(
    cast_effect_table(fit, newdata = dat, shifts = c(-1, 0, 1), verbose = FALSE),
    "non-zero"
  )
})

test_that("cast_necessity scores an out-of-sample knockout cost per driver", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 300)
  nec <- cast_necessity(dat, variables = c("x1", "x2", "x4"), k = 3,
                        num_trees = 60, seed = 27, verbose = FALSE)
  expect_s3_class(nec, "cast_necessity")
  expect_setequal(nec$necessity$variable, c("x1", "x2", "x4"))
  expect_true(all(c("mean_dAUC", "pct_folds_positive", "necessary") %in%
                    names(nec$necessity)))
  expect_identical(nec$k, 3L)
  expect_identical(dim(nec$fold_dauc), c(3L, 3L))
  # The true drivers must cost more to drop than the pure noise predictor.
  d <- stats::setNames(nec$necessity$mean_dAUC, nec$necessity$variable)
  expect_gt(d[["x1"]], d[["x4"]])
})

test_that("cast_necessity reads the predictor set off a screen and needs at least two", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 200)
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  nec <- cast_necessity(dat, screen = screen, k = 2, num_trees = 40,
                        seed = 28, verbose = FALSE)
  expect_setequal(nec$necessity$variable, c("x1", "x2"))
  expect_error(
    cast_necessity(dat, variables = "x1", k = 2, verbose = FALSE),
    "at least two"
  )
})
