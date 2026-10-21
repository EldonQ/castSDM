# Elastic-net screen: collinearity filter, ranking, force, coverage ----------

skip_screen_deps <- function() {
  skip_if_not_installed("glmnet")
}

feature_null_data <- function(n = 300, r = 0.98, seed = 21) {
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

test_that("screen drops a near-duplicate and keeps the true drivers", {
  skip_screen_deps()
  dat <- feature_null_data()
  scr <- cast_select(dat, ncov = 2L, maxncov = 2L, seed = 22, verbose = FALSE)
  expect_s3_class(scr, "cast_select")
  expect_identical(scr$method, "screen")
  expect_equal(sum(c("x1", "x3") %in% scr$selected), 1L)
  expect_true("x2" %in% scr$selected)
  expect_false(any(c("x4", "x5") %in% scr$selected))
  expect_true(all(scr$selected %in% scr$diagnostics$filtered))
  expect_true(all(scr$scores$selected_reason[scr$scores$selected] %in%
                    c("screen", "prespecified", "ensured")))
})

test_that("screen is deterministic under a seed", {
  skip_screen_deps()
  dat <- feature_null_data(n = 200)
  a <- cast_select(dat, seed = 99, verbose = FALSE)
  b <- cast_select(dat, seed = 99, verbose = FALSE)
  expect_identical(a$selected, b$selected)
})

test_that("var_id keeps one covariate per construct", {
  skip_screen_deps()
  set.seed(23)
  n <- 300
  base <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.3 * base)),
    scale100 = base, scale500 = 0.9 * base + 0.1 * rnorm(n),
    scale1000 = 0.85 * base + 0.2 * rnorm(n), noise = rnorm(n)
  )
  ids <- c(scale100 = "construct", scale500 = "construct",
           scale1000 = "construct", noise = "noise")
  scr <- cast_select(dat, var_id = ids, seed = 24, verbose = FALSE)
  expect_equal(sum(c("scale100", "scale500", "scale1000") %in% scr$selected), 1L)
})

test_that("ncov truncates the elastic-net screen", {
  skip_screen_deps()
  dat <- feature_null_data()
  scr2 <- cast_select(dat, ncov = 2L, seed = 25, verbose = FALSE)
  expect_identical(scr2$diagnostics$ncov, 2L)
  expect_lte(length(scr2$selected), 2L)
  n_pres <- sum(dat$presence)
  scr_default <- cast_select(dat, seed = 25, verbose = FALSE)
  expect_identical(scr_default$diagnostics$ncov,
                   as.integer(min(12L, ceiling(log2(n_pres)))))
})

test_that("keep forces prespecified predictors through the screen", {
  skip_screen_deps()
  dat <- feature_null_data(n = 200)
  scr <- suppressWarnings(cast_select(dat, keep = "x5", seed = 26, verbose = FALSE))
  expect_true("x5" %in% scr$selected)
  row <- scr$scores[scr$scores$variable == "x5", ]
  expect_true(row$kept_by_design)
  expect_identical(row$selected_reason, "prespecified")
  reported <- cast_importance(scr)$effects
  expect_true(all(reported$kept_by_design[reported$variable == "x5"]))
})

test_that("ensure tops up a declared coverage class", {
  skip_screen_deps()
  dat <- feature_null_data()
  flags <- c(x1 = FALSE, x2 = FALSE, x3 = FALSE, x4 = TRUE, x5 = FALSE)
  scr <- cast_select(dat, ensure = list(climate = flags),
                     ncov = 2L, maxncov = 2L, seed = 27, verbose = FALSE)
  expect_true("x4" %in% scr$selected)
  expect_identical(scr$scores$selected_reason[scr$scores$variable == "x4"],
                   "ensured")
  expect_true("x4" %in% scr$diagnostics$ensured)
})

test_that("scores and diagnostics carry the screen schema", {
  skip_screen_deps()
  dat <- feature_null_data(n = 250)
  s <- cast_select(dat, keep = "x1", seed = 28, verbose = FALSE)
  expect_true(all(c("step_added", "loss_gain", "selected_reason",
                    "kept_by_design", "selected") %in% names(s$scores)))
  d <- s$diagnostics
  expect_true(all(c("keep", "corcut", "ncov", "maxncov", "ensured",
                    "filtered", "lambda", "coefficients", "n_pres") %in%
                    names(d)))
  expect_identical(d$keep, "x1")
  expect_true(length(d$filtered) >= 1L)
  expect_true(all(s$selected %in% union(d$filtered, d$keep)))
  expect_s3_class(cast_importance(s), "cast_importance")
  expect_true(all(c("step_added", "loss_gain") %in%
                    names(cast_importance(s)$effects)))
})

test_that("invalid selection arguments fail explicitly", {
  skip_screen_deps()
  dat <- feature_null_data(n = 150)
  expect_error(cast_select(dat, method = "typo", verbose = FALSE), "screen")
  expect_error(cast_select(dat, corcut = 2, verbose = FALSE), "corcut")
  expect_error(cast_select(dat, ncov = 0, verbose = FALSE), "ncov")
  expect_error(cast_select(dat, maxncov = 0, verbose = FALSE), "maxncov")
  expect_error(cast_select(dat, seed = "x", verbose = FALSE), "seed")
})

test_that("screen plots as a ggplot", {
  skip_screen_deps()
  skip_if_not_installed("ggplot2")
  dat <- feature_null_data(n = 250)
  s <- cast_select(dat, keep = "x1", seed = 30, verbose = FALSE)
  expect_s3_class(plot(cast_importance(s)), "ggplot")
  expect_s3_class(plot(s), "ggplot")
})

test_that("a null signal still returns a valid, documented screen", {
  skip_screen_deps()
  set.seed(31)
  n <- 300
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, 0.3),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  s <- suppressWarnings(cast_select(dat, seed = 32, verbose = FALSE))
  expect_s3_class(s, "cast_select")
  expect_true(all(s$selected %in% s$diagnostics$filtered))
  expect_true(all(s$scores$selected_reason %in%
                    c("screen", "prespecified", "ensured", "excluded")))
  expect_false(any(c("fallback", "p_value", "null_threshold",
                     "interventional_effect") %in% names(s$scores)))
})
