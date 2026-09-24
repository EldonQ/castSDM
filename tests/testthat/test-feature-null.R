# Stage 2 forward selection: admission by inner-CV loss improvement ----------

test_that("forward search recovers true parents and drops a redundant proxy", {
  skip_if_not_installed("ranger")
  set.seed(21)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- 0.98 * x1 + sqrt(1 - 0.98^2) * rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - 1.2 * x2)),
    x1 = x1, x2 = x2, x3 = x3, x4 = rnorm(n)
  )
  s <- cast_select(dat, num_trees = 150, seed = 22, verbose = FALSE)
  expect_s3_class(s, "cast_select")
  expect_identical(s$diagnostics$status, "selected")
  # The true driver is retained ...
  expect_true("x2" %in% s$selected)
  # ... exactly one of the near-duplicate pair survives (thinning or search) ...
  expect_equal(sum(c("x1", "x3") %in% s$selected), 1L)
  # ... and the pure noise predictor is not admitted.
  expect_false("x4" %in% s$selected)
  # Reasons use the forward vocabulary only.
  expect_true(all(s$scores$selected_reason[s$scores$selected] %in%
                    c("forward", "prespecified")))
  expect_true(nrow(s$diagnostics$path) >= 1L)
})

test_that("a null signal returns an empty set, not a fallback set", {
  skip_if_not_installed("ranger")
  set.seed(23)
  n <- 300
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, 0.3),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  expect_warning(
    s <- cast_select(dat, num_trees = 60, seed = 24, verbose = FALSE),
    "empty set"
  )
  expect_length(s$selected, 0L)
  expect_identical(s$diagnostics$status, "empty_selection")
  expect_true(all(s$scores$selected_reason == "excluded"))
  expect_false(any(c("fallback", "p_value", "null_threshold",
                     "interventional_effect") %in% names(s$scores)))
})

test_that("keep forces prespecified predictors at step 0", {
  skip_if_not_installed("ranger")
  set.seed(25)
  n <- 200
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.2 * rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x5 = rnorm(n)
  )
  s <- suppressWarnings(cast_select(dat, keep = "x5", num_trees = 40,
                                    seed = 26, verbose = FALSE))
  expect_true("x5" %in% s$selected)
  row <- s$scores[s$scores$variable == "x5", ]
  expect_identical(row$step_added, 0L)
  expect_identical(row$selected_reason, "prespecified")
  expect_true(row$kept_by_design)
  reported <- cast_importance(s)$effects
  expect_true(all(reported$kept_by_design[reported$variable == "x5"]))
})

test_that("scores and diagnostics carry the forward path schema", {
  skip_if_not_installed("ranger")
  set.seed(27)
  n <- 250
  x1 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.4 * x1)),
    x1 = x1, x2 = rnorm(n), x3 = rnorm(n)
  )
  # keep forces a non-empty path so the schema is exercised deterministically.
  s <- cast_select(dat, keep = "x1", num_trees = 60, seed = 28,
                   verbose = FALSE)
  expect_true(all(c("step_added", "loss_gain", "selected_reason",
                    "kept_by_design", "selected") %in% names(s$scores)))
  d <- s$diagnostics
  expect_true(all(c("path", "metric", "inner_method", "tolerance",
                    "n_folds", "null_loss", "final_loss", "n_fits",
                    "status") %in% names(d)))
  expect_true(d$n_fits >= 1L)
  expect_true(d$status %in% c("selected", "empty_selection"))
  expect_s3_class(cast_importance(s), "cast_importance")
  expect_true(all(c("step_added", "loss_gain") %in%
                    names(cast_importance(s)$effects)))
})

test_that("metric = 'auc' runs the same forward search on rank loss", {
  skip_if_not_installed("ranger")
  set.seed(29)
  n <- 250
  x1 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.4 * x1)),
    x1 = x1, x2 = rnorm(n), x3 = rnorm(n)
  )
  s <- cast_select(dat, metric = "auc", num_trees = 60, seed = 30,
                   verbose = FALSE)
  expect_identical(s$diagnostics$metric, "auc")
  expect_true("x1" %in% s$selected)
})

test_that("invalid stage-2 arguments fail explicitly", {
  set.seed(31)
  dat <- data.frame(presence = rep(0:1, 60),
                    x1 = rnorm(120), x2 = rnorm(120), x3 = rnorm(120))
  expect_error(cast_select(dat, tolerance = -1, verbose = FALSE), "tolerance")
  expect_error(cast_select(dat, n_folds = 1L, verbose = FALSE), "n_folds")
  expect_error(cast_select(dat, metric = "typo", verbose = FALSE), "metric")
  expect_error(cast_select(dat, num_trees = 0L, verbose = FALSE), "num_trees")
})

test_that("a huge tolerance stops the search with an empty set", {
  skip_if_not_installed("ranger")
  set.seed(32)
  n <- 250
  x1 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.4 * x1)),
    x1 = x1, x2 = rnorm(n), x3 = rnorm(n)
  )
  # No count cap exists: stopping is driven by the tolerance alone.
  expect_warning(
    s <- cast_select(dat, tolerance = 10, num_trees = 60, seed = 33,
                     verbose = FALSE),
    "empty set"
  )
  expect_length(s$selected, 0L)
})

test_that("forward path plots as a ggplot", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("ggplot2")
  set.seed(34)
  n <- 250
  x1 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.4 * x1)),
    x1 = x1, x2 = rnorm(n), x3 = rnorm(n)
  )
  s <- cast_select(dat, num_trees = 60, seed = 35, verbose = FALSE)
  expect_s3_class(plot(cast_importance(s)), "ggplot")
  expect_s3_class(plot(s), "ggplot")
})
