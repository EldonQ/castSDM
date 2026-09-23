# Two-stage selection and the shift effect table ------------------------------

# x3 is a near-copy of x1 with no effect of its own: stage 1 must drop it.
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

test_that("stage 2 thresholds against the conditional null, not zero", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data()
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 19,
                     seed = 23, verbose = FALSE)
  expect_true(all(is.finite(scr$diagnostics$null_threshold)))
  expect_equal(scr$diagnostics$null_quantile, 0.95)
  expect_equal(scr$diagnostics$n_perm, 19L)
  expect_match(scr$diagnostics$null_method, "within-stratum")
  # Pure noise predictors have effects above zero but not above the null.
  noise <- scr$scores[scr$scores$variable %in% c("x4", "x5"), ]
  expect_true(all(noise$p_value > 0.05 | !noise$selected))
  # The retained set must be a subset of the stage-1 survivors.
  expect_true(all(scr$selected %in% scr$diagnostics$stage1_kept))
  # Single-statistic screen: no second attribution column.
  expect_false("perm_importance" %in% names(scr$scores))
})

test_that("stage 1 ranks a U-shaped driver above noise (poly2 signal)", {
  skip_if_not_installed("ranger")
  set.seed(31)
  n <- 400
  xu <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(-1 + 2.2 * xu^2)),
    xu = xu, noise = rnorm(n), x3 = rnorm(n), x4 = rnorm(n), x5 = rnorm(n)
  )
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 9,
                     seed = 32, verbose = FALSE)
  rnk <- stats::setNames(scr$scores$stage1_rank, scr$scores$variable)
  expect_lt(rnk[["xu"]], rnk[["noise"]])
})

test_that("ncov caps the retained set and records the reason", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 400)
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 9,
                     ncov = 1L, seed = 33, verbose = FALSE)
  expect_lte(length(scr$selected), 1L)
  expect_true(scr$diagnostics$ncov == 1L)
  expect_true(all(scr$scores$selected_reason[scr$scores$selected] %in%
                    c("null+top-ncov", "fallback-top-ncov")))
})

test_that("method = 'full' keeps every predictor and skips both stages", {
  dat <- make_collinear_data(n = 120)
  scr <- cast_select(dat, method = "full", verbose = FALSE)
  expect_setequal(scr$selected, c("x1", "x2", "x3", "x4", "x5"))
  expect_identical(scr$method, "full")
})

test_that("retired or invalid selection requests fail without running a replacement", {
  dat <- make_collinear_data(n = 150)
  expect_error(cast_select(dat, alpha = 0.01, min_vars = 3L), "unused argument")
  expect_error(cast_select(dat, environment = "x1"), "unused argument")
  for (method in list("cpi", "dml", "rf", "tramicp", "typo", "two", NULL, character(),
                      NA_character_, 1, c("full", "two_stage"))) {
    expect_error(cast_select(dat, method = method), "method.*must be one of")
  }
})

test_that("cast_effect_table reports one raw shift with magnitude and sign", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 200)
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 60,
                  seed = 25, verbose = FALSE)
  eff <- cast_effect_table(fit, newdata = dat, shift = list(x1 = 1, x2 = 1),
                           verbose = FALSE)
  expect_s3_class(eff, "cast_effect_table")
  expect_setequal(eff$driver, c("x1", "x2"))
  expect_true(all(c("shift_raw", "mean_abs_dHSS", "mean_signed_dHSS",
                    "n", "n_supported", "support", "masked") %in% names(eff)))
  expect_equal(eff$shift_raw, c(1, 1))
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
    cast_effect_table(fit, newdata = dat, shift = 0, verbose = FALSE),
    "non-zero"
  )
})

test_that("removed knockout and scenario products forward to the shift products", {
  expect_error(cast_necessity(data.frame()), "removed in 0.12.0")
  expect_error(cast_sensitivity(), "removed in 0.12.0")
})

test_that("prespecified predictors survive thinning, null threshold and cap", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 150)
  required <- c("x1", "x3")
  expect_warning(scr <- cast_select(dat, keep = required, ncov = 2L,
    num_trees = 20L, n_perm = 1L, seed = 41, verbose = FALSE), "No predictor exceeded")
  expect_identical(scr$selected, required)
  rows <- scr$scores$variable %in% required
  expect_true(all(scr$scores$kept_by_design[rows]))
  expect_false(any(scr$scores$collinear_thinned[rows]))
  expect_false(any(scr$scores$passed_null[rows]))
  expect_false(any(scr$scores$fallback[rows]))
  expect_true(all(scr$scores$selected_reason[rows] == "prespecified"))
  expect_identical(scr$diagnostics$keep, required)
  expect_true(all(required %in% scr$diagnostics$stage1_kept))
  reported <- cast_importance(scr)$effects
  expect_true(all(reported$kept_by_design[reported$variable %in% required]))
  expect_true(all(reported$selected_reason[reported$variable %in% required] == "prespecified"))
  expect_false("perm_importance" %in% names(reported))
})

test_that("prespecified retention leaves only the remaining cap for optional predictors", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 150)
  scr <- suppressWarnings(cast_select(dat, keep = "x5", ncov = 2L,
    num_trees = 20L, n_perm = 1L, seed = 42, verbose = FALSE))
  expect_length(scr$selected, 2L)
  expect_true("x5" %in% scr$selected)
  expect_equal(sum(scr$scores$fallback), 1L)
  full <- cast_select(dat, method = "full", keep = "x5", verbose = FALSE)
  expect_equal(full$scores$variable[full$scores$kept_by_design], "x5")
  expect_setequal(full$selected, get_env_vars(dat))
})

test_that("invalid or unaccommodated prespecified sets fail explicitly", {
  dat <- make_collinear_data(n = 120)
  for (required in list("unknown", "lon", "presence", c("x1", "x1"), NA_character_, 1)) {
    expect_error(cast_select(dat, keep = required, verbose = FALSE), "keep")
  }
  expect_error(cast_select(dat, keep = c("x1", "x3"), ncov = 1L,
    verbose = FALSE), "increase the cap")
  dat$x5 <- 0
  expect_error(cast_select(dat, keep = "x5", verbose = FALSE), "vary in the training data")
})

test_that("nested CV retains the same prespecified variables within each training fold", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 200)
  required <- c("x1", "x3")
  cv <- suppressWarnings(cast_cv(dat, k = 2L, models = "rf", rf_ntree = 20L,
    select_args = list(keep = required, ncov = 2L, num_trees = 20L, n_perm = 1L),
    seed = 43, verbose = FALSE))
  expect_length(cv$screens, 2L)
  expect_true(all(vapply(cv$screens, function(s) {
    !is.null(s) && identical(s$diagnostics$keep, required) && all(required %in% s$selected)
  }, logical(1))))
})

test_that("the high-level pipeline forwards prespecified retention to CV", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 160)
  required <- c("x1", "x3")
  result <- suppressWarnings(cast(dat, models = "rf", do_predict = FALSE,
    do_cv = TRUE, cv_k = 2L, select_num_trees = 20L, select_n_perm = 1L,
    select_ncov = 2L, select_keep = required, seed = 44, verbose = FALSE))
  expect_identical(result$screen$selected, required)
  expect_false(is.null(result$cv))
  expect_true(all(vapply(result$cv$screens, function(s) {
    !is.null(s) && identical(s$diagnostics$keep, required) && all(required %in% s$selected)
  }, logical(1))))
})
