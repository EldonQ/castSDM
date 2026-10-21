# Elastic-net screen and the shift effect table -------------------------------

skip_screen_deps <- function() {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mgcv")
}

# x3 is a near-copy of x1 with no effect of its own: the filter must drop it.
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

test_that("the collinearity filter drops a near-duplicate predictor", {
  skip_screen_deps()
  dat <- make_collinear_data()
  scr <- cast_select(dat, seed = 22, verbose = FALSE)
  expect_s3_class(scr, "cast_select")
  expect_identical(scr$method, "screen")
  # One of the two collinear partners must be thinned, never both kept.
  expect_true(any(scr$scores$collinear_thinned[scr$scores$variable %in% c("x1", "x3")]))
  expect_false(all(c("x1", "x3") %in% scr$diagnostics$filtered))
})

test_that("the screen selects from the filtered set and records ranks", {
  skip_screen_deps()
  dat <- make_collinear_data()
  scr <- cast_select(dat, seed = 23, verbose = FALSE)
  d <- scr$diagnostics
  expect_true(all(scr$selected %in% union(d$filtered, d$keep)))
  expect_true(all(c("filtered", "lambda", "coefficients", "n_pres",
                    "final_ranks") %in% names(d)))
  expect_true(nrow(d$final_ranks) >= 1L)
  # No null distribution, no p-values, no count cap leftovers.
  expect_false(any(c("p_value", "null_threshold", "interventional_effect",
                     "passed_null", "fallback") %in% names(scr$scores)))
  # Screened admissions carry an elastic-net rank; prespecified entries are
  # marked by design.
  adm <- scr$scores$selected_reason[scr$scores$selected] == "screen"
  expect_true(all(is.finite(scr$scores$step_added[scr$scores$selected][adm])))
})

test_that("the filter ranks a U-shaped driver above noise (poly2 signal)", {
  skip_screen_deps()
  set.seed(31)
  n <- 400
  xu <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(-1 + 2.2 * xu^2)),
    xu = xu, noise = rnorm(n), x3 = rnorm(n), x4 = rnorm(n), x5 = rnorm(n)
  )
  scr <- cast_select(dat, seed = 32, verbose = FALSE)
  p <- stats::setNames(scr$scores$stage1_p, scr$scores$variable)
  expect_true(is.finite(p[["xu"]]) && is.finite(p[["noise"]]))
  expect_lt(p[["xu"]], p[["noise"]])
})

test_that("method = 'full' keeps every predictor and skips both stages", {
  dat <- make_collinear_data(n = 120)
  scr <- cast_select(dat, method = "full", verbose = FALSE)
  expect_setequal(scr$selected, c("x1", "x2", "x3", "x4", "x5"))
  expect_identical(scr$method, "full")
  expect_true(all(is.na(scr$scores$step_added)))
})

test_that("retired or invalid selection requests fail without running a replacement", {
  dat <- make_collinear_data(n = 150)
  expect_error(cast_select(dat, alpha = 0.01, min_vars = 3L), "unused argument")
  expect_error(cast_select(dat, environment = "x1"), "unused argument")
  expect_error(cast_select(dat, n_perm = 19L), "unused argument")
  expect_error(cast_select(dat, shift_size = 1), "unused argument")
  expect_error(cast_select(dat, num_trees = 60L), "unused argument")
  for (method in list("cpi", "dml", "rf", "typo", "two_stage", "two", NULL,
                      character(), NA_character_, 1, c("full", "screen"))) {
    expect_error(cast_select(dat, method = method), "method.*must be")
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

test_that("prespecified predictors survive thinning and are marked by design", {
  skip_screen_deps()
  dat <- make_collinear_data(n = 150)
  required <- c("x1", "x3")
  scr <- suppressWarnings(cast_select(dat, keep = required,
    seed = 41, verbose = FALSE))
  expect_true(all(required %in% scr$selected))
  rows <- scr$scores$variable %in% required
  expect_true(all(scr$scores$kept_by_design[rows]))
  expect_false(any(scr$scores$collinear_thinned[rows]))
  expect_true(all(scr$scores$selected_reason[rows] == "prespecified"))
  expect_identical(scr$diagnostics$keep, required)
  expect_true(all(scr$scores$step_added[rows] == 0L))
  reported <- cast_importance(scr)$effects
  expect_true(all(reported$kept_by_design[reported$variable %in% required]))
  expect_true(all(reported$selected_reason[reported$variable %in% required] == "prespecified"))
  expect_false("perm_importance" %in% names(reported))
})

test_that("prespecified retention coexists with screened admissions", {
  skip_screen_deps()
  dat <- make_collinear_data(n = 150)
  scr <- suppressWarnings(cast_select(dat, keep = "x5",
    seed = 42, verbose = FALSE))
  expect_true("x5" %in% scr$selected)
  expect_true(scr$scores$kept_by_design[scr$scores$variable == "x5"])
  full <- cast_select(dat, method = "full", keep = "x5", verbose = FALSE)
  expect_equal(full$scores$variable[full$scores$kept_by_design], "x5")
  expect_setequal(full$selected, get_env_vars(dat))
})

test_that("invalid or unaccommodated prespecified sets fail explicitly", {
  dat <- make_collinear_data(n = 120)
  for (required in list("unknown", "lon", "presence", c("x1", "x1"), NA_character_, 1)) {
    expect_error(cast_select(dat, keep = required, verbose = FALSE), "keep")
  }
  dat$x5 <- 0
  expect_error(cast_select(dat, keep = "x5", verbose = FALSE), "vary in the training data")
})

test_that("nested CV retains the same prespecified variables within each training fold", {
  skip_screen_deps()
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 200)
  required <- c("x1", "x3")
  cv <- suppressWarnings(cast_cv(dat, k = 2L, models = "rf", rf_ntree = 20L,
    select_args = list(keep = required),
    seed = 43, verbose = FALSE))
  expect_length(cv$screens, 2L)
  expect_true(all(vapply(cv$screens, function(s) {
    !is.null(s) && identical(s$diagnostics$keep, required) && all(required %in% s$selected)
  }, logical(1))))
})

test_that("stage functions forward prespecified retention to CV folds", {
  skip_screen_deps()
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 160)
  required <- c("x1", "x3")
  split <- cast_prepare(dat, seed = 44, verbose = FALSE)
  screen <- suppressWarnings(cast_select(
    split$train, keep = required, seed = 44, verbose = FALSE))
  cv <- suppressWarnings(cast_cv(
    split$train, select_args = list(keep = required),
    k = 2L, models = "rf", rf_ntree = 20L, seed = 44, verbose = FALSE))
  expect_true(all(required %in% screen$selected))
  expect_false(is.null(cv))
  expect_true(all(vapply(cv$screens, function(s) {
    !is.null(s) && identical(s$diagnostics$keep, required) && all(required %in% s$selected)
  }, logical(1))))
})

test_that("explicit predictor sets preserve membership and caller order", {
  dat <- make_collinear_data(n = 120)
  requested <- c("x3", "x1", "x2")
  scr <- cast_select(dat, method = "full", env_vars = requested, verbose = FALSE)
  expect_identical(scr$selected, requested)
  expect_identical(scr$scores$variable, requested)
  # A real numeric predictor may deliberately override a reserved metadata name.
  dat$fraction <- dat$x4
  requested <- c("fraction", "x1", "x2")
  expect_identical(cast_select(dat, method = "full", env_vars = requested,
                               verbose = FALSE)$selected, requested)
})

test_that("invalid explicit predictor sets fail before selection", {
  dat <- make_collinear_data(n = 120)
  for (bad in list(c("x1", "x2", "missing"), c("x1", "x1", "x2"),
                   c("x1", "x2", "lon"), c("x1", "x2", "presence"),
                   c("x1", "x2", NA_character_), 1:3)) {
    expect_error(cast_select(dat, method = "full", env_vars = bad,
                             verbose = FALSE), "env_vars")
  }
  dat$label <- rep("river", nrow(dat))
  expect_error(cast_select(dat, method = "full",
                           env_vars = c("x1", "x2", "label"), verbose = FALSE),
                "Non-numeric predictor")
})

test_that("fold-constant predictors are dropped and metadata follow survivors", {
  skip_screen_deps()
  dat <- make_collinear_data(n = 120)
  dat$x5 <- 0
  requested <- c("x3", "x1", "x2", "x5")
  ids <- c(x5 = "id5", x2 = "id2", x1 = "id1", x3 = "id3")
  categories <- c(x1 = "climate", x2 = "terrain", x5 = "other", x3 = "water")
  flags <- c(x2 = FALSE, x5 = FALSE, x3 = FALSE, x1 = TRUE)
  expect_warning(scr <- cast_select(dat, env_vars = requested,
    var_id = ids, category = categories, ensure = list(climate = flags),
    seed = 45, verbose = FALSE), "lack supported variation")
  expect_false("x5" %in% scr$selected)
  expect_false("x5" %in% scr$scores$variable)
  expect_true("x1" %in% scr$selected)
  # x1 is both the declared coverage member and a screened survivor; the
  # reason records how it entered, so either admission path is valid.
  expect_true(scr$scores$selected_reason[scr$scores$variable == "x1"] %in%
                c("screen", "ensured"))
  expect_error(suppressWarnings(cast_select(dat, env_vars = requested,
    keep = "x5", verbose = FALSE)), "vary in the training data")
  expect_error(cast_select(dat, method = "full", env_vars = c("x3", "x1", "x2"),
    var_id = c(x3 = "id3", x1 = "id1"), verbose = FALSE), "var_id.*cover")
  expect_error(cast_select(dat, method = "full", env_vars = c("x3", "x1", "x2"),
    category = c(x3 = "water", x1 = "climate"), verbose = FALSE), "category.*cover")
})
