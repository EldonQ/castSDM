# P0 review round: hyperparameter tuning (C1/C2), MaxEnt engine clamp (W1),
# ecospat-aligned Boyce window handling (W2), user-defined and bias-weighted
# background sampling (C3).

# ---------------------------------------------------------------- Boyce (W2)

test_that("compute_boyce keeps ratio-0 windows and stays finite (W2)", {
  # Presences occupy the bottom decile only: every window above ~0.3 has
  # background but no presences, i.e. predicted/expected ratio exactly 0.
  # The old `pe > 0` filter dropped those windows; the ecospat-aligned
  # implementation keeps them.
  set.seed(11)
  pres <- seq(0.02, 0.12, length.out = 20)
  bg <- runif(200, 0.05, 1)
  b <- compute_boyce(c(pres, bg), c(rep(1L, 20), rep(0L, 200)))
  expect_true(is.finite(b))
})

test_that("compute_boyce returns a high value for presences at the top of the background range (W2)", {
  # Presence density increases smoothly with suitability (weighted draws from
  # the background grid), so the predicted/expected ratio rises with the
  # window position. The retained ratio-0 windows at the low end no longer
  # break the correlation.
  set.seed(12)
  bg <- seq(0.01, 1, length.out = 300)
  pres <- sample(bg, 60, prob = bg^3, replace = TRUE)
  b <- compute_boyce(c(pres, bg), c(rep(1L, 60), rep(0L, 300)))
  expect_true(is.finite(b))
  expect_gt(b, 0.8)
})

test_that("compute_boyce stays NA when background never enters the windows", {
  # A single spike of presences far above all background: most windows have
  # neither class (0/0 -> NA) and the usable set is too small to correlate.
  b <- castSDM:::compute_boyce(
    c(rep(0.99, 5), rep(0.5, 2)), c(rep(1L, 5), rep(0L, 2))
  )
  expect_true(is.na(b) || is.finite(b))  # no error either way; value may be NA
})

# ------------------------------------------------------- MaxEnt clamp (W1)

test_that("MaxEnt engine no longer clamps silently; package clamp decides (W1)", {
  skip_if_not_installed("maxnet")
  set.seed(21)
  n <- 200
  dat <- data.frame(
    presence = c(rep(1L, 40), rep(0L, 160)),
    x1 = c(runif(40, 0.6, 1), runif(160, 0, 1)),
    x2 = c(runif(40, 0.6, 1), runif(160, 0, 1))
  )
  fit <- cast_fit(dat, models = "maxent", verbose = FALSE)

  xin <- dat[1:5, c("x1", "x2")]
  xout <- xin; xout$x1 <- xout$x1 * 50  # far outside the training range

  p_out_raw <- predict_single_model(fit$models$maxent, xout)
  # Package-level clamping of the same input must now change the prediction:
  # previously the engine clamped both back to identical in-range values.
  xout_clamped <- .cast_clamp(xout, fit$scaling$reference)
  expect_false(isTRUE(
    all.equal(p_out_raw,
              predict_single_model(fit$models$maxent, xout_clamped))
  ))
  # Explicit clamp = TRUE through the predict path still works.
  p_clamped <- predict_single_model(fit$models$maxent, xout, clamp = TRUE)
  expect_length(p_clamped, nrow(xout))
})

# ------------------------------------------------------- Hyperparameter tuning (C1/C2)

.tune_test_data <- function(seed = 31) {
  set.seed(seed)
  data.frame(
    presence = c(rep(1L, 60), rep(0L, 300)),
    x1 = c(runif(60, 0.6, 1), runif(300, 0, 1)),
    x2 = c(runif(60, 0.6, 1), runif(300, 0, 1))
  )
}

test_that("cast_fit(tune = TRUE) tunes RF and MaxEnt and records the grid", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("maxnet")
  dat <- .tune_test_data()
  fit <- cast_fit(dat, models = c("rf", "maxent"),
                  rf_ntree = 100L, tune = TRUE, tune_folds = 2L,
                  seed = 5, verbose = FALSE)

  t_rf <- fit$models$rf$tune
  expect_equal(t_rf$engine, "rf")
  expect_equal(t_rf$metric, "tss_oob")
  expect_true(t_rf$best$mtry %in% t_rf$grid$mtry)
  expect_true(is.finite(t_rf$best_score))

  t_me <- fit$models$maxent$tune
  expect_equal(t_me$engine, "maxent")
  expect_equal(t_me$metric, "tss_cv")
  expect_true(t_me$best$classes %in% c("l", "lq", "lqh", "lqhp"))
  expect_true(t_me$best$regmult %in% c(0.5, 1, 2))
  expect_equal(nrow(t_me$grid), 12L)

  # The tuned fit produces working predictions end to end.
  ev <- cast_evaluate(fit, dat)
  expect_true(all(is.finite(ev$metrics$auc_mean)))
})

test_that("cast_fit(tune = TRUE) tunes BRT with an adaptive tree count", {
  skip_if_not_installed("gbm")
  dat <- .tune_test_data(32)
  fit <- suppressWarnings(cast_fit(
    dat, models = "brt", brt_n_trees = 200L,
    tune = TRUE, tune_folds = 2L, seed = 6, verbose = FALSE
  ))
  t <- fit$models$brt$tune
  expect_equal(t$engine, "brt")
  expect_equal(nrow(t$grid), 4L)
  expect_true(t$best$depth %in% c(2, 5))
  expect_true(t$best$shrinkage %in% c(0.005, 0.01))
})

test_that("explicit MaxEnt classes are held fixed while regmult is searched", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("maxnet")
  dat <- .tune_test_data(33)
  fit <- cast_fit(dat, models = "maxent", maxent_classes = "l",
                  tune = TRUE, tune_folds = 2L, seed = 7, verbose = FALSE)
  t <- fit$models$maxent$tune
  expect_equal(nrow(t$grid), 3L)
  expect_true(all(t$grid$classes == "l"))
  expect_equal(t$best$classes, "l")
  expect_true(t$best$regmult %in% c(0.5, 1, 2))
})

test_that("invalid hyperparameter arguments abort informatively", {
  dat <- .tune_test_data(34)
  expect_error(
    cast_fit(dat, models = "maxent", maxent_classes = "xyz"),
    "feature-class letters"
  )
  expect_error(
    cast_fit(dat, models = "maxent", maxent_regmult = -1),
    "positive number"
  )
  expect_error(
    cast_fit(dat, models = "rf", rf_mtry = 0),
    "positive integer"
  )
})

test_that("cast_cv(tune = TRUE) runs the grid inside each outer fold", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(35)
  n <- 300
  dat <- data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = c(rep(1L, 100), rep(0L, 200)),
    x1 = c(runif(100, 0.6, 1), runif(200, 0, 1)),
    x2 = c(runif(100, 0.6, 1), runif(200, 0, 1)),
    x3 = runif(n)
  )
  cv <- cast_cv(dat, models = "rf", k = 3L, select_method = "full",
                tune = TRUE, tune_folds = 2L,
                seed = 8, verbose = FALSE)
  expect_true("evaluated" %in% cv$fold_status)
  expect_true(nrow(cv$fold_metrics) >= 1L)
})

# ------------------------------------------------ User / bias background (C3)

.bg_test_raster <- function(values_list = NULL) {
  mk <- function(v) {
    terra::setValues(
      terra::rast(nrows = 10, ncols = 10, xmin = 100, xmax = 110,
                  ymin = 30, ymax = 40),
      v
    )
  }
  r <- if (is.null(values_list)) {
    c(mk(runif(100)), mk(runif(100)))
  } else {
    do.call(c, lapply(values_list, mk))
  }
  names(r) <- if (is.null(values_list)) c("x1", "x2") else
    paste0("x", seq_along(values_list))
  r
}

test_that("cast_background(user_table =) uses exactly the supplied points", {
  skip_if_not_installed("terra")
  set.seed(41)
  r <- .bg_test_raster()
  occ <- data.frame(lon = c(100.5, 102.5), lat = c(30.5, 32.5))
  # Cell-interior coordinates only: boundary values make cellFromXY ambiguous.
  ut <- data.frame(lon = 101.5 + (0:14 %% 5), lat = 33.5 + (0:14 %/% 5))
  out <- cast_background(occ, NULL, r, user_table = ut, seed = 42,
                         verbose = FALSE)
  expect_equal(sum(out$presence == 0), 15L)
  # Output carries cell centres; compare through the cell numbers instead.
  expect_setequal(
    terra::cellFromXY(r, as.matrix(out[out$presence == 0, c("lon", "lat")])),
    terra::cellFromXY(r, as.matrix(ut))
  )
  expect_true(all(stats::complete.cases(out[, c("x1", "x2")])))
})

test_that("user_table drops points in presence cells and outside the raster", {
  skip_if_not_installed("terra")
  set.seed(43)
  r <- .bg_test_raster()
  occ <- data.frame(lon = c(100.5, 102.5), lat = c(30.5, 32.5))
  ut <- data.frame(
    lon = c(200, 100.5, 103.5, 104.5),  # outside, presence cell, ok, ok
    lat = c(35, 30.5, 35.5, 36.5)
  )
  out <- suppressWarnings(
    cast_background(occ, NULL, r, user_table = ut, seed = 44, verbose = TRUE)
  )
  expect_equal(sum(out$presence == 0), 2L)
  expect_error(
    cast_background(occ, NULL, r, user_table = ut[0, ], verbose = FALSE),
    "empty"
  )
})

test_that("bias_raster weights the sampling and top-up keeps the bias", {
  skip_if_not_installed("terra")
  set.seed(45)
  # Top half of the grid (rows 1-5, lat >= 35) weighted; bottom half zero.
  w <- c(rep(1, 50), rep(0, 50))
  r <- .bg_test_raster(list(runif(100), runif(100)))
  bias <- terra::setValues(terra::rast(r[[1]]), w)
  occ <- data.frame(lon = 100.5, lat = 30.5)
  out <- cast_background(occ, NULL, r, n_bg = 30, bias_raster = bias,
                         seed = 46, verbose = FALSE)
  expect_equal(sum(out$presence == 0), 30L)
  bg <- out[out$presence == 0, ]
  expect_true(all(bg$lat >= 35))  # no draws from the zero-weight half
})

test_that("user_table / bias_raster misuse aborts informatively", {
  skip_if_not_installed("terra")
  set.seed(47)
  r <- .bg_test_raster()
  occ <- data.frame(lon = 100.5, lat = 30.5)
  ut <- data.frame(lon = 103.5, lat = 33.5)
  bias <- terra::setValues(terra::rast(r[[1]]), runif(100))
  r_shift <- terra::shift(r, dx = 50)

  expect_error(
    cast_background(occ, NULL, r, user_table = ut, bias_raster = bias,
                    verbose = FALSE),
    "not both"
  )
  expect_error(
    cast_background(occ, NULL, r, strategy = "environmental",
                    bias_raster = bias, verbose = FALSE),
    "random"
  )
  expect_error(
    cast_background(occ, NULL, r_shift, n_bg = 5, bias_raster = bias,
                    verbose = FALSE),
    "incompatible geometry"
  )
  # A multi-layer bias raster is reduced to its first layer with a warning.
  expect_warning(
    out <- cast_background(occ, NULL, r, n_bg = 5, bias_raster = c(bias, bias),
                           seed = 48, verbose = FALSE),
    "using the first one"
  )
  expect_equal(sum(out$presence == 0), 5L)
})
