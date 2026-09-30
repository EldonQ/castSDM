# A2 deep-review round: S1 threshold-free/kappa metrics, S7 environmental
# blocking, S3 area of applicability, S6 grouped/decay/permutation ensembles.

# ---- Shared fixtures --------------------------------------------------------

.make_a2_data <- function(seed = 301, n = 240) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  data.frame(
    lon = runif(n, 0, 1), lat = runif(n, 0, 1),
    presence = rbinom(n, 1, plogis(1.5 * x1)),
    x1 = x1, x2 = x2, x3 = x3
  )
}

.make_a2_two_model_fixture <- function(seed = 303, n = 200) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 0, 1), lat = runif(n, 0, 1),
    presence = rbinom(n, 1, plogis(1.2 * x1 - x2)),
    x1 = x1, x2 = x2
  )
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 40, seed = seed + 1L, verbose = FALSE)
  # Hand-built cv metrics (same pattern as test-ensemble-raster-renorm.R):
  # rf scores mean(0.8, 0.6, 0.5) = 0.6333; gam mean(0.7, 0.55, 0.45) = 0.5667.
  cv <- new_cast_cv(
    metrics = data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.90, 0.85), auc_sd = c(0.1, 0.1),
      auc_n_folds = c(3, 3),
      tss_mean = c(0.60, 0.55), tss_sd = c(0.1, 0.1),
      tss_n_folds = c(3, 3),
      cbi_mean = c(0.50, 0.45), cbi_sd = c(0.1, 0.1),
      cbi_n_folds = c(3, 3),
      n_folds = c(3, 3), n_selected_mean = c(2, 2)
    ),
    fold_metrics = data.frame(), folds = integer(n), k = 3L,
    block_method = "grid", thresholds = c(rf = 0.5, gam = 0.5)
  )
  list(fit = fit, cv = cv, dat = dat)
}

# ---- S1: kappa / omission@E / MPA -------------------------------------------

test_that("compute_kappa reproduces hand-computed confusion tables", {
  obs <- c(1L, 1L, 0L, 0L)
  pred <- c(0.9, 0.8, 0.2, 0.1)
  expect_equal(castSDM:::compute_kappa(obs, pred, 0.5), 1)
  obs2 <- c(1L, 0L, 1L, 0L)
  pred2 <- c(0.9, 0.8, 0.2, 0.1)
  expect_equal(castSDM:::compute_kappa(obs2, pred2, 0.5), 0)
  expect_true(is.na(castSDM:::compute_kappa(rep(1L, 4), pred, 0.5)))
  expect_true(is.na(castSDM:::compute_kappa(obs, pred, NA_real_)))
  # NAs in pred are dropped, not propagated
  expect_equal(castSDM:::compute_kappa(obs, c(NA, pred[-1]), 0.5),
               castSDM:::compute_kappa(obs[-1], pred[-1], 0.5))
})

test_that("compute_omission_at and compute_mpa match their definitions", {
  pred <- seq(0.1, 1, by = 0.1)
  obs <- c(1L, rep(0L, 9))
  # e = 0.1 -> cut = type-1 q90 of all predictions = 0.9; presence 0.1 < 0.9
  expect_equal(castSDM:::compute_omission_at(pred, obs, e = 0.1), 1)
  # e = 0.9 -> cut = 0.1; presence 0.1 is not below the cut
  expect_equal(castSDM:::compute_omission_at(pred, obs, e = 0.9), 0)
  expect_true(is.na(castSDM:::compute_omission_at(c(.5, .6), c(0L, 0L))))

  obs_m <- c(0L, 0L, 0L, 0L, rep(1L, 6))
  # keep = 0.9 -> cut = type-1 q10 of presence predictions = 0.5;
  # fraction of ALL predictions >= 0.5 is 6/10
  expect_equal(castSDM:::compute_mpa(pred, obs_m, keep = 0.9), 0.6)
  expect_true(is.na(castSDM:::compute_mpa(pred, c(rep(1L, 4), rep(0L, 6)))))
})

test_that("evaluate_model_full reports the new threshold-free metrics", {
  set.seed(11)
  pred <- plogis(rnorm(200))
  obs <- rbinom(200, 1, pred)
  met <- castSDM:::evaluate_model_full(pred, obs)
  expect_true(all(c("kappa", "omission_5", "omission_10", "mpa")
                  %in% names(met)))
  expect_true(is.finite(met[["kappa"]]))
  expect_true(met[["omission_5"]] >= 0 && met[["omission_5"]] <= 1)
  # empty input falls back to the NA vector of the same length
  met_na <- castSDM:::evaluate_model_full(c(NA, NA), c(1L, 0L))
  expect_true(all(is.na(met_na)))
  expect_identical(names(met_na), names(met))
})

# ---- S7: environmental blocking ---------------------------------------------

test_that("make_spatial_folds supports env blocking and guards degeneracy", {
  set.seed(21)
  n <- 60
  env <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  f <- castSDM:::make_spatial_folds(runif(n), runif(n), k = 3,
                                    method = "env", env = env)
  expect_setequal(unique(f), 1:3)
  expect_length(f, n)
  # Two identical columns carry no environmental variance
  expect_error(
    castSDM:::make_spatial_folds(runif(n), runif(n), k = 3, method = "env",
                                 env = cbind(x1 = rep(1, n), x2 = rep(2, n))),
    "non-zero variance"
  )
})

test_that("cast_cv(block_method = 'env') runs and reports the method", {
  skip_if_not_installed("ranger")
  dat <- .make_a2_data(seed = 41, n = 180)
  cv <- cast_cv(dat, select_method = "full", k = 3, models = "rf",
                rf_ntree = 40, block_method = "env", seed = 42,
                verbose = FALSE)
  expect_identical(cv$block_method, "env")
  expect_true(nrow(cv$fold_metrics) > 0)
})

test_that("cast_prepare labels the env split differently", {
  dat <- .make_a2_data(seed = 43, n = 180)
  prep <- cast_prepare(dat, train_fraction = 0.7, split = "spatial",
                       block_method = "env", seed = 44, verbose = FALSE)
  expect_match(prep$split, "environmental block")
  prep_sp <- cast_prepare(dat, train_fraction = 0.7, split = "spatial",
                          block_method = "grid", seed = 45, verbose = FALSE)
  expect_match(prep_sp$split, "spatial block")
})

# ---- S3: area of applicability ----------------------------------------------

test_that("AOA space standardisation and VI weighting are exact", {
  set.seed(51)
  X <- data.frame(x1 = rnorm(200, 10, 3), x2 = rnorm(200, -5, 0.5))
  prep <- castSDM:::.cast_aoa_prep(X, X[1:5, ])
  expect_equal(colMeans(prep$train), c(x1 = 0, x2 = 0), tolerance = 1e-10)
  expect_equal(apply(prep$train, 2, sd), c(x1 = 1, x2 = 1), tolerance = 1e-10)
  # identical rows project onto each other exactly
  expect_equal(prep$new[1, ], prep$train[1, ], tolerance = 1e-10)

  prep_w <- castSDM:::.cast_aoa_prep(X, X[1:5, ], vi = c(x1 = 3, x2 = 1))
  expect_true(prep_w$weighted)
  # columns scaled by sqrt(w / sum(w)) = sqrt(0.75), sqrt(0.25)
  expect_equal(apply(prep_w$train, 2, sd), c(x1 = sqrt(0.75), x2 = 0.5),
               tolerance = 1e-10)
  expect_error(
    castSDM:::.cast_aoa_prep(X, X, vi = c(x1 = 1)),
    "missing importance"
  )
  # NA rows are median-imputed, not propagated
  X_na <- X[1:5, ]
  X_na[1, "x1"] <- NA
  prep_na <- castSDM:::.cast_aoa_prep(X, X_na)
  expect_true(all(is.finite(prep_na$new)))
})

test_that("AOA DI separates coincident from distant points", {
  set.seed(52)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  prep <- castSDM:::.cast_aoa_prep(X, X)
  di <- castSDM:::.cast_aoa_di(prep$train, prep$new)
  # each training point's nearest neighbour is itself (distance 0)
  expect_equal(di, rep(0, 50), tolerance = 1e-12)
  far <- data.frame(x1 = rep(100, 3), x2 = rep(-100, 3))
  prep_far <- castSDM:::.cast_aoa_prep(X, far)
  di_far <- castSDM:::.cast_aoa_di(prep$train, prep_far$new)
  expect_true(all(di_far > 50))
})

test_that("cast_cv(aoa = TRUE) calibrates a threshold from held-out folds", {
  skip_if_not_installed("ranger")
  dat <- .make_a2_data(seed = 61, n = 240)
  cv <- cast_cv(dat, select_method = "full", k = 3, models = "rf",
                rf_ntree = 40, aoa = TRUE, seed = 62, verbose = FALSE)
  expect_true(isTRUE(cv$aoa$enabled))
  expect_identical(cv$aoa$method, "cv_heldout_di_p95")
  expect_true(is.finite(cv$aoa$threshold) && cv$aoa$threshold >= 0)
  expect_length(cv$aoa$di_oof, nrow(dat))
  expect_true(all(cv$aoa$di_oof[is.finite(cv$aoa$di_oof)] >= 0))
  expect_true(cv$aoa$n_di > 0)
  expect_false(cv$aoa$weighted)
  # plain cv carries no calibration
  cv0 <- cast_cv(dat, select_method = "full", k = 3, models = "rf",
                 rf_ntree = 40, seed = 63, verbose = FALSE)
  expect_null(cv0$aoa)
  expect_true(all(c("kappa", "omission_5", "omission_10", "mpa")
                  %in% names(cv0$fold_metrics)))
  expect_true("kappa_mean" %in% names(cv0$metrics))
})

test_that("cast_predict flags the area of applicability", {
  skip_if_not_installed("ranger")
  dat <- .make_a2_data(seed = 61, n = 240)
  cv <- cast_cv(dat, select_method = "full", k = 3, models = "rf",
                rf_ntree = 40, aoa = TRUE, seed = 62, verbose = FALSE)
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 64, verbose = FALSE)

  inside <- data.frame(lon = runif(20), lat = runif(20),
                       x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
  pr <- cast_predict(fit, inside, aoa = TRUE, cv = cv, extrapolation = FALSE)
  expect_true(all(c("aoa", "aoa_di") %in% names(pr$predictions)))
  expect_type(pr$predictions$aoa, "logical")
  expect_true(all(pr$predictions$aoa_di >= 0, na.rm = TRUE))
  # deep extrapolation must fall outside the calibrated applicability
  far <- data.frame(lon = c(0.5, 0.5), lat = c(0.5, 0.5),
                    x1 = c(50, -50), x2 = c(0, 0), x3 = c(0, 0))
  pr_far <- cast_predict(fit, far, aoa = TRUE, cv = cv,
                         extrapolation = FALSE)
  expect_true(all(!pr_far$predictions$aoa))
  expect_true(all(pr_far$predictions$aoa_di > cv$aoa$threshold))

  # vi-mismatch warning: threshold calibrated uniform, DI weighted
  expect_warning(
    cast_predict(fit, inside, aoa = TRUE, cv = cv,
                 vi = c(x1 = 2, x2 = 1, x3 = 1), extrapolation = FALSE),
    "uniform"
  )
  # aoa = TRUE without a calibrated cv errors
  expect_error(
    cast_predict(fit, inside, aoa = TRUE, extrapolation = FALSE),
    "cast_cv"
  )
})

test_that("cast_ensemble_raster writes AOA layers and rejects bad calibrations", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("terra")
  dat <- .make_a2_data(seed = 71, n = 240)
  cv <- cast_cv(dat, select_method = "full", k = 3, models = "rf",
                rf_ntree = 40, aoa = TRUE, seed = 72, verbose = FALSE)
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 73, verbose = FALSE)

  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, nlyrs = 3)
  terra::values(r) <- cbind(rnorm(100), rnorm(100), rnorm(100))
  names(r) <- c("x1", "x2", "x3")
  td <- tempfile("aoarast")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  res <- cast_ensemble_raster(fit, cv, r, output_dir = td,
                              extrapolation = FALSE, aoa_cv = cv,
                              verbose = FALSE)
  expect_true(file.exists(res$aoa_path))
  expect_true(file.exists(res$aoa_di_path))
  aoa_vals <- terra::values(terra::rast(res$aoa_path))
  expect_true(all(na.omit(aoa_vals) %in% c(0, 1)))
  di_vals <- terra::values(terra::rast(res$aoa_di_path))
  expect_true(all(na.omit(di_vals) >= 0))

  # a cv without calibration is rejected
  cv_bad <- new_cast_cv(
    metrics = data.frame(model = "rf", auc_mean = 0.9, auc_sd = 0.1,
                         auc_n_folds = 3, tss_mean = 0.6, tss_sd = 0.1,
                         tss_n_folds = 3, cbi_mean = 0.5, cbi_sd = 0.1,
                         cbi_n_folds = 3, n_folds = 3, n_selected_mean = 2),
    fold_metrics = data.frame(), folds = integer(1), k = 3L,
    block_method = "grid", thresholds = c(rf = 0.5)
  )
  expect_error(
    cast_ensemble_raster(fit, cv_bad, r, output_dir = td, aoa_cv = cv_bad,
                         verbose = FALSE),
    "cast_cv"
  )
})

# ---- S6: decay weighting, grouped ensembles, permutation importance ---------

test_that("decay power weighting behaves as documented", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("mgcv")
  fx <- .make_a2_two_model_fixture()

  expect_error(
    cast_ensemble(fx$fit, fx$cv, fx$dat, decay = -1),
    "decay"
  )
  w1 <- cast_ensemble(fx$fit, fx$cv, fx$dat, decay = 1)$weights
  w0 <- cast_ensemble(fx$fit, fx$cv, fx$dat, decay = 0)$weights
  w3 <- cast_ensemble(fx$fit, fx$cv, fx$dat, decay = 3)$weights
  # decay = 0 degenerates to equal weights (as.numeric strips weight attrs)
  expect_equal(as.numeric(w0), c(0.5, 0.5))
  # higher scores get more weight; decay > 1 sharpens further
  expect_gt(w1[["rf"]], w1[["gam"]])
  expect_gt(w3[["rf"]], w1[["rf"]])
  # raw scores are untouched by the exponent
  ens3 <- cast_ensemble(fx$fit, fx$cv, fx$dat, decay = 3)
  expect_equal(unname(ens3$model_scores),
               unname(cast_ensemble(fx$fit, fx$cv, fx$dat)$model_scores),
               tolerance = 1e-12)
  expect_identical(ens3$decay, 3)
})

test_that("cast_ensemble_by builds one ensemble per group", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("mgcv")
  fx <- .make_a2_two_model_fixture()

  gb <- cast_ensemble_by(fx$fit, fx$cv, fx$dat,
                         groups = list(trees = "rf",
                                       both = c("rf", "gam")))
  expect_s3_class(gb, "cast_ensemble_grouped")
  expect_s3_class(gb$ensembles$trees, "cast_ensemble")
  expect_s3_class(gb$ensembles$both, "cast_ensemble")
  for (nm in c("trees", "both")) {
    expect_true(all(c(paste0("hss_", nm), paste0("hss_sd_", nm),
                      paste0("binary_", nm)) %in% names(gb$predictions)))
  }
  expect_identical(nrow(gb$predictions), nrow(fx$dat))
  # group ensembles equal standalone ensembles on the same model subset
  expect_equal(gb$ensembles$trees$weights,
               cast_ensemble(fx$fit, fx$cv, fx$dat, models = "rf")$weights)

  expect_error(
    cast_ensemble_by(fx$fit, fx$cv, fx$dat,
                     groups = list(c("rf", "gam"))),
    "named"
  )
  expect_error(
    cast_ensemble_by(fx$fit, fx$cv, fx$dat, groups = list(a = "bogus")),
    "available"
  )
})

test_that("cast_ensemble_importance ranks the informative predictor first", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("mgcv")
  # presence depends on x1 only; x3 is pure noise
  set.seed(83)
  n <- 200
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, plogis(1.8 * rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  dat$presence <- with(dat, rbinom(n, 1, plogis(1.8 * x1)))
  screen <- new_cast_select(c("x1", "x3"),
                            data.frame(variable = c("x1", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 40, seed = 84, verbose = FALSE)
  cv <- new_cast_cv(
    metrics = data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.9, 0.85), auc_sd = c(0.1, 0.1), auc_n_folds = c(3, 3),
      tss_mean = c(0.6, 0.55), tss_sd = c(0.1, 0.1), tss_n_folds = c(3, 3),
      cbi_mean = c(0.5, 0.45), cbi_sd = c(0.1, 0.1), cbi_n_folds = c(3, 3),
      n_folds = c(3, 3), n_selected_mean = c(2, 2)
    ),
    fold_metrics = data.frame(), folds = integer(n), k = 3L,
    block_method = "grid", thresholds = c(rf = 0.5, gam = 0.5)
  )

  imp <- cast_ensemble_importance(fit, cv, dat, n_perm = 6, seed = 85,
                                  verbose = FALSE)
  expect_s3_class(imp, "cast_ensemble_importance")
  expect_setequal(names(imp$importance),
                  c("variable", "importance_mean", "importance_sd"))
  expect_setequal(imp$importance$variable, c("x1", "x3"))
  expect_gt(imp$importance$importance_mean[
    imp$importance$variable == "x1"],
    imp$importance$importance_mean[
      imp$importance$variable == "x3"])
  expect_identical(imp$n_perm, 6L)
  expect_identical(imp$metric, "pearson")
  expect_length(imp$baseline, n)

  expect_error(
    cast_ensemble_importance(fit, cv, dat[1:2, ], n_perm = 2),
    "at least 3 rows"
  )
  expect_error(
    cast_ensemble_importance(fit, cv, dat, n_perm = 0),
    "positive integer"
  )
})
