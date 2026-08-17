# Raster ensemble: weight renormalisation, MESS layer, clamp ---------------

.make_ens_fixture <- function(seed = 201) {
  set.seed(seed)
  n <- 300
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 0, 1), lat = runif(n, 0, 1),
    presence = rbinom(n, 1, plogis(1.2 * x1 - x2)),
    x1 = x1, x2 = x2, x3 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 40, seed = seed + 1L, verbose = FALSE)
  # Both composite scores >= 0.5, so the initial weights are non-trivial
  # (rf ~ 0.53, gam ~ 0.47) and a dropped gam must trigger renormalisation.
  cv <- new_cast_cv(
    metrics = data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.90, 0.85), auc_sd = c(0.1, 0.1),
      tss_mean = c(0.60, 0.55), tss_sd = c(0.1, 0.1),
      cbi_mean = c(0.50, 0.45), cbi_sd = c(0.1, 0.1),
      n_folds = c(3, 3), n_selected_mean = c(3, 3)
    ),
    fold_metrics = data.frame(), folds = integer(n), k = 3L,
    block_method = "grid", thresholds = c(rf = 0.5, gam = 0.5)
  )
  list(fit = fit, cv = cv)
}

.make_ens_raster <- function(seed = 203, widen = FALSE) {
  set.seed(seed)
  scale_f <- if (widen) 4 else 1   # widen pushes cells off the training range
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, nlyrs = 3)
  terra::values(r) <- cbind(rnorm(400, sd = scale_f),
                            rnorm(400, sd = scale_f),
                            rnorm(400, sd = scale_f))
  names(r) <- c("x1", "x2", "x3")
  r
}

test_that("cast_ensemble_raster renormalises weights after excluding a model", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("terra")
  fx <- .make_ens_fixture()
  fit <- fx$fit
  fit$models$gam$model <- NULL  # gam predictions become all-NA -> excluded

  r <- .make_ens_raster()
  td <- tempfile("ensrast")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  expect_warning(
    res <- cast_ensemble_raster(fit, fx$cv, r, output_dir = td,
                                extrapolation = FALSE, verbose = FALSE),
    "non-finite"
  )
  hss <- terra::values(terra::rast(res$hss_path), mat = FALSE)

  # Reference: rf-only predictions over the same cells. Without
  # renormalisation the raster would hold weights["rf"] * HSS_rf (~0.53x).
  grid <- terra::as.data.frame(r, na.rm = FALSE)
  rf_ref <- cast_predict(fit, grid, models = "rf",
                         extrapolation = FALSE)$predictions$HSS_rf
  expect_equal(unname(hss), unname(rf_ref), tolerance = 1e-5)

  # Only one contributing model -> cross-model SD layer is all NA.
  sdv <- terra::values(terra::rast(res$hss_sd_path), mat = FALSE)
  expect_true(all(is.na(sdv)))
})

test_that("cast_ensemble_raster writes a MESS layer and clamp leaves MESS untouched", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("terra")
  fx <- .make_ens_fixture(seed = 211)
  r <- .make_ens_raster(seed = 213, widen = TRUE)

  td <- tempfile("ensmess")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res <- cast_ensemble_raster(fx$fit, fx$cv, r, output_dir = td,
                              verbose = FALSE)
  expect_true(!is.null(res$mess_path))
  expect_true(file.exists(res$mess_path))
  mess <- terra::values(terra::rast(res$mess_path), mat = FALSE)
  expect_true(all(is.finite(mess)))
  expect_true(any(mess < 0))   # widened raster leaves the training envelope

  # clamp = TRUE changes predictions but never the MESS layer (MESS is
  # computed on the unclamped input).
  td2 <- tempfile("ensclamp")
  on.exit(unlink(td2, recursive = TRUE), add = TRUE)
  res2 <- cast_ensemble_raster(fx$fit, fx$cv, r, output_dir = td2,
                               clamp = TRUE, verbose = FALSE)
  mess2 <- terra::values(terra::rast(res2$mess_path), mat = FALSE)
  expect_equal(unname(mess2), unname(mess), tolerance = 1e-5)

  # extrapolation = FALSE -> no MESS layer, NULL path.
  td3 <- tempfile("ensnomess")
  on.exit(unlink(td3, recursive = TRUE), add = TRUE)
  res3 <- cast_ensemble_raster(fx$fit, fx$cv, r, output_dir = td3,
                               extrapolation = FALSE, verbose = FALSE)
  expect_null(res3$mess_path)
})

test_that("cast_ensemble_raster skip check also covers hss_sd and mess", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("terra")
  fx <- .make_ens_fixture(seed = 221)
  r <- .make_ens_raster(seed = 223)

  td <- tempfile("ensskip")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res1 <- cast_ensemble_raster(fx$fit, fx$cv, r, output_dir = td,
                               verbose = FALSE)
  expect_true(!is.null(res1$weights))

  # All outputs present -> skipped, weights NULL.
  res2 <- cast_ensemble_raster(fx$fit, fx$cv, r, output_dir = td,
                               verbose = FALSE)
  expect_null(res2$weights)

  # hss_sd removed -> must recompute instead of skipping.
  unlink(res1$hss_sd_path)
  res3 <- cast_ensemble_raster(fx$fit, fx$cv, r, output_dir = td,
                               verbose = FALSE)
  expect_true(!is.null(res3$weights))
  expect_true(file.exists(res3$hss_sd_path))

  # mess removed -> must recompute instead of skipping.
  unlink(res1$mess_path)
  res4 <- cast_ensemble_raster(fx$fit, fx$cv, r, output_dir = td,
                               verbose = FALSE)
  expect_true(!is.null(res4$weights))
  expect_true(file.exists(res4$mess_path))
})
