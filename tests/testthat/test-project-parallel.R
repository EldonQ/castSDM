test_that("parallel future scenarios preserve raster values and scenario order", {
  skip_if_not_installed("terra")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("ranger")
  set.seed(19)
  d <- data.frame(lon = runif(80), lat = runif(80),
    x = runif(80, -1, 1), presence = rep(c(0L, 1L), 40))
  fit <- cast_fit(d, models = "rf", rf_ntree = 30L, verbose = FALSE)
  cv <- new_cast_cv(metrics = data.frame(model = "rf", auc_mean = .8,
    tss_mean = .6, n_folds = 2), fold_metrics = data.frame(),
    folds = rep(1:2, 40), k = 2L, block_method = "grid", thresholds = c(rf = .5),
    oof = data.frame(obs = d$presence, HSS_rf = seq(.05, .95, length.out = 80)),
    pred_folds = c(rf = 2L))
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  terra::values(r) <- seq(-.8, .8, length.out = 9); names(r) <- "x"
  futures <- list(a = r + .1, b = r + .2)
  names(futures$a) <- names(futures$b) <- "x"
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  a <- cast_project_raster(fit, cv, r, futures, tempfile(), method = "equal", verbose = FALSE)
  b <- cast_project_raster(fit, cv, r, futures, tempfile(), method = "equal",
    parallel = TRUE, verbose = FALSE)
  expect_identical(names(a$future), names(b$future))
  expect_equal(a$stats, b$stats)
  for (s in names(futures)) expect_equal(
    terra::values(terra::rast(a$future[[s]]$hss_path)),
    terra::values(terra::rast(b$future[[s]]$hss_path)))
})
