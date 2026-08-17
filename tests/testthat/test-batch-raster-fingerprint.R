# C3/M13 regressions: SpatRaster-safe batch cache signatures.

test_that(".cast_cfg_fingerprint replaces SpatRaster with a stable fingerprint", {
  skip_if_not_installed("terra")
  mk_lyr <- function(v) {
    terra::setValues(
      terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5),
      v
    )
  }
  r <- c(mk_lyr(runif(25)), mk_lyr(runif(25)))
  names(r) <- c("x1", "x2")

  cfg <- list(
    raster_stack = r,
    future_rasters = list(ssp1 = r),
    raster_mask = NULL,
    response = "presence"
  )
  fp <- .cast_cfg_fingerprint(cfg)
  expect_false(inherits(fp$raster_stack, "SpatRaster"))
  expect_false(inherits(fp$future_rasters$ssp1, "SpatRaster"))
  expect_identical(fp$response, "presence")

  # Digesting the fingerprint is safe (that is the whole point of C3) and
  # deterministic for unchanged rasters.
  sig1 <- .cast_digest(list(cfg = fp))
  sig2 <- .cast_digest(list(cfg = .cast_cfg_fingerprint(cfg)))
  expect_identical(sig1, sig2)

  # ... but changes when the raster values change.
  cfg_alt <- cfg
  cfg_alt$raster_stack <- r * 2
  sig3 <- .cast_digest(list(cfg = .cast_cfg_fingerprint(cfg_alt)))
  expect_false(identical(sig1, sig3))
})

test_that("cast_batch runs with an in-memory SpatRaster raster_stack (C3)", {
  skip_if_not_installed("terra")
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(101)
  n <- 120
  dat <- data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = rbinom(n, 1, 0.5),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  r <- c(
    terra::setValues(terra::rast(nrows = 10, ncols = 10, xmin = 100,
                                 xmax = 110, ymin = 30, ymax = 40), runif(100)),
    terra::setValues(terra::rast(nrows = 10, ncols = 10, xmin = 100,
                                 xmax = 110, ymin = 30, ymax = 40), runif(100)),
    terra::setValues(terra::rast(nrows = 10, ncols = 10, xmin = 100,
                                 xmax = 110, ymin = 30, ymax = 40), runif(100))
  )
  names(r) <- c("x1", "x2", "x3")

  td <- tempfile("batchrast")
  # Reaching the expectation at all is the regression: before the fix the
  # step-cache digest of the SpatRaster cfg segfaulted the session.
  batch <- cast_batch(
    list(sp_a = dat),
    models = "rf",
    do_cv = TRUE, cv_k = 2L,
    select_method = "rf", select_min_vars = 2, select_num_trees = 40,
    min_occ = 20L,
    raster_stack = r,
    parallel = FALSE, seed = 102, output_dir = td, verbose = FALSE
  )
  expect_true(inherits(batch$results$sp_a, "cast_result"))
  expect_true(
    file.exists(file.path(td, "sp_a", "rasters", "current_hss_ensemble.tif"))
  )
  unlink(td, recursive = TRUE)
})

test_that("predict-step cache invalidates when env_data changes (M13)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(103)
  n <- 120
  dat <- data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = rbinom(n, 1, 0.5),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  env1 <- data.frame(
    lon = runif(40, 100, 110), lat = runif(40, 30, 40),
    x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40)
  )
  env2 <- env1
  env2$x1 <- env2$x1 + 10

  run_batch <- function(env) {
    cast_batch(
      list(sp_a = dat), env_data = env,
      models = "rf", do_cv = FALSE, do_ensemble = FALSE,
      select_method = "rf", select_min_vars = 2, select_num_trees = 40,
      parallel = FALSE, seed = 104, output_dir = td, verbose = FALSE
    )
  }
  td <- tempfile("batchenv")
  run_batch(env1)
  ckpt <- file.path(td, "sp_a", ".steps", "predict.rds")
  sig1 <- readRDS(ckpt)$params$signature

  run_batch(env2)
  sig2 <- readRDS(ckpt)$params$signature
  expect_false(identical(sig1, sig2))

  # Same inputs again: the signature is deterministic (cache hit, file
  # untouched, same signature readable).
  run_batch(env1)
  sig3 <- readRDS(ckpt)$params$signature
  expect_identical(sig1, sig3)
  unlink(td, recursive = TRUE)
})
