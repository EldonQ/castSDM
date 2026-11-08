# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
test_that("worker identity checks reject real differences with specific diagnostics", {
  expected <- list(version = "0.28.1", path = "/library/castSDM", R = "4.4.3", code = "same code")
  expect_true(.cast_check_compare_worker(expected, expected, 1L))
  for (field in names(expected)) {
    actual <- expected; actual[[field]] <- "different"
    expect_error(.cast_check_compare_worker(expected, actual, 2L),
      if (field == "code") "core scenario code" else field)
  }
})

test_that("installed workers ignore startup profiles and pin the loaded installation", {
  skip_if_not(file.exists(file.path(getNamespaceInfo(asNamespace("castSDM"), "path"), "Meta", "package.rds")),
    "Installed namespace bootstrap is exercised in R CMD check")
  root <- tempfile("worker-startup-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  profile <- file.path(root, "profile.R")
  writeLines("options(castSDM.worker_profile_marker = TRUE)", profile)
  withr::local_envvar(R_PROFILE_USER = profile)
  # The worker's startup library deliberately resolves castSDM to a different
  # installation. Namespace-carrying source references used to load it before
  # the initializer could set .libPaths(). A malformed dummy cannot test this.
  installed <- getNamespaceInfo(asNamespace("castSDM"), "path")
  expect_true(file.copy(installed, root, recursive = TRUE))
  previous <- .libPaths(); on.exit(.libPaths(previous), add = TRUE)
  withr::local_envvar(R_LIBS = paste(c(root, previous), collapse = .Platform$path.sep),
    R_LIBS_USER = root)
  .libPaths(c(root, previous))
  cl <- .cast_start_compare_cluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  probe <- function() list(profile = getOption("castSDM.worker_profile_marker"),
    path = normalizePath(getNamespaceInfo("castSDM", "path"), winslash = "/"))
  environment(probe) <- baseenv()
  actual <- parallel::clusterCall(cl, probe)
  expected_path <- normalizePath(getNamespaceInfo("castSDM", "path"), winslash = "/")
  expect_true(all(vapply(actual, function(x) is.null(x$profile), logical(1))))
  expect_true(all(vapply(actual, function(x) identical(x$path, expected_path), logical(1))))
})

test_that("additive GAM term reuse equals direct clamped calibrated prediction", {
  skip_if_not_installed("mgcv")
  set.seed(276)
  d <- data.frame(lon = runif(140), lat = runif(140), x = rnorm(140), z = rnorm(140))
  d$presence <- rbinom(nrow(d), 1, plogis(d$x - d$z))
  fit <- cast_fit(d, models = "gam", prevalence_target = .4, verbose = FALSE)
  a <- d[1:17, c("x", "z")]; b <- a; b$x <- b$x + 4; b$z <- b$z - .3
  for (clamp in c(FALSE, TRUE)) {
    ctx <- .cast_prediction_context(fit, c(gam = 1), a, b, clamp)
    for (switched in list(character(), "x", "z", c("x", "z"))) {
      x <- a; x[switched] <- b[switched]
      if (clamp) x <- .cast_clamp(x, fit$scaling$reference)
      direct <- predict_single_model(fit$models$gam, x)
      expect_equal(as.numeric(ctx$predict(switched)), as.numeric(direct), tolerance = 1e-12)
    }
  }
  expect_named(.cast_gam_terms(fit$models$gam), c("x", "z"))
})

test_that("callback fingerprints survive JIT compilation and follow captured values", {
  multiplier <- .2
  callback <- function(model, newdata) plogis(multiplier * newdata$x)
  before <- .cast_callback_identity(callback)
  for (i in 1:8) callback(NULL, data.frame(x = 1))
  expect_identical(.cast_callback_identity(callback), before)
  multiplier <- .4
  expect_false(identical(.cast_callback_identity(callback), before))
})

test_that("threshold sensitivity uses saved predictions with inclusive boundaries", {
  a <- data.frame(x = c(-1, 0, 1), z = 0)
  fit <- cast_as_fit(list(one = 1), a,
    function(model, newdata) stats::plogis(newdata$x + newdata$z))
  b <- a; b$x <- b$x + .5
  out <- cast_compare(fit, a, b, list(first = "x", second = "z"), threshold = .5,
    support_probs = NULL, verbose = FALSE)
  sensitivity <- summary(out, thresholds = c(.5, .7))
  expect_equal(sensitivity[sensitivity$threshold == .5, names(out$binary)], out$binary,
    ignore_attr = TRUE)
  expect_true(all(sensitivity$n_current_range[sensitivity$threshold == .5] == 2))
  # Corrupt the callback after comparison: sensitivity must not call it.
  fit$models$one$predict_function <- function(...) stop("no new prediction")
  expect_equal(summary(out, thresholds = .5)$n_current_range, rep(2, 8))
})

test_that("unusable OOF data never invent an ensemble threshold", {
  a <- data.frame(x = 1:12, z = 0)
  fit <- cast_as_fit(list(one = 1), a, function(model, newdata) plogis(newdata$x))
  cv <- structure(list(oof = data.frame(obs = rep(0:1, 6), HSS_one = NA_real_),
    thresholds = c(one = .5)), class = "cast_cv")
  out <- cast_compare(fit, a, a, list(first = "x", second = "z"), cv = cv,
    weights = c(one = 1), verbose = FALSE)
  expect_null(out$threshold)
  expect_null(out$binary)
  expect_identical(out$threshold_source, "not_supplied")
})

test_that("projected thinning preserves the inclusive distance and newest records", {
  d <- data.frame(lon = c(0, 1000, 1100, NA), lat = 0, year = c(2000, 2001, 2010, 2020))
  out <- cast_thin(d, 1, projected = TRUE, order = "year_desc", year_col = "year", verbose = FALSE)
  expect_equal(out$lon, c(0, 1100))
  expect_equal(nrow(cast_thin(d[1:2, ], 1, projected = TRUE, seed = 1, verbose = FALSE)), 2L)
})

test_that("streaming raster centroids agree with direct weighted coordinates", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 400, ncols = 300, xmin = 0, xmax = 300, ymin = 0, ymax = 400)
  h <- r; b <- r
  terra::values(h) <- seq_len(terra::ncell(r)) / terra::ncell(r)
  terra::values(b) <- as.integer(seq_len(terra::ncell(r)) %% 3 == 0)
  keep <- which(terra::values(b, mat = FALSE) == 1)
  xy <- terra::xyFromCell(r, keep); w <- terra::values(h, mat = FALSE)[keep]
  expected <- colSums(xy * w)/sum(w)
  expect_equal(unname(unlist(.weighted_centroid_raster(h, b))), unname(expected), tolerance = 1e-12)
})

test_that("installed scenario workers equal sequential execution", {
  skip_if_not(file.exists(file.path(getNamespaceInfo(asNamespace("castSDM"), "path"), "Meta", "package.rds")),
    "PSOCK uses the installed package; exercised in R CMD check")
  a <- data.frame(x = seq(-1, 1, length.out = 12), z = rep(0:1, 6))
  coefficient <- .7
  callback <- function(model, newdata) stats::plogis(coefficient * newdata$x + newdata$z)
  fit <- cast_as_fit(list(one = 1), a, callback)
  futures <- setNames(lapply(seq_len(9L), function(i) transform(a, x = x + i*.1)),
    paste0("future_", seq_len(9L)))
  one <- cast_compare(fit, a, futures, list(first = "x", second = "z"), verbose = FALSE)
  two <- cast_compare(fit, a, futures, list(first = "x", second = "z"),
    resources = list(cores = 3, workers = 3), verbose = FALSE)
  expect_equal(two$summary, one$summary)
  scientific <- function(results) lapply(results, function(x) { x$timings <- NULL; x })
  expect_equal(scientific(two$results), scientific(one$results))
  skip_if_not_installed("terra")
  raster <- terra::rast(nrows = 3, ncols = 4, nlyrs = 2)
  terra::values(raster) <- as.matrix(a); names(raster) <- names(a)
  future_rasters <- lapply(futures, function(x) {
    r <- terra::rast(raster); terra::values(r) <- as.matrix(x); r
  })
  mask <- raster[[1]]; terra::values(mask) <- c(0, rep(1, 11))
  directory <- tempfile("compare-workers-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  expected <- cast_compare(fit, a, futures, list(first = "x", second = "z"),
    mask = c(0, rep(1, 11)), threshold = .5, verbose = FALSE)
  actual <- cast_compare(fit, raster, future_rasters, list(first = "x", second = "z"),
    mask = mask, threshold = .5, output_dir = directory,
    resources = list(cores = 3, workers = 3, memory_mb = 96), verbose = FALSE)
  expect_equal(actual$summary, expected$summary, tolerance = 1e-12)
  expect_equal(actual$binary, expected$binary)
  expect_true(all(vapply(actual$results, function(x) file.exists(x$raster_file), logical(1))))
})
