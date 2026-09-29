# P1 review round, W4: replicated pseudo-absence sets (cast_background n_rep)
# and repeated spatial blocking (cast_cv n_repeat).

.w4_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 100, xmax = 110,
                   ymin = 30, ymax = 40)
  s <- c(terra::setValues(r, runif(100)), terra::setValues(r, runif(100)))
  names(s) <- c("x1", "x2")
  s
}

test_that("cast_background(n_rep > 1) returns independent replicate sets", {
  stk <- .w4_raster()
  set.seed(81)
  occ <- data.frame(lon = runif(12, 100, 110), lat = runif(12, 30, 40))

  one <- cast_background(occ, raster_stack = stk, n_bg = 20, seed = 4,
                         verbose = FALSE)
  expect_s3_class(one, "data.frame")   # backwards compatible single set

  reps <- cast_background(occ, raster_stack = stk, n_bg = 20, n_rep = 3,
                          seed = 4, verbose = FALSE)
  expect_s3_class(reps, "cast_background_reps")
  expect_equal(length(reps), 3L)
  expect_setequal(names(reps), c("rep1", "rep2", "rep3"))
  # Same structure, different draws.
  for (r in reps) {
    expect_s3_class(r, "data.frame")
    expect_equal(names(r), names(one))
    expect_equal(sum(r$presence == 0L), sum(one$presence == 0L))
  }
  expect_false(isTRUE(all.equal(reps$rep1, reps$rep2)))
  # Reproducible: the same seed regenerates identical replicate sets.
  again <- cast_background(occ, raster_stack = stk, n_bg = 20, n_rep = 3,
                           seed = 4, verbose = FALSE)
  expect_identical(reps, again)

  expect_error(
    cast_background(occ, raster_stack = stk, n_rep = 0, verbose = FALSE),
    "integer >= 1"
  )
})

test_that("cast_cv(n_repeat > 1) aggregates across blocking repeats", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(82)
  n <- 300
  dat <- data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = c(rep(1L, 100), rep(0L, 200)),
    x1 = c(runif(100, 0.6, 1), runif(200, 0, 1)),
    x2 = c(runif(100, 0.6, 1), runif(200, 0, 1)),
    x3 = runif(n)
  )
  cv <- cast_cv(dat, models = "rf", k = 3L, select_method = "full",
                n_repeat = 3L, seed = 9, verbose = FALSE)
  expect_equal(cv$n_repeat, 3L)
  # Fold metrics carry a replicate column spanning all repeats.
  expect_true("replicate" %in% names(cv$fold_metrics))
  expect_setequal(unique(cv$fold_metrics$replicate), 1:3)
  # Aggregated metrics exist and the between-repeat sd is reported.
  expect_true(nrow(cv$metrics) >= 1L)
  expect_true("tss_mean" %in% names(cv$metrics))
  # Downstream-compatible fields still present: first-repeat oof/thresholds,
  # k-length fold_status, pooled selection frequencies in [0, 1].
  expect_length(cv$fold_status, 3L)
  expect_false(is.null(cv$oof))
  expect_false(is.null(cv$thresholds))
  expect_true(all(cv$selection_freq$freq >= 0 & cv$selection_freq$freq <= 1))
})

test_that("cast_cv(n_repeat = 1) keeps the legacy fold_metrics format", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(83)
  n <- 240
  dat <- data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = c(rep(1L, 80), rep(0L, 160)),
    x1 = c(runif(80, 0.6, 1), runif(160, 0, 1)),
    x2 = c(runif(80, 0.6, 1), runif(160, 0, 1)),
    x3 = runif(n)
  )
  cv <- cast_cv(dat, models = "rf", k = 3L, select_method = "full",
                seed = 13, verbose = FALSE)
  expect_equal(cv$n_repeat, 1L)
  expect_false("repeat" %in% names(cv$fold_metrics))
  expect_true(nrow(cv$metrics) >= 1L)

  expect_error(
    cast_cv(dat, models = "rf", k = 3L, n_repeat = 0, verbose = FALSE),
    "integer >= 1"
  )
})
