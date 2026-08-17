# M19 regressions: mask/stack geometry guard, all-layer NA validity,
# duplicate-free top-up, and preservation of extra occurrence columns.

.test_bg_raster <- function(na_cells = integer(0)) {
  mk_lyr <- function(v) {
    terra::setValues(
      terra::rast(nrows = 10, ncols = 10, xmin = 100, xmax = 110,
                  ymin = 30, ymax = 40),
      v
    )
  }
  x2v <- runif(100)
  x2v[na_cells] <- NA
  r <- c(mk_lyr(runif(100)), mk_lyr(x2v))
  names(r) <- c("x1", "x2")
  r
}

test_that("cast_background aborts on mask/stack geometry mismatch (M19)", {
  skip_if_not_installed("terra")
  r <- .test_bg_raster()
  mask_r <- !is.na(r[[1]])
  mask_r[mask_r == 0] <- NA
  sa <- structure(list(mask = mask_r), class = "cast_study_area")
  r_shift <- terra::shift(r, dx = 50)
  occ <- data.frame(lon = runif(5, 100, 110), lat = runif(5, 30, 40))
  expect_error(
    cast_background(occ, sa, r_shift, n_bg = 20, verbose = FALSE),
    "incompatible geometry"
  )
})

test_that("cast_background never samples cells that are NA in any layer (M19)", {
  skip_if_not_installed("terra")
  set.seed(31)
  # NA cells in x2 only: first-layer checks would miss them.
  r <- .test_bg_raster(na_cells = 1:20)
  occ <- data.frame(
    lon = 100.5 + 0:4, lat = rep(30.5, 5), source = letters[1:5]
  )
  out <- cast_background(occ, NULL, r, n_bg = 30, seed = 32, verbose = FALSE)

  expect_equal(sum(out$presence == 0), 30)
  env_cols <- c("x1", "x2")
  expect_true(all(stats::complete.cases(out[, env_cols])))

  # No background cell is sampled twice (top-up exclusion works).
  bg <- out[out$presence == 0, ]
  expect_equal(anyDuplicated(bg[, c("lon", "lat")]), 0L)

  # Presence cells are excluded from the background sample.
  pres_cells <- paste(out$lon[out$presence == 1], out$lat[out$presence == 1])
  bg_cells <- paste(bg$lon, bg$lat)
  expect_false(any(bg_cells %in% pres_cells))
})

test_that("cast_background preserves extra occurrence columns (M19)", {
  skip_if_not_installed("terra")
  set.seed(33)
  r <- .test_bg_raster()
  occ <- data.frame(
    lon = 100.5 + 0:4, lat = rep(30.5, 5), source = letters[1:5]
  )
  out <- cast_background(occ, NULL, r, n_bg = 20, seed = 34, verbose = FALSE)
  expect_true("source" %in% names(out))
  expect_identical(out$source[out$presence == 1], letters[1:5])
  expect_true(all(is.na(out$source[out$presence == 0])))
})
