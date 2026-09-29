# P1 review round, W5: SRE (surface range envelope) and disk pseudo-absence
# strategies (Barbet-Massin et al. 2012). Both restrict the candidate pool;
# the NA top-up loop must stay inside the same pool.

.w5_test_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 100, xmax = 110,
                   ymin = 30, ymax = 40)
  col_vals <- seq(0.1, 1, length.out = 10)
  row_vals <- seq(0.1, 1, length.out = 10)
  # x1 increases with column, x2 with row (cell order is row-major).
  x1 <- terra::setValues(r, rep(col_vals, times = 10))
  x2 <- terra::setValues(r, rep(row_vals, each = 10))
  stk <- c(x1, x2)
  names(stk) <- c("x1", "x2")
  stk
}

test_that("SRE keeps background inside the presence environmental envelope", {
  stk <- .w5_test_raster()
  vals1 <- as.numeric(terra::values(stk[[1]]))
  vals2 <- as.numeric(terra::values(stk[[2]]))
  set.seed(51)
  pres_cells <- sample(100L, 12L)
  xy <- terra::xyFromCell(stk, pres_cells)
  occ <- data.frame(lon = xy[, 1], lat = xy[, 2])

  bg <- cast_background(occ, raster_stack = stk, n_bg = 30,
                        strategy = "sre", sre_quantile = 0,
                        seed = 5, verbose = FALSE)
  bg0 <- bg[bg$presence == 0L, ]
  expect_gte(nrow(bg0), 1L)
  # Full min-max envelope (q = 0): every background point lies inside the
  # presence range of BOTH variables.
  expect_true(all(bg0$x1 >= min(vals1[pres_cells]) &
                  bg0$x1 <= max(vals1[pres_cells])))
  expect_true(all(bg0$x2 >= min(vals2[pres_cells]) &
                  bg0$x2 <= max(vals2[pres_cells])))
  # And no background point sits on a presence cell.
  bg_cells <- terra::cellFromXY(stk, as.matrix(bg0[, c("lon", "lat")]))
  expect_length(intersect(bg_cells, pres_cells), 0L)
})

test_that("disk keeps background inside the distance band", {
  stk <- .w5_test_raster()
  set.seed(52)
  pres_cells <- sample(100L, 3L)
  xy <- terra::xyFromCell(stk, pres_cells)
  occ <- data.frame(lon = xy[, 1], lat = xy[, 2])

  bg <- cast_background(occ, raster_stack = stk, n_bg = 25,
                        strategy = "disk", disk_min = 1.5, disk_max = 2.5,
                        seed = 6, verbose = FALSE)
  bg0 <- bg[bg$presence == 0L, ]
  expect_gte(nrow(bg0), 1L)
  dmin <- vapply(seq_len(nrow(bg0)), function(i) {
    min(sqrt((bg0$lon[i] - occ$lon)^2 + (bg0$lat[i] - occ$lat)^2))
  }, numeric(1))
  expect_true(all(dmin >= 1.5 - 1e-9 & dmin <= 2.5 + 1e-9))
})

test_that("disk defaults resolve; invalid bands abort", {
  stk <- .w5_test_raster()
  xy <- terra::xyFromCell(stk, c(11L, 88L))
  occ <- data.frame(lon = xy[, 1], lat = xy[, 2])
  # Defaults (disk_min = 0, disk_max = half the shortest extent side = 5).
  bg <- cast_background(occ, raster_stack = stk, n_bg = 20,
                        strategy = "disk", seed = 7, verbose = FALSE)
  expect_gte(sum(bg$presence == 0L), 1L)

  expect_error(
    cast_background(occ, raster_stack = stk, strategy = "disk",
                    disk_min = 3, disk_max = 1, verbose = FALSE),
    "disk_min < disk_max"
  )
  expect_error(
    cast_background(occ, raster_stack = stk, strategy = "sre",
                    sre_quantile = 0.6, verbose = FALSE),
    "\\[0, 0.5\\)"
  )
  # Band beyond the raster extent: empty candidate pool is an error, not a
  # silent fall-back to random.
  expect_error(
    cast_background(occ, raster_stack = stk, n_bg = 10, strategy = "disk",
                    disk_min = 100, disk_max = 200, verbose = FALSE),
    "No candidate background cells"
  )
})

test_that("bias_raster still combines with random only", {
  stk <- .w5_test_raster()
  xy <- terra::xyFromCell(stk, c(11L, 88L))
  occ <- data.frame(lon = xy[, 1], lat = xy[, 2])
  bias <- terra::setValues(terra::rast(stk[[1]]), runif(100))
  expect_error(
    cast_background(occ, raster_stack = stk, strategy = "sre",
                    bias_raster = bias, verbose = FALSE),
    "random"
  )
  expect_error(
    cast_background(occ, raster_stack = stk, strategy = "disk",
                    bias_raster = bias, verbose = FALSE),
    "random"
  )
})
