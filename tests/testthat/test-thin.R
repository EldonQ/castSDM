# P1 review round, W6: cast_thin() — greedy Haversine distance thinning of
# occurrence records with a seeded random visiting order.

test_that("cast_thin removes records closer than min_dist_km", {
  set.seed(61)
  # 0.5 degrees of latitude apart: ~55.6 km great-circle.
  occ <- data.frame(
    lon = c(10, 10, 20), lat = c(50, 50.5, -20), id = 1:3
  )
  tight <- cast_thin(occ, min_dist_km = 60, verbose = FALSE)
  expect_equal(nrow(tight), 2L)
  expect_true(3L %in% tight$id)          # the far record always survives
  expect_equal(sort(tight$id), c(1L, 3L))

  loose <- cast_thin(occ, min_dist_km = 40, verbose = FALSE)
  expect_equal(nrow(loose), 3L)          # 55.6 km apart -> all retained
})

test_that("cast_thin is deterministic under a seed", {
  set.seed(62)
  occ <- data.frame(
    lon = runif(60, 100, 105), lat = runif(60, 30, 35)
  )
  a <- cast_thin(occ, min_dist_km = 30, seed = 11, verbose = FALSE)
  b <- cast_thin(occ, min_dist_km = 30, seed = 11, verbose = FALSE)
  expect_identical(a, b)
  expect_true(nrow(a) < nrow(occ))
  # Retained records really are pairwise >= min_dist apart (Haversine).
  lon_r <- a$lon * pi / 180; lat_r <- a$lat * pi / 180
  R <- 6371.0088
  ok <- TRUE
  for (i in seq_len(nrow(a) - 1L)) {
    dlat <- lat_r[(i + 1L):nrow(a)] - lat_r[i]
    dlon <- lon_r[(i + 1L):nrow(a)] - lon_r[i]
    h <- sin(dlat / 2)^2 + cos(lat_r[i]) * cos(lat_r[(i + 1L):nrow(a)]) *
      sin(dlon / 2)^2
    d <- 2 * R * asin(pmin(1, sqrt(h)))
    ok <- ok && all(d >= 30 - 1e-6)
  }
  expect_true(ok)
})

test_that("cast_thin keeps extra columns and original row order", {
  occ <- data.frame(
    lon = c(0, 0, 10), lat = c(0, 0.01, 10),
    species = "x", date = c("2020-01-01", "2020-02-01", "2020-03-01")
  )
  out <- cast_thin(occ, min_dist_km = 5, seed = 3, verbose = FALSE)
  expect_setequal(names(out), names(occ))
  expect_true(all(out$lon <= 10 & out$lon >= 0))
  expect_equal(out$lon, sort(out$lon))   # original order preserved
})

test_that("cast_thin validates inputs and handles NA coordinates", {
  occ <- data.frame(lon = c(0, 0.5, NA), lat = c(0, 0.5, 10))
  expect_warning(
    out <- cast_thin(occ, min_dist_km = 100, seed = 1, verbose = TRUE),
    "missing coordinates"
  )
  expect_equal(nrow(out), 1L)

  expect_error(cast_thin(data.frame(lon = 1, lat = 1), min_dist_km = 0),
               "positive number")
  expect_error(cast_thin(data.frame(lon = 1, lat = 1), min_dist_km = -1),
               "positive number")
  expect_error(cast_thin(data.frame(lon = 1, lat = 1, x = 2),
                         lon_col = "xx"), "missing coordinate")
  # Projected coordinates are rejected rather than silently mis-thinned.
  proj <- data.frame(lon = c(500000, 501000), lat = c(4000000, 4001000))
  expect_error(cast_thin(proj, min_dist_km = 2), "decimal degrees")
  # All-NA coordinates cannot be thinned.
  expect_error(
    cast_thin(data.frame(lon = c(NA, NA), lat = c(NA, NA)), min_dist_km = 1),
    "no finite coordinates"
  )
})
