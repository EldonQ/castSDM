# Tile-based prediction: NA write-back and future plan hygiene -------------

.make_tiled_fixture <- function(seed = 301, n = 150) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 0, 1), lat = runif(n, 0, 1),
    presence = rbinom(n, 1, plogis(1.2 * x1 - x2)),
    x1 = x1, x2 = x2, x3 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
           seed = seed + 1L, verbose = FALSE)
}

.make_tiled_raster <- function(seed = 303, na_cells = integer(0)) {
  set.seed(seed)
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 1,
                   ymin = 0, ymax = 1, nlyrs = 3)
  v1 <- rnorm(400); v1[na_cells] <- NA_real_
  terra::values(r) <- cbind(v1, rnorm(400), rnorm(400))
  names(r) <- c("x1", "x2", "x3")
  r
}

test_that("cast_predict_tiled writes NA back over cells with NA covariates", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("terra")
  fit <- .make_tiled_fixture()
  na_cells <- 1:25
  r <- .make_tiled_raster(na_cells = na_cells)

  td <- tempfile("tiledna")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  out <- cast_predict_tiled(fit, r, output_dir = td, tile_size = 10L,
                            verbose = FALSE)
  vals <- terra::values(terra::rast(out$rasters[["rf"]]), mat = FALSE)

  # NA-covariate cells stay NA (no median-imputed fabricated values)...
  expect_true(all(is.na(vals[na_cells])))
  # ...while every fully-observed cell carries a finite prediction.
  expect_true(all(is.finite(vals[-na_cells])))
})
