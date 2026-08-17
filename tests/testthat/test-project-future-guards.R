# M16 regressions: cast_project() grid alignment, per-scenario failure
# handling, and .save_prediction_tif() resolution inference.

.make_fit_cv <- function(dat) {
  screen <- new_cast_select(
    c("x1", "x2"), data.frame(variable = c("x1", "x2")), method = "manual"
  )
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 7, verbose = FALSE)
  cv <- new_cast_cv(
    metrics = data.frame(
      model = "rf", auc_mean = 0.8, auc_sd = 0.1, tss_mean = 0.5,
      tss_sd = 0.1, cbi_mean = 0.4, cbi_sd = 0.1, n_folds = 3,
      n_selected_mean = 2
    ),
    fold_metrics = data.frame(), folds = rep(1:3, length.out = nrow(dat)),
    k = 3L, block_method = "grid", thresholds = c(rf = 0.5)
  )
  list(fit = fit, cv = cv)
}

.make_grid <- function(seed = 21, n = 40) {
  set.seed(seed)
  data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = rbinom(n, 1, plogis(rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n)
  )
}

test_that("cast_project aborts on a row-count mismatch (M16)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .make_grid()
  fc <- .make_fit_cv(dat)
  cur <- dat[, c("lon", "lat", "x1", "x2")]
  expect_error(
    cast_project(fc$fit, fc$cv, cur, list(s1 = cur[1:20, ])),
    "one row per current-grid row"
  )
})

test_that("cast_project re-aligns a re-ordered future grid on lon/lat (M16)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .make_grid()
  fc <- .make_fit_cv(dat)
  cur <- dat[, c("lon", "lat", "x1", "x2")]
  fut <- cur
  fut$x1 <- fut$x1 + 1

  p_ref <- cast_project(fc$fit, fc$cv, cur, list(s1 = fut))
  set.seed(99)
  ord <- sample.int(nrow(fut))
  p_ord <- cast_project(fc$fit, fc$cv, cur, list(s1 = fut[ord, ]))

  expect_equal(p_ref$changes$s1, p_ord$changes$s1)
  expect_equal(p_ref$stats, p_ord$stats)
})

test_that("cast_project records an NA stats row for a failed scenario (M16)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .make_grid()
  fc <- .make_fit_cv(dat)
  cur <- dat[, c("lon", "lat", "x1", "x2")]
  fut_ok <- cur
  fut_ok$x1 <- fut_ok$x1 + 1
  fut_bad <- cur
  fut_bad$x1 <- NULL   # missing predictor -> ensemble fails for this scenario

  expect_warning(
    p <- cast_project(fc$fit, fc$cv, cur,
                      list(bad = fut_bad, good = fut_ok)),
    "failed"
  )
  expect_true("good" %in% names(p$future))
  expect_false(is.null(p$changes$good))
  bad_row <- p$stats[p$stats$scenario == "bad", ]
  expect_true(is.na(bad_row$n_gain))
  expect_true(is.na(bad_row$pct_change))
  good_row <- p$stats[p$stats$scenario == "good", ]
  expect_false(is.na(good_row$n_gain))
})

test_that(".cast_change_classes guards against length mismatch (M16)", {
  expect_error(.cast_change_classes(c(0, 1), c(0, 1, 0)),
               "differ in length")
  expect_identical(.cast_change_classes(c(0, 1, 1, 0, NA), c(1, 0, 1, 0, 0)),
                   c("gain", "loss", "stable_present", "stable_absent", NA))
})

test_that(".save_prediction_tif infers resolution from the grid (M16)", {
  skip_if_not_installed("terra")
  lons <- seq(100, 101, by = 0.1)
  lats <- seq(30, 31, by = 0.1)
  df <- data.frame(
    lon = rep(lons, times = length(lats)),
    lat = rep(lats, each = length(lons)),
    v = runif(length(lons) * length(lats))
  )
  path <- tempfile(fileext = ".tif")
  .save_prediction_tif(df, "lon", "lat", "v", path)
  expect_true(file.exists(path))
  r <- terra::rast(path)
  expect_equal(terra::res(r), c(0.1, 0.1), tolerance = 1e-6)
  expect_true(grepl("4326", terra::crs(r)))
  got <- as.numeric(terra::extract(r, cbind(df$lon[1], df$lat[1]))[1])
  expect_equal(got, df$v[1])
  unlink(path)
})

test_that("cast_project_raster skips a failed scenario and continues (M16)", {
  skip_if_not_installed("terra")
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .make_grid()
  fc <- .make_fit_cv(dat)

  mk_rast <- function(with_x2 = TRUE, jitter = 0) {
    mk_lyr <- function(v) {
      terra::setValues(
        terra::rast(nrows = 8, ncols = 8, xmin = 100, xmax = 110,
                    ymin = 30, ymax = 40),
        v
      )
    }
    r <- mk_lyr(runif(64) + jitter)
    if (with_x2) r <- c(r, mk_lyr(runif(64)))
    names(r) <- if (with_x2) c("x1", "x2") else "x1"
    r
  }
  r_cur <- mk_rast()
  r_ok  <- mk_rast(jitter = 0.5)
  r_bad <- mk_rast(with_x2 = FALSE)

  td <- tempfile("projrast")
  expect_warning(
    out <- cast_project_raster(
      fc$fit, fc$cv, r_cur,
      list(ok = r_ok, bad = r_bad),
      output_dir = td, verbose = FALSE
    ),
    "failed"
  )
  expect_true(file.exists(file.path(td, "rasters", "ok_change_class.tif")))
  st <- out$stats
  expect_true(is.na(st$n_gain[st$scenario == "bad"]))
  expect_false(is.na(st$n_gain[st$scenario == "ok"]))

  # clamp passthrough reservation: accepted either way, warns only when the
  # installed cast_ensemble_raster() lacks clamp support.
  has_clamp <- "clamp" %in% names(formals(cast_ensemble_raster))
  if (!has_clamp) {
    expect_warning(
      cast_project_raster(fc$fit, fc$cv, r_cur, list(ok = r_ok),
                          output_dir = tempfile("projrast2"),
                          clamp = TRUE, verbose = FALSE),
      "clamp"
    )
  }
  unlink(td, recursive = TRUE)
})
