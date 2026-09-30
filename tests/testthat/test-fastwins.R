# P2 fast-win round: W9 biomod2 change coding, S2 range-change stats,
# S4 masking-rate report, S8 environmental binning methods, S9 buffer
# units, S12 out-of-raster occurrence warning.

# ---- Shared fixtures --------------------------------------------------------

.mk_fit_cv <- function(dat) {
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

.mk_grid <- function(seed = 21, n = 40) {
  set.seed(seed)
  data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = rbinom(n, 1, plogis(rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n)
  )
}

.mk_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 100, xmax = 110,
                   ymin = 30, ymax = 40)
  col_vals <- seq(0.1, 1, length.out = 10)
  row_vals <- seq(0.1, 1, length.out = 10)
  x1 <- terra::setValues(r, rep(col_vals, times = 10))
  x2 <- terra::setValues(r, rep(row_vals, each = 10))
  stk <- c(x1, x2)
  names(stk) <- c("x1", "x2")
  stk
}

# ---- W9: change coding ------------------------------------------------------

test_that(".cast_change_codes exposes both coding conventions", {
  cc <- .cast_change_codes("cast")
  bc <- .cast_change_codes("biomod2")
  expect_identical(unlist(cc, use.names = TRUE),
                   c(gain = 1L, loss = -1L, stable_present = 2L,
                     stable_absent = 0L))
  expect_identical(unlist(bc, use.names = TRUE),
                   c(gain = 1L, loss = -2L, stable_present = -1L,
                     stable_absent = 0L))
  expect_error(.cast_change_codes("nope"))
})

test_that("cast_project stats gain range columns and biomod2 percentages", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .mk_grid()
  fc <- .mk_fit_cv(dat)
  cur <- dat[, c("lon", "lat", "x1", "x2")]
  fut <- cur
  fut$x1 <- fut$x1 + 1

  p <- cast_project(fc$fit, fc$cv, cur, list(s1 = fut))
  expect_true(all(c("n_current_range", "n_future_range", "pct_loss", "pct_gain")
                  %in% names(p$stats)))
  st <- p$stats[p$stats$scenario == "s1", ]
  expect_identical(st$n_current_range, st$n_loss + st$n_stable_present)
  expect_identical(st$n_future_range, st$n_gain + st$n_stable_present)
  if (st$n_current_range > 0) {
    expect_equal(st$pct_loss,
                 round(100 * st$n_loss / st$n_current_range, 2))
  }
  if (st$n_future_range > 0) {
    expect_equal(st$pct_gain,
                 round(100 * st$n_gain / st$n_future_range, 2))
  }
})

test_that("failed scenarios carry NA in the new stats columns (S2)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .mk_grid()
  fc <- .mk_fit_cv(dat)
  cur <- dat[, c("lon", "lat", "x1", "x2")]
  fut_ok <- cur
  fut_ok$x1 <- fut_ok$x1 + 1
  fut_bad <- cur
  fut_bad$x1 <- NULL

  expect_warning(
    p <- cast_project(fc$fit, fc$cv, cur,
                      list(bad = fut_bad, good = fut_ok)),
    "failed"
  )
  bad <- p$stats[p$stats$scenario == "bad", ]
  expect_true(is.na(bad$n_current_range))
  expect_true(is.na(bad$pct_gain))
})

test_that("cast_project honours coding in the saved change GeoTIFF (W9)", {
  skip_if_not_installed("terra")
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .mk_grid()
  fc <- .mk_fit_cv(dat)
  cur <- dat[, c("lon", "lat", "x1", "x2")]
  fut <- cur
  fut$x1 <- fut$x1 + 1

  td <- tempfile("w9proj")
  p <- cast_project(fc$fit, fc$cv, cur, list(s1 = fut), save_dir = td,
                    coding = "biomod2")
  r <- terra::rast(file.path(td, "s1_change.tif"))
  v <- sort(unique(terra::values(r)[, 1]))
  v <- v[is.finite(v)]
  expect_true(all(v %in% c(1, -2, -1, 0)))
  ch <- p$changes$s1
  n_sp <- sum(ch$change == "stable_present", na.rm = TRUE)
  skip_if(n_sp == 0)
  expect_true(-1 %in% v)   # biomod2 encodes stable_present as -1
  expect_false(2 %in% v)   # cast coding would mark it as 2
  unlink(td, recursive = TRUE)
})

test_that("cast_project_raster counts classes under the chosen coding (W9/S2)", {
  skip_if_not_installed("terra")
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- .mk_grid()
  fc <- .mk_fit_cv(dat)

  mk_rast <- function(jitter = 0) {
    base <- terra::rast(nrows = 8, ncols = 8, xmin = 100, xmax = 110,
                        ymin = 30, ymax = 40)
    r <- c(terra::setValues(base, runif(64) + jitter),
           terra::setValues(base, runif(64)))
    names(r) <- c("x1", "x2")
    r
  }
  td <- tempfile("w9rast")
  out <- cast_project_raster(fc$fit, fc$cv, mk_rast(),
                             list(ok = mk_rast(0.5)),
                             output_dir = td, coding = "biomod2",
                             verbose = FALSE)
  st <- out$stats[out$stats$scenario == "ok", ]
  expect_identical(st$n_current_range, st$n_loss + st$n_stable_present)
  expect_identical(st$n_future_range, st$n_gain + st$n_stable_present)

  r <- terra::rast(file.path(td, "rasters", "ok_change_class.tif"))
  v <- sort(unique(terra::values(r)[, 1]))
  v <- v[is.finite(v)]
  expect_true(all(v %in% c(1, -2, -1, 0)))
  expect_false(2 %in% v)
  tab <- terra::freq(r)
  get_n <- function(val) {
    hit <- tab$count[tab$value == val]
    if (length(hit)) sum(hit) else 0L
  }
  expect_identical(as.integer(get_n(-1L)), as.integer(st$n_stable_present))
  unlink(td, recursive = TRUE)
})

# ---- S4: masking-rate report ------------------------------------------------

test_that("print.cast_effect_table reports the masking rate (S4)", {
  df <- data.frame(
    driver = c("a", "b"), shift_raw = 1,
    mean_abs_dHSS = c(0.2, NA_real_),
    mean_signed_dHSS = c(-0.2, NA_real_),
    n = 100L, n_supported = c(80L, 10L),
    support = c(0.8, 0.1), masked = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  attr(df, "min_support") <- 0.5
  attr(df, "intervention") <- "do(driver += shift_raw)"
  class(df) <- c("cast_effect_table", "data.frame")
  # cli_text() signals a message condition (stderr): capture that stream.
  out <- capture.output(print(df), type = "message")
  expect_true(any(grepl("Masked 1/2", out)))
  expect_true(any(grepl("0.100", out)))
})

# ---- S8: environmental binning methods --------------------------------------

test_that("environmental background supports kmeans binning (S8)", {
  skip_if_not_installed("terra")
  stk <- .mk_raster()
  set.seed(53)
  xy <- terra::xyFromCell(stk, sample(100L, 8L))
  occ <- data.frame(lon = xy[, 1], lat = xy[, 2])

  bg <- cast_background(occ, raster_stack = stk, n_bg = 30,
                        strategy = "environmental",
                        bin_method = "kmeans", bin_k = 5,
                        seed = 8, verbose = FALSE)
  bg0 <- bg[bg$presence == 0L, ]
  expect_gte(nrow(bg0), 1L)
  expect_false(anyNA(bg0$x1))

  # pca2d stays the default and honours bin_k without error.
  bg2 <- cast_background(occ, raster_stack = stk, n_bg = 30,
                         strategy = "environmental", bin_k = 4,
                         seed = 8, verbose = FALSE)
  expect_gte(sum(bg2$presence == 0L), 1L)

  # Validation.
  expect_error(
    cast_background(occ, raster_stack = stk, strategy = "environmental",
                    bin_k = 1, verbose = FALSE),
    ">= 2"
  )
  expect_error(
    cast_background(occ, raster_stack = stk, strategy = "environmental",
                    bin_method = "nope", verbose = FALSE),
    "should be one of"
  )
})

# ---- S9: buffer units -------------------------------------------------------

test_that(".cast_buffer_train_idx supports km exclusion bands (S9)", {
  # Test row 1 at (0, 0). Row 2 sits 0.5 deg of latitude away (~55.6 km,
  # KEPT at a 50 km band), row 3 3 deg of longitude (~334 km, kept),
  # row 4 ~49.7 km diagonal (excluded: closer than the band).
  lon <- c(0, 0, 3, 0.2)
  lat <- c(0, 0.5, 0, 0.4)
  idx <- .cast_buffer_train_idx(lon, lat, 1L, buffer = 50,
                                buffer_unit = "km")
  expect_true(2 %in% idx)
  expect_false(4 %in% idx)
  expect_true(3 %in% idx)
  # The degree path is unchanged: a 0.45 deg band keeps row 2 (0.5 deg)
  # and drops row 4 (0.447 deg).
  idx_deg <- .cast_buffer_train_idx(lon, lat, 1L, buffer = 0.45,
                                    buffer_unit = "deg")
  expect_true(2 %in% idx_deg)
  expect_false(4 %in% idx_deg)
})

test_that("cast_cv rejects projected coordinates with buffer_unit km (S9)", {
  set.seed(1)
  d <- data.frame(lon = runif(30, 0, 4e5), lat = runif(30, 2.9e6, 3.3e6),
                  presence = rep(c(1L, 0L), 15L), x1 = rnorm(30))
  expect_error(
    cast_cv(d, select_method = NULL, k = 2, models = "rf",
            buffer = 1000, buffer_unit = "km", verbose = FALSE),
    "decimal-degree"
  )
})

test_that("cast_cv runs with a km buffer on lon/lat data (S9)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(11)
  n <- 240
  d <- data.frame(
    lon = runif(n, 100, 110), lat = runif(n, 30, 40),
    presence = c(rep(1L, 80), rep(0L, 160)),
    x1 = c(runif(80, 0.6, 1), runif(160, 0, 1)),
    x2 = c(runif(80, 0.6, 1), runif(160, 0, 1)),
    x3 = runif(n)
  )
  cv <- cast_cv(d, models = "rf", k = 3L, select_method = "full",
                rf_ntree = 30, buffer = 30, buffer_unit = "km",
                seed = 3, verbose = FALSE)
  expect_s3_class(cv, "cast_cv")
  expect_true(nrow(cv$metrics) >= 1L)
})

# ---- S12: out-of-raster occurrence warning ----------------------------------

test_that("cell_thin warns on occurrences outside the raster (S12)", {
  skip_if_not_installed("terra")
  stk <- .mk_raster()
  xy <- terra::xyFromCell(stk, c(11L, 55L))
  occ <- rbind(
    data.frame(lon = xy[, 1], lat = xy[, 2]),
    data.frame(lon = 200, lat = 10),      # outside the raster extent
    data.frame(lon = 105.5, lat = 35.5)   # inside, kept
  )
  expect_warning(
    bg <- cast_background(occ, raster_stack = stk, n_bg = 10,
                          strategy = "random", cell_thin = TRUE,
                          seed = 9, verbose = FALSE),
    "outside"
  )
  expect_false(any(abs(bg$lon - 200) < 1e-8))
})
