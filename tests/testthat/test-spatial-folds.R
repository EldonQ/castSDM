# Contiguous spatial blocking (C2) and the buffer exclusion band -------------

# Same-fold rate among horizontally/vertically adjacent grid cells.
adjacent_same_fold <- function(lon, lat, k, method, seed = 1) {
  folds <- make_spatial_folds(lon, lat, k, method = method, seed = seed)
  side <- ceiling(sqrt(k * 2))
  xb <- castSDM:::.cast_bin(lon, side)
  yb <- castSDM:::.cast_bin(lat, side)
  cells <- unique(data.frame(xb = xb, yb = yb, fold = folds))
  same <- 0L; total <- 0L
  for (i in seq_len(nrow(cells))) {
    for (j in seq_len(nrow(cells))) {
      if (j <= i) next
      is_adj <- (cells$xb[i] == cells$xb[j] && abs(cells$yb[i] - cells$yb[j]) == 1) ||
        (cells$yb[i] == cells$yb[j] && abs(cells$xb[i] - cells$xb[j]) == 1)
      if (is_adj) {
        total <- total + 1L
        same <- same + (cells$fold[i] == cells$fold[j])
      }
    }
  }
  same / total
}

test_that("grid folds are spatially contiguous, unlike legacy grid_random", {
  g <- expand.grid(x = seq_len(20), y = seq_len(20))
  rate_new <- adjacent_same_fold(g$x, g$y, k = 5L, method = "grid")
  rate_old <- adjacent_same_fold(g$x, g$y, k = 5L, method = "grid_random")
  expect_gt(rate_new, rate_old)
  expect_gte(rate_new, 0.4)  # contiguous blocks share most cell borders
  expect_lte(rate_old, 0.2)  # legacy packing interleaves folds
})

test_that("grid_random preserves the legacy count-balanced packing", {
  g <- expand.grid(x = seq_len(12), y = seq_len(12))
  f1 <- make_spatial_folds(g$x, g$y, 4L, method = "grid_random")
  f2 <- make_spatial_folds(g$x, g$y, 4L, method = "grid_random")
  expect_identical(f1, f2)  # deterministic
  sizes <- as.numeric(table(f1))
  expect_equal(sum(sizes), nrow(g))
  expect_lte(max(sizes) - min(sizes), 16L)  # balanced up to one grid cell
})

test_that("buffer drops training rows inside the exclusion band", {
  lon <- c(0, 1, 2, 5, 6, 7); lat <- rep(0, 6)
  test_idx <- c(1L, 2L)  # test points at 0 and 1
  tr <- castSDM:::.cast_buffer_train_idx(lon, lat, test_idx, buffer = 1.5)
  expect_equal(tr, c(4L, 5L, 6L))  # the point at 2 is within 1.5 of point 1
  expect_equal(castSDM:::.cast_buffer_train_idx(lon, lat, test_idx, 0),
               c(3L, 4L, 5L, 6L))
})

test_that("cast_cv accepts contiguous blocking and a buffer", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(11)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - x2)),
    x1 = x1, x2 = x2, x3 = x3
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  cv <- suppressWarnings(cast_cv(
    dat, screen = screen, select_method = NULL, k = 3L, models = "rf",
    block_method = "grid", buffer = 1, rf_ntree = 50L, seed = 5,
    verbose = FALSE
  ))
  expect_s3_class(cv, "cast_cv")
  expect_identical(cv$block_method, "grid")
})

test_that("cast_cv warns about selection leakage with a fixed screen", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(12)
  n <- 150
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, plogis(x1)),
    x1 = x1, x2 = x2, x3 = x3
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  expect_warning(
    cast_cv(dat, screen = screen, select_method = NULL, k = 2L,
            models = "rf", rf_ntree = 50L, seed = 5, verbose = FALSE),
    "leak"
  )
})
