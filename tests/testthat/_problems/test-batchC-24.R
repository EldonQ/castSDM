# Extracted from test-batchC.R:24

# test -------------------------------------------------------------------------
skip_if_not_installed("terra")
skip_if_not_installed("ranger")
skip_if_not_installed("pROC")
set.seed(60)
n <- 120
occ <- data.frame(
    lon = runif(20, 100, 110), lat = runif(20, 30, 40)
  )
r <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 99, xmax = 111, ymin = 29, ymax = 41
  )
r$x1 <- terra::setValues(r, runif(100))
r$x2 <- terra::setValues(r, runif(100))
res <- cast_rep(
    occ, r, n_reps = 2, models = "rf",
    select_method = "rf", select_min_vars = 2,
    do_cv = FALSE, refit_full = FALSE,
    min_bg = 60, seed = 61, verbose = FALSE
  )
