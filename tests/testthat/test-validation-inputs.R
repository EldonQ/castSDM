# Input validation (M5/F11), method=NULL (F14), empty consensus --------------

make_val_data <- function(n = 120, seed = 21) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, plogis(x1)),
    x1 = x1, x2 = x2, x3 = x3
  )
}

test_that("validate_species_data rejects non-binary responses, warns on NA", {
  dat <- make_val_data()
  dat$presence[1] <- 2
  expect_error(validate_species_data(dat), "binary")

  dat <- make_val_data()
  dat$presence[3] <- NA
  expect_warning(validate_species_data(dat), "missing")
})

test_that("cast_fit rejects non-numeric predictors and non-binary responses", {
  skip_if_not_installed("ranger")
  dat <- make_val_data()
  dat$x2 <- factor(dat$x2 > 0)
  # a screen can force a non-numeric column into the fit path
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  expect_error(
    cast_fit(dat, screen = screen, models = "rf", rf_ntree = 20L,
             verbose = FALSE),
    "Non-numeric"
  )

  dat <- make_val_data()
  dat$presence <- dat$presence + 1L  # values in {1, 2}
  expect_error(
    cast_fit(dat, models = "rf", rf_ntree = 20L, verbose = FALSE),
    "binary"
  )
})

test_that("cast_fit coerces a logical response instead of failing downstream", {
  skip_if_not_installed("ranger")
  dat <- make_val_data()
  dat$presence <- as.logical(dat$presence)
  fit <- cast_fit(dat, models = "rf", rf_ntree = 20L, seed = 1,
                  verbose = FALSE)
  expect_s3_class(fit, "cast_fit")
  pred <- predict_single_model(fit$models$rf,
                               as.data.frame(dat[, c("x1", "x2", "x3")]))
  expect_true(all(is.finite(pred)))
})

test_that("cast_evaluate rejects non-numeric predictors in test_data", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_val_data()
  fit <- cast_fit(dat, models = "rf", rf_ntree = 20L, seed = 1,
                  verbose = FALSE)
  test_dat <- dat
  test_dat$x3 <- as.character(test_dat$x3)
  expect_error(cast_evaluate(fit, test_dat), "Non-numeric")
})

test_that("cast_select(method = NULL) errors instead of silently using cpi", {
  dat <- make_val_data()
  expect_error(cast_select(dat, method = NULL, verbose = FALSE), "method")
})

test_that("empty consensus warns and cast_fit aborts on an empty screen", {
  skip_if_not_installed("ranger")
  dat <- make_val_data()
  cv <- new_cast_cv(
    metrics = data.frame(), fold_metrics = data.frame(),
    folds = rep(1:2, length.out = nrow(dat)), k = 2L, block_method = "grid",
    thresholds = c(rf = 0.5),
    selection_freq = data.frame(variable = c("x1", "x2"), freq = c(0.2, 0.1),
                                stringsAsFactors = FALSE)
  )
  expect_warning(cons <- cast_consensus(cv, threshold = 0.5), "Empty consensus")
  expect_length(cons$selected, 0L)
  expect_error(
    cast_fit(dat, screen = cons, models = "rf", rf_ntree = 20L,
             verbose = FALSE),
    "empty"
  )
})
