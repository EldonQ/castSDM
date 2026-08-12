# Pipeline-stage smoke tests ------------------------------------------------

make_syn_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - 1.2 * x2)),
    x1 = x1, x2 = x2, x3 = x3
  )
}

make_syn_grid <- function(dat, n = 200) {
  grid <- data.frame(
    lon = runif(n, min(dat$lon), max(dat$lon)),
    lat = runif(n, min(dat$lat), max(dat$lat))
  )
  for (v in c("x1", "x2", "x3")) grid[[v]] <- rnorm(n)
  grid
}

test_that("cast_prepare splits spatially and reports the strategy", {
  dat <- make_syn_data()
  sp <- cast_prepare(dat, train_fraction = 0.7, seed = 5, verbose = FALSE)
  expect_equal(nrow(sp$train) + nrow(sp$test), nrow(dat))
  expect_true(nzchar(sp$split))
  expect_setequal(sp$env_vars, c("x1", "x2", "x3"))
})

test_that("cast_fit/evaluate/predict run on the rf backend", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_syn_data()
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 50,
                  seed = 6, verbose = FALSE)
  expect_s3_class(fit, "cast_fit")
  expect_true("rf" %in% names(fit$models))

  ev <- cast_evaluate(fit, dat)
  expect_true(all(c("auc_mean", "tss_mean", "cbi_mean") %in% names(ev$metrics)))

  grid <- make_syn_grid(dat)
  pred <- cast_predict(fit, grid)
  expect_true("HSS_rf" %in% names(pred$predictions))
  expect_true(all(c("mess", "extrapolating") %in% names(pred$predictions)))
})

test_that("cast_ensemble produces HSS and binary columns", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_syn_data()
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 50, seed = 7, verbose = FALSE)
  cv <- new_cast_cv(
    metrics = data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.8, 0.7), auc_sd = c(0.1, 0.1),
      tss_mean = c(0.5, 0.4), tss_sd = c(0.1, 0.1),
      cbi_mean = c(0.4, 0.3), cbi_sd = c(0.1, 0.1),
      n_folds = c(3, 3), n_selected_mean = c(3, 3)
    ),
    fold_metrics = data.frame(), folds = integer(nrow(dat)), k = 3L,
    block_method = "grid",
    thresholds = c(rf = 0.5, gam = 0.5)
  )
  ens <- cast_ensemble(fit, cv, make_syn_grid(dat), method = "weighted")
  expect_true(all(c("lon", "lat", "hss_ensemble", "binary_ensemble") %in%
                  names(ens$predictions)))
  expect_equal(sum(ens$weights), 1)
})

test_that("cast_project computes change classes and stats", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_syn_data()
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 50,
                  seed = 8, verbose = FALSE)
  cv <- new_cast_cv(
    metrics = data.frame(
      model = "rf",
      auc_mean = 0.8, auc_sd = 0.1, tss_mean = 0.5, tss_sd = 0.1,
      cbi_mean = 0.4, cbi_sd = 0.1, n_folds = 3, n_selected_mean = 3
    ),
    fold_metrics = data.frame(), folds = integer(nrow(dat)), k = 3L,
    block_method = "grid", thresholds = c(rf = 0.5)
  )
  current_grid <- make_syn_grid(dat)
  future_grid <- make_syn_grid(dat)
  future_grid$x1 <- future_grid$x1 + 2  # strong shift
  proj <- cast_project(fit, cv, current_grid,
                       future_envs = list(ssp585_2050 = future_grid))
  expect_s3_class(proj, "cast_project")
  expect_true(all(c("gain", "loss", "stable_present", "stable_absent") %in%
                  proj$changes[[1]]$change))
  expect_true(all(c("n_gain", "n_loss", "pct_change") %in% names(proj$stats)))
})

test_that("cast_esm fits bivariate sub-models for rare species", {
  set.seed(9)
  n <- 80
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = c(rep(1, 12), rep(0, n - 12)),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n)
  )
  esm <- cast_esm(dat, top_k = 4, base_algo = "glm", seed = 10,
                  verbose = FALSE)
  expect_s3_class(esm, "cast_esm")
  expect_equal(ncol(esm$pairs), choose(4, 2))
  pred <- predict_cast_esm(esm, dat)
  expect_true(all(is.finite(pred)))
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("cast_vif removes collinear columns", {
  skip_if_not_installed("car")
  set.seed(12)
  n <- 60
  a <- rnorm(n); b <- a + rnorm(n, 0, 0.001)
  dat <- data.frame(a = a, b = b, c = rnorm(n), d = rnorm(n))
  out <- cast_vif(dat, threshold = 10, verbose = FALSE)
  expect_true(length(out$selected) < ncol(dat))
  expect_false(all(c("a", "b") %in% out$selected))
})

test_that("cast() refits final models on the full data set", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_syn_data(n = 200)
  grid <- make_syn_grid(dat)
  res <- cast(dat, env_data = grid, models = "rf",
              select_method = "rf", select_min_vars = 2,
              do_cv = TRUE, cv_k = 2, refit_full = TRUE,
              seed = 45, verbose = FALSE)
  expect_s3_class(res, "cast_result")
  expect_s3_class(res$fit_full, "cast_fit")
  expect_true("HSS_rf" %in% names(res$predict$predictions))
})

test_that("cast() reports when the ensemble is skipped (no CV)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_syn_data(n = 150)
  grid <- make_syn_grid(dat)
  msg <- character(0)
  withCallingHandlers(
    res <- cast(dat, env_data = grid, models = "rf",
                select_method = "rf", select_min_vars = 2,
                do_cv = FALSE, refit_full = FALSE,
                seed = 46, verbose = TRUE),
    message = function(m) { msg <<- c(msg, conditionMessage(m));
                            invokeRestart("muffleMessage") }
  )
  expect_null(res$ensemble)
  expect_true(any(grepl("Ensemble skipped", msg)))
})
