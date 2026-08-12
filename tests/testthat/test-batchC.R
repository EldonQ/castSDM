# Tests for the 0.7.0 replicate / uncertainty / ODMAP / min_occ features ----

test_that("cast_rep aggregates replicate metrics and selection frequency", {
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
  r$x3 <- terra::setValues(r, runif(100))

  res <- cast_rep(
    occ, r, n_reps = 2, models = "rf",
    select_method = "rf", select_min_vars = 2,
    do_cv = FALSE, refit_full = FALSE,
    min_bg = 60, seed = 61, verbose = FALSE
  )
  expect_s3_class(res, "cast_rep")
  expect_equal(res$n_reps, 2)
  expect_true(all(c("auc_mean", "auc_sd", "tss_mean", "tss_sd") %in%
                  names(res$metrics)))
  expect_true(all(c("variable", "freq") %in% names(res$selection_freq)))
  expect_true(all(res$selection_freq$freq >= 0 & res$selection_freq$freq <= 1))
})

test_that("ensemble predictions carry a cross-model SD column", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(62)
  n <- 120
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 40, seed = 63, verbose = FALSE)
  cv <- new_cast_cv(
    metrics = data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.8, 0.7), auc_sd = c(0.1, 0.1),
      tss_mean = c(0.5, 0.4), tss_sd = c(0.1, 0.1),
      cbi_mean = c(0.4, 0.3), cbi_sd = c(0.1, 0.1),
      n_folds = c(3, 3), n_selected_mean = c(3, 3)
    ),
    fold_metrics = data.frame(), folds = integer(n), k = 3L,
    block_method = "grid", thresholds = c(rf = 0.5, gam = 0.5)
  )
  ens <- cast_ensemble(fit, cv, dat, method = "weighted")
  expect_true("hss_sd" %in% names(ens$predictions))
  expect_true(all(is.finite(ens$predictions$hss_sd)))
})

test_that("counterfactual carries cross-model delta SD", {
  skip_if_not_installed("ranger")
  set.seed(64)
  n <- 150
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, plogis(rnorm(n))), x1 = rnorm(n), x2 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 40, seed = 65, verbose = FALSE)
  cf <- cast_counterfactual(fit, dat, variable = "x1", shift = 1,
                            shift_type = "sd")
  expect_true("delta_sd" %in% names(cf$predictions))
  expect_true(is.finite(cf$summary$mean_abs_delta_sd))
})

test_that("cast_report_odmap writes the five ODMAP sections", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(66)
  n <- 150
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, plogis(rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  res <- cast(dat, models = "rf", select_method = "rf", select_min_vars = 2,
              do_cv = FALSE, refit_full = FALSE, seed = 67, verbose = FALSE)
  path <- tempfile(fileext = ".md")
  out <- cast_report_odmap(res, path = path,
                           meta = list(taxon = "Testus testus"), verbose = FALSE)
  expect_identical(out, path)
  txt <- paste(readLines(path), collapse = "\n")
  for (sec in c("## 1. Overview", "## 2. Data", "## 3. Model",
                "## 4. Assessment", "## 5. Prediction")) {
    expect_true(grepl(sec, txt, fixed = TRUE), info = sec)
  }
  expect_true(grepl("Testus testus", txt, fixed = TRUE))
})

test_that("batch routes rare species to ESM and skips below esm_min", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(68)
  n <- 120
  dat <- data.frame(
    lon = runif(n, 100, 120), lat = runif(n, 28, 40),
    presence = rbinom(n, 1, plogis(rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  rare <- dat
  rare$presence <- 0
  rare$presence[seq_len(8)] <- 1          # 8 presences -> ESM route
  tiny <- rare
  tiny$presence[seq_len(6)] <- 0          # 2 presences -> skipped

  td <- tempfile("batchc")
  batch <- cast_batch(
    list(rare_sp = rare, tiny_sp = tiny, ok_sp = dat),
    models = "rf",
    do_cv = FALSE,
    select_method = "rf", select_min_vars = 2, select_num_trees = 40,
    min_occ = 20L, esm_min = 5L,
    parallel = FALSE, seed = 69, output_dir = td, verbose = FALSE
  )
  expect_true(inherits(batch$results$rare_sp, "cast_result"))
  expect_true(isTRUE(batch$results$rare_sp$esm_used))
  expect_true("esm" %in% names(batch$results$rare_sp$fit$models))
  expect_null(batch$results$tiny_sp)
  expect_true(inherits(batch$results$ok_sp, "cast_result"))
  unlink(td, recursive = TRUE)
})
