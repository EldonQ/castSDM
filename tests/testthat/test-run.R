# cast_run(): the single production pipeline ----------------------------------

run_smoke_data <- function(n = 200, seed = 41) {
  set.seed(seed)
  x1 <- rnorm(n)
  data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.4 * x1)),
    x1 = x1, x2 = rnorm(n), x3 = rnorm(n)
  )
}

run_smoke_protocol <- function(seed = 7L, ...) {
  cast_protocol(
    models = "rf", select_method = "full",
    cv_k = 3L, block_method = "grid", aoa = FALSE,
    rf_ntree = 30L, prevalence_target = NULL, tune = FALSE,
    seed = seed, ...
  )
}

test_that("cast_protocol validates its fields and prints effective values", {
  p <- cast_protocol()
  expect_s3_class(p, "cast_protocol")
  msgs <- utils::capture.output(print(p), type = "message")
  expect_true(any(grepl("effective", msgs, ignore.case = TRUE)))
  expect_error(cast_protocol(models = "nope"), "engines")
  expect_error(cast_protocol(cv_k = 1L), "cv_k")
  expect_error(cast_protocol(train_fraction = 1), "train_fraction")
  expect_error(cast_protocol(prevalence_target = 2), "prevalence_target")
  expect_error(cast_protocol(select_method = "typo"), "select_method")
  expect_error(cast_protocol(block_method = "cluster"), "block_method")
  expect_error(cast_protocol(ensemble_method = "consensus"), "ensemble_method")
  expect_error(cast_protocol(seed = "x"), "seed")
  expect_error(.cast_validate_protocol(list(a = 1)),
               "missing field")
})

test_that("cast_run follows the evaluation protocol on the df path", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- run_smoke_data()
  res <- cast_run(dat, env_data = dat[1:10, , drop = FALSE],
                  protocol = run_smoke_protocol(), verbose = FALSE)
  expect_s3_class(res, "cast_result")
  expect_s3_class(res$protocol, "cast_protocol")
  # One-shot hold-out evaluation at the out-of-fold thresholds.
  expect_identical(res$eval$threshold_source, "out_of_fold_cv")
  expect_equal(res$eval$thresholds, res$cv$thresholds)
  # Selection / CV / thresholds were estimated on the training split only.
  expect_identical(nrow(res$screen$scores) > 0, TRUE)
  expect_true(length(res$cv$thresholds) == 1L)
  # The published fit uses the full data set.
  expect_s3_class(res$fit_full, "cast_fit")
  # Ensemble threshold is optimised on the out-of-fold surface it applies to.
  expect_true(is.finite(res$ensemble$threshold))
  expect_identical(res$ensemble$method, "weighted")
  expect_true(all(c("HSS_rf") %in% names(res$predict$predictions)))
  expect_true(all(c("hss_ensemble", "binary_ensemble") %in%
                    names(res$ensemble$predictions)))
})

test_that("cast_run writes stage artifacts, receipt manifest and report", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- run_smoke_data(n = 160, seed = 42)
  td <- file.path(tempdir(), paste0("castrun-", as.integer(runif(1, 1, 1e9))))
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res <- cast_run(dat, protocol = run_smoke_protocol(seed = 8L),
                  output_dir = td, verbose = FALSE)
  expect_true(file.exists(file.path(td, "manifest.rds")))
  expect_true(file.exists(file.path(td, "odmap_report.md")))
  expect_true(file.exists(file.path(td, "tables", "cv_metrics.csv")))
  expect_true(file.exists(file.path(td, "tables", "holdout_metrics.csv")))
  expect_true(file.exists(file.path(td, "tables", "stage_receipts.csv")))
  expect_true(file.exists(file.path(td, "stages", "01_prepare.rds")))
  expect_true(file.exists(file.path(td, "stages", "06_refit.rds")))
  expect_identical(res$report_path, file.path(td, "odmap_report.md"))
  manifest <- readRDS(file.path(td, "manifest.rds"))
  expect_s3_class(manifest$protocol, "cast_protocol")
  expect_true(length(manifest$receipts) >= 6L)
  stages <- vapply(manifest$receipts, `[[`, character(1), "stage")
  expect_true(all(c("01_prepare", "02_select", "03_fit", "04_cv",
                    "05_evaluate", "06_refit") %in% stages))
  # Receipts carry content digests for every recorded file.
  for (r in manifest$receipts) {
    expect_true(nrow(r$outputs) >= 1L)
    expect_false(anyNA(r$outputs$md5))
    expect_true(all(file.exists(r$outputs$file)))
  }
  receipts_tab <- utils::read.csv(file.path(td, "tables", "stage_receipts.csv"))
  expect_true(nrow(receipts_tab) >= 6L)
})

test_that("cast_run restores the caller's RNG state", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- run_smoke_data(n = 120, seed = 43)
  set.seed(99)
  state <- get(".Random.seed", envir = globalenv())
  invisible(cast_run(dat, protocol = run_smoke_protocol(seed = 3L),
                     verbose = FALSE))
  expect_identical(get(".Random.seed", envir = globalenv()), state)
})

test_that("cast_run forwards selection metadata into the nested CV folds", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mgcv")
  dat <- run_smoke_data(n = 180, seed = 44)
  protocol <- cast_protocol(
    models = "rf", select_method = "screen", select_keep = "x2",
    cv_k = 2L, block_method = "grid", aoa = FALSE, rf_ntree = 30L,
    prevalence_target = NULL, seed = 5L
  )
  res <- cast_run(dat, protocol = protocol, verbose = FALSE)
  expect_true("x2" %in% res$screen$selected)
  expect_true(all(vapply(res$cv$screens, function(s) {
    !is.null(s) && "x2" %in% s$selected
  }, logical(1))))
})

test_that("cast_run rejects raster projections without an output_dir", {
  skip_if_not_installed("terra")
  dat <- run_smoke_data(n = 60, seed = 45)
  r <- terra::rast(nrows = 6, ncols = 6, nlyrs = 3)
  names(r) <- c("x1", "x2", "x3")
  terra::values(r) <- matrix(rnorm(108), ncol = 3)
  expect_error(
    cast_run(dat, env_raster = r, protocol = run_smoke_protocol(),
             verbose = FALSE),
    "output_dir"
  )
})
