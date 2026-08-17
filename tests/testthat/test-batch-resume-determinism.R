# M12/M14 regressions: name-derived per-species seeds make resumed batches
# identical to uninterrupted runs; unreadable result files are re-run.

test_that(".cast_species_seed is name-derived and stable", {
  expect_identical(.cast_species_seed(42, "ab"), .cast_species_seed(42, "ab"))
  expect_false(identical(.cast_species_seed(42, "ab"),
                         .cast_species_seed(42, "ba")))
  expect_null(.cast_species_seed(NULL, "ab"))
})

.test_species <- function() {
  set.seed(111)
  mk <- function(off) {
    n <- 100
    data.frame(
      lon = runif(n, 100, 110), lat = runif(n, 30, 40),
      presence = rbinom(n, 1, plogis(rnorm(n))),
      x1 = rnorm(n) + off, x2 = rnorm(n), x3 = rnorm(n)
    )
  }
  list(alpha_sp = mk(0), beta_sp = mk(1), gamma_sp = mk(2))
}

.test_batch_args <- list(
  models = "rf", do_cv = FALSE, do_predict = FALSE, do_ensemble = FALSE,
  refit_full = FALSE,
  select_method = "rf", select_min_vars = 2, select_num_trees = 40,
  min_occ = 5L, esm_min = 3L,
  parallel = FALSE, seed = 42, verbose = FALSE
)

test_that("resumed batch matches an uninterrupted run (M12)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  sp <- .test_species()

  td1 <- tempfile("b1")
  td2 <- tempfile("b2")
  # Split run: first two species, then resume the remainder.
  do.call(cast_batch, c(list(species_list = sp[1:2], output_dir = td1),
                        .test_batch_args))
  do.call(cast_batch_resume, c(list(output_dir = td1, species_list = sp),
                               .test_batch_args))
  # Uninterrupted reference run.
  do.call(cast_batch, c(list(species_list = sp, output_dir = td2),
                        .test_batch_args))

  for (s in names(sp)) {
    r1 <- readRDS(file.path(td1, s, "cast_result.rds"))
    r2 <- readRDS(file.path(td2, s, "cast_result.rds"))
    expect_identical(r1$screen$selected, r2$screen$selected)
    expect_equal(r1$eval$metrics, r2$eval$metrics)
  }
  unlink(c(td1, td2), recursive = TRUE)
})

test_that("resume re-runs a species whose cast_result.rds is unreadable (M14)", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  sp <- .test_species()["alpha_sp"]

  td <- tempfile("b3")
  dir.create(file.path(td, "alpha_sp"), recursive = TRUE)
  rds <- file.path(td, "alpha_sp", "cast_result.rds")
  writeLines("truncated garbage", rds)

  expect_warning(
    do.call(cast_batch_resume,
            c(list(output_dir = td, species_list = sp), .test_batch_args)),
    "unreadable"
  )
  expect_true(inherits(readRDS(rds), "cast_result"))
  unlink(td, recursive = TRUE)
})

test_that("print.cast_batch shows the assembled metric columns (M17)", {
  b <- new_cast_batch(
    species_metrics = data.frame(
      species = "sp1", model = "rf", auc = 0.81, tss = 0.52, cbi = 0.43
    ),
    species = "sp1", models = "rf"
  )
  txt <- paste(capture.output(print(b)), collapse = "\n")
  expect_true(grepl("auc", txt, fixed = TRUE))
  expect_true(grepl("0.81", txt, fixed = TRUE))
})
