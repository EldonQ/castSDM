comparison_fit <- function() {
  reference <- data.frame(x = seq(-1, 1, length.out = 8), z = rep(c(.2, .8), 4), static = 1)
  cast_as_fit(list(a = c(.4, .2), b = c(.8, -.1)), reference,
    function(model, newdata) stats::plogis(model[1] * newdata$x + model[2] * newdata$z))
}

test_that("paired comparisons agree with hand-computed model contrasts", {
  fit <- comparison_fit()
  a <- fit$scaling$reference
  b <- a; b$x <- b$x + .2; b$z <- b$z - .1
  out <- cast_compare(fit, a, b, groups = list(climate = "x", land = "z"),
    weights = c(a = .25, b = .75), threshold = .5, verbose = FALSE)
  ca <- a; ca$x <- b$x
  lu <- a; lu$z <- b$z
  predicted <- lapply(list(a, b, ca, lu), function(x) {
    vapply(fit$models, function(m) m$predict_function(m$model, x), numeric(nrow(x)))
  })
  effects <- lapply(1:2, function(m) .cast_factorial_values(do.call(cbind,
    lapply(predicted, function(x) x[, m]))))
  expect_equal(out$results$future$data$joint_change, effects[[1]][, 3] * .25 + effects[[2]][, 3] * .75)
  expect_equal(out$results$future$data$sd_joint_change, apply(cbind(effects[[1]][, 3], effects[[2]][, 3]), 1, stats::sd))
  expect_equal(out$summary$mean, unname(colMeans(effects[[1]] * .25 + effects[[2]] * .75)))
  expect_equal(sum(out$weights), 1)
  expect_equal(out$threshold_source, "supplied")
})

test_that("single row and single model do not invent disagreement", {
  fit <- comparison_fit()
  a <- fit$scaling$reference[1, , drop = FALSE]
  b <- a; b$x <- b$x + .1
  out <- cast_compare(fit, a, b, list(first = "x", second = "z"),
    weights = c(a = 1, b = 0), verbose = FALSE)
  expect_true(all(is.na(out$summary$mean_model_sd)))
  expect_null(out$binary)
  expect_equal(nrow(out$results$future$data), 1L)
})

test_that("common domain rejects missing current and any requested future", {
  fit <- comparison_fit(); a <- fit$scaling$reference
  b <- a; b$x[2] <- NA_real_
  c <- a; c$z[3] <- Inf
  a$x[1] <- NA_real_
  out <- cast_compare(fit, a, list(b = b, c = c), list(one = "x", two = "z"), verbose = FALSE)
  expect_true(all(out$summary$n_valid == 5))
  expect_true(all(is.na(out$results$b$data$joint_change[1:3])))
  per <- cast_compare(fit, a, list(b = b, c = c), list(one = "x", two = "z"),
    domain = "scenario", verbose = FALSE)
  expect_true(all(per$summary$n_valid == 6))
})

test_that("group and static declarations are validated", {
  fit <- comparison_fit(); a <- fit$scaling$reference
  expect_error(cast_compare(fit, a, a, list(one = "x", two = "x")), "disjoint")
  b <- a; b$static <- 2
  expect_error(cast_compare(fit, a, b, list(one = "x", two = "z"), verbose = FALSE), "Static predictor")
  expect_error(cast_compare(fit, a, a["x"], list(one = "x", two = "z")), "both changing")
})

test_that("active prediction failures cannot change paired weights or domains", {
  a <- data.frame(x = -1:1, z = 0)
  fit <- cast_as_fit(list(good = 1, bad = 0), a, function(model, newdata)
    ifelse(model == 0 & newdata$x > 0, NA_real_, plogis(newdata$x)))
  expect_error(cast_compare(fit, a, a, list(one = "x", two = "z"), verbose = FALSE),
    "non-finite predictions")
  out <- cast_compare(fit, a, a, list(one = "x", two = "z"),
    weights = c(good = 1, bad = 0), verbose = FALSE)
  expect_true(all(out$summary$n_valid == 3))
})

test_that("table receipts verify models and artifact contents", {
  fit <- comparison_fit(); a <- fit$scaling$reference
  outdir <- tempfile("compare-")
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  run <- function(f = fit) cast_compare(f, a, a, list(one = "x", two = "z"),
    output_dir = outdir, verbose = FALSE)
  one <- run(); two <- run()
  expect_true(all(two$events$cache_hit))
  fit$models$a$model[1] <- 3
  three <- run(fit)
  expect_false(three$events$cache_hit[three$events$stage == "current_prediction"])
  expect_false(three$events$cache_hit[three$events$stage == "scenario"])
  cast_compare(fit, a, a, list(one = "x", two = "z"), threshold = .5,
    output_dir = outdir, verbose = FALSE)
  expect_true(file.exists(file.path(outdir, "tables", "binary.csv")))
  run()
  expect_false(file.exists(file.path(outdir, "tables", "binary.csv")))
})

test_that("raster and table paths share values and counts", {
  skip_if_not_installed("terra")
  fit <- comparison_fit(); a <- fit$scaling$reference
  b <- a; b$x <- b$x + .1; b$z <- b$z - .1
  r <- terra::rast(nrows = 2, ncols = 4, nlyrs = 3)
  terra::values(r) <- as.matrix(a); names(r) <- names(a)
  s <- terra::rast(r); terra::values(s) <- as.matrix(b)
  outdir <- tempfile("compare-raster-")
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  one <- cast_compare(fit, a, b, list(one = "x", two = "z"), threshold = .5, verbose = FALSE)
  two <- cast_compare(fit, r, s, list(one = "x", two = "z"), threshold = .5,
    output_dir = outdir, verbose = FALSE)
  expect_equal(two$summary, one$summary, tolerance = 1e-12)
  expect_equal(two$binary, one$binary)
  three <- cast_compare(fit, r, s, list(one = "x", two = "z"), threshold = .5,
    output_dir = outdir, verbose = FALSE)
  expect_true(all(three$events$cache_hit))
  artifact <- three$results$future$raster_file
  writeBin(as.raw(1:10), artifact)
  expect_error(cast_compare(fit, r, s, list(one = "x", two = "z"), threshold = .5,
    output_dir = outdir, verbose = FALSE), "Stale/unverified")
  repaired <- cast_compare(fit, r, s, list(one = "x", two = "z"), threshold = .5,
    output_dir = outdir, overwrite = TRUE, verbose = FALSE)
  expect_equal(repaired$summary, one$summary, tolerance = 1e-12)
})
