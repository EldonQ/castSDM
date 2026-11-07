# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
comparison_fit <- function() {
  reference <- data.frame(x = seq(-1, 1, length.out = 8), z = rep(c(.2, .8), 4), static = 1)
  cast_as_fit(list(a = c(.4, .2), b = c(.8, -.1)), reference,
    function(model, newdata) stats::plogis(model[1] * newdata$x + model[2] * newdata$z))
}

test_that("unchanged corners reuse exact predictions at row level", {
  n_predicted <- 0L
  a <- data.frame(x = 1:4, z = 1:4)
  fit <- cast_as_fit(list(a = 1), a, function(model, newdata) {
    n_predicted <<- n_predicted + nrow(newdata)
    plogis(.1*newdata$x + .2*newdata$z + .03*newdata$x*newdata$z)
  })
  b <- a; b$x[c(2, 4)] <- b$x[c(2, 4)] + .1
  b$z[c(3, 4)] <- b$z[c(3, 4)] + .1
  base <- matrix(plogis(.1*a$x + .2*a$z + .03*a$x*a$z), ncol = 1,
    dimnames = list(NULL, "a"))
  block <- .cast_compare_block(fit, a, b, list(one = "x", two = "z"),
    character(), c(a = 1), base, NULL, FALSE)
  expect_equal(n_predicted, 5L) # joint: 3; partials: one row each
  c <- a; c$x <- b$x; d <- a; d$z <- b$z
  expected <- vapply(list(a, b, c, d), function(x)
    plogis(.1*x$x + .2*x$z + .03*x$x*x$z), numeric(4))
  expect_equal(unname(block$faces), unname(expected))
})

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
  expect_true(all(is.na(out$agreement$mean_sign_agreement)))
  expect_null(out$binary)
  expect_equal(nrow(out$results$future$data), 1L)
})

test_that("common support uses all futures and shares the same baseline counts", {
  fit <- comparison_fit(); a <- fit$scaling$reference
  b <- c <- a
  b$x[2] <- 5; c$x[3] <- -5
  out <- cast_compare(fit, a, list(b = b, c = c), list(one = "x", two = "z"),
    support_probs = c(0, 1), threshold = .5, thresholds = c(.3, .7), verbose = FALSE)
  expect_equal(out$summary$n_supported, rep(7, 16))
  expect_equal(out$summary_common_support$n_valid, rep(6, 16))
  expect_equal(out$results$b$data$common_support, c(1, 0, 0, 1, 1, 1, 1, 1))
  expect_equal(out$binary_common_support$n_current_range[1:8],
    out$binary_common_support$n_current_range[9:16])
  sensitivity <- summary(out, thresholds = c(.3, .5, .7))
  expect_setequal(sensitivity$domain, c("full", "supported", "common_supported"))
  idx <- sensitivity$scenario == "b" & sensitivity$threshold == .3 &
    sensitivity$domain == "common_supported" & sensitivity$component == "joint_change"
  keep <- out$results$b$data$common_support == 1
  expect_equal(sensitivity$estimate[idx], sum(out$results$b$data$joint[keep] >= .3) -
    sum(out$results$b$data$current[keep] >= .3))
})

test_that("paired opposing responses report model sign conflict", {
  a <- data.frame(x = seq(-1, 1, length.out = 6), z = 0)
  fit <- cast_as_fit(list(up = 1, down = -1), a,
    function(model, newdata) plogis(model * newdata$x))
  b <- a; b$x <- b$x + .1
  out <- cast_compare(fit, a, b, list(one = "x", two = "z"),
    support_probs = NULL, verbose = FALSE)
  joint <- out$agreement[out$agreement$component == "joint_change", ]
  expect_equal(joint$mixed_sign_share, 1)
  expect_equal(joint$mean_sign_agreement, .5)
  expect_equal(out$results$future$data$sign_conflict_joint, rep(1, 6))
  expect_null(out$summary_common_support)
})

test_that("extra thresholds reuse predictions and persistent raster scans agree", {
  skip_if_not_installed("terra")
  fit <- comparison_fit(); a <- fit$scaling$reference
  b <- a; b$x <- b$x + .1
  r <- terra::rast(nrows = 2, ncols = 4, nlyrs = 3)
  terra::values(r) <- as.matrix(a); names(r) <- names(a)
  s <- terra::rast(r); terra::values(s) <- as.matrix(b)
  td <- tempfile("threshold-cache-")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  one <- cast_compare(fit, r, s, list(one = "x", two = "z"),
    threshold = .5, thresholds = .3, output_dir = td, verbose = FALSE)
  path <- one$results$future$raster_file
  before <- tools::md5sum(path)
  two <- cast_compare(fit, r, s, list(one = "x", two = "z"),
    threshold = .5, thresholds = c(.3, .7), output_dir = td, verbose = FALSE)
  expect_true(all(two$events$cache_hit))
  expect_identical(tools::md5sum(path), before)
  table <- cast_compare(fit, a, b, list(one = "x", two = "z"),
    threshold = .5, thresholds = c(.3, .7), verbose = FALSE)
  expect_equal(summary(two, thresholds = c(.7, .3)), summary(table, thresholds = c(.7, .3)))
  expect_equal(summary(one, thresholds = .6), summary(table, thresholds = .6))
  three <- cast_compare(fit, r, s, list(one = "x", two = "z"),
    threshold = .5, thresholds = c(.3, .7), output_dir = td, verbose = FALSE)
  expect_true(all(three$events$cache_hit))
  expect_equal(two$threshold_sensitivity, three$threshold_sensitivity)
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
