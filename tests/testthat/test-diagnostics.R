# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
test_that("average precision handles constants, ties and known rankings", {
  expect_equal(compute_average_precision(c(1, 0, 0, 0), rep(.2, 4)), .25)
  expect_equal(compute_average_precision(c(1, 0, 1, 0), c(.9, .8, .7, .6)), 5/6)
  y <- c(1, 0, 1, 0); p <- c(.9, .9, .1, .1)
  expect_equal(compute_average_precision(y, p), .5)
  expect_identical(compute_average_precision(y, p), compute_average_precision(y[c(2,1,4,3)], p))
  expect_true(is.na(compute_average_precision(c(.5, 0, 1), c(.1, .2, .3))))
  expect_true(is.na(compute_average_precision(rep(0, 3), 1:3)))
  skip_if_not_installed("pROC")
  row <- .cast_diagnostic_metrics(rep(.2, 4), c(1,0,0,0), NA_real_, "constant", "baseline", "none")
  expect_equal(row$average_precision, row$sampled_prevalence)
  expect_true(is.na(row$pr_auc))
})

test_that("native and shifted outputs share exactly one engine prediction", {
  calls <- 0L
  fit <- cast_as_fit(list(a = 1), data.frame(x = 1:4), function(model, newdata) {
    calls <<- calls + 1L
    rep(.2, nrow(newdata))
  })
  info <- fit$models$a; info$offset <- 2
  pair <- predict_single_model(info, data.frame(x = 1:4), output = "both")
  expect_equal(calls, 1L)
  expect_equal(pair$native, rep(.2, 4))
  expect_equal(pair$score, rep(plogis(qlogis(.2) + 2), 4))
})

test_that("spatial run diagnoses a frozen ensemble without requesting maps", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(717)
  n <- 180
  dat <- data.frame(lon = runif(n, 100, 120), lat = runif(n, 20, 40),
    x = rnorm(n), z = rnorm(n))
  dat$presence <- rbinom(n, 1, plogis(dat$x))
  p <- cast_protocol(models = "rf", select_method = "full", block_method = "grid",
    cv_k = 3, rf_ntree = 25, aoa = FALSE, ensemble_min_score = 0,
    prevalence_target = .5)
  run <- cast_run(dat, protocol = p, verbose = FALSE)
  expect_null(run$ensemble)
  expect_true(run$cv$clamp)
  expect_true(run$eval$clamp)
  d <- run$diagnostics
  expect_s3_class(d, "cast_diagnostics")
  expect_true("ensemble" %in% d$holdout$model)
  expect_equal(sum(d$folds$n), nrow(run$cv$oof))
  expect_equal(d$folds$n, d$folds$n_presence + d$folds$n_background)
  expect_equal(run$cv$native_oof$obs, run$cv$oof$obs)
  for (i in 1:3) {
    idx <- run$cv$folds == i
    b <- run$cv$offsets$offset[run$cv$offsets$fold == i]
    expect_equal(.cast_apply_offset(run$cv$native_oof$HSS_rf[idx], b),
      run$cv$oof$HSS_rf[idx])
  }
  ens <- d$holdout[d$holdout$model == "ensemble", ]
  expect_equal(ens$threshold, run$cv$thresholds[["rf"]])
  obs <- run$eval$predictions$obs; pred <- run$eval$predictions$HSS_rf
  expect_equal(ens$brier, mean((pred-obs)^2))
  expect_equal(ens$omission_at_threshold, mean(pred[obs == 1] < ens$threshold))
  q <- unname(quantile(run$cv$oof$HSS_rf[run$cv$oof$obs == 1], .05, type = 1))
  expect_equal(ens$training_oof_presence05_threshold, q)
  expect_equal(ens$omission_at_presence05, mean(pred[obs == 1] < q))
  changed <- run
  changed$eval$predictions$obs <- 1L-obs
  again <- cast_diagnose(changed)
  expect_identical(again$ensemble_policy, d$ensemble_policy)
  expect_equal(again$holdout$threshold, d$holdout$threshold)
  expect_identical(again$score_views_oof, d$score_views_oof)
  expect_false(isTRUE(all.equal(again$holdout$brier, d$holdout$brier)))
  expect_true(all(c("declared", "native_same_fit") %in% d$score_views_oof$view))
})

test_that("evaluation clamps uniformly and never selects on hold-out labels", {
  skip_if_not_installed("pROC")
  ref <- data.frame(x = -1:1)
  fit <- cast_as_fit(list(a = 1), ref,
    function(model, newdata) plogis(newdata$x))
  test <- data.frame(x = c(-5, 5), presence = 0:1)
  clamped <- cast_evaluate(fit, test, thresholds = c(a = .5), clamp = TRUE)
  expect_equal(clamped$predictions$HSS_a, plogis(c(-1, 1)))
  raw <- cast_evaluate(fit, test, thresholds = c(a = .5), clamp = FALSE)
  expect_equal(raw$predictions$HSS_a, plogis(c(-5, 5)))
  expect_warning(unassessed <- cast_evaluate(fit, test), "binary metrics are unassessed")
  expect_true(is.na(unassessed$metrics$tss_mean))
  expect_identical(unassessed$threshold_source, "unavailable")
})

test_that("diagnostic ensemble does not renormalize failed model rows", {
  pred <- data.frame(HSS_a = c(.2, .8), HSS_b = c(NA, .4))
  policy <- list(weights = c(a = .25, b = .75), method = "weighted")
  expect_equal(.cast_diagnostic_combine(pred, policy), c(NA, .5))
})

test_that("repeated spatial CV records every scheduled fold in run artifacts", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(2718)
  dat <- data.frame(lon = runif(100, 90, 120), lat = runif(100, 20, 40),
    presence = rep(0:1, 50), x = rnorm(100), z = rnorm(100))
  td <- tempfile("repeat-diagnostics-")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  p <- cast_protocol(models = "rf", select_method = "full", block_method = "grid",
    cv_k = 2, n_repeat = 2, rf_ntree = 20, aoa = FALSE, ensemble_min_score = 0)
  run <- cast_run(dat, protocol = p, output_dir = td, verbose = FALSE)
  status <- read.csv(file.path(td, "tables/cv_fold_status.csv"))
  expect_equal(nrow(status), 4L)
  expect_equal(status$replicate, c(1L, 1L, 2L, 2L))
  expect_identical(status$status, run$cv$repeat_status$status)
  expect_equal(nrow(run$cv$oof), nrow(run$cv$fold_info))
})
