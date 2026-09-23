# Two-stage selection, effect table, and the necessity audit ----------------

# x3 is a near-copy of x1 with no effect of its own: stage 1 must drop it,
# and the necessity audit must find x1 substitutable when x3 is kept.
make_collinear_data <- function(n = 300, r = 0.98, seed = 21) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- r * x1 + sqrt(1 - r^2) * rnorm(n)
  data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - 1.2 * x2)),
    x1 = x1, x2 = x2, x3 = x3, x4 = rnorm(n), x5 = rnorm(n)
  )
}

test_that("stage 1 drops a near-duplicate predictor", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data()
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 9,
                     seed = 22, verbose = FALSE)
  expect_s3_class(scr, "cast_select")
  expect_identical(scr$method, "two_stage")
  # One of the two collinear partners must be thinned, never both kept.
  expect_true(any(scr$scores$collinear_thinned[scr$scores$variable %in% c("x1", "x3")]))
  expect_false(all(c("x1", "x3") %in% scr$diagnostics$stage1_kept))
})

test_that("stage 2 thresholds against the permutation null, not zero", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data()
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 19,
                     seed = 23, verbose = FALSE)
  expect_true(all(is.finite(scr$diagnostics$null_threshold)))
  expect_equal(scr$diagnostics$null_quantile, 0.95)
  expect_equal(scr$diagnostics$n_perm, 19L)
  # Pure noise predictors have importance above zero but not above the null.
  noise <- scr$scores[scr$scores$variable %in% c("x4", "x5"), ]
  expect_true(all(noise$p_value > 0.05 | !noise$selected))
  # The retained set must be a subset of the stage-1 survivors.
  expect_true(all(scr$selected %in% scr$diagnostics$stage1_kept))
})

test_that("stage 1 ranks a U-shaped driver above noise (poly2 signal)", {
  skip_if_not_installed("ranger")
  set.seed(31)
  n <- 400
  xu <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(-1 + 2.2 * xu^2)),
    xu = xu, noise = rnorm(n), x3 = rnorm(n), x4 = rnorm(n), x5 = rnorm(n)
  )
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 9,
                     seed = 32, verbose = FALSE)
  rnk <- stats::setNames(scr$scores$stage1_rank, scr$scores$variable)
  expect_lt(rnk[["xu"]], rnk[["noise"]])
})

test_that("ncov caps the retained set and records the reason", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 400)
  scr <- cast_select(dat, method = "two_stage", num_trees = 60, n_perm = 9,
                     ncov = 1L, seed = 33, verbose = FALSE)
  expect_lte(length(scr$selected), 1L)
  expect_true(scr$diagnostics$ncov == 1L)
  expect_true(all(scr$scores$selected_reason[scr$scores$selected] %in%
                    c("null+top-ncov", "fallback-top-ncov")))
})

test_that("method = 'full' keeps every predictor and skips both stages", {
  dat <- make_collinear_data(n = 120)
  scr <- cast_select(dat, method = "full", verbose = FALSE)
  expect_setequal(scr$selected, c("x1", "x2", "x3", "x4", "x5"))
  expect_identical(scr$method, "full")
})

test_that("retired or invalid selection requests fail without running a replacement", {
  dat <- make_collinear_data(n = 150)
  expect_error(cast_select(dat, alpha = 0.01, min_vars = 3L), "unused arguments")
  for (method in list("cpi", "dml", "rf", "typo", "two", NULL, character(),
                      NA_character_, 1, c("full", "two_stage"))) {
    expect_error(cast_select(dat, method = method), "method.*must be one of")
  }
})

test_that("cast_effect_table averages a symmetric shift set and splits magnitude from sign", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 200)
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 60,
                  seed = 25, verbose = FALSE)
  eff <- cast_effect_table(fit, newdata = dat, verbose = FALSE)
  expect_s3_class(eff, "cast_effect_table")
  expect_setequal(eff$driver, c("x1", "x2"))
  expect_true(all(c("mean_abs_dHSS", "mean_signed_dHSS") %in% names(eff)))
  # The frozen shift set is +/- 1, 2 SD: four interventions per driver.
  expect_true(all(eff$n_shifts == 4L))
  # Magnitude is non-negative and never smaller than |direction|.
  expect_true(all(eff$mean_abs_dHSS >= 0))
  expect_true(all(eff$mean_abs_dHSS >= abs(eff$mean_signed_dHSS)))
  # x1 raises suitability, x2 lowers it (the data-generating signs).
  expect_gt(eff$mean_signed_dHSS[eff$driver == "x1"], 0)
  expect_lt(eff$mean_signed_dHSS[eff$driver == "x2"], 0)
})

test_that("cast_effect_table rejects a zero shift", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 120)
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 26, verbose = FALSE)
  expect_error(
    cast_effect_table(fit, newdata = dat, shifts = c(-1, 0, 1), verbose = FALSE),
    "non-zero"
  )
})

test_that("cast_necessity scores an out-of-sample knockout cost per driver", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 300)
  nec <- cast_necessity(dat, variables = c("x1", "x2", "x4"), k = 3,
                        num_trees = 60, seed = 27, verbose = FALSE)
  expect_s3_class(nec, "cast_necessity")
  expect_setequal(nec$necessity$variable, c("x1", "x2", "x4"))
  expect_true(all(c("mean_dAUC", "pct_folds_positive", "necessary") %in%
                    names(nec$necessity)))
  expect_identical(nec$k, 3L)
  expect_identical(dim(nec$fold_dauc), c(3L, 3L))
  # The true drivers must cost more to drop than the pure noise predictor.
  d <- stats::setNames(nec$necessity$mean_dAUC, nec$necessity$variable)
  expect_gt(d[["x1"]], d[["x4"]])
})

test_that("cast_necessity honours injected folds and rejects a wrong-length vector", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 320)
  my_folds <- rep(1:4, length.out = nrow(dat))
  nec <- cast_necessity(dat, variables = c("x1", "x2", "x4"), folds = my_folds,
                        num_trees = 60, seed = 29, verbose = FALSE)
  expect_identical(nec$k, 4L)
  expect_identical(ncol(nec$fold_dauc), 4L)
  expect_identical(nec$block_method, "custom")
  expect_identical(nec$folds, as.integer(my_folds))
  expect_error(
    cast_necessity(dat, variables = c("x1", "x2"), folds = 1:10, verbose = FALSE),
    "one entry per row"
  )
})

test_that("cast_necessity reads the predictor set off a screen and needs at least two", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 200)
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  nec <- cast_necessity(dat, screen = screen, k = 2, num_trees = 40,
                        seed = 28, verbose = FALSE)
  expect_setequal(nec$necessity$variable, c("x1", "x2"))
  expect_error(
    cast_necessity(dat, variables = "x1", k = 2, verbose = FALSE),
    "at least two"
  )
})

test_that("prespecified predictors survive thinning, null threshold and cap", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 150)
  required <- c("x1", "x3")
  expect_warning(scr <- cast_select(dat, keep = required, ncov = 2L,
    num_trees = 20L, n_perm = 1L, seed = 41, verbose = FALSE), "No predictor exceeded")
  expect_identical(scr$selected, required)
  rows <- scr$scores$variable %in% required
  expect_true(all(scr$scores$kept_by_design[rows]))
  expect_false(any(scr$scores$collinear_thinned[rows]))
  expect_false(any(scr$scores$passed_null[rows]))
  expect_false(any(scr$scores$fallback[rows]))
  expect_true(all(scr$scores$selected_reason[rows] == "prespecified"))
  expect_identical(scr$diagnostics$keep, required)
  expect_true(all(required %in% scr$diagnostics$stage1_kept))
  reported <- cast_importance(scr)$effects
  expect_true(all(reported$kept_by_design[reported$variable %in% required]))
  expect_true(all(reported$selected_reason[reported$variable %in% required] == "prespecified"))
})

test_that("prespecified retention leaves only the remaining cap for optional predictors", {
  skip_if_not_installed("ranger")
  dat <- make_collinear_data(n = 150)
  scr <- suppressWarnings(cast_select(dat, keep = "x5", ncov = 2L,
    num_trees = 20L, n_perm = 1L, seed = 42, verbose = FALSE))
  expect_length(scr$selected, 2L)
  expect_true("x5" %in% scr$selected)
  expect_equal(sum(scr$scores$fallback), 1L)
  full <- cast_select(dat, method = "full", keep = "x5", verbose = FALSE)
  expect_equal(full$scores$variable[full$scores$kept_by_design], "x5")
  expect_setequal(full$selected, get_env_vars(dat))
})

test_that("invalid or unaccommodated prespecified sets fail explicitly", {
  dat <- make_collinear_data(n = 120)
  for (required in list("unknown", "lon", "presence", c("x1", "x1"), NA_character_, 1)) {
    expect_error(cast_select(dat, keep = required, verbose = FALSE), "keep")
  }
  expect_error(cast_select(dat, keep = c("x1", "x3"), ncov = 1L,
    verbose = FALSE), "increase the cap")
  dat$x5 <- 0
  expect_error(cast_select(dat, keep = "x5", verbose = FALSE), "vary in the training data")
})

test_that("nested CV retains the same prespecified variables within each training fold", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 200)
  required <- c("x1", "x3")
  cv <- suppressWarnings(cast_cv(dat, k = 2L, models = "rf", rf_ntree = 20L,
    select_args = list(keep = required, ncov = 2L, num_trees = 20L, n_perm = 1L),
    seed = 43, verbose = FALSE))
  expect_length(cv$screens, 2L)
  expect_true(all(vapply(cv$screens, function(s) {
    !is.null(s) && identical(s$diagnostics$keep, required) && all(required %in% s$selected)
  }, logical(1))))
})

test_that("the high-level pipeline forwards prespecified retention to CV", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  dat <- make_collinear_data(n = 160)
  required <- c("x1", "x3")
  result <- suppressWarnings(cast(dat, models = "rf", do_predict = FALSE,
    do_cv = TRUE, cv_k = 2L, select_num_trees = 20L, select_n_perm = 1L,
    select_ncov = 2L, select_keep = required, seed = 44, verbose = FALSE))
  expect_identical(result$screen$selected, required)
  expect_false(is.null(result$cv))
  expect_true(all(vapply(result$cv$screens, function(s) {
    !is.null(s) && identical(s$diagnostics$keep, required) && all(required %in% s$selected)
  }, logical(1))))
})

# Real tramicp result layout, including the literal Empty sentinel.
make_tramicp_result <- function(p = c(0.01, 0.8, 0.02, 0.9)) {
  sets <- list("Empty", "X1", "X2", c("X1", "X2"))
  structure(list(candidate_causal_predictors = "X1",
    set_pvals = stats::setNames(p, c("Empty", "X1", "X2", "X1+X2")),
    predictor_pvals = c(X1 = 0.02, X2 = 0.8),
    tests = lapply(seq_along(sets), function(i)
      structure(list(set = sets[[i]], test = list(p.value = p[i])),
                class = "dICPtest"))), class = "dICP")
}

make_tramicp_data <- function() {
  data.frame(presence = rep(0:1, 4), x1 = seq_len(8),
             x2 = c(2, 4, 1, 3, 8, 6, 7, 5), campaign = rep(1:2, each = 4))
}

test_that("tramicp intersects all accepted subsets and distinguishes empty statuses", {
  backend <- make_tramicp_result()
  local_mocked_bindings(.cast_run_tramicp = function(...) backend, .package = "castSDM")
  cases <- list(
    selected = c(0.01, 0.8, 0.02, 0.9),
    empty_set_accepted = c(0.8, 0.8, 0.02, 0.9),
    empty_intersection = c(0.01, 0.8, 0.8, 0.9),
    no_accepted_sets = c(0.01, 0.02, 0.03, 0.05))
  for (status in names(cases)) {
    backend <- make_tramicp_result(cases[[status]])
    scr <- cast_select(make_tramicp_data(), method = "tramicp",
                       environment = "campaign", verbose = FALSE,
                       n_perm = 0, shift_size = -1, max_rows = 0, maxncov = 0)
    expect_s3_class(scr, "cast_select")
    expect_identical(scr$method, "tramicp")
    expect_identical(scr$diagnostics$status, status)
    expect_identical(scr$selected, if (status == "selected") "x1" else character(0))
    expect_identical(scr$diagnostics$intersection, scr$selected)
    expect_equal(nrow(scr$diagnostics$tested_sets), 4L)
    expect_identical(scr$diagnostics$tested_sets$set[[1]], character(0))
    expect_identical(scr$diagnostics$backend_result, backend)
    expect_identical(names(scr$scores), c("variable", "selected", "selected_reason"))
    expect_true(all(scr$scores$selected_reason == ifelse(scr$scores$selected,
      "invariant_intersection", "not_identified")))
    expect_match(scr$diagnostics$interpretation, "neither adjustment sufficiency")
  }
})

test_that("tramicp uses strict alpha acceptance and retains full intersections", {
  backend <- make_tramicp_result(c(0.05, 0.05, 0.05, 0.051))
  local_mocked_bindings(.cast_run_tramicp = function(...) backend, .package = "castSDM")
  scr <- cast_select(make_tramicp_data(), method = "tramicp", environment = "campaign",
                     icp_alpha = 0.05, maxncov = 1L)
  expect_identical(scr$selected, c("x1", "x2"))
  expect_identical(scr$diagnostics$tested_sets$accepted, c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(scr$diagnostics$accepted_sets, list(c("x1", "x2")))
})

test_that("tramicp excludes numeric environments and preserves names and aligned rows", {
  dat <- make_tramicp_data()
  names(dat)[2:3] <- c("a + `b`", "Empty")
  dat[[2]][2] <- NA_real_
  dat$presence[3] <- NA_real_
  dat$campaign[4] <- NA_real_
  dat$lon <- NA_real_ # Coordinates do not participate in complete-row filtering.
  seen <- NULL
  local_mocked_bindings(.cast_run_tramicp = function(data, formula, icp_alpha, verbose) {
    seen <<- data
    expect_identical(all.vars(formula), c("Y", "X1", "X2"))
    expect_identical(icp_alpha, 0.1)
    expect_false(verbose)
    make_tramicp_result()
  }, .package = "castSDM")
  scr <- cast_select(dat, method = "tramicp", environment = "campaign",
                     icp_alpha = 0.1, verbose = FALSE)
  rows <- c(1L, 5L, 6L, 7L, 8L)
  expect_identical(names(seen), c("X1", "X2", "Y", "E"))
  expect_identical(seen$X1, dat[[2]][rows])
  expect_identical(seen$Y, as.numeric(dat$presence[rows]))
  expect_identical(seen$E, factor(dat$campaign[rows]))
  expect_identical(scr$diagnostics$row_ids, rows)
  expect_identical(scr$diagnostics$excluded_row_ids, 2:4)
  expect_identical(scr$diagnostics$environment_counts, table(seen$E, dnn = "campaign"))
  expect_identical(scr$selected, "a + `b`")
  expect_identical(scr$scores$variable, c("a + `b`", "Empty"))
  expect_identical(scr$diagnostics$tested_sets$set[[4]], c("a + `b`", "Empty"))
})

test_that("tramicp requires an explicit valid environment column", {
  dat <- make_tramicp_data()
  for (env in list(NULL, NA_character_, "missing", "presence", 1, dat$campaign,
                   c("campaign", "x1"))) {
    expect_error(cast_select(dat, method = "tramicp", environment = env), "environment")
  }
  for (method in c("full", "two_stage")) {
    expect_error(cast_select(dat, method = method, environment = "campaign"), "only supported")
  }
  dat$campaign <- factor(rep("a", 8), levels = c("a", "b"))
  expect_error(cast_select(dat, method = "tramicp", environment = "campaign"),
               "at least two environment levels")
  dat <- make_tramicp_data()
  dat$x1[5:8] <- NA_real_
  expect_error(cast_select(dat, method = "tramicp", environment = "campaign"),
               "after shared complete-row filtering")
  dat <- make_tramicp_data()
  dat$campaign[1] <- Inf
  expect_error(cast_select(dat, method = "tramicp", environment = "campaign"), "finite categorical")
})

test_that("tramicp guards computation rather than selecting top K or imposing keep", {
  dat <- make_tramicp_data()
  local_mocked_bindings(.cast_run_tramicp = function(...) stop("backend should not run"),
                        .package = "castSDM")
  expect_error(cast_select(dat, method = "tramicp", environment = "campaign",
    icp_max_predictors = 1L), "no top-K truncation")
  for (budget in list(0, -1, 1.5, Inf, NA_real_, c(1, 2), "2")) {
    expect_error(cast_select(dat, method = "tramicp", environment = "campaign",
      icp_max_predictors = budget), "icp_max_predictors")
  }
  for (alpha in list(0, 1, Inf, NA_real_, c(0.1, 0.2), "0.05")) {
    expect_error(cast_select(dat, method = "tramicp", environment = "campaign",
      icp_alpha = alpha), "icp_alpha")
  }
  expect_error(cast_select(dat, method = "tramicp", environment = "campaign",
    keep = "x1"), "intersection must not be altered")
  for (cap in list(1L, Inf)) {
    expect_error(cast_select(dat, method = "tramicp", environment = "campaign",
      ncov = cap), "intersection must not be capped")
  }
  # Even constant candidates count: they must not be silently prescreened.
  dat$x2 <- 1
  expect_error(cast_select(dat, method = "tramicp", environment = "campaign",
    icp_max_predictors = 1L), "Full subset search")
})

test_that("tramicp validates predictor types, names, response and finite data", {
  dat <- make_tramicp_data()
  for (x in list(letters[1:8], factor(rep(1:2, 4)), rep(Inf, 8))) {
    bad <- dat
    bad$x1 <- x
    expect_error(cast_select(bad, method = "tramicp", environment = "campaign"),
                 "numeric predictor|finite predictors")
  }
  for (nms in list(c("presence", "x1", "x1", "campaign"),
                   c("presence", "", "x2", "campaign"))) {
    bad <- dat
    names(bad) <- nms
    expect_error(cast_select(bad, method = "tramicp", environment = "campaign"), "column names")
  }
  bad <- dat
  bad$presence <- 2
  expect_error(cast_select(bad, method = "tramicp", environment = "campaign"), "binary 0/1")
  bad$presence <- 1
  expect_error(cast_select(bad, method = "tramicp", environment = "campaign"), "two response classes")
  bad <- dat
  bad$x1[bad$presence == 0] <- NA_real_
  expect_error(cast_select(bad, method = "tramicp", environment = "campaign"), "two response classes")
})

test_that("tramicp supports a single predictor without two-stage validation", {
  dat <- make_tramicp_data()[, c("presence", "x1", "campaign")]
  backend <- make_tramicp_result()
  backend$tests <- backend$tests[1:2]
  backend$set_pvals <- backend$set_pvals[1:2]
  local_mocked_bindings(.cast_run_tramicp = function(...) backend, .package = "castSDM")
  scr <- cast_select(dat, method = "tramicp", environment = "campaign")
  expect_identical(scr$selected, "x1")
  expect_equal(nrow(scr$diagnostics$tested_sets), 2L)
})

test_that("tramicp rejects failed and incomplete backend results without fallback", {
  backend <- make_tramicp_result()
  local_mocked_bindings(.cast_run_tramicp = function(...) backend, .package = "castSDM")
  bad_results <- list()
  for (p in list(NA_real_, NaN, Inf, -0.1, 1.1, NULL)) {
    bad <- make_tramicp_result()
    bad$tests[[2]]$test$p.value <- p
    bad_results[[length(bad_results) + 1L]] <- bad
  }
  bad <- make_tramicp_result()
  bad$tests <- bad$tests[-1]
  bad_results <- c(bad_results, list(bad))
  bad <- make_tramicp_result()
  bad$tests[[2]] <- bad$tests[[1]]
  bad_results <- c(bad_results, list(bad))
  bad <- make_tramicp_result()
  bad$tests[[2]]$set <- "unknown"
  bad_results <- c(bad_results, list(bad))
  bad <- make_tramicp_result()
  bad$set_pvals[2] <- NA_real_
  bad_results <- c(bad_results, list(bad))
  bad <- make_tramicp_result()
  bad$set_pvals[2] <- 0.3
  bad_results <- c(bad_results, list(bad))
  for (bad in bad_results) {
    backend <- bad
    expect_error(cast_select(make_tramicp_data(), method = "tramicp",
      environment = "campaign"), "failed, nonfinite or incomplete subset tests")
  }
})

test_that("tramicp backend exceptions are explicit", {
  local_mocked_bindings(.cast_run_tramicp = function(...) stop("test backend failure"),
                        .package = "castSDM")
  expect_error(cast_select(make_tramicp_data(), method = "tramicp",
    environment = "campaign"), "tramicp backend failed: test backend failure")
})

test_that("tramicp integration agrees with the direct backend, not truth recovery", {
  skip_if_not_installed("tramicp")
  set.seed(611)
  dat <- data.frame(presence = rep(0:1, 60), x1 = rnorm(120),
                    campaign = rep(c("a", "b"), each = 60))
  # Mutable backend options must not enable crossfitting or early stopping.
  withr::local_options(crossfit = TRUE, stop_if_empty_set_invariant = TRUE)
  controls <- tramicp::dicp_controls(type = "residual", test = "gcm.test",
    alpha = 0.05, residuals = getFromNamespace("residuals.binglm", "tramicp"),
    crossfit = FALSE, stop_if_empty_set_invariant = FALSE)
  direct_data <- data.frame(X1 = dat$x1, Y = as.numeric(dat$presence), E = factor(dat$campaign))
  set.seed(612)
  direct <- tramicp::glmICP(Y ~ X1, data = direct_data, env = ~ E,
    family = stats::binomial(link = "logit"), type = "residual", test = "gcm.test",
    controls = controls, alpha = 0.05, greedy = FALSE, max_size = NULL,
    mandatory = NULL, verbose = FALSE)
  scr <- cast_select(dat, method = "tramicp", environment = "campaign",
                     seed = 612, verbose = FALSE)
  expect_equal(scr$diagnostics$backend_result$set_pvals, direct$set_pvals)
  accepted <- lapply(direct$tests, function(x) {
    if (identical(x$set, "Empty")) character(0) else x$set
  })[direct$set_pvals > 0.05]
  expected <- if (length(accepted)) Reduce(intersect, accepted) else character(0)
  expect_identical(scr$selected, if (length(expected)) "x1" else character(0))
  expect_equal(nrow(scr$diagnostics$tested_sets), 2L)
  expect_false(scr$diagnostics$backend_result$controls$crossfit)
  expect_false(scr$diagnostics$backend_result$controls$stop_if_empty_set_invariant)
  expect_true(getOption("crossfit"))
  expect_true(getOption("stop_if_empty_set_invariant"))
})
