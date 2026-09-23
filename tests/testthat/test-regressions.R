# Regression tests for defects fixed in 0.7.0 -------------------------------

test_that("single-predictor MaxEnt preserves default features and background augmentation", {
  skip_if_not_installed("maxnet")
  set.seed(19)
  for (n_pres in c(8L, 12L, 20L)) {
    x <- c(rnorm(n_pres, 1), rnorm(150))
    x[2] <- x[1]
    x[n_pres + 1L] <- x[3]
    dat <- data.frame(lon = seq_along(x), lat = 0,
                      presence = c(rep(1L, n_pres), rep(0L, 150)), x = x)
    fit <- cast_fit(dat, models = "maxent", verbose = FALSE)
    expect_identical(fit$env_vars, "x")
    m <- fit$models$maxent$model
    f <- maxnet::maxnet.formula(dat$presence, dat["x"])
    if (n_pres < 10L) f <- stats::update(f, ~ . + 1)
    reference <- maxnet::maxnet(dat$presence,
      data.frame(x = x, unused = 0), f = f)
    expect_equal(m$betas, reference$betas)
    expect_equal(m$alpha, reference$alpha)
    expect_equal(m$entropy, reference$entropy)
    expect_false("(Intercept)" %in% names(m$betas))
    p <- predict_single_model(fit$models$maxent, dat["x"])
    expect_true(all(is.finite(p) & p >= 0 & p <= 1))
    expect_equal(p, as.numeric(predict(reference, dat["x"], type = "logistic")))
    expect_equal(length(predict_single_model(fit$models$maxent, dat[1, "x", drop = FALSE])), 1L)
  }
})

test_that("spatial CV reports fold fitting failures instead of silently discarding them", {
  skip_if_not_installed("pROC")
  testthat::local_mocked_bindings(
    make_spatial_folds = function(...) rep(rep(1:2, each = 2L), 15L),
    cast_fit = function(...) stop("deliberate engine failure")
  )
  dat <- data.frame(lon = seq_len(60), lat = 0, presence = rep(0:1, 30),
                    x = seq_len(60), x2 = sin(seq_len(60)), x3 = cos(seq_len(60)))
  warnings <- character()
  expect_error(withCallingHandlers(
    cast_cv(dat, select_method = "full", k = 2L, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), "All spatial CV folds failed")
  expect_equal(warnings, paste0("Model fitting failed in fold ", 1:2,
                                ": deliberate engine failure"))
})

test_that("change classes mark non-binary cells as NA, never empty string", {
  ch <- .cast_change_classes(
    cur_bin = c(1, 0, NA, 1, 0),
    fut_bin = c(1, NA, 1, NA, 0)
  )
  expect_identical(ch, c("stable_present", NA, NA, NA, "stable_absent"))
  expect_false(any(ch == "", na.rm = TRUE))
})

test_that("cast_vif rejects non-finite predictors with a targeted error", {
  skip_if_not_installed("car")
  dat <- data.frame(
    a = rnorm(50), b = rnorm(50), c = rnorm(50),
    d = c(1, rep(Inf, 49))
  )
  expect_error(cast_vif(dat, threshold = 10, verbose = FALSE),
               "Non-finite")
})

test_that("cast_predict reports missing predictors by name", {
  skip_if_not_installed("ranger")
  set.seed(1)
  n <- 100
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 2, verbose = FALSE)
  bad_grid <- data.frame(x1 = rnorm(5))
  expect_error(cast_predict(fit, bad_grid), "missing fitted predictor")
})



test_that("ensemble excludes models with non-finite predictions and warns", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(3)
  n <- 120
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, plogis(rnorm(n))),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = c("rf", "gam"),
                  rf_ntree = 40, seed = 4, verbose = FALSE)
  # Sabotage the gam model so its predictions are NA.
  fit$models$gam$model <- NULL

  cv <- new_cast_cv(
    metrics = data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.8, 0.7), auc_sd = c(0.1, 0.1),
      tss_mean = c(0.5, 0.4), tss_sd = c(0.1, 0.1),
      cbi_mean = c(0.4, 0.3), cbi_sd = c(0.1, 0.1),
      n_folds = c(3, 3), n_selected_mean = c(3, 3)
    ),
    fold_metrics = data.frame(), folds = integer(n), k = 3L,
    block_method = "grid",
    thresholds = c(rf = 0.5, gam = 0.5)
  )

  expect_warning(
    expect_warning(ens <- cast_ensemble(fit, cv, dat, method = "weighted"),
                   "non-finite"),
    "no usable out-of-fold")
  expect_true(all(is.finite(ens$predictions$hss_ensemble)))
  # rf keeps the full renormalised weight.
  expect_equal(unname(ens$weights["rf"]), 1)
})

test_that("cast_effect_map drives every engine, not only the matrix-tolerant one", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("gbm")
  skip_if_not_installed("terra")
  set.seed(5)
  n <- 240
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - x2)), x1 = x1, x2 = x2
  )
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  # brt and gam predict via model.frame(), which rejects a matrix; rf does not.
  fit <- cast_fit(dat, screen = screen, models = c("rf", "brt"),
                  rf_ntree = 40, seed = 6, verbose = FALSE)

  stack <- terra::rast(nrows = 12, ncols = 12, nlyrs = 2,
                       xmin = 70, xmax = 130, ymin = 20, ymax = 50)
  names(stack) <- c("x1", "x2")
  terra::values(stack) <- cbind(rnorm(terra::ncell(stack)),
                               rnorm(terra::ncell(stack)))

  em <- cast_effect_map(fit, stack, drivers = "x1", block_rows = 5L,
                        verbose = FALSE)
  expect_setequal(names(em), c("dHSS_x1", "absdHSS_x1", "support_x1"))
  vals <- terra::values(em[["dHSS_x1"]], mat = FALSE)
  expect_true(any(is.finite(vals)))
  # x1 raises suitability, so the sign-aligned mean effect must be positive.
  expect_gt(mean(vals, na.rm = TRUE), 0)
})

test_that("CV selector arguments cannot replace the outer training data", {
  dat <- data.frame(lon = 1:12, lat = 0, presence = rep(0:1, 6), x = 1:12)
  for (arg in c("data", "response", "method")) {
    expect_error(cast_cv(dat, select_args = stats::setNames(list(dat), arg)),
                 "cannot override")
  }
  for (args in list(list(1), list(a = 1, a = 2), "bad")) {
    expect_error(cast_cv(dat, select_args = args), "uniquely named list")
  }
})

test_that("CV passes buffered training rows to the selector and keeps statuses", {
  skip_if_not_installed("pROC")
  dat <- data.frame(lon = 1:120, lat = 0, presence = rep(0:1, 60),
                    x = sin(1:120))
  folds <- rep(1:3, each = 40)
  seen <- list()
  testthat::local_mocked_bindings(
    make_spatial_folds = function(...) folds,
    cast_select = function(data, response, method, ...) {
      i <- length(seen) + 1L
      seen[[i]] <<- data
      expect_identical(method, "two_stage")
      if (i == 2L) stop("test backend failure")
      selected <- if (i == 1L) character(0) else "x"
      new_cast_select(selected, data.frame(variable = "x"), method = method,
        diagnostics = list(status = if (i == 1L) "empty_selection" else "selected"))
    },
    cast_fit = function(...) list(env_vars = "x", scaling = list(impute = c(x = 0)),
                                  models = list(gam = list(model = TRUE))),
    predict_single_model = function(info, x) plogis(x$x),
    evaluate_model_full = function(...) c(auc = .7, tss = .4, cbi = .3)
  )
  expect_warning(cv <- cast_cv(dat, k = 3, models = "gam", buffer = 1.5,
    select_method = "two_stage", verbose = FALSE), "test backend failure")
  for (i in 1:3) {
    idx <- .cast_buffer_train_idx(dat$lon, dat$lat, which(folds == i), 1.5)
    expect_identical(seen[[i]], dat[idx, , drop = FALSE])
    expect_false(any(seen[[i]]$lon %in% dat$lon[folds == i]))
  }
  expect_length(cv$screens, 3L)
  expect_length(cv$selections, 3L)
  expect_identical(cv$screens[[1]]$diagnostics$status, "empty_selection")
  expect_null(cv$screens[[2]])
  expect_identical(cv$screens[[3]]$selected, "x")
  expect_identical(cv$fold_status, c("empty_selection", "selection_error", "evaluated"))
  expect_equal(cv$selection_freq$freq, 1 / 3)
  expect_true(all(is.na(cv$oof$HSS_gam[folds != 3L])))
  expect_identical(cv$fold_metrics$fold, 3L)
})

test_that("CV all-empty errors retain the selection diagnostics", {
  skip_if_not_installed("pROC")
  dat <- data.frame(lon = 1:40, lat = 0, presence = rep(0:1, 20), x = 1:40)
  testthat::local_mocked_bindings(
    make_spatial_folds = function(...) rep(1:2, each = 20),
    cast_select = function(...) new_cast_select(character(0), data.frame(variable = "x"),
      method = "two_stage", diagnostics = list(status = "empty_selection"))
  )
  err <- tryCatch(cast_cv(dat, k = 2L, verbose = FALSE), error = identity)
  expect_s3_class(err, "cast_cv_no_evaluable_folds")
  expect_identical(err$fold_status, rep("empty_selection", 2))
  expect_length(err$screens, 2)
})

test_that("the high-level pipeline forwards selection settings", {
  dat <- data.frame(lon = 1:100, lat = 0, presence = rep(0:1, 50),
                    x = sin(1:100))
  seen_select <- NULL
  seen_cv <- NULL
  testthat::local_mocked_bindings(
    cast_select = function(data, ..., num_trees, n_perm, ncov, keep) {
      seen_select <<- list(data = data, num_trees = num_trees,
                           n_perm = n_perm, ncov = ncov, keep = keep)
      new_cast_select("x", data.frame(variable = "x"), method = "two_stage")
    },
    cast_fit = function(...) list(),
    cast_cv = function(data, ..., select_args) {
      seen_cv <<- list(data = data, args = select_args)
      NULL
    },
    cast_evaluate = function(...) list()
  )
  cast(dat, select_num_trees = 111L, select_n_perm = 7L, select_ncov = 3L,
       select_keep = "x", do_predict = FALSE, verbose = FALSE, seed = 8)
  expect_lt(nrow(seen_select$data), nrow(dat))
  expect_identical(seen_select$num_trees, 111L)
  expect_identical(seen_select$n_perm, 7L)
  expect_identical(seen_select$ncov, 3L)
  expect_identical(seen_select$keep, "x")
  expect_identical(seen_cv$args$num_trees, 111L)
  expect_identical(seen_cv$args$n_perm, 7L)
  expect_identical(seen_cv$args$ncov, 3L)
  expect_identical(seen_cv$args$keep, "x")
})

test_that("the pipeline falls back to hold-out evaluation when spatial CV fails", {
  dat <- data.frame(lon = 1:100, lat = 0, presence = rep(0:1, 50),
                    x = sin(1:100))
  evaluated <- FALSE
  local_mocked_bindings(
    cast_select = function(...) new_cast_select("x", data.frame(variable = "x")),
    cast_fit = function(...) list(),
    cast_cv = function(...) stop("deliberate CV failure"),
    cast_evaluate = function(...) { evaluated <<- TRUE; list() }
  )
  expect_warning(result <- cast(dat, do_predict = FALSE, verbose = FALSE, seed = 8),
                 "Falling back to hold-out")
  expect_null(result$cv)
  expect_true(evaluated)
})

