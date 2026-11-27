test_that("the fixture retains real sample design and native grids", {
  x <- real_case(); grids <- real_grids()
  expect_equal(nrow(x$data), 11817L)
  expect_equal(sum(x$data$presence), 1817)
  expect_equal(sum(x$data$presence == 0), 10000)
  expect_identical(sort(x$source_row), seq_len(nrow(x$data)))
  expect_equal(sum(x$partition == "holdout"), 3542L)
  expect_true(all(is.na(x$fold[x$partition == "holdout"])))
  expect_setequal(unique(x$fold[x$partition == "train"]), 1:5)
  expect_identical(names(grids$current), x$variables)
  expect_true(terra::compareGeom(grids$current, grids$future))
  expect_equal(terra::ncell(grids$current), 25L)
})

test_that("ordinary, accelerated, tabular and scenario predictions share one clamp", {
  x <- real_case(); fit <- real_fit()
  old <- options(castSDM.predict_threads = 1L); on.exit(options(old))
  at <- x$partition == "holdout"
  current <- x$data[at, x$variables]; future <- x$future[at, ]
  weights <- setNames(rep(1/length(fit$models), length(fit$models)), names(fit$models))
  ref <- fit$scaling$reference
  outside <- Reduce(`|`, lapply(x$variables, function(v)
    future[[v]] < min(ref[[v]]) | future[[v]] > max(ref[[v]])))
  expect_true(any(outside))
  gam <- fit$models$gam$model
  selected <- gam$smooth[[1L]]$label
  regular <- castSDM:::.cast_predict_gam(gam, future, type = "terms", terms = selected)
  direct <- castSDM:::.cast_selected_gam_terms(gam, future, selected)
  expect_equal(colnames(direct), colnames(regular))
  expect_equal(as.numeric(direct), as.numeric(regular), tolerance = 1e-12)
  all_terms <- castSDM:::.cast_predict_gam(gam, future, type = "terms")
  rebuilt <- castSDM:::.cast_selected_gam_terms(gam, future)
  expect_identical(colnames(rebuilt), colnames(all_terms))
  expect_identical(attr(rebuilt, "constant"), attr(all_terms, "constant"))
  expect_equal(as.numeric(rebuilt), as.numeric(all_terms), tolerance = 1e-12)
  for (clamp in c(FALSE, TRUE)) {
    context <- castSDM:::.cast_prediction_context(fit, weights, current, future, clamp)
    scenario <- context$predict(x$variables)
    public <- cast_predict(fit, future, clamp = clamp, extrapolation = FALSE)$predictions
    matrix <- castSDM:::.cast_predict_matrix(fit, future, names(weights), clamp = clamp)
    expect_equal(unname(matrix), unname(scenario), tolerance = 1e-12)
    for (model in names(weights)) {
      ordinary <- castSDM:::predict_single_model(fit$models[[model]], future,
        clamp = clamp, reference = ref)
      accelerated <- castSDM:::.cast_response_predictor(fit$models[[model]], ref, clamp)
      expect_equal(accelerated(future), ordinary, tolerance = 1e-12)
      expect_equal(unname(scenario[, model]), ordinary, tolerance = 1e-12)
      expect_equal(public[[paste0("HSS_", model)]], ordinary, tolerance = 1e-12)
    }
  }
})

test_that("fixed profiles clamp against the full fit rather than profile medians", {
  x <- real_case(); fit <- real_fit()
  profile <- cast_response_curves(fit, variables = "bio08", type = "fixed")
  reference <- fit$scaling$reference
  for (model in names(fit$models)) {
    curve <- profile[profile$model == model, ]
    data <- reference[rep(1L, nrow(curve)), , drop = FALSE]
    for (v in names(data)) data[[v]] <- median(reference[[v]])
    data$bio08 <- curve$value
    expected <- castSDM:::predict_single_model(fit$models[[model]], data,
      clamp = TRUE, reference = reference)
    expect_equal(curve$prediction, expected, tolerance = 1e-12)
  }
})
