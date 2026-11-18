test_that("VIF screens actual training predictors without altering observations", {
  x <- real_case(); data <- x$data[x$partition == "train", c("presence", x$variables)]
  original <- data
  screen <- cast_select(data, env_vars = x$variables, verbose = FALSE)
  expect_identical(data, original)
  expect_true(all(screen$selected %in% x$variables))
  expect_true(length(screen$selected) > 0L)
  selected <- screen$scores[screen$scores$selected, ]
  expect_true(all(selected$final_vif <= 10))
  expect_identical(cast_select(data, env_vars = x$variables, verbose = FALSE)$selected,
    screen$selected)
})

test_that("scoped fit identities reuse unchanged fits and detect changed real references", {
  fit <- real_fit(); x <- real_case()
  old <- options(castSDM.fit_identity_cache = new.env(parent = emptyenv()))
  on.exit(options(old))
  a <- castSDM:::.cast_fit_identity(fit)
  expect_identical(castSDM:::.cast_fit_identity(fit), a)
  fit$scaling$reference <- x$data[, x$variables]
  b <- castSDM:::.cast_fit_identity(fit)
  expect_false(identical(a, b))
  expect_identical(castSDM:::.cast_fit_identity(fit), b)
})

test_that("cached actual data are verified before deserialization", {
  x <- real_case(); path <- tempfile(fileext = ".rds")
  key <- castSDM:::.cast_digest(x$data)
  castSDM:::.cast_cache_write(x$data, path, key)
  expect_identical(castSDM:::.cast_cache_read(path, key), x$data)
  expect_null(castSDM:::.cast_cache_read(path, paste0(key, "_different")))
  receipt <- readRDS(path); receipt$digest <- "damaged receipt"
  saveRDS(receipt, path)
  expect_null(castSDM:::.cast_cache_read(path, key))
})

test_that("GAM identity excludes temporary formula frames but keeps actual model values", {
  fit <- real_fit()
  frame <- new.env(parent = asNamespace("mgcv"))
  frame$execution_note <- "before save"
  model <- fit$models$gam$model
  # Change runtime metadata only; actual observations, formula and fit remain.
  environment(model$formula) <- frame
  fit$models$gam$model <- model
  key <- castSDM:::.cast_fit_identity(fit)
  frame$execution_note <- "after save"
  expect_identical(castSDM:::.cast_fit_identity(fit), key)
  restored <- unserialize(serialize(fit, NULL, version = 3L))
  expect_identical(castSDM:::.cast_fit_identity(restored), key)
  expect_identical(environment(real_fit()$models$gam$model$dinfo$gp$pf), asNamespace("mgcv"))
  actual <- real_case()$data[, fit$env_vars, drop = FALSE]
  expect_equal(castSDM:::predict_single_model(restored$models$gam, actual),
    castSDM:::predict_single_model(fit$models$gam, actual), tolerance = 1e-12)
  fit$scaling$reference <- actual
  expect_false(identical(castSDM:::.cast_fit_identity(fit), key))
})
