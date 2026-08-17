# Replicate aggregation: per-model columns and replicate-wise filtering ------

# Minimal stand-ins for cast()/cast_background() results so the aggregation
# logic in cast_rep() can be tested deterministically and fast.
.fake_cast_factory <- function(env_data, drop_predict_at = integer(0)) {
  calls <- 0L
  function(species_data, env_data = NULL, models, ...) {
    calls <<- calls + 1L
    metrics <- data.frame(
      model = c("rf", "gam"),
      auc_mean = c(0.8, 0.7), tss_mean = c(0.5, 0.4),
      cbi_mean = c(0.4, 0.3)
    )
    screen <- list(scores = data.frame(variable = "x1"), selected = "x1")
    predict <- NULL
    if (!is.null(env_data) && !calls %in% drop_predict_at) {
      predict <- list(predictions = data.frame(
        lon = env_data$lon, lat = env_data$lat,
        HSS_rf  = rep(0.1 * calls, nrow(env_data)),
        HSS_gam = rep(0.2 * calls, nrow(env_data))
      ))
    }
    list(eval = list(metrics = metrics), screen = screen, predict = predict)
  }
}

test_that("cast_rep aggregates per model and filters replicate-wise", {
  skip_if_not_installed("terra")
  set.seed(321)
  env <- data.frame(lon = runif(8, 100, 101), lat = runif(8, 30, 31),
                    x1 = runif(8))
  occ <- data.frame(lon = runif(10, 100, 101), lat = runif(10, 30, 31))

  fake_cast <- .fake_cast_factory(env, drop_predict_at = 2L)
  res <- with_mocked_bindings(
    cast_rep(occ, raster_stack = NULL, env_data = env, n_reps = 3,
             models = c("rf", "gam"), seed = 1, verbose = FALSE),
    cast = fake_cast,
    cast_background = function(occurrences, ...) occurrences
  )
  expect_s3_class(res, "cast_rep")
  expect_equal(res$n_reps, 3)

  # Per-model aggregation columns, not a single first-model pair.
  expect_true(all(c("hss_mean_rf", "hss_sd_rf", "hss_mean_gam", "hss_sd_gam")
                  %in% names(res$prediction)))
  expect_false(any(c("hss_mean", "hss_sd") %in% names(res$prediction)))

  # Replicate 2 has no prediction surface: filtered replicate-wise, so the
  # means come from replicates 1 and 3 only (not NULL for the whole grid).
  expect_equal(res$prediction$hss_mean_rf, rep(mean(c(0.1, 0.3)), 8))
  expect_equal(res$prediction$hss_mean_gam, rep(mean(c(0.2, 0.6)), 8))
  expect_equal(res$prediction$hss_sd_rf, rep(stats::sd(c(0.1, 0.3)), 8))
})

test_that("cast_rep leaves prediction NULL without env_data", {
  skip_if_not_installed("terra")
  occ <- data.frame(lon = runif(10, 100, 101), lat = runif(10, 30, 31))
  fake_cast <- .fake_cast_factory(NULL)
  res <- with_mocked_bindings(
    cast_rep(occ, raster_stack = NULL, env_data = NULL, n_reps = 2,
             models = c("rf", "gam"), seed = 1, verbose = FALSE),
    cast = fake_cast,
    cast_background = function(occurrences, ...) occurrences
  )
  expect_null(res$prediction)
})

test_that("plot.cast_rep spans negative metrics and uses per-model columns", {
  skip_if_not_installed("ggplot2")
  set.seed(325)
  x <- structure(list(
    reps = list(),
    metrics = data.frame(
      model = "rf",
      auc_mean = 0.8, auc_sd = 0.1,
      tss_mean = -0.1, tss_sd = 0.2,   # negative TSS must stay visible
      cbi_mean = 0.4, cbi_sd = 0.1, n_reps = 3
    ),
    selection_freq = data.frame(variable = "x1", freq = 1),
    prediction = NULL,
    n_reps = 3
  ), class = "cast_rep")
  p1 <- plot.cast_rep(x)
  expect_identical(p1$coordinates$limits$y, c(-1, 1))

  # With per-model prediction columns the right panel maps the first pair.
  x$prediction <- data.frame(
    lon = 1:5, lat = 1:5,
    hss_mean_rf = runif(5), hss_sd_rf = runif(5, 0, 0.1),
    hss_mean_gam = runif(5), hss_sd_gam = runif(5, 0, 0.1)
  )
  expect_no_error(p <- plot.cast_rep(x))
  expect_true(inherits(p, c("patchwork", "ggplot")))
})
