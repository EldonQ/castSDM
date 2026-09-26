test_that("evaluation plots retain negative, unbounded and unavailable scores", {
  skip_if_not_installed("ggplot2")
  x <- structure(list(metrics = data.frame(model = c("rf", "gam"),
    auc_mean = c(0.8, 0.7), tss_mean = c(0.2, NA_real_),
    cbi_mean = c(-0.4, NA_real_), logloss_mean = c(2.3, Inf))), class = "cast_eval")
  p <- plot(x, metrics = c("auc", "tss", "cbi", "logloss"))
  expect_equal(nrow(p$data), 8L)
  b <- ggplot2::ggplot_build(p)
  expect_true(any(b$data[[2]]$x == -0.4, na.rm = TRUE))
  expect_true(any(b$data[[2]]$x == 2.3, na.rm = TRUE))
  expect_equal(sum(b$data[[3]]$label == "n.a."), 3L)
  x$metrics[c("auc_mean", "tss_mean", "cbi_mean")] <- NA_real_
  p <- plot(x)
  expect_equal(nrow(p$data), 6L)
  expect_equal(sum(ggplot2::ggplot_build(p)$data[[3]]$label == "n.a."), 6L)
})

test_that("CV plots retain every planned fold and metric-specific availability", {
  skip_if_not_installed("ggplot2")
  x <- structure(list(k = 5L, block_method = "grid", thresholds = c(rf = 0.5, gam = 0.5),
    fold_status = c("evaluated", "evaluated", "evaluated", "single_class", "evaluated"),
    fold_metrics = data.frame(model = rep(c("rf", "gam"), each = 4),
      fold = rep(c(1L, 2L, 3L, 5L), 2),
      auc = c(0.6, 0.7, 0.8, 0.9, 0.5, 0.6, 0.7, 0.8),
      cbi = c(-0.4, NA, NA, NA, rep(NA_real_, 4)))), class = "cast_cv")
  p <- plot(x, metric = "cbi")
  expect_equal(nrow(p$data), 10L)
  expect_equal(levels(p$data$fold), as.character(1:5))
  expect_equal(sum(is.finite(p$data$cbi)), 1L)
  expect_false(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1))))
  b <- ggplot2::ggplot_build(p)
  expect_true(any(b$data[[1]]$y == -0.4, na.rm = TRUE))
  expect_equal(sum(b$data[[2]]$label == "n.a."), 9L)
  expect_match(p$scales$get_scales("x")$labels[4], "single class")
  expect_match(p$scales$get_scales("colour")$labels[1], "1/5")
  x$fold_metrics$cbi <- NA_real_
  expect_equal(sum(ggplot2::ggplot_build(plot(x, metric = "cbi"))$data[[2]]$label == "n.a."), 10L)
  x$fold_metrics <- x$fold_metrics[x$fold_metrics$model == "rf", ]
  expect_equal(nrow(plot(x)$data), 10L)
})
