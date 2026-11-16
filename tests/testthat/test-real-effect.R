test_that("shift interfaces require explicit changes and reuse paired predictions", {
  x <- real_case(); fit <- real_fit()
  current <- x$data[x$partition == "holdout", x$variables]
  future <- x$future[x$partition == "holdout", ]
  bounds <- castSDM:::.cast_support_bounds(fit$scaling$reference, c(.01, .99))
  # Select an existing pair supported in both backgrounds. Its actual bio08
  # change supplies the explicit shift; no environmental value is invented.
  supported <- castSDM:::.cast_attribution_inside(current, bounds) &
    castSDM:::.cast_attribution_inside(future, bounds) & current$bio08 != future$bio08
  i <- which(supported)[1L]
  expect_true(is.finite(i))
  a <- current[i, , drop = FALSE]; b <- a; b$bio08 <- future$bio08[i]
  shift <- future$bio08[i]-current$bio08[i]
  expect_error(cast_effect_table(fit, a, drivers = "bio08", verbose = FALSE), "explicit")
  expect_error(cast_effect_table(fit, a, drivers = c("bio02", "bio08"),
    shift = shift, verbose = FALSE), "named shifts")
  for (clamp in c(FALSE, TRUE)) {
    table <- cast_effect_table(fit, a, drivers = "bio08", shift = shift,
      clamp = clamp, verbose = FALSE)
    paired <- cast_attribution(fit, a, b, drivers = "bio08", clamp = clamp, verbose = FALSE)
    expect_equal(table$mean_signed_dHSS, paired$data$total_change, tolerance = 1e-12)
    expect_equal(table$n_supported, 1)
  }
})

test_that("unlimited exploration retains every real CV location", {
  x <- real_case(); grids <- real_grids()
  # Saved projected coordinates remain unchanged; extraction can legitimately
  # return NA where this small native crop does not contain those locations.
  info <- data.frame(x$data[x$partition == "train", c("lon", "lat")],
    fold = x$fold[x$partition == "train"])
  project <- list(run = list(cv = list(fold_info = info)))
  sample <- castSDM:::.cast_explore_sample(project, grids$current, grids$future, NULL)
  expect_identical(sample$row, seq_len(nrow(info)))
  expect_equal(nrow(sample$current), nrow(info))
  expect_equal(sample$xy, as.matrix(info[, c("lon", "lat")]))
})
