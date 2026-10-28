test_that("factorial attribution returns the 2x2 contrasts", {
  cur <- terra::rast(matrix(0.2, 4, 4))
  clim <- terra::rast(matrix(0.6, 4, 4))
  lulc <- terra::rast(matrix(0.7, 4, 4))
  joint <- terra::rast(matrix(0.9, 4, 4))
  support <- terra::rast(matrix(1, 4, 4))
  sd <- terra::rast(matrix(0.01, 4, 4))

  fa <- cast_factorial_attribution(
    cur, joint, clim, lulc,
    threshold = 0.5, support = support,
    current_sd = sd, joint_sd = sd, climate_sd = sd, lulc_sd = sd,
    scenario = "test", ssp = "ssp585", period = "2071_2100")

  expect_s3_class(fa, "cast_factorial")
  expect_equal(fa$summary$mean[fa$summary$component == "climate_main"], 0.4)
  expect_equal(fa$summary$mean[fa$summary$component == "lulc_main"], 0.5)
  expect_equal(fa$summary$mean[fa$summary$component == "joint_change"], 0.7)
  expect_equal(fa$summary$mean[fa$summary$component == "interaction"], -0.2)
  expect_equal(fa$summary$mean[fa$summary$component == "climate_conditional"], 0.2)
  expect_equal(fa$summary$mean[fa$summary$component == "lulc_conditional"], 0.3)
  expect_true(all(fa$summary$support_share == 1))
  expect_true(all(is.finite(fa$summary$mean_surface_sd)))
  expect_equal(fa$binary$estimate[fa$binary$component == "joint_change"], 16)
  expect_equal(fa$binary$estimate[fa$binary$component == "interaction"], -16)
  expect_true(any(grepl("model_based_causal_scenario_attribution",
                        fa$contract$value)))
})

test_that("factorial attribution writes a contrast stack", {
  cur <- terra::rast(matrix(0.2, 4, 4))
  clim <- terra::rast(matrix(0.3, 4, 4))
  lulc <- terra::rast(matrix(0.4, 4, 4))
  joint <- terra::rast(matrix(0.5, 4, 4))
  out <- tempfile(fileext = ".tif")
  fa <- cast_factorial_attribution(
    cur, joint, clim, lulc, filename = out, overwrite = TRUE)
  expect_true(file.exists(out))
  r <- terra::rast(out)
  expect_true(all(c("climate_main", "lulc_main", "joint_change",
                    "interaction", "climate_conditional",
                    "lulc_conditional") %in% names(r)))
})

test_that("factorial attribution rejects misaligned and partial inputs", {
  cur <- terra::rast(matrix(0.2, 4, 4))
  clim <- terra::rast(matrix(0.3, 4, 4))
  lulc <- terra::rast(matrix(0.4, 4, 4))
  joint <- terra::rast(matrix(0.5, 4, 4))
  bad <- terra::rast(matrix(0.5, 3, 3))
  expect_error(cast_factorial_attribution(cur, bad, clim, lulc),
               "not aligned")
  expect_error(cast_factorial_attribution(
    cur, joint, clim, lulc, current_binary = terra::rast(matrix(1, 4, 4))),
    "all four binary")
})
