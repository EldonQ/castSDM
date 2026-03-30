test_that("cast_prepare works with ovis_ammon data", {
  data(ovis_ammon, package = "castSDM")
  split <- cast_prepare(ovis_ammon, train_fraction = 0.7, seed = 42)
  expect_true(is.list(split))
  expect_true(nrow(split$train) > 0)
  expect_true(nrow(split$test) > 0)
  expect_true(length(split$env_vars) >= 11)
  expect_equal(nrow(split$train) + nrow(split$test), nrow(ovis_ammon))
})

test_that("cast_vif works with numeric env data", {
  data(ovis_ammon, package = "castSDM")
  env_cols <- c(
    "bio02", "bio15", "bio19", "elevation",
    "aridityindexthornthwaite", "maxtempcoldest",
    "tri", "topowet", "nontree", "etccdi_cwd", "landcover_igbp"
  )
  result <- cast_vif(ovis_ammon, threshold = 10, exclude = c(
    "HID", "lon", "lat", "species", "sid", "family",
    "category", "presence", "fraction"
  ))
  expect_true(is.list(result))
  expect_true(length(result$selected) > 0)
  expect_true(is.data.frame(result$vif_log))
})

test_that("example datasets load correctly", {
  data(ovis_ammon, package = "castSDM")
  data(china_env_grid, package = "castSDM")
  expect_equal(nrow(ovis_ammon), 3648)
  expect_equal(nrow(china_env_grid), 3648)
  expect_true("presence" %in% names(ovis_ammon))
  expect_true("lon" %in% names(china_env_grid))
})
