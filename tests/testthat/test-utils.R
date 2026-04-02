test_that("validate_species_data rejects non-data.frame", {
  expect_error(
    validate_species_data(list(a = 1)),
    class = "rlang_error"
  )
})

test_that("validate_species_data detects missing columns", {
  df <- data.frame(lon = 1, lat = 2)
  expect_error(
    validate_species_data(df),
    "presence"
  )
})

test_that("validate_species_data passes with correct input", {
  df <- data.frame(lon = 1, lat = 2, presence = 1, bio1 = 15)
  expect_invisible(validate_species_data(df))
})

test_that("get_env_vars extracts correct columns", {
  df <- data.frame(
    lon = c(1, 2, 3), lat = c(2, 3, 4), presence = c(1, 0, 1),
    bio1 = c(15.0, 18.0, 12.0), bio12 = c(1200, 800, 1500),
    name = c("sp1", "sp2", "sp3")
  )
  result <- get_env_vars(df)
  expect_equal(result, c("bio1", "bio12"))
})

test_that("normalize01 handles constant input", {
  expect_equal(normalize01(rep(5, 3)), rep(0.5, 3))
})

test_that("normalize01 scales to [0, 1]", {
  result <- normalize01(c(0, 5, 10))
  expect_equal(result, c(0, 0.5, 1))
})
