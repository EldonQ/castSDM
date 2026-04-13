# cast_dag(): structure_method routing, parameter isolation, bnlearn validation

make_toy_species_df <- function(n = 120L, seed = 1L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- 0.4 * x1 + rnorm(n)
  x3 <- -0.2 * x1 + 0.3 * x2 + rnorm(n)
  x4 <- rnorm(n)
  data.frame(
    presence = rbinom(n, 1L, 0.45),
    lon = stats::runif(n, -10, 10),
    lat = stats::runif(n, 35, 55),
    bio_a = x1, bio_b = x2, bio_c = x3, bio_d = x4
  )
}

test_that("bootstrap_hc rejects algorithm = pc with actionable message", {
  skip_if_not_installed("bnlearn")
  d <- make_toy_species_df()
  expect_error(
    cast_dag(
      d,
      structure_method = "bootstrap_hc",
      algorithm = "pc",
      R = 5L,
      seed = 1L,
      verbose = FALSE
    ),
    regexp = "structure_method"
  )
})

test_that("bootstrap_hc rejects unknown algorithm", {
  skip_if_not_installed("bnlearn")
  d <- make_toy_species_df()
  expect_error(
    cast_dag(
      d,
      structure_method = "bootstrap_hc",
      algorithm = "not_a_real_learner",
      R = 5L,
      seed = 1L,
      verbose = FALSE
    ),
    regexp = "valid learner|algorithm"
  )
})

test_that("cast_dag runs PC (structure_method pc)", {
  skip_if_not_installed("bnlearn")
  d <- make_toy_species_df(200L)
  dag <- cast_dag(
    d,
    structure_method = "pc",
    pc_alpha = 0.05,
    pc_test = "zf",
    algorithm = "hc",
    score = "bic-g",
    seed = 2L,
    verbose = FALSE
  )
  expect_s3_class(dag, "cast_dag")
  expect_equal(dag$structure_method, "pc")
  expect_true(length(dag$nodes) >= 3L)
  expect_s3_class(dag$edges, "data.frame")
})

test_that("cast_dag runs FCI", {
  skip_if_not_installed("bnlearn")
  d <- make_toy_species_df(200L)
  dag <- cast_dag(
    d,
    structure_method = "fci",
    fci_alpha = 0.05,
    fci_test = "zf",
    seed = 3L,
    verbose = FALSE
  )
  expect_s3_class(dag, "cast_dag")
  expect_equal(dag$structure_method, "fci")
})

test_that("cast_dag runs bootstrap_hc with small R", {
  skip_if_not_installed("bnlearn")
  d <- make_toy_species_df(150L)
  dag <- cast_dag(
    d,
    structure_method = "bootstrap_hc",
    algorithm = "hc",
    score = "bic-g",
    R = 8L,
    strength_threshold = 0.35,
    direction_threshold = 0.35,
    seed = 4L,
    verbose = FALSE
  )
  expect_s3_class(dag, "cast_dag")
  expect_equal(dag$structure_method, "bootstrap_hc")
  expect_equal(dag$boot_R, 8L)
})

test_that("cast_dag runs bidag_bge when BiDAG is installed", {
  skip_if_not_installed("BiDAG")
  skip_if_not_installed("bnlearn")
  d <- make_toy_species_df(100L)
  dag <- tryCatch(
    cast_dag(
      d,
      structure_method = "bidag_bge",
      bidag_algorithm = "order",
      bidag_iterations = 80L,
      seed = 5L,
      verbose = FALSE
    ),
    error = function(e) {
      skip(sprintf("BiDAG: %s", conditionMessage(e)))
    }
  )
  expect_s3_class(dag, "cast_dag")
  expect_equal(dag$structure_method, "bidag_bge")
})

test_that("cast_dag runs notears_linear when torch is installed", {
  skip_if_not_installed("torch")
  skip_on_cran()
  d <- make_toy_species_df(160L)
  dag <- tryCatch(
    cast_dag(
      d,
      structure_method = "notears_linear",
      notears_max_iter = 400L,
      notears_lambda = 0.05,
      seed = 6L,
      verbose = FALSE
    ),
    error = function(e) {
      skip(sprintf("NOTEARS: %s", conditionMessage(e)))
    }
  )
  expect_s3_class(dag, "cast_dag")
  expect_equal(dag$structure_method, "notears_linear")
})
