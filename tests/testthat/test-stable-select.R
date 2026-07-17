test_that("domains are reproducible and class-complete", {
  set.seed(10)
  n <- 320
  dat <- data.frame(
    lon = rep(1:4, each = n / 4) + rnorm(n, sd = 0.1),
    lat = rep(1:4, each = n / 4) + rnorm(n, sd = 0.1),
    presence = rep(rep(c(0, 1), each = n / 8), 4),
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  a <- cast_domains(dat, k = 4, seed = 5, min_n = 30, min_class = 5)
  b <- cast_domains(dat, k = 4, seed = 5, min_n = 30, min_class = 5)
  expect_identical(a, b)
  expect_gte(nlevels(a), 2)
  expect_equal(length(a), n)
})

test_that("stable selector returns auditable diagnostics", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ranger")
  set.seed(11)
  n <- 360
  domain <- factor(rep(1:4, each = n / 4))
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  noise <- replicate(8, rnorm(n))
  y <- rbinom(n, 1, plogis(1.8 * x1 - 1.2 * x2))
  dat <- data.frame(
    lon = rep(1:4, each = n / 4) + rnorm(n),
    lat = rep(1:4, each = n / 4) + rnorm(n),
    presence = y, x1 = x1, x2 = x2, noise
  )
  out <- cast_select(
    dat, method = "stable", domains = domain, max_vars = 6,
    min_vars = 2, num_trees = 60, seed = 12, verbose = FALSE
  )
  expect_s3_class(out, "cast_select")
  expect_gte(length(out$selected), 2)
  expect_true(all(out$selected %in% names(dat)))
  expect_named(
    out$diagnostics,
    c("domains", "n_domains", "candidate_variables", "p_domain",
      "p_interaction", "invariance_pass", "domain_logloss",
      "search_path", "use_invariance")
  )
  expect_false("causal_role" %in% names(out$scores))
})
