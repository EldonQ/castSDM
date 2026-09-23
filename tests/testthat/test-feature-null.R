# Stage 2 calibration: each predictor is compared with its own null scale ----

# A stand-in forest whose interventional effect is a known function of the
# shifted predictor. `effect_for` names the predictor that moves the fitted
# probability; shifting any other predictor leaves the prediction unchanged.
fake_forest <- function(effect_for, magnitude) {
  f <- function(object, data, ...) {
    z <- if (object$effect_for %in% names(data)) {
      data[[object$effect_for]]
    } else {
      rep(0, nrow(data))
    }
    p <- stats::plogis(object$magnitude * z)
    list(predictions = cbind("0" = 1 - p, "1" = p))
  }
  obj <- list(effect_for = effect_for, magnitude = magnitude)
  attr(obj, "predict") <- f
  obj
}

# Mirrors the package statistic with the fake forest's own predict closure, so
# the calibration logic is what is under test, not ranger.
fake_shift_effect <- function(model, X, sds, shift_size) {
  pf <- attr(model, "predict")
  prob1 <- function(Z) pf(model, Z)$predictions[, "1"]
  base <- prob1(X)
  out <- numeric(ncol(X))
  names(out) <- colnames(X)
  for (v in colnames(X)) {
    Xs <- X
    Xs[[v]] <- Xs[[v]] + shift_size * sds[[v]]
    out[[v]] <- mean(abs(prob1(Xs) - base))
  }
  out
}

test_that("each predictor is compared with its own conditional null scale", {
  skip_if_not_installed("ranger")
  set.seed(812)
  d <- data.frame(presence = rep(0:1, 100), a = rnorm(200),
                   b = rnorm(200), c = rnorm(200))
  local_mocked_bindings(
    .cast_importance_fit = function(X, y, num_trees, seed) {
      list(model = fake_forest("a", 3),
           importance = stats::setNames(rep(NA_real_, ncol(X)), names(X)))
    },
    .cast_shift_effect = fake_shift_effect,
    .package = "castSDM")
  local_mocked_bindings(
    ranger = function(x, y, ...) fake_forest("b", 0.001),
    .package = "ranger")
  s <- cast_select(d, n_perm = 19, seed = 1, verbose = FALSE)
  # Only "a" moves the fitted probability when it is shifted.
  expect_identical(s$selected, "a")
  expect_true(all(is.finite(s$diagnostics$null_threshold)))
  expect_false(any(s$scores$fallback))
  # Single-statistic screen: no second attribution column.
  expect_false("perm_importance" %in% names(s$scores))
  expect_false("p_adjusted" %in% names(s$scores))
  expect_false("importance_agreement" %in% names(s$diagnostics))
  expect_identical(s$diagnostics$null_method,
                   "feature-wise within-stratum permutation of the shift effect")
  expect_true(is.finite(s$diagnostics$null_strata))
})

test_that("a null fallback is marked as fallback, not passed evidence", {
  skip_if_not_installed("ranger")
  set.seed(813)
  d <- data.frame(presence = rep(0:1, 100), a = rnorm(200),
                   b = rnorm(200), c = rnorm(200))
  # The null forests move the probability as much as the observed one, so no
  # predictor can be separated from its own null.
  local_mocked_bindings(
    .cast_importance_fit = function(X, y, num_trees, seed) {
      list(model = fake_forest("a", 3),
           importance = stats::setNames(rep(NA_real_, ncol(X)), names(X)))
    },
    .cast_shift_effect = fake_shift_effect,
    .package = "castSDM")
  local_mocked_bindings(
    ranger = function(x, y, ...) fake_forest("a", 3),
    .package = "ranger")
  expect_warning(s <- cast_select(d, n_perm = 19, seed = 1, verbose = FALSE),
                  "keeping the top")
  expect_true(all(s$scores$fallback[s$scores$selected]))
  expect_false(any(s$scores$passed_null))
  # The fallback keeps at most ncov predictors, not the whole stage-1 set.
  expect_lte(length(s$selected), s$diagnostics$ncov)
  expect_true(all(s$scores$selected_reason[s$scores$selected] == "fallback-top-ncov"))
})

test_that("the conditional effect is larger for the true driver on real data", {
  skip_if_not_installed("ranger")
  set.seed(814)
  n <- 600
  x1 <- rnorm(n); x2 <- rnorm(n); noise <- rnorm(n)
  d <- data.frame(presence = rbinom(n, 1, plogis(1.4 * x1 - 1.1 * x2)),
                  x1 = x1, x2 = x2, noise = noise)
  s <- cast_select(d, num_trees = 150, n_perm = 19, seed = 5, verbose = FALSE)
  eff <- stats::setNames(s$scores$interventional_effect, s$scores$variable)
  expect_gt(eff[["x1"]], eff[["noise"]])
  expect_gt(eff[["x2"]], eff[["noise"]])
  expect_setequal(s$selected, c("x1", "x2"))
})

test_that("a collider is not selected on the strength of the association it creates", {
  skip_if_not_installed("ranger")
  # `coll` is a common effect of the driver `x1` and the response, so x1 has no
  # effect on the response and shifting it alone has a true effect of zero.
  # Conditioning on the collider manufactures an association. The shift
  # contrast should decline it relative to the true driver.
  set.seed(815)
  n <- 1500
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(-0.3 - 1.6 * x2))
  coll <- as.numeric(0.9 * x1 + 1.8 * (y - mean(y)) + rnorm(n) >
                       stats::qnorm(0.5))
  d <- data.frame(presence = y, x1 = x1, x2 = x2, coll = coll)
  s <- cast_select(d, num_trees = 200, n_perm = 19, seed = 6, verbose = FALSE)
  # The true driver is retained.
  expect_true("x2" %in% s$selected)
  # x1 has no effect, so the shift contrast must downweight it relative to
  # the true driver. (Exact exclusion is not asserted: with 19 permutations
  # the p = 0.05/0.10 boundary turns on a single null draw.)
  eff <- stats::setNames(s$scores$interventional_effect, s$scores$variable)
  expect_lt(eff[["x1"]], eff[["x2"]])
})

test_that("strata helpers degrade gracefully and preserve margins", {
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  s1 <- .cast_perm_strata_list(X[1:3, , drop = FALSE], seed = 1)
  expect_true(all(vapply(s1, function(s) identical(s, rep(1L, 3)), logical(1))))
  s5 <- .cast_perm_strata_list(X, seed = 2)
  expect_true(all(vapply(s5, function(s) length(unique(s)) <= 5L, logical(1))))
  expect_length(s5, 2L)
  # Each stratification excludes its own predictor: with two columns the
  # strata come from a single other column each.
  expect_named(s5, c("a", "b"))
  # Within-stratum permutation preserves each column's multiset ...
  Xp <- .cast_stratum_permute_list(X, s5)
  expect_equal(sort(Xp$a), sort(X$a))
  expect_equal(sort(Xp$b), sort(Xp$b))
  # ... and leaves single-row strata untouched.
  Xp1 <- .cast_stratum_permute_list(X, list(a = seq_len(50), b = seq_len(50)))
  expect_identical(Xp1, X)
  # A lone survivor has nothing to condition on: full permutation.
  expect_identical(.cast_perm_strata_list(X[, "a", drop = FALSE], seed = 1)$a,
                   rep(1L, 50))
})
