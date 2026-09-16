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

test_that("each predictor is compared with its own null scale", {
  skip_if_not_installed("ranger")
  set.seed(812)
  d <- data.frame(presence = rep(0:1, 100), a = rnorm(200),
                   b = rnorm(200), c = rnorm(200))
  local_mocked_bindings(
    .cast_importance_fit = function(X, y, num_trees, seed) {
      list(model = fake_forest("a", 3),
           importance = c(a = 1, b = 0.02, c = 0))
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
  # The permutation-importance diagnostic is still reported alongside.
  expect_true("perm_importance" %in% names(s$scores))
  expect_true(any(is.finite(s$scores$perm_importance)))
  expect_true(is.finite(s$diagnostics$importance_agreement))
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
           importance = stats::setNames(rep(0, ncol(X)), names(X)))
    },
    .cast_shift_effect = fake_shift_effect,
    .package = "castSDM")
  local_mocked_bindings(
    ranger = function(x, y, ...) fake_forest("a", 3),
    .package = "ranger")
  expect_warning(s <- cast_select(d, n_perm = 19, seed = 1, verbose = FALSE),
                  "keeping the stage-1 set")
  expect_true(all(s$scores$fallback))
  expect_false(any(s$scores$passed_null))
})

test_that("the interventional effect is larger for the true driver on real data", {
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
  expect_true(is.finite(s$diagnostics$importance_agreement))
})

test_that("a collider is not selected on the strength of the association it creates", {
  skip_if_not_installed("ranger")
  # `coll` is a common effect of the driver `x1` and the response, so x1 has no
  # effect on the response and shifting it alone has a true effect of zero.
  # Conditioning on the collider manufactures an association, which is what
  # permutation importance rewards and what the shift contrast should decline.
  # This is the small-scale form of the benchmark's collider scenario.
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
  # x1 has no effect, so the shift contrast must not clear its own null.
  expect_false("x1" %in% s$selected)
})

