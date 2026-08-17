# Vectorised MESS vs the original double-loop implementation ---------------
#
# `.cast_mess_loop()` below is a verbatim copy of the pre-vectorisation
# implementation (per-cell sum(ref < pi) + per-cell branch loop). It serves
# as the oracle: the vectorised `.cast_mess()` in R/predict-spatial.R must
# reproduce it branch for branch, including degenerate ranges and NA
# semantics.

.cast_mess_loop <- function(reference, newdata) {
  vars <- intersect(names(reference), names(newdata))
  if (!length(vars)) return(rep(NA_real_, nrow(newdata)))
  sim <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(vars))
  for (j in seq_along(vars)) {
    v <- vars[j]
    ref <- reference[[v]][is.finite(reference[[v]])]
    p <- as.numeric(newdata[[v]])
    if (!length(ref)) next
    mn <- min(ref); mx <- max(ref); rng <- mx - mn
    n <- length(ref)
    f <- vapply(p, function(pi) {
      if (!is.finite(pi)) return(NA_real_)
      100 * sum(ref < pi) / n
    }, numeric(1))
    s <- numeric(length(p))
    for (i in seq_along(p)) {
      pi <- p[i]; fi <- f[i]
      if (!is.finite(pi)) { s[i] <- NA_real_; next }
      if (rng <= 0) {
        s[i] <- if (pi == mn) 100 else -100
      } else if (fi == 0) {
        s[i] <- (pi - mn) / rng * 100
      } else if (fi <= 50) {
        s[i] <- 2 * fi
      } else if (fi < 100) {
        s[i] <- 2 * (100 - fi)
      } else {
        s[i] <- (mx - pi) / rng * 100
      }
    }
    sim[, j] <- s
  }
  apply(sim, 1L, function(r) if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE))
}

test_that("vectorised .cast_mess matches the loop oracle on random data", {
  for (seed in 101:105) {
    set.seed(seed)
    n_ref <- sample(20:80, 1)
    n_new <- sample(50:200, 1)
    reference <- data.frame(
      a = rnorm(n_ref),
      b = runif(n_ref, -2, 2),
      c = rnorm(n_ref, 5, 3)
    )
    newdata <- data.frame(
      a = rnorm(n_new),
      b = runif(n_new, -4, 4),   # reaches outside the reference range
      c = rnorm(n_new, 5, 6)
    )
    newdata$a[sample(n_new, 3)] <- NA
    reference$b[sample(n_ref, 2)] <- NA
    expect_identical(.cast_mess(reference, newdata),
                     .cast_mess_loop(reference, newdata))
  }
})

test_that("vectorised .cast_mess handles ties and strict-below quantiles", {
  # f must count reference values strictly below p: pi == ref value must NOT
  # count that value (left.open findInterval vs sum(ref < pi)).
  reference <- data.frame(a = c(1, 2, 3, 4, 5))
  newdata <- data.frame(a = c(3, 5, 1, 0, 6))
  expect_identical(.cast_mess(reference, newdata),
                   .cast_mess_loop(reference, newdata))
  out <- .cast_mess(reference, newdata)
  expect_equal(out[1], 80)            # f = 40 -> 2 * 40
  expect_equal(out[3], 0)             # pi == mn, f = 0 -> (1-1)/4 * 100
  expect_true(out[4] < 0)             # below reference min
  expect_true(out[5] < 0)             # above reference max
})

test_that("vectorised .cast_mess keeps degenerate-range and NA semantics", {
  reference <- data.frame(a = rep(3, 10), b = 1:10)  # range(a) == 0
  newdata <- data.frame(
    a = c(3, 4, NA, 2),
    b = c(5, 0, 11, NA)
  )
  expect_identical(.cast_mess(reference, newdata),
                   .cast_mess_loop(reference, newdata))
  out <- .cast_mess(reference, newdata)
  expect_equal(out[1], 80)            # a == mn -> 100; b=5 -> f=40 -> 2*40; min
  expect_equal(out[2], -100)          # a off the degenerate range
  expect_true(out[3] < 0)             # a NA, b above max -> negative
  expect_equal(out[4], -100)          # a off range, b NA -> min over finite

  # All-NA rows stay NA; non-intersecting variables are ignored.
  expect_identical(.cast_mess(reference, data.frame(a = NA, b = NA)),
                   NA_real_)
  expect_identical(.cast_mess(data.frame(z = 1:5), data.frame(a = 1:3)),
                   rep(NA_real_, 3))
})

test_that("cast_ensemble retains MESS extrapolation flags", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(110)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - 1.2 * x2)),
    x1 = x1, x2 = x2
  )
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 111, verbose = FALSE)
  cv <- new_cast_cv(
    metrics = data.frame(
      model = "rf",
      auc_mean = 0.8, auc_sd = 0.1, tss_mean = 0.5, tss_sd = 0.1,
      cbi_mean = 0.4, cbi_sd = 0.1, n_folds = 3, n_selected_mean = 2
    ),
    fold_metrics = data.frame(), folds = integer(n), k = 3L,
    block_method = "grid", thresholds = c(rf = 0.5)
  )
  grid <- data.frame(
    lon = runif(50, 70, 130), lat = runif(50, 20, 50),
    x1 = rnorm(50), x2 = rnorm(50)
  )
  ens <- cast_ensemble(fit, cv, grid, method = "weighted")
  expect_true(all(c("mess", "extrapolating") %in% names(ens$predictions)))
  expect_identical(ens$predictions$extrapolating, ens$predictions$mess < 0)
})

test_that("cast_predict aborts on factor/character predictors", {
  skip_if_not_installed("ranger")
  set.seed(112)
  n <- 100
  dat <- data.frame(
    lon = runif(n), lat = runif(n),
    presence = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n)
  )
  screen <- new_cast_select(c("x1", "x2"),
                            data.frame(variable = c("x1", "x2")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 40,
                  seed = 113, verbose = FALSE)
  bad <- data.frame(x1 = rnorm(5), x2 = factor(letters[1:5]))
  expect_error(cast_predict(fit, bad), "[Nn]on-numeric")
  bad2 <- data.frame(x1 = rnorm(5), x2 = letters[1:5])
  expect_error(cast_predict(fit, bad2), "[Nn]on-numeric")
})
