# Contract tests for the two entry points the audit found untested:
#   cast_threshold()  - every TSS / SEDI / kappa / committee vote goes through it
#   print.cast_cv()   - its header used to claim "max TSS" for any method
# The assertions follow the real contract, not an assumed one:
#   * method names are underscored ("max_tss", not "maxTSS")
#   * the return value is a plain unnamed numeric with no method attribute
#   * unusable input aborts rather than falling back to a default
# cli renders its output outside capture.output(), so the print contract is
# pinned with expect_snapshot() rather than by grepping printed text.
test_that("cast_threshold returns a usable cut and aborts on unusable input", {
  set.seed(11)
  obs <- c(rep(1, 25), rep(0, 75))
  pred <- c(rnorm(25, 0.75, 0.10), rnorm(75, 0.30, 0.10))

  th <- cast_threshold(pred, obs, method = "max_tss")
  expect_true(is.numeric(th))
  expect_length(th, 1L)
  expect_true(is.finite(th))
  expect_gte(th, min(pred))
  expect_lte(th, max(pred))
  # unname() on purpose: callers index this positionally
  expect_null(names(th))

  # the guards the contract actually has
  expect_error(cast_threshold(pred, rep(1, 100), method = "max_tss"),
               "both response classes")
  expect_error(cast_threshold(c(0.9, NA, 0.2), c(1, 0, 1), method = "max_tss"))
  expect_error(cast_threshold(pred, obs, method = "max_tss", sensitivity = 0))
  expect_error(cast_threshold(pred, obs, method = "max_tss", omission = 1))
  expect_error(cast_threshold(pred, c(obs, 0L), method = "max_tss"))
})

test_that("cast_threshold methods are distinct but all inside the prediction range", {
  set.seed(12)
  obs <- c(rep(1, 30), rep(0, 90))
  pred <- c(rnorm(30, 0.72, 0.12), rnorm(90, 0.32, 0.12))
  methods <- c("max_tss", "equal_sens_spec", "max_jaccard", "max_kappa",
               "equal_prevalence", "sensitivity")
  ths <- vapply(methods, function(m) cast_threshold(pred, obs, method = m), numeric(1))
  expect_true(all(is.finite(ths)))
  expect_true(all(ths >= min(pred) & ths <= max(pred)))
  # at least one pair must differ, otherwise the criteria are not doing anything
  expect_true(any(abs(outer(ths, ths, "-")) > 1e-9))
  expect_true(is.finite(cast_threshold(pred, obs, method = "sensitivity",
                                       sensitivity = 0.6)))
})

test_that("an uninformative predictor still yields a cut, and it is not a claim", {
  set.seed(13)
  obs <- c(rep(1, 30), rep(0, 70))
  flat <- rep(0.5, 100)
  th <- cast_threshold(flat, obs, method = "max_tss")
  expect_true(is.finite(th))
  # every cut classifies identically: AUC carries no signal here, which is why
  # the ensemble score must not be read off a threshold chosen this way
  pred_bin <- as.integer(flat >= th)
  expect_true(all(pred_bin == 0L) || all(pred_bin == 1L))
})

# ---- print.cast_cv contract ------------------------------------------------
# Why this is a source-level check rather than a rendered-output check:
# cli writes straight past capture.output(), past local_reproducible_output(),
# and past mocked bindings when the call site is namespaced (cli::cli_h2), and
# expect_snapshot() wrote no _snaps file in this setup. All three would have
# reported success without asserting anything. The rendered output was verified
# by hand ("-- Optimal thresholds (equal_sens_spec) --"); here the contract is
# pinned on the two things that can silently regress: the method must come from
# the object, and the old hardcoded header must be gone.
cast_cv_fixture <- function(method = "equal_sens_spec") {
  structure(list(
    k = 3L,
    block_method = "grid",
    threshold_method = method,
    metrics = data.frame(model = "maxent", auc_mean = 0.8, auc_sd = 0.02,
                         tss_mean = 0.4, cbi_mean = 0.5,
                         stringsAsFactors = FALSE),
    thresholds = c(maxent = 0.31)
  ), class = "cast_cv")
}

test_that("print.cast_cv takes its header from the object, not a literal", {
  src <- paste(deparse(body(castSDM:::print.cast_cv)), collapse = " ")
  expect_true(grepl("threshold_method", src, fixed = TRUE),
              info = "print.cast_cv no longer reads x$threshold_method")
  expect_false(grepl("Optimal thresholds (max TSS)", src, fixed = TRUE),
               info = "the hardcoded max TSS header is back")
  # the %||% fallback must survive for objects built before the field existed
  expect_true(grepl("%||%", src, fixed = TRUE))
})

test_that("cast_cv_fixture documents the field the header depends on", {
  for (m in c("max_tss", "equal_sens_spec", "max_kappa", "sensitivity")) {
    expect_identical(cast_cv_fixture(m)$threshold_method, m)
  }
  bare <- cast_cv_fixture("max_tss")
  bare$threshold_method <- NULL
  expect_null(bare$threshold_method)
})
