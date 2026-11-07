# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
test_that("training stages resume across output directories and invalidate scientific changes", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  set.seed(835)
  d <- data.frame(lon = runif(120, 70, 130), lat = runif(120, 20, 50),
    x = rnorm(120), z = rnorm(120))
  d$presence <- rbinom(120, 1, plogis(d$x))
  p <- cast_protocol(models = "rf", select_method = "full", rf_ntree = 15,
    cv_k = 2, block_method = "grid", aoa = FALSE, seed = 8)
  root <- tempfile("training-cache-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  run <- function(where, protocol = p) cast_run(d, protocol = protocol,
    output_dir = file.path(root, where), cache_dir = file.path(root, "shared"), verbose = FALSE)
  first <- run("a"); second <- run("b")
  hits <- function(x) vapply(x$receipts, function(r) isTRUE(r$cache_hit), logical(1))
  expect_true(all(hits(second)), info = paste("Missed stages:",
    paste(vapply(second$receipts[!hits(second)], `[[`, character(1), "stage"), collapse = ", ")))
  expect_equal(second$cv$oof, first$cv$oof)
  p$rf_ntree <- 20
  third <- run("c", p)
  expect_equal(hits(third)[1:2], c(TRUE, TRUE))
  expect_false(any(hits(third)[3:6]))
  split <- cast_prepare(d, block_method = "grid", seed = 8, verbose = FALSE)
  expect_identical(.cast_verify_prepared(d, split), split)
  split$train$x[1] <- split$train$x[1] + 1
  expect_error(.cast_verify_prepared(d, split), "do not match")
})

test_that("source metadata does not invalidate serialized scientific values", {
  source <- srcfilecopy("source-example.R", "presence ~ x")
  alias <- srcfilealias("source-example.R", source)
  call <- quote(presence ~ x)
  attr(call, "srcref") <- srcref(alias, c(1L, 1L, 1L, 12L, 1L, 12L))
  value <- list(call = call, coefficient = .7, data = data.frame(x = 1:3))
  key <- .cast_digest(value)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  .cast_cache_write(value, path, key)
  restored <- .cast_cache_read(path, key)
  expect_false(is.null(restored))
  expect_equal(restored$data, value$data)
  expect_identical(.cast_digest(restored), key)
  source$lines <- "unrelated source display bookkeeping"
  expect_identical(.cast_digest(value), key)
  value$coefficient <- .8
  expect_false(identical(.cast_digest(value), key))
  value$coefficient <- .7
  value$call[[3L]] <- quote(z)
  expect_false(identical(.cast_digest(value), key))
})
