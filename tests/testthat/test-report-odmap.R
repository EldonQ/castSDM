# ODMAP reporting ------------------------------------------------------------

make_odmap_result <- function(n = 300, seed = 11) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  dat <- data.frame(
    lon = runif(n, 70, 130), lat = runif(n, 20, 50),
    presence = rbinom(n, 1, plogis(1.5 * x1 - 1.2 * x2)),
    x1 = x1, x2 = x2, x3 = x3
  )
  screen <- new_cast_select(c("x1", "x2", "x3"),
                            data.frame(variable = c("x1", "x2", "x3")),
                            method = "manual")
  fit <- cast_fit(dat, screen = screen, models = "rf", rf_ntree = 50,
                  seed = 6, verbose = FALSE)
  new_cast_result(screen, fit, cast_evaluate(fit, dat))
}

test_that("cast_report_odmap writes all five ODMAP sections", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  res <- make_odmap_result()
  path <- tempfile(fileext = ".md")
  on.exit(unlink(path), add = TRUE)
  expect_equal(cast_report_odmap(res, path = path, verbose = FALSE), path)
  txt <- readLines(path)
  expect_true(all(c("## 1. Overview", "## 2. Data", "## 3. Model",
                    "## 4. Assessment", "## 5. Prediction") %in% txt))
  expect_true(any(grepl("x1", txt, fixed = TRUE)))
})

test_that("cast_report_odmap uses meta and flags unanswerable items", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")
  res <- make_odmap_result()
  path <- tempfile(fileext = ".md")
  on.exit(unlink(path), add = TRUE)
  cast_report_odmap(res, path = path, verbose = FALSE,
                    meta = list(taxon = "Panthera uncia",
                                purpose = "Range shift under SSP scenarios"))
  txt <- paste(readLines(path), collapse = "\n")
  expect_match(txt, "Panthera uncia", fixed = TRUE)
  expect_match(txt, "Range shift under SSP scenarios", fixed = TRUE)
  # Author-supplied items stay explicitly open rather than being invented.
  expect_match(txt, "**Authors**: [report]", fixed = TRUE)
  # This result carries no CV and no prediction, which must be stated, not
  # silently rendered as an empty field.
  expect_match(txt, "[report] (no spatial CV recorded in this object)",
               fixed = TRUE)
  expect_match(txt, "[report] (no spatial prediction recorded in this object)",
               fixed = TRUE)
})

test_that("cast_report_odmap rejects a non-cast_result", {
  expect_error(cast_report_odmap(list(a = 1), path = tempfile()),
               "cast_result")
})
