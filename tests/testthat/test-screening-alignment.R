test_that("old explanation and role helper APIs are not exported", {
  exports <- getNamespaceExports("castSDM")
  expect_false(any(c(
    paste0("cast_", "sh", "ap_xgb"),
    paste0("cast_", "sh", "ap_fit"),
    paste0("cast_", "explain"),
    paste0("cast_", "roles"),
    paste0("cast_", "backdoor"),
    paste0("cast_", "report")
  ) %in% exports))
})

test_that("response_markov_blanket prefers Stage 1 MB metadata", {
  dag <- new_cast_dag(
    edges = data.frame(
      from = c("bio01", "bio02"),
      to = c("presence", "bio03"),
      strength = 1,
      direction = 1,
      stringsAsFactors = FALSE
    ),
    nodes = c("bio01", "bio02", "bio03", "presence"),
    boot_R = NA_integer_,
    strength_threshold = 0.7,
    direction_threshold = 0.6,
    score = "mb_first:fast.iamb",
    structure_method = "mb_first",
    response_node = "presence",
    metadata = list(
      mb_vars_stage1 = c("bio01", "bio03"),
      response_as_sink = TRUE
    )
  )

  mb <- castSDM:::response_markov_blanket(
    dag,
    "presence",
    env_vars = c("bio01", "bio02", "bio03")
  )

  expect_setequal(mb$all, c("bio01", "bio03"))
  expect_setequal(mb$parents, c("bio01", "bio03"))
  expect_length(mb$children, 0)
  expect_length(mb$co_parents, 0)
})

test_that("cast_select emits current screening role names", {
  skip_if_not_installed("ranger")
  set.seed(1)
  n <- 80
  dat <- data.frame(
    lon = runif(n),
    lat = runif(n),
    presence = rbinom(n, 1, 0.5),
    bio01 = rnorm(n),
    bio02 = rnorm(n),
    bio03 = rnorm(n),
    bio04 = rnorm(n),
    bio05 = rnorm(n),
    stringsAsFactors = FALSE
  )
  dag <- new_cast_dag(
    edges = data.frame(
      from = c("bio01", "bio02"),
      to = c("presence", "presence"),
      strength = 1,
      direction = 1,
      stringsAsFactors = FALSE
    ),
    nodes = c("bio01", "bio02", "bio03", "bio04", "bio05", "presence"),
    boot_R = NA_integer_,
    strength_threshold = 0.7,
    direction_threshold = 0.6,
    score = "mb_first:fast.iamb",
    structure_method = "mb_first",
    response_node = "presence",
    metadata = list(
      mb_vars_stage1 = c("bio01", "bio02"),
      response_as_sink = TRUE
    )
  )

  screen <- cast_select(
    dag,
    dat,
    min_vars = 2,
    min_fraction = 0,
    num_trees = 20,
    seed = 1,
    verbose = FALSE
  )

  expect_true(all(screen$roles$role %in% c(
    "mb_direct", "mb_associated", "importance_added", "importance_screened"
  )))
  expect_true(all(c("bio01", "bio02") %in% screen$selected))
})

test_that("cast_batch exposes prediction and ensemble controls", {
  f <- names(formals(cast_batch))
  expect_true(all(c("do_predict", "do_ensemble", "ensemble_method") %in% f))
})
