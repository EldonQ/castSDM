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
  expect_true("causal_role" %in% names(screen$roles))
  expect_true(all(screen$roles$causal_role %in% c(
    "causal_core", "causal_adjuster", "predictive_rescue", "unstable_rejected"
  )))
})

test_that("response-focused DAG plot keeps isolated response node", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("igraph")

  dag <- new_cast_dag(
    edges = data.frame(
      from = "bio01", to = "bio02", strength = 1, direction = 1,
      stringsAsFactors = FALSE
    ),
    nodes = c("bio01", "bio02", "presence"),
    boot_R = NA_integer_,
    strength_threshold = 0.7,
    direction_threshold = 0.6,
    score = "mb_first:fast.iamb",
    structure_method = "mb_first",
    response_node = "presence",
    metadata = list(response_as_sink = TRUE)
  )
  p <- plot(dag)
  built <- ggplot2::ggplot_build(p)
  labels <- unlist(lapply(built$data, function(x) x$label), use.names = FALSE)
  expect_true("presence" %in% labels)
})

test_that("cast_refute_screen returns diagnostic tables", {
  skip_if_not_installed("ranger")
  set.seed(2)
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
    dag, dat, min_vars = 2, min_fraction = 0, num_trees = 20,
    seed = 2, verbose = FALSE
  )
  refute <- cast_refute_screen(
    dag, screen, dat, reps = 2, num_trees = 20,
    seed = 2, verbose = FALSE
  )
  expect_s3_class(refute, "cast_refute")
  expect_true(all(c("test", "overlap_fraction") %in% names(refute$tests)))
  expect_true(all(c("variable", "subset_frequency") %in% names(refute$summary)))
})

test_that("cast_batch exposes prediction and ensemble controls", {
  f <- names(formals(cast_batch))
  expect_true(all(c("do_predict", "do_ensemble", "ensemble_method") %in% f))
})

test_that("plot.cast_batch handles more than ten species", {
  skip_if_not_installed("ggplot2")
  sm <- expand.grid(
    species = paste0("sp", seq_len(20)),
    model = c("rf", "brt"),
    fold = seq_len(2),
    stringsAsFactors = FALSE
  )
  sm$auc <- stats::runif(nrow(sm), 0.6, 0.9)
  sm$tss <- stats::runif(nrow(sm), 0.2, 0.7)
  sm$cbi <- stats::runif(nrow(sm), 0.4, 0.9)
  batch <- new_cast_batch(
    species_metrics = sm,
    species = unique(sm$species),
    models = c("rf", "brt"),
    results = NULL
  )
  p <- plot(batch)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})
