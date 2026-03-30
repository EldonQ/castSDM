test_that("cast_dag constructor creates correct class", {
  dag <- new_cast_dag(
    edges = data.frame(from = "bio1", to = "presence", strength = 0.8,
                       direction = 0.9),
    nodes = c("bio1", "presence"),
    boot_R = 100L,
    strength_threshold = 0.7,
    direction_threshold = 0.6
  )
  expect_s3_class(dag, "cast_dag")
  expect_equal(nrow(dag$edges), 1)
  expect_equal(dag$boot_R, 100L)
})

test_that("cast_ate constructor creates correct class", {
  ate <- new_cast_ate(
    estimates = data.frame(variable = "bio1", coef = 0.5, se = 0.1,
                           p_value = 0.01, significant = TRUE),
    K = 2L
  )
  expect_s3_class(ate, "cast_ate")
  expect_equal(ate$K, 2L)
})

test_that("cast_result constructor creates correct class", {
  result <- new_cast_result(
    dag = new_cast_dag(data.frame(), character(0), 100L, 0.7, 0.6),
    ate = new_cast_ate(data.frame(), 2L),
    screen = new_cast_screen(character(0), data.frame(), c(w_dag = 0.3, w_ate = 0.3, w_imp = 0.4)),
    roles = new_cast_roles(data.frame()),
    fit = new_cast_fit(list(), character(0), character(0), list(means = numeric(0), sds = numeric(0))),
    eval = new_cast_eval(data.frame())
  )
  expect_s3_class(result, "cast_result")
  expect_null(result$predict)
  expect_null(result$cate)
})
