# castSDM single-species example: Ovis ammon
#
# This example demonstrates the current response-focused MB screening workflow.
# It intentionally uses screen$roles directly; removed role helpers are not
# part of the active package API.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools or install castSDM before running this example.")
}

pkg_root <- normalizePath(file.path("E:/Package/cast"), winslash = "/", mustWork = FALSE)
if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  devtools::load_all(pkg_root)
} else {
  library(castSDM)
}

data(Ovis_ammon)
data(china_env_grid)

split <- cast_prepare(
  Ovis_ammon,
  train_fraction = 0.7,
  seed = 42,
  verbose = TRUE
)

dag <- cast_dag(
  split$train,
  structure_method = "mb_first",
  include_response = TRUE,
  response_as_sink = TRUE,
  seed = 42,
  verbose = TRUE
)

screen <- cast_select(
  dag,
  split$train,
  min_vars = 5,
  min_fraction = 0.3,
  seed = 42,
  verbose = TRUE
)

print(screen)
print(screen$roles)

fit <- cast_fit(
  split$train,
  screen = screen,
  dag = dag,
  models = c("rf", "brt", "maxent", "gam"),
  seed = 42,
  verbose = TRUE
)

eval <- cast_evaluate(fit, split$test)
print(eval)

cv <- cast_cv(
  Ovis_ammon,
  screen = screen,
  dag = dag,
  models = c("rf", "brt"),
  k = 5,
  seed = 42,
  verbose = TRUE
)

pred <- cast_predict(fit, china_env_grid)
ens <- cast_ensemble(fit, cv, china_env_grid, method = "weighted")

print(pred)
print(ens)
