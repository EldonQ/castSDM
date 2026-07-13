# castSDM single-species example: Ovis ammon
#
# This example demonstrates the role-constrained causal core.

library(castSDM)

data(Ovis_ammon)
data(china_env_grid)

split <- cast_prepare(
  Ovis_ammon,
  train_fraction = 0.7,
  seed = 42,
  verbose = TRUE
)

screen <- cast_select(
  data = split$train,
  min_vars = 5,
  cor_threshold = 0.8,
  seed = 42,
  verbose = TRUE
)

print(screen)
print(screen$roles)

fit <- cast_fit(
  split$train,
  screen = screen,
  models = c("rf", "brt", "maxent", "gam"),
  seed = 42,
  verbose = TRUE
)

eval <- cast_evaluate(fit, split$test)
print(eval)

cv <- cast_cv(
  Ovis_ammon,
  screen = screen,
  models = c("rf", "brt"),
  k = 5,
  seed = 42,
  verbose = TRUE
)

pred <- cast_predict(fit, china_env_grid)
ens <- cast_ensemble(fit, cv, china_env_grid, method = "weighted")

print(pred)
print(ens)
