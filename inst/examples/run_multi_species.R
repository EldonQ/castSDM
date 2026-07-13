# castSDM multi-species example
#
# This example uses the current v0.4 role-constrained selector.

library(castSDM)

data(Ovis_ammon)
data(china_env_grid)

species_list <- list(
  Ovis_ammon = Ovis_ammon
)

result <- cast_batch(
  species_list = species_list,
  env_data = china_env_grid,
  models = c("rf", "brt", "maxent", "gam"),
  output_dir = "castSDM_multi_species_example",
  parallel = FALSE,
  seed = 42,
  select_method = "causal_prior_rf",
  select_min_vars = 5,
  select_min_fraction = 0,
  select_cor_threshold = 0.8,
  do_predict = TRUE,
  do_ensemble = TRUE,
  ensemble_method = "weighted",
  do_cate = FALSE
)

print(result)
