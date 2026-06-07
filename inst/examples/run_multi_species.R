# castSDM multi-species example
#
# This example uses bundled package data and the current v0.3 API:
# response-focused MB screening and no removed role helpers.

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
  dag_structure_method = "mb_first",
  dag_include_response = TRUE,
  dag_response_as_sink = TRUE,
  select_min_vars = 5,
  select_min_fraction = 0.3,
  do_predict = TRUE,
  do_ensemble = TRUE,
  ensemble_method = "weighted",
  do_cate = FALSE
)

print(result)
