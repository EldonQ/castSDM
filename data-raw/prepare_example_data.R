## prepare_example_data.R
## ============================================================================
## Convert raw CSV example data to .rda for castSDM package
##
## Data: Ovis ammon (Argali sheep) occurrence + environment on ISEA3H Res9
##       hexagonal grid, China region. VIF-screened to 11 env variables.
##
## Source: EcoISEA3H project (causal SDM framework)
## ============================================================================

# -- Species occurrence data --
ovis_ammon <- read.csv(
  "inst/extdata/CAST_Ovis_ammon_Res9_screened.csv",
  stringsAsFactors = FALSE
)

# -- Full environmental grid (for spatial prediction) --
china_env_grid <- read.csv(
  "inst/extdata/China_EnvData_Res9_Screened.csv",
  stringsAsFactors = FALSE
)

# Save as .rda (compressed)
usethis::use_data(ovis_ammon, overwrite = TRUE, compress = "xz")
usethis::use_data(china_env_grid, overwrite = TRUE, compress = "xz")
