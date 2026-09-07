## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)


## ----data---------------------------------------------------------------------
library(castSDM)

data_file <- system.file("extdata",
  "CAST_Alces_alces_Res9_screened.csv", package = "castSDM")
dat <- read.csv(data_file)

# Keep the workflow light for this vignette: subsample rows.
set.seed(1)
dat <- dat[c(
  sample(which(dat$presence == 1), 150),
  sample(which(dat$presence == 0), 150)
), ]

head(dat[, c("lon", "lat", "presence", "bio01", "bio12", "elevation")])


## ----pipeline-----------------------------------------------------------------
result <- cast(
  dat,
  models = c("rf", "gam"),
  select_method = "cpi",
  select_num_trees = 100,
  select_min_vars = 3,      # small example: floor the screen at 3 predictors
  do_cv = TRUE,
  cv_k = 3,
  seed = 42,
  verbose = FALSE
)

summary(result)
plot(result$screen)


## ----importance---------------------------------------------------------------
eff <- cast_importance(result)               # tidy CPI table + CIs
plot(eff)                                    # coefficient (forest) plot


## ----sensitivity--------------------------------------------------------------
# Shift the first predictor retained by the screen (guaranteed to be fitted).
cf <- cast_sensitivity(result$fit, newdata = dat,
                       variable = result$fit$env_vars[1], shift = 1,
                       shift_type = "sd")
plot(cf, basemap = "none")


## ----predict------------------------------------------------------------------
# Build a small prediction grid from the data range.
grid <- expand.grid(
  lon = seq(min(dat$lon), max(dat$lon), length.out = 20),
  lat = seq(min(dat$lat), max(dat$lat), length.out = 20)
)
env_vars <- result$fit$env_vars
for (v in env_vars) {
  grid[[v]] <- median(dat[[v]], na.rm = TRUE)
}

pred <- cast_predict(result$fit, grid)
plot(pred, basemap = "none")

