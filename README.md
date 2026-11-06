# castSDM

**Auditable scenario comparisons for species distribution models.**

`cast_compare()` constructs current, first-group-only, second-group-only and
joint predictions from a frozen fit. Fixed positive model weights and an
explicit finite domain make the comparisons auditable. Results contain
continuous contrasts, optional thresholded range counts, training-support
diagnostics and paired model SD.

## Installation

```r
install.packages("remotes")
remotes::install_github("EldonQ/castSDM")
install.packages(c("ranger", "maxnet", "mgcv", "glmnet", "pROC", "blockCV", "FNN", "sf", "terra"))
```

## A complete small comparison

This example supplies a fitted GLM through the public external-model adapter.
External models supply their callback and training reference; evaluation
scores are never invented for them.

```r
library(castSDM)
set.seed(42)
current <- data.frame(temperature = rnorm(120), land = runif(120))
observations <- transform(current,
  presence = rbinom(120, 1, plogis(temperature - land)))
model <- glm(presence ~ temperature + land, observations, family = binomial())
fit <- cast_as_fit(list(glm = model), current,
  function(model, newdata) predict(model, newdata, type = "response"))
future <- transform(current, temperature = temperature + .5, land = land + .1)
comparison <- cast_compare(fit, current, future,
  groups = list(climate = "temperature", land_use = "land"))
summary(comparison)
plot(comparison)
```

Without usable cross-validation or a supplied threshold there is no binary
range estimate. `summary(comparison, thresholds = c(.3, .5, .7))` calculates
sensitivity from saved predictions without fitting or predicting again.

## Fit, validate, compare

`cast_run()` handles screening, spatial validation, independent hold-out
evaluation and full-data refitting. Scientific settings live in one
`cast_protocol()`. Screening and tuning repeat inside outer CV folds; hold-out
data do not determine ensemble weights or thresholds.

```r
# species_data: lon, lat, presence (0/1), and numeric predictors.
protocol <- cast_protocol(models = c("rf", "maxent", "gam"), seed = 42)
run <- cast_run(species_data, protocol = protocol,
  output_dir = "results/training", cache_dir = "cache/training")
comparison <- cast_compare(run, current_raster, future_rasters,
  groups = list(climate = climate_variables, land_use = land_use_variables),
  mask = analysis_mask, output_dir = "results/comparison",
  resources = list(cores = 8, workers = 2, memory_mb = 4096))
```

Tables and aligned `terra::SpatRaster` inputs share a numerical kernel.
Unlisted fitted predictors are recorded as static. Their future layers may be
omitted; supplied static values must equal current values. Inputs must already
share grid, extent and CRS. No automatic resampling occurs.

The [single-species example](examples/single_species/README.md) supplies
occurrence QC, domain masks, exact future-file mapping and saved-output plots.
Training and comparison use public package functions.

## Interpretation

Let `C` be current, `A` first-group-only, `B` second-group-only and `J` joint.
Legacy `climate` and `lulc` prefixes denote the first and second groups even
when those groups have different names.

| Field | Definition |
|---|---|
| `climate_main` | `A - C`, baseline simple effect |
| `lulc_main` | `B - C`, baseline simple effect |
| `joint_change` | `J - C` |
| `interaction` | `J - A - B + C` |
| `climate_conditional` | `J - B` |
| `lulc_conditional` | `J - A` |
| `climate_average` | `((A - C) + (J - B))/2` |
| `lulc_average` | `((B - C) + (J - A))/2` |

These are **model-conditional scenario contrasts** on a fitted suitability
surface. They do not establish ecological causation, occurrence probabilities
or realised distributions. Historical `*_main` fields mean baseline simple
effects; the averaged factorial effects have explicit `*_average` fields.

Each cell uses all positive-weight models in every corner. Model SD is the
unweighted sample SD of within-model contrasts among active models; it is
`NA` with one model and is not a confidence interval. Training support uses
marginal quantiles on unclamped inputs and cannot establish joint support.
Full-domain and supported summaries remain separate. Requested clamping
applies uniformly to all engines.

Binary predictions use `suitability >= threshold`. Binary summaries count
valid grid cells or table rows, not square kilometres. Area interpretation
requires an appropriate equal-area grid or a separate area-weighted analysis.

`cast_factorial_attribution()` accepts four already-created single-layer
surfaces or equal-length numeric vectors. Ensemble surfaces cannot recover
paired model disagreement; optional surface-SD propagation is labelled as a
separate diagnostic. `cast_attribution()` provides optional driver-only and
retained-driver comparisons. Individual effects need not add to total change.

## Execution and reproducibility

Raster processing skips empty masked blocks and publishes artifacts atomically.
Current predictions are shared across futures. Additive univariate logit GAMs
reuse exact terms for partial scenarios; other structures use ordinary
prediction. A total core and memory budget bounds scenario workers and RF
threads. Workers require an installed matching package.

Training-stage cache keys depend on data, scientific settings and relevant
code, independently of output paths and future scenarios. Comparison receipts
also verify model/input content and output hashes before reuse. OS file locks
prevent simultaneous writes. Resume in the same directory; choose a new
comparison directory or `overwrite = TRUE` after changing scientific inputs.
Caching avoids repeated work and does not replace scientific validation.

See `vignette("castSDM")` and `vignette("castSDM-advanced")` for executable
examples and the detailed contracts.

## License

GPL (>= 3).
