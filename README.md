# castSDM

`castSDM` is an end-to-end species distribution modelling toolkit for causal
variable selection, nested spatial validation, standard SDM ensembles, raster
prediction, and future projection.

Its distinguishing method is `cast_select(method = "cpi")` (the default),
a conditional predictive impact (CPI) selector. For each candidate
predictor it replaces the variable with a knockoff given every other
predictor and measures the loss of predictive accuracy, then keeps the
predictors whose conditional contribution survives Benjamini-Hochberg FDR
control. An optional double machine-learning selector
(`method = "dml"`) reports signed Neyman-orthogonal partial-linear effects
with confidence intervals via [cast_effect()], and
[cast_counterfactual()] maps single-driver what-if shifts on the current
climate.

The only ecological choice is the FDR level (`alpha`). Cross-fitting folds and
the random-forest nuisance learner are method defaults, so the selector avoids
the many hand-tuned knobs of traditional screening. Selection reports adjusted
effect sizes with confidence intervals; it targets interpretable driver
analysis but does not claim automatic causal discovery from observational data.

## Core workflow

```r
library(castSDM)

result <- cast(
  species_data,
  env_data = prediction_grid,
  select_method = "cpi",
  models = c("rf", "brt", "maxent", "gam"),
  do_cv = TRUE,
  seed = 42
)

summary(result)
plot(result$screen)          # conditional impacts, FDR-adjusted
```

Signed effect reporting with the DML engine:

```r
screen_dml <- cast_select(species_data, method = "dml")
eff <- cast_effect(screen_dml)   # tidy DML effect table + CIs
plot(eff)                        # coefficient (forest) plot
```

The pipeline runs:

```text
prepare -> DML selection -> fit -> nested spatial CV
        -> evaluate -> predict -> ensemble -> projection
```

Selection is repeated inside every outer spatial training fold. Held-out folds
never influence variable choice or tuning.

## Causal interpretation

Two functions turn a fitted workflow into causal-flavoured evidence without a
new estimation engine:
```r
eff <- cast_effect(result)               # tidy DML effect table + CIs
plot(eff)                                # coefficient (forest) plot

cf <- cast_counterfactual(               # g-computation what-if on current
  result$fit, newdata = current_grid,    #   climate: shift one driver, hold
  variable = "bio1", shift = 1           #   the rest fixed, map the change
)
plot(cf, basemap = "china")              # diverging map of change in HSS
```

## Main functions

| Stage | Functions |
|---|---|
| Study design | `cast_study_area()`, `cast_background()` |
| Preparation | `cast_prepare()`, `get_env_vars()`, `cast_vif()` |
| Causal selection | `cast_select()` |
| Modelling | `cast_fit()`, `cast_esm()` |
| Validation | `cast_cv()`, `cast_evaluate()` |
| Prediction | `cast_predict()`, `cast_predict_tiled()` |
| Ensemble/projection | `cast_ensemble()`, `cast_project()` |
| Causal interpretation | `cast_effect()`, `cast_counterfactual()` |
| Batch/resume | `cast_batch()`, `cast_batch_resume()` |

## Model backends

- Random Forest: `ranger`
- BRT: `gbm`
- MaxEnt: `maxnet`
- GAM: `mgcv`
- DML selector: `DoubleML` + `mlr3` / `mlr3learners` (random-forest nuisance)

Optional functionality is activated only when its suggested package is
installed. Package installation is never attempted during a model run.

## Installation

```r
pak::pak("EldonQ/castSDM")
```

For the full workflow:

```r
install.packages(c(
  "DoubleML", "mlr3", "mlr3learners", "ranger", "gbm", "maxnet", "mgcv",
  "pROC", "ggplot2", "sf", "terra", "future", "future.apply", "yaml"
))
```

## Interpretation

- Selected variables are FDR-significant partial-effect drivers, adjusted for
  the remaining predictors.
- The DML effect is an association purged of the measured confounders; it is
  not proof of a manipulable causal mechanism.
- Counterfactual maps are purely interpretive what-if summaries on the current
  climate; they do not extrapolate to future scenarios.
- Future projections assume that the learned response relationship remains
  applicable under the projected environment.

## License

GPL (>= 3)
