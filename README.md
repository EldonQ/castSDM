# castSDM

`castSDM` is an end-to-end species distribution modelling toolkit for stable
variable selection, nested spatial validation, standard SDM ensembles, raster
prediction, and future projection.

Its distinguishing method is `cast_select(method = "stable")`. The selector:

1. builds data-derived environments from the training data;
2. creates a small response-aware predictor shortlist;
3. tests whether the species response is conditionally invariant across those
   environments;
4. removes variables while preserving invariance and leave-domain-out
   predictive performance.

This causal-inspired criterion targets transferable prediction. It does not
claim that observational occurrence data identify true ecological causes.

## Core workflow

```r
library(castSDM)

result <- cast(
  species_data,
  env_data = prediction_grid,
  select_method = "stable",
  models = c("rf", "brt", "maxent", "gam"),
  do_cv = TRUE,
  seed = 42
)

summary(result)
plot(result$screen)
```

The pipeline runs:

```text
prepare -> stable selection -> fit -> nested spatial CV
        -> evaluate -> predict -> ensemble -> projection
```

Selection is repeated inside every outer spatial training fold. Held-out folds
never influence variable choice or tuning.

## Main functions

| Stage | Functions |
|---|---|
| Study design | `cast_study_area()`, `cast_background()` |
| Preparation | `cast_prepare()`, `get_env_vars()`, `cast_vif()` |
| Stable selection | `cast_domains()`, `cast_select()` |
| Modelling | `cast_fit()`, `cast_esm()` |
| Validation | `cast_cv()`, `cast_evaluate()` |
| Prediction | `cast_predict()`, `cast_predict_tiled()` |
| Ensemble/projection | `cast_ensemble()`, `cast_project()` |
| Batch/resume | `cast_batch()`, `cast_batch_resume()` |

## Model backends

- Random Forest: `ranger`
- BRT: `gbm`
- MaxEnt: `maxnet`
- GAM: `mgcv`
- stable shortlisting: `glmnet` + `ranger`

Optional functionality is activated only when its suggested package is
installed. Package installation is never attempted during a model run.

## Installation

```r
pak::pak("EldonQ/castSDM")
```

For the full workflow:

```r
install.packages(c(
  "glmnet", "ranger", "gbm", "maxnet", "mgcv", "pROC",
  "ggplot2", "sf", "terra", "future", "future.apply", "yaml"
))
```

## Interpretation

- Selected variables are stable predictive-driver candidates.
- A passed invariance test is evidence of transferability, not proof of a
  manipulable causal effect.
- Future projections assume that the learned response relationship remains
  applicable under the projected environment.

## License

GPL (>= 3)
