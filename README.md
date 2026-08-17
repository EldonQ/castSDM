# castSDM

`castSDM` solves one problem in species distribution modelling: **separating
true ecological drivers from collinear bystanders**. Standard variable
selection — correlation filters, stepwise VIF, univariate tests, permutation
importance, or the embedded `covsel` screen of N-SDM — ranks predictors by
*marginal* association or pairwise collinearity. In the highly collinear
environmental stacks typical of SDMs (bioclimatic layers often share R² > 0.98
of their variance), a marginal screen cannot tell a driver from a collinear
proxy: the proxy scores high on marginal association and is retained even
though it carries no information once the driver is in the model. That
misattributes ecological drivers and degrades transferability under climate
change.

Its distinguishing method is `cast_select(method = "cpi")` (the default),
a conditional predictive impact (CPI) selector that asks the conditional
question instead. For each candidate predictor it replaces the variable with
a Gaussian knockoff given every other predictor and measures the loss of
predictive accuracy, then keeps the predictors whose conditional contribution
survives Benjamini-Hochberg FDR control. An optional double machine-learning
selector (`method = "dml"`) reports signed Neyman-orthogonal partial-linear
effects with confidence intervals via [cast_effect()], and
[cast_counterfactual()] maps single-driver what-if shifts on the current
climate.

The screen is embedded in a complete workflow — nested spatial
cross-validation, ensemble fitting, projection, and counterfactual mapping —
so the selected set can be used end to end, and every retention decision is
recorded in an auditable, predictor-level table. The only ecological choice is
the FDR level (`alpha`). Cross-fitting folds and the random-forest nuisance
learner are method defaults, so the selector avoids the many hand-tuned knobs
of traditional screening. castSDM targets interpretable, conditional driver
selection but does not claim automatic causal discovery from observational
data.

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

cf <- cast_counterfactual(result$fit_full %||% result$fit, # g-computation what-if
  newdata = current_grid,    #   climate: shift one driver, hold
  variable = "bio1", shift = 1   #   the rest fixed, map the change
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
| Replicates & uncertainty | `cast_rep()`, `hss_sd` layers, `delta_sd` |
| Reporting | `cast_report_odmap()` |
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

- Selected variables are FDR-significant conditional drivers, adjusted for
  the remaining predictors; scores flag `fallback` (kept via the `min_vars`
  floor) and `forced` (kept via `force_include`) retentions.
- Replicate runs (`cast_rep()`) quantify the sensitivity of metrics,
  selection, and maps to the arbitrary choice of background points; ensemble
  predictions carry a cross-model `hss_sd` uncertainty layer.
- Batch runs route species with few presences to the ESM pipeline
  automatically (`min_occ` / `esm_min` gates) and flag the result with
  `esm_used = TRUE`.
- The DML effect is an association purged of the measured confounders; it is
  not proof of a manipulable causal mechanism.
- Counterfactual maps are purely interpretive what-if summaries on the current
  climate; they do not extrapolate to future scenarios.
- Future projections assume that the learned response relationship remains
  applicable under the projected environment.
- `cast_report_odmap()` renders the analysis settings as an ODMAP-aligned
  report (Zurell et al. 2020).

## License

GPL (>= 3)
