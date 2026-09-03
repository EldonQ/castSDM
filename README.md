# castSDM

`castSDM` solves one problem in species distribution modelling: **causal effect analysis**: for every environmental driver, an interventional
effect size and a 1 km causal effect map that show where habitat increases or
decreases when that driver changes and all others are held fixed. Standard variable
selection 鈥?correlation filters, stepwise VIF, univariate tests, permutation
importance, or the embedded `covsel` screen of N-SDM 鈥?ranks predictors by
*marginal* association or pairwise collinearity. In the highly collinear
environmental stacks typical of SDMs (bioclimatic layers often share R虏 > 0.98
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
survives Benjamini-Hochberg FDR control. [cast_importance()] turns that screen
into a tidy per-predictor table with confidence intervals, and
[cast_sensitivity()] maps single-driver what-if shifts on the current climate.

The screen is embedded in a complete workflow 鈥?nested spatial
cross-validation, ensemble fitting, projection, and sensitivity mapping 鈥?so
the selected set can be used end to end, and every retention decision is
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

Conditional importance and sensitivity reporting:

```r
eff <- cast_importance(result)                # tidy CPI table + CIs
plot(eff)                                    # coefficient (forest) plot

cf <- cast_sensitivity(result$fit,           # what-if: shift one driver,
  newdata = current_grid,                    #   hold the rest fixed,
  variable = "bio1", shift = 1)              #   map the change in HSS
plot(cf, basemap = "china")                  # diverging map of change
```

The pipeline runs:

```text
prepare -> conditional (CPI) selection -> fit -> nested spatial CV
        -> evaluate -> predict -> ensemble -> projection
```

Selection is repeated inside every outer spatial training fold. Held-out folds
never influence variable choice or tuning.

## Main functions

| Stage | Functions |
|---|---|
| Study design | `cast_study_area()`, `cast_background()` |
| Preparation | `cast_prepare()`, `get_env_vars()`, `cast_vif()` |
| Conditional selection | `cast_select()` |
| Modelling | `cast_fit()` |
| Validation | `cast_cv()`, `cast_evaluate()` |
| Prediction | `cast_predict()`, `cast_predict_tiled()` |
| Ensemble/projection | `cast_ensemble()`, `cast_project()` |
| Conditional interpretation | `cast_importance()`, `cast_sensitivity()` |
| Reporting | `cast_report_odmap()` |

## Model backends

- Random Forest: `ranger`
- BRT: `gbm`
- MaxEnt: `maxnet`
- GAM: `mgcv`

Optional functionality is activated only when its suggested package is
installed. Package installation is never attempted during a model run.

## Installation

```r
pak::pak("EldonQ/castSDM")
```

For the full workflow:

```r
install.packages(c(
  "mlr3", "mlr3learners", "ranger", "gbm", "maxnet", "mgcv",
  "pROC", "ggplot2", "sf", "terra", "future", "future.apply"
))
```

## Interpretation

- Selected variables are FDR-significant conditional drivers, adjusted for
  the remaining predictors; scores flag `fallback` (kept via the `min_vars`
  floor) and `forced` (kept via `force_include`) retentions.
- Ensemble predictions carry a cross-model `hss_sd` uncertainty layer.
- The CPI estimate is an association purged of the measured predictors; it is
  not proof of a manipulable causal mechanism.
- Sensitivity maps are purely interpretive what-if summaries on the current
  climate; they do not extrapolate to future scenarios.
- Future projections assume that the learned response relationship remains
  applicable under the projected environment.
- `cast_report_odmap()` renders the analysis settings as an ODMAP-aligned
  report (Zurell et al. 2020).

## License

GPL (>= 3)
