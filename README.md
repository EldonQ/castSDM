# castSDM

`castSDM` is a species distribution modelling package built around one
question that conventional SDM pipelines leave unanswered: **which
environmental drivers can the data actually support attributing habitat
change to, and where?**

Two products carry that:

- **Interventional effect analysis.** For every driver, `cast_effect_table()`
  and `cast_effect_map()` shift the driver while holding every other driver
  fixed (a do-intervention, evaluated by g-computation over the fitted
  ensemble) and report the resulting change in habitat suitability, as a
  magnitude and a direction, per driver and per 1 km cell.
- **A necessity audit.** `cast_necessity()` refits the model without each
  driver inside every spatial training fold and measures the held-out AUC it
  costs. Sensitivity without necessity does not identify a driver: a
  predictor the model responds to strongly can still be freely replaceable by
  a collinear partner. Read the pair; where they disagree, attribution to
  that driver is not identified and should be reported as such.

Variable selection is a separate concern, handled by `cast_select()` and
scoped to what selection can honestly deliver: **parsimony and projection
robustness, not causal attribution.** It is a two-stage, literature-standard
procedure:

1. **Collinearity thinning.** Rank predictors by absolute marginal
   association, then greedily keep predictors whose pairwise correlation with
   all kept predictors is at most 0.7 (Dormann et al. 2013).
2. **Importance filter against a permutation null.** Fit a probability random
   forest on the survivors, refit it on permuted responses to build the null
   distribution of permutation importance, and keep predictors above the 95th
   percentile of that null (Altmann et al. 2010). Permutation importance
   under the null is centred on zero, not bounded by it, so a bare
   `importance > 0` rule retains roughly half of all uninformative
   predictors; the null quantile is the calibrated replacement.

Selection is re-run inside every outer spatial training fold, so held-out
folds never influence variable choice or tuning.

## Core workflow

```r
library(castSDM)

result <- cast(
  species_data,
  env_data = prediction_grid,
  models = c("rf", "brt", "maxent", "gam"),
  do_cv = TRUE,
  seed = 42
)

summary(result)
plot(result$screen)          # stage-1 status and stage-2 permutation null
```

Attribution: the effect/necessity pair.

```r
eff <- cast_effect_table(result$fit, newdata = current_grid)
eff                                  # rank by mean_abs_dHSS, read the sign

nec <- cast_necessity(species_data, screen = result$screen, k = 5)
nec                                  # mean_dAUC lost by dropping each driver

emap <- cast_effect_map(result$fit, current_stack)   # dHSS_* / absdHSS_*
```

Single-driver what-if maps on the current climate:

```r
cf <- cast_sensitivity(result$fit,
  newdata = current_grid,
  variable = "bio1", shift = 1)
plot(cf, basemap = "china")
```

The pipeline runs:

```text
prepare -> two-stage selection -> fit -> nested spatial CV -> evaluate
        -> predict -> ensemble -> projection -> attribution audit
```

## Main functions

| Stage | Functions |
|---|---|
| Study design | `cast_study_area()`, `cast_background()` |
| Preparation | `cast_prepare()`, `get_env_vars()`, `cast_vif()` |
| Variable selection | `cast_select()`, `cast_importance()` |
| Modelling | `cast_fit()` |
| Validation | `cast_cv()`, `cast_evaluate()` |
| Prediction | `cast_predict()`, `cast_predict_tiled()` |
| Ensemble/projection | `cast_ensemble()`, `cast_project()` |
| Attribution | `cast_effect_table()`, `cast_effect_map()`, `cast_necessity()`, `cast_sensitivity()` |
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
  "ranger", "gbm", "maxnet", "mgcv",
  "pROC", "ggplot2", "sf", "terra", "future", "future.apply"
))
```

## Interpretation

- Selected variables are a parsimonious, projection-robust predictor set.
  They are not a list of causes; use the effect/necessity pair for that.
- Effect tables and maps are model-based interventional estimates. They
  assume no unobserved confounding and are not proof of a manipulable
  mechanism.
- A large `mean_abs_dHSS` with `mean_dAUC` near zero means the driver is
  substitutable: report it as unidentified rather than as a driver.
- Ensemble predictions carry a cross-model `hss_sd` uncertainty layer.
- Sensitivity maps are interpretive what-if summaries on the current
  climate; they do not extrapolate to future scenarios.
- Future projections assume the learned response relationship remains
  applicable under the projected environment.
- `cast_report_odmap()` renders the analysis settings as an ODMAP-aligned
  report (Zurell et al. 2020).

## References

Altmann, A., Tolosi, L., Sander, O. & Lengauer, T. (2010). Permutation
importance: a corrected feature importance measure. *Bioinformatics* 26:
1340-1347.

Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal with
it in ecological studies. *Ecography* 36: 27-46.

## License

GPL (>= 3)
