# castSDM

`castSDM` is a species distribution modelling package built around one
question: **if a driver moved, what would happen to predicted suitability, and
where is that answer actually supported by data?**

Two products carry that:

- **Shift effects.** `cast_effect_table()` and `cast_effect_map()`
  apply a single raw-unit shift to one driver while holding every other
  driver at its observed value (for example `bio01 + 2` degrees), and
  report the resulting change in predicted suitability as a magnitude and a
  direction, per driver and per input raster cell.
- **Hard range-masking.** Rows or cells outside the training quantile box —
  at baseline or after shifting — are masked, not ranked. Masked drivers
  report `NA` estimates. High coverage does not establish joint positivity
  or causal identification.

## Variable selection selects on the conditional effect

`cast_select()` is a two-stage screen:

1. **Collinearity thinning.** Rank predictors by a univariate
   quadratic-logistic signal, then greedily keep predictors whose pairwise
   correlation with all kept predictors is at most 0.7 (Dormann et al. 2013).
2. **Conditional effect against a conditional null, capped by sample
   size.** Fit a probability random forest on the survivors and measure how
   far each predictor moves the fitted probability **when it is shifted
   while every other predictor stays at its observed value**. Recompute the
   same statistic on forests refitted to within-stratum permuted predictors
   (strata group rows with similar values on the other survivors, one
   stratification per predictor) to build a feature-wise null, keep
   predictors with Monte Carlo tail probability <= 0.05, and cap the set at
   `ceiling(log2(n_presence))` predictors (at most 12), ordered by the
   effect.

Shuffling only among comparable rows keeps the null inside the observed joint
support, unlike marginal permutation, which is governed by the model's
extrapolation behaviour (Hooker, Mentch & Zhou 2021). The reported
p-values are exploratory scores, not confirmatory tests with guaranteed error
control: the response-dependent stage-1 filter is not rerun under the null.
This screen does not identify a causal adjustment set.

Selection is re-run inside every outer spatial training fold, so held-out folds
never influence variable choice or tuning.

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
plot(result$screen)          # conditional effect vs the within-stratum null
```

Attribution (give ecologically meaningful raw-unit shifts where possible):

```r
eff <- cast_effect_table(result$fit, newdata = current_grid,
                         shift = list(bio01 = 2))
eff                                  # rank supported drivers; masked read NA

emap <- cast_effect_map(result$fit, current_stack,
                        shift = list(bio01 = 2))   # dHSS_* / absdHSS_* / support_*
```

The pipeline runs:

```text
prepare -> two-stage selection -> fit -> nested spatial CV -> evaluate
        -> predict -> ensemble -> projection -> shift attribution
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
| Shift effects | `cast_effect_table()`, `cast_effect_map()` |
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

- Effect tables and maps are **model-based shift
  contrasts** (g-computation / standardization) for one stated raw-unit
  shift. A causal interpretation needs consistency of that shift, a
  justified adjustment set, no uncontrolled confounding, joint
  positivity, and adequate response and observation models. These assumptions
  are not established by the software. The estimator is not doubly robust.
- Masking reports the fraction of evaluated rows or cells whose baseline
  and shifted predictor vectors both fall inside the training quantile box;
  masked units carry `NA` estimates instead of extrapolated rankings.
  This range diagnostic cannot detect holes inside the box or establish joint
  positivity. `min_support` (default 0.5) is a descriptive coverage rule,
  not causal estimability.
- Holding other predictors fixed does not guarantee a feasible intervention.
  In particular, mediators, colliders and deterministically linked predictors
  need scientific treatment rather than automatic inclusion or selection.
- On presence/background data the scale is **relative suitability, not
  occurrence probability**.
- `cast_select(keep = c("exposure", "confounder"))` protects a prespecified
  exposure/adjustment set from thinning, null screening and the predictor cap;
  `ncov` must accommodate that set. `cast(select_keep = ...)` also protects it
  inside nested CV. Optional candidates must still be scientifically admissible.
  `kept_by_design` is not statistical evidence, and screening does not identify
  a sufficient adjustment set.
- Selected variables are ranked by conditional effect.
  Selection is a model-based screen, not a list of causes.
- Future projections assume the learned response relationship remains
  applicable under the projected environment.
- `cast_report_odmap()` renders the analysis settings as an ODMAP-aligned
  report (Zurell et al. 2020).

## References

Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal with
it in ecological studies. *Ecography* 36: 27-46.

Hooker, G., Mentch, L. & Zhou, S. (2021). Unrestricted permutation forces
extrapolation: variable importance requires at least one more model, or there
is no free variable importance. *Statistics and Computing*
31: 82. DOI: 10.1007/s11222-021-10057-z.

## License

GPL (>= 3)
