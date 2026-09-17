# castSDM

`castSDM` is a species distribution modelling package built around one
question: **if a driver moved, what would happen to predicted suitability, and
where is that answer actually supported by data?**

Three products carry that:

- **Interventional effects.** `cast_effect_table()` and `cast_effect_map()`
  shift one driver while holding every other driver at its observed value, and
  report the resulting change in predicted suitability as a magnitude and a
  direction, per driver and per input raster cell.
- **Response shape.** `cast_dose_response()` sweeps the size of the shift, so a
  saturating or threshold response stays visible instead of collapsing to one
  number.
- **A positivity diagnostic.** `cast_effect_support()` reports the fraction of
  observed predictor vectors that are still inside the training support after
  the shift. A shift answered mostly by extrapolation is reported as such,
  never silently answered anyway.

## Variable selection selects on the interventional effect

`cast_select()` is a two-stage screen:

1. **Collinearity thinning.** Rank predictors by a univariate
   quadratic-logistic signal, then greedily keep predictors whose pairwise
   correlation with all kept predictors is at most 0.7 (Dormann et al. 2013).
2. **Interventional effect against a permutation null, capped by sample
   size.** Fit a probability random forest on the survivors and measure how
   far each predictor moves the fitted probability **when it is shifted
   while every other predictor stays at its observed value**. Recompute the
   same statistic on forests refitted to a permuted response to build a
   feature-wise null, keep predictors with Monte Carlo tail probability
   <= 0.05 (Altmann et al. 2010), and cap the set at
   `ceiling(log2(n_presence))` predictors (at most 12), ordered by the
   effect.

Why the second stage is a shift rather than a permutation: permuting a
predictor breaks its correlation with every other predictor, so for collinear
predictors the permuted rows leave the observed data support and the score is
governed by the model's extrapolation behaviour rather than by the predictor's
influence (Hooker, Mentch & Zhou 2021). Shifting one predictor and holding the
rest at their observed values keeps the contrast inside the support and answers
a question that has a definition.

The same forest's permutation importance is still reported, as a diagnostic.
`screen$diagnostics$importance_agreement` gives the Spearman agreement between
the two rankings; large disagreement marks predictors the forest leans on but
barely responds to when changed — the signature of a collinear stand-in.

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
plot(result$screen)          # interventional effect vs the permuted-response null
```

Attribution:

```r
eff <- cast_effect_table(result$fit, newdata = current_grid)
eff                                  # rank by mean_abs_dHSS, read the sign

emap <- cast_effect_map(result$fit, current_stack)   # dHSS_* / absdHSS_*

dr <- cast_dose_response(result$fit, "bio1", shift = seq(-3, 3, by = 0.5))
plot(dr)                             # response shape, hollow points = off support

sp <- cast_effect_support(result$fit, shift = c(1, 2, 3))
plot(sp)                             # positivity by driver and shift size

nec <- cast_necessity(species_data, screen = result$screen, k = 5)
nec                                  # mean_dAUC lost by dropping each driver
```

The pipeline runs:

```text
prepare -> two-stage selection -> fit -> nested spatial CV -> evaluate
        -> predict -> ensemble -> projection -> interventional attribution
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
| Interventional effects | `cast_effect_table()`, `cast_effect_map()`, `cast_dose_response()`, `cast_effect_support()`, `cast_sensitivity()`, `cast_necessity()` |
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

- Effect tables, maps and curves are **model-based interventional
  contrasts** (g-computation / standardization). They assume consistency, no
  unobserved confounding given the adjustment set, and positivity; the support
  column reports the third. They are not doubly robust and not TMLE.
- The `support` column is a marginal-quantile (hyper-rectangle) approximation
  to the joint support. It flags gross off-support shifts, not subtler holes
  inside the observed box.
- On presence/background data the scale is **relative suitability, not
  occurrence probability**.
- Selected variables are a parsimonious set ranked by interventional effect.
  Selection is a model-based screen, not a list of causes.
- A large `mean_abs_dHSS` with small `mean_dAUC` shows model response with
  limited RF knockout cost. It does not prove substitutability or its cause.
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

Hooker, G., Mentch, L. & Zhou, S. (2021). Unrestricted permutation forces
extrapolation: variable importance requires at least one more model, or there
is no free variable importance. *WIREs Data Mining and Knowledge Discovery*
11: e1421.

## License

GPL (>= 3)
