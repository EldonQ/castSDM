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

## Variable selection is a spatial forward search

`cast_select()` is a two-stage screen:

1. **Collinearity thinning.** Rank predictors by a univariate
   quadratic-logistic signal, then greedily keep predictors whose pairwise
   correlation with all kept predictors is at most 0.7 (Dormann et al. 2013).
2. **Forward selection on inner spatial-CV loss.** Starting from prespecified
   predictors (if any), repeatedly admit the survivor that most improves the
   held-out loss of a probability random forest — spatial folds when
   coordinates exist — and stop when no candidate passes a paired 2-SE
   improvement guard. A redundant proxy adds nothing once its parents are in,
   so it is never admitted. There is no predictor-count cap and no fallback
   set: an empty selection honestly means nothing beat the intercept-only
   model (Meyer et al. 2018, 2019).

Stopping is performance-driven, not count-driven. The screen serves
parsimony for interpretation and projection; it does not identify a causal
adjustment set.

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
plot(result$screen)          # forward-selection path: admission steps and gains
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
  exposure/adjustment set from thinning and from the forward search: kept
  predictors enter at step 0 no matter what the data say.
  `cast(select_keep = ...)` also protects it inside nested CV. Optional
  candidates must still be scientifically admissible. `kept_by_design` is not
  statistical evidence, and screening does not identify a sufficient
  adjustment set.
- Selected variables are ordered by forward-admission step.
  Selection is a model-based screen, not a list of causes.
- Future projections assume the learned response relationship remains
  applicable under the projected environment.
- `cast_report_odmap()` renders the analysis settings as an ODMAP-aligned
  report (Zurell et al. 2020).

## References

Dormann, C. F. et al. (2013). Collinearity: a review of methods to deal with
it in ecological studies. *Ecography* 36: 27-46.

Meyer, H. et al. (2018). Improving performance of spatio-temporal machine
learning models using forward feature selection and target-oriented
validation. *Environmental Modelling & Software* 101: 1-9.
DOI: 10.1016/j.envsoft.2017.12.001.

Meyer, H. et al. (2019). Importance of spatial predictor variable selection
in machine learning applications — moving from data reproduction to spatial
prediction. *Ecological Modelling* 411: 108815.
DOI: 10.1016/j.ecolmodel.2019.108815.

## License

GPL (>= 3)
