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
- **A range-extrapolation diagnostic.** `cast_effect_support()` reports the
  fraction of evaluated rows inside the training quantile box both before and
  after shifting. Effects are still returned outside this box. High coverage
  does not establish joint positivity or causal identification.

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

Permutation can create unobserved predictor combinations, making importance
sensitive to model extrapolation (Hooker, Mentch & Zhou 2021). Additive shifts
can also leave joint support. Their purpose here is to evaluate a specified
change in a predictor, not to guarantee support, causality or superiority.
The response-dependent stage-1 filter is not rerun under the null; the reported
p-values are exploratory scores, not confirmatory tests with guaranteed error
control. This screen does not identify a causal adjustment set.

The same forest's permutation importance is reported as a diagnostic.
`screen$diagnostics$importance_agreement` gives the Spearman agreement between
the two rankings. Disagreement does not establish that a predictor is a proxy.

Selection is re-run inside every outer spatial training fold, so held-out folds
never influence variable choice or tuning.

## Optional invariant selection and scenario contrasts

`method = "tramicp"` uses the optional `tramicp` package's binary-logistic
invariant causal prediction. Supply an environment **column name**, based on
study design rather than arbitrary partitions or CV folds. The column is never
a predictor. Correct model specification, valid environments, adequate measured
causes/confounders and the method's sampling assumptions remain necessary;
spatial blocking alone does not establish them.

```r
screen <- cast_select(
  training_data, method = "tramicp", environment = "environment_group", seed = 42
)
print(screen)
screen$diagnostics

cv <- cast_cv(
  training_data, select_method = "tramicp",
  select_args = list(environment = "environment_group"),
  models = "gam", k = 5, seed = 42
)

if (!length(screen$selected)) stop("No invariant predictors identified; inspect diagnostics.")
fit <- cast_fit(training_data, screen = screen, models = "gam")
scenario <- cast_sensitivity(
  fit, prediction_grid, variable = "bio01",
  shift = 1, shift_type = "raw", model = "gam", backend = "marginaleffects"
)
scenario$predictions
plot(scenario, basemap = "none")
```

`environment_group` is a placeholder for user-supplied scientific group labels;
no such column is generated automatically. The example's `bio01` must be a
retained predictor in your data. Choose the scenario variable and its raw-unit
change for ecological reasons, after checking its measurement units.
The high-level `cast()` accepts `select_method = "tramicp"` and
`select_environment = "environment_group"` with the same semantics.

Selection tests all predictor subsets, with no response-based prescreen or
fallback ranking. `icp_max_predictors = 12` is a computational guard, not a
selection cap. `ncov` and `keep` are not supported for this method. Empty
intersections, accepted empty sets, and no accepted sets are reported separately;
an empty screen stops fitting instead of silently using other predictors.
These results do not prove absence of ecological causes or provide a sufficient
causal adjustment set. Failed/empty CV folds remain in the fold diagnostics and
selection-frequency denominator; they do not contribute predictive metrics.
When all invariant-selection CV folds are unevaluable, `cast()` stops with the
same diagnostic error instead of substituting hold-out evaluation. Retired
selection arguments and unknown method names raise errors, not replacement runs.

The optional `marginaleffects` backend currently supports one explicitly chosen
GAM only. It returns baseline, shifted prediction and their rowwise difference
without standard errors. `delta_sd` is not a confidence interval. Native
multi-model contrasts remain available. Both are model-based scenario responses,
not automatically causal effects or calibrated occurrence probabilities.
Install optional backends yourself when needed; model runs never install them.

For Windows R 4.4, run the following in R to install the optional backends and
required dependencies from CRAN binaries, without compiling source packages:

```r
install.packages(
  c("tramicp", "marginaleffects"),
  repos = "https://cloud.r-project.org", type = "binary"
)
```

To verify the local source checkout, run this separately after installation
(adjust the checkout path if needed):

```r
stopifnot(
  requireNamespace("tramicp", quietly = TRUE),
  requireNamespace("marginaleffects", quietly = TRUE),
  utils::packageVersion("tramicp") >= "0.1-0",
  utils::packageVersion("marginaleffects") >= "0.31.0"
)
results <- testthat::test_local(
  "E:/Package/cast",
  filter = "selection-and-attribution|effect-consistency|regressions",
  reporter = "summary", stop_on_failure = TRUE
)
checks <- as.data.frame(results)
stopifnot(!any(checks$skipped), !any(checks$failed), !any(checks$error))
```

The dependency checks deliberately stop rather than silently skipping the real
backend test. These tests check software agreement and training-fold isolation,
not ecological cause recovery. They do not run a species workflow or modify
raw data or previous validation outputs.

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

emap <- cast_effect_map(result$fit, current_stack)   # dHSS_* / absdHSS_* / support_*

dr <- cast_dose_response(result$fit, "bio1", shift = seq(-3, 3, by = 0.5))
plot(dr)                             # hollow points = box coverage below 0.5

sp <- cast_effect_support(result$fit, shift = c(1, 2, 3))
plot(sp)                             # quantile-box coverage by driver and shift

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
  contrasts** (g-computation / standardization). A causal interpretation needs
  consistency, a justified adjustment set, no uncontrolled confounding, joint
  positivity, and adequate response and observation models. These assumptions
  are not established by the software. The estimator is not doubly robust.
- The `support` column reports the fraction of evaluated rows whose baseline
  and shifted predictor vectors both fall inside the training quantile box.
  This range diagnostic cannot detect holes inside the box or establish joint
  positivity. `range_supported` marks box coverage of at least 0.5, not causal
  estimability. Contrasts are not filtered by this threshold.
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
- The `support_<driver>` raster averages coverage across shifts at each cell;
  table support is the minimum across shifts of spatial coverage. Neither is
  a conditional-support test. Low-coverage effects are still returned.
- Selected optional variables are ranked by interventional effect.
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
is no free variable importance. *Statistics and Computing*
31: 82. DOI: 10.1007/s11222-021-10057-z.

## License

GPL (>= 3)
