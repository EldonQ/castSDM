# castSDM

**Factorial causal scenario attribution for species distribution models.**

castSDM separates two questions that are often conflated: how much do climate
and land use each contribute to a projected change, and how do individual
drivers explain the joint future? The first is answered by a 2 x 2 factorial
design; the second by the within-scenario driver-only/retained decomposition.

```r
# Four aligned suitability surfaces from one fitted SDM and frozen weights.
fa <- cast_factorial_attribution(
  current = current_hss,
  joint   = joint_hss,
  climate = climate_only_hss,
  lulc    = lulc_only_hss,
  threshold = 0.43
)
fa$summary          # climate main, land-use main, joint, interaction
fa$binary           # thresholded range counts and percentage changes
fa$contract         # claim level and explicit assumptions
plot(fa)            # contrast layers
```

The factorial contrasts are:

- `climate_main = climate - current`;
- `lulc_main = lulc - current`;
- `joint_change = joint - current`;
- `interaction = joint - climate - lulc + current`;
- `climate_conditional = joint - lulc` and
  `lulc_conditional = joint - climate`.

The other scenario component is fixed by the design, so the factorial layer
does **not** require a DAG or a user-supplied adjustment set. It is a
model-based causal scenario attribution under explicit assumptions: the
scenario components are well-defined manipulations, the fitted response is
valid for those scenarios, no omitted component changes the contrast, and the
counterfactual values are supported. It is not an identified real-world causal
effect from observational data. Continuous and binary summaries, support,
optional model-SD layers and a claim contract are returned together.

```r
# Optional within-scenario explanation, using one fitted model.
fit <- cast_as_fit(list(gam = fitted_gam), training_predictors,
  function(model, newdata) predict(model, newdata, type = "response"))
attribution <- cast_attribution(fit, current, future)
plot(attribution, contrast = "retained")
```

`effect_<driver>` changes only that predictor to its actual future values.
`retained_<driver>` compares the full future with a future where that predictor
stays current. Their difference exposes context-dependent responses; individual
effects are not additive. `total_change` is the complete future/current contrast.
The remainder also includes changes in any drivers omitted from the comparison.

Tables and aligned rasters use the same numerical implementation. Raster outputs
are written in bounded blocks; plots use the package's existing cartographic
assets. Raw comparisons and within-range summaries are reported separately.
Cross-model standard deviations describe model disagreement, not confidence
intervals. Presence/background predictions are relative suitability; the
factorial contrasts are causal scenario attributions under assumptions, not
identified ecological causal effects or realised colonisation.

Existing fitting, spatial validation, projection and single-shift effect APIs
remain available. Mature engine packages perform model fitting.
`cast_as_fit()` accepts an explicit prediction callback; `cast_cv(fold_ids = ...)`
accepts spatial groups designed with tools such as blockCV. External fits need
an explicit training reference and comparable suitability scales; they do not
automatically inherit valid evaluation scores.

The production pipeline is `cast_run()`, driven by one `cast_protocol()`:
the independent elastic-net variable screen, RF / MaxEnt / GAM engines,
nested spatial cross-validation on blockCV blocks, a one-shot hold-out
evaluation at out-of-fold thresholds, full-data refit and ensemble.
`cast_run(output_dir = ...)` writes every stage to disk with content-digest
receipts in `manifest.rds`, together with an ODMAP-aligned report.

Typical case reports use four main figures: model evaluation, habitat change,
factorial causal attribution (climate main, land-use main, interaction), and
support/model disagreement. Within-scenario driver comparisons and detailed
diagnostics remain available without filling the default report.

## Installation

```r
pak::pak("EldonQ/castSDM")
```

For the full feature set:

```r
install.packages(c(
  "ranger", "maxnet", "mgcv", "glmnet",
  "pROC", "ggplot2", "sf", "terra", "future", "future.apply"
))
```

Two vignettes ship with the package: the getting-started walkthrough of the
default pipeline, and *Beyond the defaults* (`vignette("castSDM-advanced")`)
covering hyperparameter tuning, prevalence weighting, background designs and
replicate sets, replicated spatial CV, response curves, and biomod2-compatible
change coding.

## Interpretation

- Driver-only maps answer: how would suitability change if this predictor took
  its future values while the others stayed current? Retained-driver maps ask
  how much the complete future differs when that predictor stays current.
- The factorial contrasts are model-based causal scenario attributions under
  the assumptions recorded in `fa$contract`: well-defined scenario components,
  a valid fitted response for those scenarios, no omitted component, supported
  counterfactuals and no interference. They are not identified real-world causal
  effects from observational data.
- Driver-only and retained-driver comparisons are within-scenario explanations
  of a fitted response. Correlated predictors, missing ecological processes and
  observation bias can change their meaning; they do not establish which
  environmental factor caused observed habitat change.
- Training-range flags show whether the compared environments fall within each
  predictor's training quantiles. They cannot detect every unsupported joint
  combination. Raw scenario contrasts and supported summaries are separate.
- Single-shift effect products use one declared raw-unit or training-SD shift.
  Their hard masking can produce `NA`; this records insufficient coverage.
- `cast_select(keep = ...)` preserves scientifically chosen predictors.
  Variable selection and embedded ranks describe a predictive screen, not evidence
  that an environmental variable causes the distribution.
- Future projections assume the learned response remains applicable. Static
  layers stay static unless supplied with future values. Thresholded habitat
  area additionally depends on the chosen operating threshold and spatial mask.
- `cast_report_odmap()` records the analysis settings in an ODMAP-aligned report
  (Zurell et al. 2020).

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
