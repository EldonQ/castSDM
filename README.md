# castSDM

**Environmental drivers and spatial attribution of habitat change.**

中文入口：[castSDM回答什么、四张结果图怎么读](tools/workflows/START_HERE_ZH.md)。

castSDM asks where an environmental change improves or reduces predicted
habitat suitability, and how individual drivers contribute to future habitat
change. Use one fitted SDM, paired current/future environments and frozen model
weights throughout each comparison.

```r
# Models can be fitted in another tool; supply its prediction method.
fit <- cast_as_fit(list(gam = fitted_gam), training_predictors,
  function(model, newdata) predict(model, newdata, type = "response"))
attribution <- cast_attribution(fit, current, future)
plot(attribution, contrast = "retained")
plot(attribution, type = "map", driver = "temperature")
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
intervals. Presence/background predictions are relative suitability; these
comparisons do not identify ecological causal effects or realized colonisation.

Existing fitting, spatial validation, projection and single-shift effect APIs
remain available. Mature engine packages perform model fitting.
`cast_as_fit()` accepts an explicit prediction callback; `cast_cv(fold_ids = ...)`
accepts spatial groups designed with tools such as blockCV. External fits need
an explicit training reference and comparable suitability scales; they do not
automatically inherit valid evaluation scores.

The [research workflows](tools/workflows/README.md) include the complete
Schizopygopsis case configuration, the independent method-validation protocol
and R figure scripts. The [validation results](tools/workflows/core_method_validation/RESULTS_20261016.md)
retain negative cases and model-dependent limitations.
The [development review](tools/workflows/POSITIONING_AND_CODE_AUDIT_20261016.md) records
the repaired code paths, remaining scientific limits and retained assets.
The maintained case report has four main figures: model evaluation, habitat
change, environmental driver comparisons, and range/model disagreement.
Detailed optional diagnostics remain available without filling the default report.
Schizopygopsis is a genus-level workflow demonstration. The separate
[species validation plan](tools/workflows/core_method_validation/ECOLOGICAL_VALIDATION_PLAN_20261016.md)
records candidate data quality and the independent ecological evidence still needed.

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

Two vignettes ship with the package: the getting-started walkthrough of the
default pipeline, and *Beyond the defaults* (`vignette("castSDM-advanced")`)
covering hyperparameter tuning, prevalence weighting, background designs and
replicate sets, replicated spatial CV, response curves, and biomod2-compatible
change coding.

## Interpretation

- Driver-only maps answer: how would suitability change if this predictor took
  its future values while the others stayed current? Retained-driver maps ask
  how much the complete future differs when that predictor stays current.
- These are comparisons of a fitted response. Correlated predictors, missing
  ecological processes and observation bias can change their meaning; the
  software does not establish which environmental factor caused observed habitat change.
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
