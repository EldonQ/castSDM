# castSDM 0.7.0 (development)

## Defect fixes and engineering hardening

* `plot.cast_select()` labels the DML axis `|DML statistic|` instead of the
  wrong `|CPI statistic|`; fold maps extend their palette beyond 10 folds.
* `cast_project()` marks cells with non-binary predictions as `NA` instead of
  an empty-string change class; summary counts ignore them explicitly.
* `cast_vif()` rejects non-finite predictors with a targeted error and fails
  loudly instead of looping when VIF cannot be computed.
* `cast_ensemble()` excludes models with non-finite predictions (with a
  warning) and renormalises weights, instead of silently zeroing them.
* Removed the dead `threshold_method` argument from `cast_ensemble()`,
  `cast_project()`, `cast_ensemble_raster()`, `cast_project_raster()`, and
  `cast_save_future_projection()`; binary thresholds always follow the
  max-TSS rule.
* `load_basemap()` restores the user's global `sf` s2 setting after reading.
* `cast_cv()` checks for `pROC` with a friendly error and isolates per-fold
  metric failures.
* `cast_select()` reports the CPI/DML 200-tree cap instead of truncating
  silently; the batch pipeline falls back to `"cpi"` (not `"dml"`).
* `cast_predict()` errors with the missing column names when a prediction
  grid lacks fitted predictors.
* CPI confidence intervals are clamped at zero (CPI is non-negative).
* Batch resource logging writes one CSV per species and merges them after
  the run, removing the parallel-append race.
* Documentation alignment: README, package overview, and NEWS now describe
  the CPI default with the DML effect-reporting layer; N-SDM is cited as
  Adde et al. (2020).

# castSDM 0.6.0

## New core

* `cast_screen_comparison()` contrasts the castSDM conditional predictive-impact
  screen against the associational baselines researchers commonly use
  (correlation filter, stepwise VIF, univariate marginal screen, and
  random-forest permutation importance) on the same shared candidate pool, with
  a `plot()` tile matrix that makes the conditional-vs-associational contrast
  explicit for the Results narrative.
* `cast_select(method = "dml")` is the new default selector. It fits a
  double machine-learning partially linear model (`DoubleML` + `mlr3`
  random-forest nuisance) per candidate predictor and keeps those whose
  Neyman-orthogonal effect survives Benjamini-Hochberg FDR control. The FDR
  level (`alpha`) is the only ecological tuning choice.
* `method = "rf"` is retained as a conventional permutation-importance
  benchmark for comparison.
* Added `cast_effect()`: a tidy DML effect table with FDR-adjusted
  significance and confidence intervals, plus a `plot()` coefficient (forest)
  plot.
* Added `cast_counterfactual()`: a g-computation what-if on the current
  climate that shifts one predictor, holds the rest fixed, and maps the change
  in habitat suitability, plus a diverging-map `plot()` method.
* `cast_cv()` re-runs selection inside every outer training fold and stores
  fold-specific selections, preventing feature-selection leakage.

## Removed

* Removed the `"stable"` selector and its invariance machinery (data-derived
  environments, binomial-GLM invariance tests, greedy minimum stable-set
  search), the `stable_no_invariance` ablation, `cast_domains()`, and the
  `glmnet` dependency. Variable choice is now the single FDR-controlled DML
  criterion.

## Existing workflow

* Retains data preparation, background sampling, RF/BRT/MaxEnt/GAM and ESM,
  tiled raster prediction, ensembles, CMIP6 projection, batch checkpointing,
  YAML configuration, and plotting.

## Plotting and stability

* China maps always show the complete, standard national frame plus the South
  China Sea inset; a range-restricted species is never cropped to its coloured
  cells. Seamless colour is achieved instead by filling the whole basemap with
  each map's null/background colour, so cells outside the accessible area
  continue the surface rather than tearing it at an artificial edge.
* Suitability maps (`plot.cast_predict()`, `plot.cast_ensemble()`) use a pale
  grey to deep teal-blue sequential ramp whose low anchor equals the basemap
  fill; range-change maps (`plot.cast_project()`) use a refined gain/loss/
  stable-present palette over a neutral `stable_absent` background; what-if maps
  (`plot.cast_counterfactual()`) diverge around a neutral grey midpoint.
* `plot.cast_select()` gains a `top` argument and now colours predictors by
  ecological class with readable labels; `plot.cast_cv()` fold maps share the
  full-China frame and carry a spatial-block narrative.
* Ensemble/counterfactual raster prediction is chunked by rows to bound peak
  memory on national 1 km grids, preventing out-of-memory segfaults.

# castSDM 0.2.0--0.5.1

Earlier experimental releases explored DAG, role-prior, CATE, RF-only, and
stable/invariance screening designs. These interfaces are retired in 0.6.0.
