# castSDM 0.6.0

## New core

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

# castSDM 0.2.0--0.5.1

Earlier experimental releases explored DAG, role-prior, CATE, RF-only, and
stable/invariance screening designs. These interfaces are retired in 0.6.0.
