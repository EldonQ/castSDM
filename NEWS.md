# castSDM 0.1.0

## New features

* Core pipeline: `cast()` for one-step causal SDM modeling.
* Modular functions: `cast_dag()`, `cast_ate()`, `cast_screen()`,
  `cast_roles()`, `cast_features()`, `cast_fit()`, `cast_evaluate()`,
  `cast_predict()`, `cast_cate()`.
* S3 classes with `print()`, `summary()`, and `plot()` methods.
* VIF-based collinearity screening via `cast_vif()`.
* CI-MLP neural network with residual blocks, focal loss, and cosine
  annealing LR schedule.
* Four DAG structure learning methods: bootstrap HC, PC, BiDAG BGe,
  and NOTEARS linear.
* E-value sensitivity analysis via `cast_ate()` (`compute_evalue = TRUE`).
* Backdoor criterion check via `cast_roles()`.
* SHAP value computation for XGBoost models via `cast_shap_xgb()`.
* MC Dropout uncertainty estimation for CI-MLP models.
* Spatial CV support in `cast_cv()`.
* Batch multi-species workflow via `cast_batch()`.
* Publication-quality spatial heatmaps (HSS & CATE) via helper functions

  in `inst/examples/cast_spatial_heatmap_helpers.R`.

## Bug fixes

* CI-MLP training loop: fixed torch dataloader iterator returning a
 symbol instead of a list on exhaustion.
* ATE multiple testing correction changed from Bonferroni to FDR (BH).
* CATE heatmap: HSS mask now applied after spatial interpolation to
  prevent nearest-neighbor infill into non-habitat areas.
* Consistency matrix: NA values displayed as "NA" on grey background
  instead of "0.000".
* `cast_dag()`: suppressed benign bnlearn v-structure warnings.
