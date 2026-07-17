# castSDM 0.6.0

## New core

* Added `cast_domains()` for data-derived spatial or environmental domains.
* Added `cast_select(method = "stable")`, the causal-inspired CAST selector.
  It combines response-aware shortlisting, binomial-GLM invariance tests, and
  a greedy minimum stable-set search.
* Added `stable_no_invariance` as an explicit ablation for validation.
* `cast_cv()` now re-runs selection inside every outer training fold and stores
  fold-specific selections, preventing feature-selection leakage.

## Simplification

* Removed DAG learning, Markov-blanket selection, causal role specifications,
  refutation utilities, CATE, and their plotting/classes/dependencies.
* Removed causal-role labels from selection outputs. Stable variables are
  candidates for transferable prediction, not automatically identified causes.
* Default selection in `cast()`, `cast_batch()`, and configuration profiles is
  now `stable` with a 12-variable shortlist ceiling.

## Existing workflow

* Retains data preparation, background sampling, RF/BRT/MaxEnt/GAM and ESM,
  tiled raster prediction, ensembles, CMIP6 projection, batch checkpointing,
  YAML configuration, and plotting.

# castSDM 0.2.0--0.5.1

Earlier experimental releases explored DAG, role-prior, CATE, and RF-only
screening designs. These interfaces are retired in 0.6.0.
