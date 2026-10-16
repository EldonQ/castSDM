# castSDM 0.19.0: scope, repairs and remaining scientific limits

Historical review. For current repair status and source comparisons, read
[POSITIONING_AND_CODE_AUDIT_20261016.md](POSITIONING_AND_CODE_AUDIT_20261016.md).
The numerical/QA claims below refer to the receipts from 2026-10-15.

## One scientific question

Which environmental changes improve or reduce predicted habitat suitability,
where do they act, and how reproducible are those spatial conclusions?

The public centre is `cast_attribution()` and its source-native plot method.
It accepts an existing `cast_fit` wrapper and paired tables or aligned rasters.
The supported engine packages still perform fitting; existing fitting and CV
APIs are retained for compatibility. This release is a focused extension and
cleanup, not a claim that every legacy modelling API has been removed.

## Findings by severity

| Severity | Location | Cause and impact | Action / remaining work |
|---|---|---|---|
| High | Case `01_config.R`: genus records, national background, 300 km buffer | The observation unit is a genus-level literature collection; land-background and buffer cells do not define fish-accessible river habitat. Strong discrimination can reflect domain separation and record bias. | Preserve the requested data and explicitly label the case. For a species claim, curate species-level records and define a basin/river sampling and projection domain. |
| High | `R/scenario-attribution.R::.cast_attribution_block()` | Single-driver replacement measures the fitted response. Correlation, unmeasured processes and observation bias prevent automatic ecological causal identification. | State both comparisons and their baselines. Report total change, driver-only and retained-driver effects, nonadditive remainder, range flags and model disagreement. |
| High | Case future alignment and ensured-climate selection | Only one of the nine selected predictors (`bio05`) has configured future values. Static covariates do not explain temporal change, and ensured climate coverage is a declared selection constraint. | Keep static-driver effects exactly zero and report the changing predictor. Broader driver attribution needs defensible future hydrology, land-use or other layers; do not invent them. |
| High | Core validation `01_attribution_validation.R` | Local effects inherit model error. Forest-seed repeatability can coexist with weak repeatability under new observations. | Complete six prespecified scenarios and 20 independent datasets each. Retain null, joint, correlated and out-of-range failures. Quantify sampling, seed and model-family differences separately. |
| Medium | `R/effect-products.R`, `R/project-future.R` | Different ensemble settings across projection and explanation would change the target response. | Single-shift effects accept named model weights; projections propagate the score exponent. Paired attribution fixes active models and weights across comparisons. |
| Medium | Case `02_pipeline.R`: thinning, mask, input manifests and CV gate | Logical thinning indices, predictor-derived mask zeros, stale caches or swallowed fold defects can silently change rows or surfaces. | Correct row indexing and audit exclusions; rasterize a constant polygon field; fingerprint source/input state; distinguish failed folds from single-class data folds and enforce coverage. |
| Medium | `R/background-sample.R` and raster projection | Scanning hundreds of nationwide layers to draw background points and transporting in-memory stacks wastes I/O and memory. | Uniform complete-case proposals without replacement; bounded raster blocks; file-backed crops; proxy transport; caller-owned future plan; one forest thread per worker. Seeded background sequences change and are documented. |
| Medium | `R/ensemble-raster.R::cast_ensemble_raster()` output buffers | The block budget limits predictor input, but full-grid output vectors remain allocated for suitability, disagreement, classes and diagnostics. Parallel worker memory therefore exceeds the advertised block budget. | File-backed inputs reduce transport costs in this release; measure process memory when choosing workers. Replace output vectors with streamed sibling GeoTIFFs in a later isolated change, with numerical and interrupted-write checks. |
| Medium | `R/effect-products.R::cast_effect_map()` output buffers | Row-block prediction still retains full-grid signed, absolute and support vectors for every driver. The legacy effect map therefore also exceeds a predictor-only block budget. | The recorded case uses disjoint row partitions with the unchanged kernel, an actual-model equivalence check and a no-resampling merge. Stream the installed function's output in a separate source change and repeat numerical validation. |
| Low | `R/project-future.R::cast_project_raster()` centroid fields | `current_centroid_lon/lat` and future counterparts carry native x/y units on projected rasters. They are metres on this Albers case, despite the field names. | Distances already use projected units; do not interpret the fields as longitude/latitude. Add explicit native-coordinate columns and CRS metadata while retaining compatibility in the next numerical revision. |
| Medium | `R/effect-products.R::cast_effect_map()` | Predicting rows that the hard range mask immediately discards wastes most computation. | Determine the existing support mask first. Predict the baseline only on the union of retained rows and each shifted surface on its supported rows; preserve unknown coverage. |
| Medium | Package cartography and case `03_figures.R` | Independent plotting implementations can drift in coordinate systems, colour meanings and handling of missing values. | Share native map helpers, preserve `inst`, reproject for display only, retain native-grid statistics, record colour clipping and every export/source. Successful export remains distinct from visual inspection. |
| Low | Active research runners | Multiple historical selection/re-export runners make the maintained path ambiguous. | Archive superseded runners before removing them from active directories. Maintain one case runner and one attribution-validation runner; preserve historical outputs as historical evidence. |

## What this release establishes

The benchmark verifies numerical implementation and measures a conditional
advantage over precisely matched global marginal PDP and uniform mean-shift
baselines. It does not benchmark all ALE/SHAP or competing tool packages.
The inspected biomod2, iml and DALEX sources inform the package boundary;
predictor replacement itself is established functionality.

The contribution is a reproducible ecological scenario workflow: actual future
environmental values, spatial driver comparisons, fixed ensemble target,
nonadditive remainder and accompanying diagnostics. Empirical ecological
mechanism and future transfer accuracy require independent observations.

Core evidence: 120 independent datasets, 492 numerical checks, zero run errors,
four source-data-linked figures; short R CMD check: zero errors and warnings.
The final documentation check reports one external-time-verification NOTE
(unable to verify current time); the preceding source check had zero notes.
Current run receipts and source manifests identify the exact evidence.

## Development safeguards

- Package cartography and citation assets are retained. The previously removed
  example CSV was restored on 2026-10-16. The vignette uses explicitly synthetic
  demonstration data; neither demonstration is method-validation evidence.
- Original working-tree changes are preserved. A pre-development archive and
  an archive of retired active runners permit recovery.
- Large occurrence inputs, model caches and rasters stay external to Git.
  Maintained research sources and compact validation evidence are published
  under `tools/workflows`, outside the installed package.
- Further slimming should follow measured duplication and adapter demand.
  Avoid another framework layer: one prediction boundary, one paired numerical
  kernel and one cartographic pathway are sufficient for the present question.
