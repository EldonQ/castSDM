# castSDM

**Environmental response maps and auditable scenario comparisons for SDMs.**

`cast_compare()` directly compares current and future environments with one
changing predictor group. Two groups additionally enable factorial contrasts. Fixed positive model weights and an
explicit finite domain make the comparisons auditable. Results contain
continuous contrasts, optional thresholded range counts, training-support
diagnostics, model sign disagreement and paired model SD. `cast_diagnose()`
connects independent ensemble evaluation to the scenario diagnostics.

## Installation

```r
install.packages("remotes")
remotes::install_github("EldonQ/castSDM")
install.packages(c("ranger", "maxnet", "mgcv", "glmnet", "pROC", "blockCV", "FNN", "sf", "terra"))
```

## Run from the interface

```r
library(castSDM)
cast_app()
```

The **Start** tab configures real data paths, occurrence records, analysis domain,
models and resources. **Check data** checks physical inputs without fitting;
**Run analysis** submits preparation, spatial validation, projections and the
complete available figure collection to an isolated process. A new workspace
uses the workflow shipped inside the installed package, with no source checkout.
Paths/settings remain explicit and replayable. The current national-library
adapter supports the documented CHELSA recipe; it is not a universal file importer.

The [CHELSA workflow](inst/workflows/single_species/README.md) uses baseline
1981–2010 and ensemble-mean climate for SSP126/370/585 in 2011–2040,
2041–2070 and 2071–2100. All nonclimate predictors stay current. Mean-climate
predictions do not quantify between-GCM uncertainty. Changed climate sources
require a fresh model run and separate output directory.

## Open a completed real case

The optional Shiny interface queries native-grid results, compares actual
scenario rasters or explicit hypothetical predictor shifts, and displays
descriptive environmental response groups. Set `case_dir` to your completed
[single-species workflow](inst/workflows/single_species/README.md) output.

```r
library(castSDM)
install.packages(c("shiny", "bslib", "leaflet", "callr", "rpart", "zip"))
cast_app(case_dir)

# Real saved locations only; no full-grid prediction or model refitting:
response <- cast_explore(case_dir, maps = FALSE)
response$groups$summary
```

Long jobs run in an isolated R process. The app preserves the complete existing
publication figures and puts new analysis artifacts in `interactive/`. Exports
include settings, source data, spatial tuning losses, native rasters and an R
replay script. Display rasters are explicitly separate from the analysis grid.

`cast_response_groups()` summarises finite model responses with a small tree,
using supplied spatial folds and one-standard-error pruning. The method follows
[published FME/cAME explanations](https://doi.org/10.32614/RJ-2024-024);
environmental groups are descriptive and need not be geographically contiguous.
The default driver overview now compares the same range-supported locations
for every factor. It shows responses in both backgrounds of the configured
comparison pair, their sign reversals and active-model sign agreement.
Predictors unchanged in the supplied scenario are labelled and not ranked.
Detailed native maps focus on one automatically selected or user-chosen
factor, avoiding national predictions for every factor on each interaction.
Tied first ranks do not trigger an arbitrary automatic choice.
SHAP contributions, CATE, ecological causal effects, actual occupancy and
probability changes are different quantities. The current presence/background
case does not identify them merely by changing the estimator or adding a GUI.
Marginal range flags cannot establish joint support or causal identification.

The implementation reuses the existing fitted models and contrast kernel.
`cast_attribution(..., drivers = character())` computes only the paired total
response. Optional named `shifts` modify future predictors inside processing
blocks without creating a full shifted raster stack. Related derived predictors
should be changed using a consistent environmental scenario.

GUI bundles also contain `view_state.rds` and `reproduce_view.R`, which restore
the comparison, selected layer, map view and point/profile query. Run the script
from the extracted bundle with the original real-case inputs available.

## Fit, validate, compare

`cast_run()` handles screening, spatial validation, independent hold-out
evaluation and full-data refitting. Scientific settings live in one
`cast_protocol()`. Screening and tuning repeat inside outer CV folds; hold-out
data do not determine ensemble weights or thresholds.

```r
# species_data: lon, lat, presence (0/1), and numeric predictors.
protocol <- cast_protocol(models = c("rf", "maxent", "gam"), seed = 42)
run <- cast_run(species_data, protocol = protocol,
  output_dir = "results/training", cache_dir = "cache/training")
comparison <- cast_compare(run, current_raster, future_rasters,
  groups = list(climate = climate_variables),
  mask = analysis_mask, output_dir = "results/comparison",
  resources = list(cores = 8, workers = 2, memory_mb = 4096))
diagnostics <- cast_diagnose(run, comparison,
  output_dir = "results/diagnostics")
diagnostics$holdout
diagnostics$folds
diagnostics$score_views_oof
```

Tables and aligned `terra::SpatRaster` inputs share a numerical kernel.
Unlisted fitted predictors are recorded as static. Their future layers may be
omitted; supplied static values must equal current values. Inputs must already
share grid, extent and CRS. No automatic resampling occurs.

The [single-species example](inst/workflows/single_species/README.md) supplies
occurrence QC, domain masks, exact future-file mapping and saved-output plots.
Training and comparison use public package functions.

## Complete publication figures

The complete numbered collection uses the established cartography, palettes,
province borders and South China Sea inset. It includes validation, habitat
change, driver explanation, support/disagreement, factorial contrasts,
response curves, ROC/PR/reliability, range classes, MESS and every future
scenario, plus the available supplementary layer/driver atlases.

```r
figures <- cast_figures(run, comparison, "results/figures",
  current = current_raster, future = reporting_future_raster, basemap = "china",
  resources = list(cores = 8, memory_mb = 4096))
# Change only rendering: no model refit, diagnostic prediction or PDP work.
figures <- cast_figures(run, comparison, "results/figures",
  basemap = "china", render_only = TRUE)
```

The first complete export computes missing frozen-model driver, PDP and MESS
products; this can be expensive. Verified numerical products are cached apart
from rendering. PNG (183 mm, 600 dpi by default), PDF, source CSVs and a figure
inventory/status manifest are written. Subsequent identical exports reuse
files only after checking data, installed code, rendering settings and file
hashes; missing or damaged files are regenerated. Missing inputs are reported explicitly.
Native-grid range classes use the frozen threshold before display-only
reprojection; original scientific rasters and counts remain unchanged.

Optional `formats = c("png", "pdf", "svg", "tiff")` adds editable SVG and
submission TIFF exports from the same cached layout. The default remains PNG/PDF.
Cairo/Arial is preferred; a system without a working Cairo text device uses
standard PDF/Helvetica and records that choice in the export manifest.

`run$diagnostics` is saved even without a map: per-engine and frozen-ensemble
hold-out scores, training-prevalence baseline, fold counts/locations/AOA,
selection frequency and first-repeat OOF native/shifted score views. The native
view uses the same fits and class weights. OOF comparisons are descriptive;
they never select a policy from hold-out labels. `prevalence_target` specifies
a weighted training mean score, not calibrated occurrence probabilities.
Brier/log loss and reliability bins concern the sampled binary response.
For presence/background data, background points are not verified absences.
The printed evaluation and diagnostic tables report non-interpolated average
precision (AP), whose constant-score baseline is the sampled positive fraction.
Historical `pr_auc` is retained as a trapezoidal area; its constant-score value
is `(1 + prevalence)/2`, so it is not compared with the AP baseline.
Legacy `omission_5/10` metrics use top-percent cutoffs of evaluation scores;
diagnostics also report omission at the frozen threshold and at training OOF
presence-score 5th/10th percentiles, with explicit threshold columns.

## Interpretation

Let `C` be current, `A` first-group-only, `B` second-group-only and `J` joint.
Legacy `climate` and `lulc` prefixes denote the first and second groups even
when those groups have different names.

| Field | Definition |
|---|---|
| `climate_main` | `A - C`, baseline simple effect |
| `lulc_main` | `B - C`, baseline simple effect |
| `joint_change` | `J - C` |
| `interaction` | `J - A - B + C` |
| `climate_conditional` | `J - B` |
| `lulc_conditional` | `J - A` |
| `climate_average` | `((A - C) + (J - B))/2` |
| `lulc_average` | `((B - C) + (J - A))/2` |

These are **model-conditional scenario contrasts** on a fitted suitability
surface. They do not establish ecological causation, occurrence probabilities
or realised distributions. Historical `*_main` fields mean baseline simple
effects; the averaged factorial effects have explicit `*_average` fields.

Each cell uses all positive-weight models in every corner. Model SD is the
unweighted sample SD of within-model contrasts among active models; it is
`NA` with one model and is not a confidence interval. Training support uses
marginal quantiles on unclamped inputs and cannot establish joint support.
Full-domain and supported summaries remain separate. Requested clamping
applies uniformly to all engines.

`summary_common_support` and `binary_common_support` intersect current and
**all** requested futures' finite inputs and training marginal bounds. These
use the same subset across scenarios. Existing `*_supported` summaries use
each scenario's subset. `agreement` quantifies positive/negative conflicts
and majority-sign agreement among active models (zero tolerance 1e-12).
Neither marginal support nor sign agreement is statistical confidence.
AOA is a separate multivariate validation diagnostic.

Binary predictions use `suitability >= threshold`. Binary summaries count
valid grid cells or table rows, not square kilometres. Area interpretation
requires an appropriate equal-area grid or a separate area-weighted analysis.

`cast_factorial_attribution()` accepts four already-created single-layer
surfaces or equal-length numeric vectors. Ensemble surfaces cannot recover
paired model disagreement; optional surface-SD propagation is labelled as a
separate diagnostic. `cast_attribution()` provides optional driver-only and
retained-driver comparisons. Individual effects need not add to total change.

## Execution and reproducibility

Raster processing skips empty masked blocks and publishes artifacts atomically.
Current predictions are shared across futures. Additive univariate logit GAMs
reuse exact terms for partial scenarios; other structures use ordinary
prediction. A total core and memory budget bounds scenario workers and RF
threads. Workers require an installed matching package.
Worker startup uses the caller's exact package library without user/site
startup profiles. Loaded identities are checked in the caller before input
hashing or raster scans; mismatch messages identify the version, path or code.
Exactly unchanged group inputs reuse the matching corner prediction per row,
after the declared clamp. Outputs use band-interleaved tiled GeoTIFFs. Saved
phase timings distinguish input reads, prediction/contrast calculations and
writes from worker/parent wall time. Budgets do not promise peak process RSS.

Training-stage cache keys depend on data, scientific settings and relevant
code, independently of output paths and future scenarios. Comparison receipts
also verify model/input content and output hashes before reuse. OS file locks
prevent simultaneous writes. Resume in the same directory; choose a new
comparison directory or `overwrite = TRUE` after changing scientific inputs.
Caching avoids repeated work and does not replace scientific validation.
Fingerprints omit R's `srcfile` display metadata so source-preserving installs
also resume across serialization; actual model/data/code changes invalidate
the receipt. Schema-1 caches are preserved and recomputed under schema 2.

See `vignette("castSDM")` and `vignette("castSDM-advanced")` for executable
examples and the detailed contracts.
Bounded performance observations are documented in [benchmarks](benchmarks/README.md);
submission evidence and limits are assessed in [MEE_READINESS.md](MEE_READINESS.md).
The [AI-assisted maintenance record](inst/AI_ASSISTANCE.md) identifies its scope.

## License

GPL (>= 3).
