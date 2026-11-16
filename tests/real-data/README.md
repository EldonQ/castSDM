# Actual-data verification

`selection.R` reads actual saved preparation rows and matching spatial folds
from a case supplied by its author. It never manufactures measurements or
labels. It checks semantic qualification, complete spatial losses, unit
conversion stability and actual engine predictions on saved validation rows.
An optional row limit bounds engineering cost; it is not a predictor-count rule.

Run with the installed package and its dependencies on the R library path:

```r
# In a shell: Rscript tests/real-data/selection.R <actual-case-directory> 3000
```

These large-case checks are opt-in and use the author's completed analysis.
Portable CI tests in `tests/testthat` use the bundled, attributed SchiSOFT/CHELSA
fixture instead. Historical simulated fixtures have been removed. A passing
engineering check is not independent ecological validation.

`prediction.R <actual-case> <QA-dir>` checks ordinary, paired and response API
predictions on every saved hold-out location, including actual future values
outside training ranges. Both clamp policies and all driver-only/driver-removed
combinations are checked against ordinary engine predictions. Matched timings
are recorded separately from full-workflow performance.

`response.R <actual-case> <QA-dir> [scenario-io-bundle]` compares accelerated ALE
with official iml and ordinary predictions on every saved reference row. No
hidden reference cap is used. The optional verified real-crop bundle written
by `scenario-io.R` also checks the eight-figure default, requested curve scope,
PNG/PDF exports, complete inventory and render-only reuse.

`stability.R <actual-case> <QA-dir>` evaluates conditional sensitivity to the
saved training versus full-data refits at every original hold-out location.
Actual future predictor values come from the declared reporting scenario.
One shared paired marginal-support subset is used for every driver and engine;
the table reports coverage, effect means, sign agreement and rank correlation.
It does not refit, subsample, estimate habitat area or prove ecological causality.


`maintenance.R <actual-case> <QA-directory>` verifies a real native-grid crop,
frozen GAM shortcut parity, controlled cancellation/resume, prediction-free
basic plots and published ALE. Keep the QA directory outside the input case.
No study-wide fitting or national scenario projection is executed.

`background.R <actual-case> <QA-directory>` compares 10,000 and the saved
background budget with fixed original training-CV folds and evaluation rows.
It refits learners on actual rows with the saved predictors: this is conditional
sensitivity, not a new variable-selection validation or a test of a new domain.
It writes discrimination/loss and finite-location response-sign comparisons;
the latter are not common-support confidence estimates.

`scenario-io.R <actual-case> <QA-directory>` checks actual native crops across
serial/PSOCK execution, verified resume, saved-model batch invariance and
paired-prediction reuse in environmental response maps. It changes no input
measurements and does not run the national research workflow.

It also verifies compact current-GAM reuse and member response summaries against
native change rasters and their frozen ensemble.

`metadata.R <converted_SchiSOFT_csv> <official_occurrence_csv>` checks original
IDs/names/coordinates and distinguishes reported sampling years from ambiguous
intervals and source-date substitutes. It changes no record coordinates or years.
`scenario-io.R` also checks equal-domain driver summaries and display-only
resampling, including native-file integrity and explicit mixed coarse support.

`display.R <actual-case> <QA-dir> [national]` checks mask-before-aggregation
against independent native-raster oracles, common-summary denominators,
per-driver double-background support, unchanged source files and shared-display
cache parity. The optional national mode checks the saved response footprint
across the full original grid, without model fitting or scenario prediction.
Comparison-cache parity remains bounded to an actual raster crop in either mode.

`temporal.R <prepare_stage.rds> <occurrence_qc.csv> <eligible_pool.csv> <QA_dir>`
checks the complete reported-year -> layer-matching -> CSV boundary using real
saved rows, original labels and original raster time slices. It also replays
legacy attributed-year input. It performs no fitting or national projection.

`interface.R <actual-case> <QA-dir>` opens an isolated shadow of actual completed artifacts, uses real installed background jobs, checks all saved reference locations and a real point query, then exports and replays the numeric request. The original case is read-only; this does not run national predictions or refit the case.
