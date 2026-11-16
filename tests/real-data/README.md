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

The input data must be obtained lawfully. This is opt-in: CI has no authorised
copy of the author's environmental/occurrence library. A passing engineering
check is not an independent ecological validation or a completed case study.
Legacy simulated unit fixtures in other test files are not run for this release.


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
