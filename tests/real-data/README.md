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
