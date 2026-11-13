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
