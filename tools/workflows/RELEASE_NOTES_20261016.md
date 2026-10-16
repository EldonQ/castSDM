# Release verification and limits

Core: 120 independent datasets, 492 passing numerical checks, zero run errors.
Case: four engines, 18 future projections and attribution stacks, 400 sampled
valid cells per scenario; maximum projection/attribution discrepancy 4.96e-8.
Four of five spatial folds were evaluable. Original full local artifacts remain.

Four core PDFs, four original case PDFs and four corrected case PDFs were
visually inspected. Corrected case PDFs are in single_species/evidence/plot_revision/figures.
Plotting-only source revisions have separately recorded baselines; computation
receipts are unmodified. The latest short R check is 0 errors, 0 warnings, 0 notes.

The core benchmark has model-dependent advantages and negative cases.
It does not establish general superiority over ALE/SHAP or ecological causality.
The genus demonstration, national terrestrial background and single changing
future predictor cannot replace independent species ecological validation.
Candidate audit/protocol are complete; ecological validation is not run.
Legacy mask/nonfinite/year-order issues remain recorded in the severity audit.
This release narrows responsibilities; it does not claim net source shrinkage.

Generated sessionInfo/ODMAP reports and preserved baseline source snapshots retain
their original whitespace for byte-level provenance. Diff whitespace checking
exempts these files; maintained executable code is checked separately.
tools/workflows is excluded from package builds. Raw occurrences, models,
large rasters and large SVG files are not included in the Git evidence bundle.
