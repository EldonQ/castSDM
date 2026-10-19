# Reproducible research workflows

These scripts mirror the maintained workflows, adapted to castSDM 0.20.0. The
completed case and core evidence were computed with their recorded 0.19.0 sources.
The [2026-10-19 raster maintenance](RASTER_REFACTOR_20261019.md) has separate
engineering evidence; it does not relabel historical experiments. The scripts
remain outside the installed package: model fitting is delegated to the engine
packages, and the package owns scenario comparisons and plotting.

## Run on the recorded workstation

Copy the contents of this directory to `E:/Package/castSDM_workflows` and keep
the package checkout at `E:/Package/cast`. Use R 4.4.3 with the packages checked
by `single_species/02_pipeline.R::required_workflow_packages()`.

```r
source("E:/Package/castSDM_workflows/single_species/00_RUN_ALL.R", encoding = "UTF-8")
```

For the independent known-response benchmark:

```sh
Rscript E:/Package/castSDM_workflows/core_method_validation/00_RUN_VALIDATION.R production
```

The case defaults to Schizopygopsis. `CAST_SPECIES` can select another configured
taxon. The data are genus-level SchiSOFT literature records, not one species.
Input paths and model settings are declared in `single_species/01_config.R`;
the 596-layer pool manifest records the raster paths and declared climate layers.
Occurrence inputs and large rasters are external data and are not redistributed
here. Current and future rasters must share the original 1 km Albers grid.

Paths preserve the executed configuration. To move to another workstation,
update `ROOT`/`PACKAGE` in the runners, `01_config.R`, and the adapter defaults.
The validation worker source path also appears in `00_RUN_VALIDATION.R`.
Sources and inputs are fingerprinted: changed sources require a fresh output
directory. Do not relabel an old result as evidence for changed code.

## Evidence and interpretation

`core_method_validation/RESEARCH_PLAN_20261015.md` defines the protocol.
`core_method_validation/RESULTS_20261016.md` states the measured advantages and failures. The complete
benchmark uses six scenarios with 20 independent datasets each, RF and GAM,
matched-model PDP and uniform-shift baselines, independent data resampling,
forest seeds, model-family disagreement, and numerical table/raster checks.

The case retains four engines, nested spatial CV, held-out evaluation,
full-data refit, all 18 configured future scenarios, attribution tables/maps,
MESS, source-native plotting and ODMAP reporting. The four default result
figures answer validation, habitat change, environmental response and range/model
disagreement. Single-shift audits and ensemble permutation are opt-in diagnostics.
National background, ensured climate coverage and genus-level records are
explicit modelling choices, not proof of ecological mechanism. All outputs
are relative suitability and fitted-model scenario comparisons.

The current core evidence is copied from `complete_production_20261016T041431Z`
into `core_method_validation/evidence`; the original paths remain in the
receipt and source manifest. `03_external_interoperability.R` actually checks
iml prediction callbacks and blockCV fold assignments. These are equivalence
checks, not evidence of superiority. Candidate-species quality and the separate
ecological validation protocol are in
`core_method_validation/ECOLOGICAL_VALIDATION_PLAN_20261016.md`.

`05_workflow_contract_comparison.R` compares the identical contrasts with an
explicit iml prediction loop. Correctly aligned outputs agree exactly; table
coordinate pairing and raster geometry rejection are checked independently.
This measures correctness on one task, not broad speed or usability superiority.

Each run exports receipts, source manifests, session information and source
data for the figures. Plotting uses R, Arial, 183 mm width, PDF/SVG and 600 dpi
PNG/TIFF. Visual inspection is recorded separately from successful export.

The full case completed as `production_20261016T041739Z_158352/Schizopygopsis`:
four engines, 18 projections and 18 attribution stacks. Compact evidence is in
`single_species/evidence`; use `plot_revision/figures` for the four corrected
main PDFs. Original exports and receipts are retained. The final plotting changes
correct unavailable-CV label overlap, square-window labels and binary legends;
they do not change numerical results. `04_REPLOT.R` provides a separate-directory
replot entry using the original computation identity.
The earlier local `06_smoke_current_api.R` is not part of the maintained mirror;
the complete production run is the accepted demonstration.

The [2026-10-16 audit](POSITIONING_AND_CODE_AUDIT_20261016.md) gives the scientific
question, direct source comparisons, severity-ranked defects and repair status.
External models can use `cast_as_fit()`; external spatial groups can use
`cast_cv(fold_ids = ...)`. These interfaces reduce duplicate fitting and blocking
machinery while preserving explicit prediction and validation responsibilities.
