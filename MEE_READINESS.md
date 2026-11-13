# castSDM: evidence and submission readiness

Engineering update: 13 November 2026 (journal guidance below was checked on
13 November 2026). This is an engineering and evidence assessment,
not an editorial decision or a promise of acceptance.

## Position and contribution

Position castSDM as an **Applications software package** for interactive,
auditable environmental response maps from frozen SDMs. `cast_app()` and
`cast_explore()` connect explicit environmental changes, native-grid queries
and descriptive environmental subgroups. `cast_compare()` supplies paired
one-group projections and optional two-group comparisons; `cast_run()` supplies optional built-in fitting and
`cast_as_fit()` accepts external engines. `cast_diagnose()` joins validation to
the interpretation of scenario comparisons.

The contribution to demonstrate is the complete comparison contract: paired
current/future predictions, optional four-corner predictions, frozen positive weights, explicit finite and common
support domains, paired model disagreement, consistent thresholds, bounded
raster execution and content-verified reuse. A 2 x 2 contrast identity, spatial
CV, RF/MaxEnt/GAM, marginal support and caching alone are not new methods.
Do not advertise a new causal estimator or calibrated occurrence probability.

Version 0.29.0 adds a local Shiny response workspace and descriptive
environmental subgroups. The grouping follows published FME/cAME concepts;
spatial-fold tuning of a tree summarises frozen-model outputs rather than
validating ecological effects. SHAP spatial explanations/environmental change
analysis already exist in [itsdm](https://doi.org/10.1111/2041-210X.14067),
and [fmeffects](https://doi.org/10.32614/RJ-2024-024) already implements finite
changes and conditional response groups. A GUI, paired subtraction or adding
these published methods does not establish methodological novelty. The
case-specific value and comparative evidence still need to be demonstrated.

Version 0.29.1 makes a narrower application question operational: for the same
supported locations and the same environmental change, does its interpretation
depend on the other predictors' background or the fitted engine? The overview
compares a shared population, excludes unchanged scenario predictors from
ranking, shows both configured-comparison backgrounds and reports sign
reversals and active-model agreement. One focused driver has native maps of
both responses and their diagnostics. These are descriptive checks of frozen
models, not an ecological causal estimator or a new uncertainty interval.

The closest comparison set must also include
[GLOSSA (MEE 2026)](https://doi.org/10.1111/2041-210x.70167), which already
provides a Shiny SDM workflow and reproducible outputs, and
[MaskSDM (MEE 2026)](https://doi.org/10.1111/2041-210x.70200), which already
studies predictor contributions with Shapley values. The proposed contribution
needs a matched real-case comparison demonstrating what the combined
background/model/support diagnostics reveal, rather than an unsupported claim
that these tools cannot perform similar analyses.

Existing tools already address important adjacent problems. [CAST](https://hannameyer.github.io/CAST/)
provides spatial model training/assessment and [AOA](https://hannameyer.github.io/CAST/reference/aoa.html).
[biomod2](https://biomodhub.github.io/biomod2/) provides multiple SDM engines,
evaluation and ensemble projection. These are documented capabilities, not an
exhaustive feature audit. A paper must compare relevant workflows fairly;
absence of a castSDM-style example on another project's front page does not
establish that the other project cannot reproduce a contrast. castSDM is a
separate project and is not the CAST package or an endorsed extension.

## Journal requirements and our own criteria

The current [MEE author guidelines](https://besjournals.onlinelibrary.wiley.com/hub/journal/2041210x/author-guidelines)
describe Applications as open-source software papers of at most 4,000 words,
including captions, statements and references. They require a worked example,
a Code Checklist, review access to anonymised code/data, appropriate archiving,
and disclosure/annotation of AI assistance. A manuscript should explain the
ecological value beyond a manual or vignette. These requirements were checked
on the review date; verify them again at submission. GUI packages must export
executable scripts for their graphical/statistical/analytical outputs. The
0.29.1 bundle includes view state and a script restoring the completed request,
selected layer, map view and point/profile query using the original inputs.

The criteria below are our assessment, not additional journal rules or
universal numerical acceptance thresholds.

| Evidence | Current scope | Remaining boundary |
|---|---|---|
| Numerical correctness | Hand-computed contrast tests, table/raster parity, paired SD, support/threshold checks, installed parallel execution | Does not establish ecological validity |
| Validation design | Nested spatial screening, training OOF thresholds/weights, frozen ensemble hold-out diagnostics, per-fold composition/AOA | A repeatedly inspected hold-out cannot validate a newly selected policy |
| Score semantics | Native/shifted scores retained from one prediction; sampled-response losses and baselines labelled | Native same-fit comparison does not test a new RF weighting or tuning policy |
| Extrapolation | Declared clamp applied in CV, evaluation and maps; raw support/AOA retained | Marginal support is not multivariate support or transferability |
| Sensitivity | All-scenario common support, model sign conflict, threshold counts with denominators | No sampling CI or significance; one GCM does not sample climate-model uncertainty |
| Software quality | Package checks, historical regression fixtures, verified receipts and bounded benchmarks; this release uses saved real-case verification | Lightweight CI does not execute historical simulated fixtures/examples; local evidence is one OS/R setup |
| Reproducibility | Source, examples, benchmark scripts and environment manifests | Deposit a reviewed release/archive DOI and authorised data before submission |

## Case-specific evidence still needed

The completed freshwater cases, including the 0.28.x case used for software
verification, are preserved. They pool Schizothoracinae
records and must not be described as a single biological species. Record years
span historical periods, and the meaning of `year` needs source verification.
GLU snapshots and climate windows differ. The domain was constructed with all
records before the hold-out. Document these choices rather than calling the
hold-out external or temporal validation.

For a defensible application paper, plan a genuinely new validation design
before policy comparison: independent source/basin or appropriately matched
temporal data, with enough presences/backgrounds in each partition. Background
years cannot be invented from occurrence years. Examine the weak spatial fold
using its actual locations, class composition and environmental coverage;
AOA alone does not explain it. Keep a used hold-out for retrospective diagnosis,
and use new validation for a changed scientific policy.

This release verifies saved real-case model outputs, real raster crops and
interactive operations. Maintenance has not run the complete research workflow
or refitted the study models. The local single-species and Shiny entries are
ready for the author to run; saved metrics retain their original provenance.

## Minimum paper evidence to assemble next

1. State one practical ecological question and show why inconsistent domains,
   weights or surface-level SD change its interpretation. Explain the user
   benefit of the comparison contract rather than claiming the arithmetic is new.
2. Use the current real records and saved model/raster outputs for reproducible
   numerical parity and software checks. Do not treat generated model responses
   as ecological causal truth or invent missing observations for validation.
   A biologically defined, independent real-data validation would strengthen
   ecological claims; the existing data cannot certify causal-effect accuracy.
3. Compare a transparent reference implementation and relevant existing tools
   on identical fits, folds, grids, masks, precision, model counts and outputs.
   Include correctness, elapsed time, memory definition and setup/serialization
   costs. Separate cold execution, verified resume and threshold-only summaries.
   The bounded scripts in `benchmarks/` do not replace a full scale benchmark.
4. Report threshold/support sensitivity, model conflict, selection stability,
   cohort/time matching and limits of the independent validation. A decline in
   the ensemble does not establish unanimous models or a causal ecological effect.
5. Verify the remote OS/R matrix, finish the official Code Checklist, obtain
   authorised data access/archiving, freeze a release and prepare an anonymised
   review bundle. Do not upload third-party data without permission.

AI-assisted development scope and reporting limits are documented in
[AI_ASSISTANCE.md](inst/AI_ASSISTANCE.md). See the local maintenance report for
the actual check counts, benchmark observations and commit identifier.


## 0.30.0 climate-only recipe and interface

The package ships its single maintained workflow under `inst/workflows/`.
The Start tab configures paths and a small scientific/resource settings list,
checks the delivered CHELSA library, and explicitly submits the full analysis.
Advanced controls are folded. Long computations run outside Shiny's main R
process, with duplicate-action protection and process-tree cancellation; see
[Posit nonblocking guidance](https://shiny.posit.co/r/articles/improve/nonblocking/)
and [input validation](https://shiny.posit.co/r/articles/improve/validation/).

The case uses CHELSA v2.1 baseline 1981–2010 and exactly the delivered
SSP126/370/585 ensemble climates. The official [CHELSA catalogue](https://www.chelsa-climate.org/datasets/chelsa_bioclim)
provides definitions; processed-library provenance is authoritative for decoded
units. A CNClim BIO23 model cannot be treated as a CHELSA fit. All nonclimate
predictors remain current; future land-use products are excluded.

The scientific question remains: **where does the same environmental change
produce different fitted suitability responses, and do those responses persist
across environmental backgrounds and algorithms?** The interface makes this
question accessible but is not itself methodological novelty. Claims require
a completed new case, comparisons with existing tools, sensitivity to validation,
background sampling and extrapolation, and useful conclusions attributable to
the diagnostic contract. This engineering release alone does not establish MEE
readiness or identify ecological causal effects from presence/background data.


### 0.30.0 checks performed on 13 November 2026

- Installed package and light `R CMD check`: 0 errors, 0 warnings, 0 notes.
  Scope: `--no-manual --no-tests --no-examples --ignore-vignettes`; historical
  simulated fixtures/examples were not executed.
- Real-library preflight: 458 current layers, 447 passing coverage/geometry,
  326 GLU or 320 HLU eligible predictors, 198 future climate layers; all nine
  exact SSP/period mappings verified and nonclimate source paths unchanged.
- Numerical regression replays the 256 saved real-case current/future scores,
  plus an actual saved grid crop. It checks paired fields, threshold counts,
  static-change rejection, multiblock reductions and installed PSOCK parity.
  Replay adapters only return exact already-recorded scores; these checks do
  not validate new CHELSA fitted engines or estimate ecological ground truth.
- Headless Chrome: setup connection, input failure/recovery, duplicate-submit
  protection, cancellation and subsequent real-input checks, desktop/mobile
  overflow. Both deployed entry points pass dry-run checks.
- No full study preparation, model fitting or national projection was run.
  New CHELSA scientific results, full figure execution and post-fit result tabs
  require the user's new completed case. `examples/validate_real_case.R` checks
  a completed real case, including each engine's paired prediction exports.
