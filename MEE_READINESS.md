# castSDM: evidence and submission readiness

Engineering update: 13 November 2026 (journal guidance below was checked on
7 November 2026). This is an engineering and evidence assessment,
not an editorial decision or a promise of acceptance.

## Position and contribution

Position castSDM as an **Applications software package** for interactive,
auditable environmental response maps from frozen SDMs. `cast_app()` and
`cast_explore()` connect explicit environmental changes, native-grid queries
and descriptive environmental subgroups. `cast_compare()` supplies paired
two-group comparisons; `cast_run()` supplies optional built-in fitting and
`cast_as_fit()` accepts external engines. `cast_diagnose()` joins validation to
the interpretation of scenario comparisons.

The contribution to demonstrate is the complete comparison contract: paired
four-corner predictions, frozen positive weights, explicit finite and common
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
on the review date; verify them again at submission.

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
