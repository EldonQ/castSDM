# Scope compared with existing open-source tools

Source inspection on 2026-10-15, with fixed GitHub revisions:

| Tool | Inspected entry | Relevant behaviour | Implication for castSDM |
|---|---|---|---|
| biomod2 | [bm_PlotResponseCurves](https://github.com/biomodhub/biomod2/blob/11a4c0fe417dc9716f8f96d0575fe41d38308ea8/R/bm_PlotResponseCurves.R) | Varies one or two predictors and fixes remaining predictors at a stated mean/median/min/max reference | Response curves are established SDM functionality; preserve interoperability instead of duplicating a modelling ecosystem |
| iml | [FeatureEffect](https://github.com/christophM/iml/blob/4818c8f81dea2e76318857fbc3d50bc78cdd22f6/R/FeatureEffect.R) | Provides PDP, ICE and ALE behind a common predictor interface | Local fitted-response contrasts are not intrinsically new; a generic predictor boundary and clear comparison target matter |
| DALEX | [predict_parts](https://github.com/ModelOriented/DALEX/blob/77ea3ad9580a3b68aaa94d8dca75d11bd527848f/R/predict_parts.R) | Delegates breakdown and SHAP methods to maintained specialist packages | Keep prediction adapters thin; do not implement a second general-purpose explanation library |

castSDM 0.19.0 centres paired current/future environmental values, driver-only
and retained-driver spatial comparisons, frozen ensemble weights, an explicit
nonadditive remainder, marginal training-range diagnostics, and source-native
maps. This is a focused ecological workflow contribution. It is not a claim
that predictor replacement, ICE or standardisation is a new estimator.

The known-response benchmark compares precisely defined global marginal PDP
and uniform-shift baselines using the same fitted models. It does not directly
benchmark the full biomod2/iml/DALEX packages, ALE, SHAP, or all existing spatial
attribution methods. The results therefore cannot establish superiority over
those packages. The source inspection informs scope and engineering design.

Existing fitting/CV APIs remain compatible in this release. Engine packages
own their fit/predict implementations; castSDM owns the ecological comparison
and reporting. Further slimming should follow verified call paths and measured
duplication, with cartographic assets and reproducible evidence preserved.
