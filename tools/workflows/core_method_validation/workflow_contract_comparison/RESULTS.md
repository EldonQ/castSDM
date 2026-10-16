# Same-task workflow comparison

A correctly aligned iml prediction loop and castSDM reproduce the same paired contrasts.
castSDM also checks raster geometry and pairs table coordinates automatically.
The naive rowwise example omits alignment; this is not a defect in iml.
The raster result agrees with the table result across block boundaries.
Five timing repetitions quantify execution cost only; no speed, usability or memory advantage is inferred.
This uses one external GLM and two changing predictors, not an ALE/SHAP or package-wide benchmark.
