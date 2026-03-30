if (file.exists("DESCRIPTION")) devtools::load_all() else library(castSDM)
data(ovis_ammon)
data(china_env_grid)

cat("=== Step 1: Prepare ===\n")
split <- cast_prepare(ovis_ammon, train_fraction = 0.7, seed = 42)
cat("Train:", nrow(split$train), "| Test:", nrow(split$test), "\n")
cat("Env vars:", length(split$env_vars), "\n")

cat("\n=== Step 2: DAG ===\n")
dag <- cast_dag(split$train, R = 30, seed = 42)
print(dag)

cat("\n=== Step 3: ATE ===\n")
ate <- cast_ate(split$train, K = 2, num_trees = 100, seed = 42)
print(ate)

cat("\n=== Step 4: Screen ===\n")
screen <- cast_screen(dag, ate, split$train, seed = 42)
print(screen)

cat("\n=== Step 5: Roles ===\n")
roles <- cast_roles(screen, dag)
print(roles)

cat("\n=== Step 6: Fit RF ===\n")
fit <- cast_fit(split$train, screen = screen, dag = dag, ate = ate,
                models = "rf", seed = 42)
print(fit)

cat("\n=== Step 7: Evaluate ===\n")
ev <- cast_evaluate(fit, split$test)
print(ev)

cat("\n=== Step 8: Predict ===\n")
pred <- cast_predict(fit, china_env_grid)
print(pred)

cat("\n=== ALL STEPS PASSED ===\n")
