# Learner-only background sensitivity on an actual saved cohort.
# Rscript background.R <actual_case> <separate_QA_directory>
# Predictors and evaluation locations/folds remain fixed. This does not
# estimate selection uncertainty or validate a different ecological domain.
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Supply the actual case and a separate QA directory.")
suppressPackageStartupMessages(library(castSDM))
case <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(qa, "/"), paste0(case, "/"))) stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
run <- readRDS(file.path(case, "training/result.rds"))
split <- readRDS(file.path(case, "training/stages/01_prepare.rds"))
comparison <- readRDS(file.path(case, "comparison/comparison.rds"))
data <- split$train
variables <- run$screen$selected
models <- run$protocol$models
folds <- run$cv$folds
stopifnot(length(folds) == nrow(data))
bundle <- readRDS(file.path(case, "case.rds"))
prepared <- data.table::fread(file.path(bundle$preparation,
  paste0(bundle$species, "_model_data.csv")), data.table = FALSE)
stopifnot(identical(as.integer(prepared$presence[split$train_index]), as.integer(data$presence)))
background <- which(prepared$presence == 0L)
budgets <- unique(c(min(10000L, length(background)), length(background)))
if (length(budgets) < 2L) stop("The saved case has too few backgrounds for this comparison.")
scores <- summaries <- effects <- list()
pool <- read.csv(file.path(case, "training_pool.csv"))
scenario <- bundle$report_scenario
table <- read.csv(file.path(case, "scenario_table.csv"))
table <- table[table$scenario == scenario, , drop = FALSE]
future <- terra::rast(table$path); names(future) <- table$layer
after <- data[, variables, drop = FALSE]
after[table$layer] <- terra::extract(future, as.matrix(data[, c("lon", "lat")]))[table$layer]
for (budget in budgets) {
  allowed <- prepared$presence[split$train_index] == 1L |
    split$train_index %in% background[seq_len(budget)]
  predicted <- matrix(NA_real_, nrow(data), length(models), dimnames = list(NULL, models))
  for (fold in sort(unique(folds))) {
    train <- folds != fold & allowed
    test <- folds == fold
    message("Actual background sensitivity: budget ", budget, "; fold ", fold,
      "; training backgrounds ", sum(train & data$presence == 0L), ".")
    fit <- cast_fit(data[train, c("presence", variables), drop = FALSE],
      models = models, gam_backend = if (is.null(run$protocol$gam_backend)) "gam" else run$protocol$gam_backend,
      rf_ntree = run$protocol$rf_ntree, num_threads = 2L,
      prevalence_target = run$protocol$prevalence_target,
      seed = run$protocol$seed, verbose = FALSE)
    prediction <- cast_predict(fit, data[test, variables, drop = FALSE],
      clamp = run$protocol$clamp)
    predicted[test, ] <- as.matrix(prediction$predictions[paste0("HSS_", models)])
  }
  stopifnot(all(is.finite(predicted)))
  scores[[as.character(budget)]] <- predicted
  summaries[[as.character(budget)]] <- do.call(rbind, lapply(models, function(model)
    cbind(background_budget = budget, n_training_background = sum(allowed & data$presence == 0L),
      castSDM:::.cast_diagnostic_metrics(predicted[, model], data$presence,
        threshold = NULL, model = model, view = "fixed-predictor sensitivity",
        threshold_source = "not_applicable"))))
  fit <- cast_fit(data[allowed, c("presence", variables), drop = FALSE], models = models,
    gam_backend = if (is.null(run$protocol$gam_backend)) "gam" else run$protocol$gam_backend, rf_ntree = run$protocol$rf_ntree,
    num_threads = 2L, prevalence_target = run$protocol$prevalence_target,
    seed = run$protocol$seed, verbose = FALSE)
  effects[[as.character(budget)]] <- cast_attribution(fit, data[, variables, drop = FALSE],
    after, drivers = table$layer, weights = comparison$weights,
    clamp = run$protocol$clamp, diagnostics = TRUE, verbose = FALSE)$data
}
metrics <- do.call(rbind, summaries)
write.csv(metrics, file.path(qa, "background_metrics.csv"), row.names = FALSE)
first <- effects[[1]]; last <- effects[[length(effects)]]
fields <- intersect(c("total_change", paste0("effect_", table$layer)), names(first))
changes <- do.call(rbind, lapply(fields, function(field) {
  a <- first[[field]]; b <- last[[field]]; ok <- is.finite(a) & is.finite(b)
  data.frame(response = field, n = sum(ok), mean_abs_difference = mean(abs(a[ok] - b[ok])),
    sign_agreement = mean(castSDM:::.cast_response_sign(a[ok]) == castSDM:::.cast_response_sign(b[ok])))
}))
write.csv(changes, file.path(qa, "background_response_sensitivity.csv"), row.names = FALSE)
saveRDS(list(budgets = budgets, metrics = metrics, responses = changes,
  interpretation = "Conditional sensitivity on saved predictors, original domain and fixed spatial folds; not independent validation or a test of the new domain."),
  file.path(qa, "background_check.rds"))
print(metrics[, c("background_budget", "model", "auc", "average_precision", "brier", "logloss")])
print(changes)
