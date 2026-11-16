# Rscript tests/real-data/prediction.R <completed_real_case> <separate_QA_dir>
# Exactness and matched execution timings on saved models and actual futures.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
library(castSDM)
case <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2L], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(tolower(qa), "/"), paste0(tolower(case), "/"))) stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
options(castSDM.predict_threads = 2L)
run <- readRDS(file.path(case, "training/result.rds"))
fit <- if (is.null(run$fit_full)) run$fit else run$fit_full
comparison <- readRDS(file.path(case, "comparison/comparison.rds"))
split <- readRDS(file.path(case, "training/stages/01_prepare.rds"))
scenario <- readRDS(file.path(case, "case.rds"))$report_scenario
mapping <- read.csv(file.path(case, "scenario_table.csv"))
drivers <- as.character(unlist(comparison$groups, use.names = FALSE))
mapping <- mapping[mapping$scenario == scenario & mapping$layer %in% drivers, ]
stopifnot(setequal(mapping$layer, drivers))
current <- split$test[, fit$env_vars, drop = FALSE]; future <- current
r <- terra::rast(mapping$path); names(r) <- mapping$layer
future[mapping$layer] <- terra::extract(r, as.matrix(split$test[, c("lon", "lat")]))
stopifnot(all(is.finite(as.matrix(current))), all(is.finite(as.matrix(future))))
reference <- fit$scaling$reference
outside <- Reduce(`|`, lapply(fit$env_vars, function(v)
  future[[v]] < min(reference[[v]]) | future[[v]] > max(reference[[v]])))
results <- timing <- list()
for (clamp in c(FALSE, TRUE)) {
  start <- proc.time()[["elapsed"]]
  context <- castSDM:::.cast_prediction_context(fit, comparison$weights, current, future, clamp)
  scenarios <- c(list(total = drivers), setNames(lapply(drivers, identity), paste0("only_", drivers)),
    setNames(lapply(drivers, function(v) setdiff(drivers, v)), paste0("without_", drivers)))
  accelerated <- lapply(scenarios, context$predict)
  elapsed_fast <- proc.time()[["elapsed"]]-start
  start <- proc.time()[["elapsed"]]
  for (name in names(scenarios)) {
    data <- current; data[scenarios[[name]]] <- future[scenarios[[name]]]
    for (model in names(fit$models)) {
      expected <- castSDM:::predict_single_model(fit$models[[model]], data,
        clamp = clamp, reference = reference)
      difference <- max(abs(expected-accelerated[[name]][, model]))
      stopifnot(is.finite(difference), difference < 1e-10)
      results[[length(results)+1L]] <- data.frame(model = model, clamp = clamp,
        combination = name, n = nrow(data), n_future_outside_training_range = sum(outside),
        max_difference = difference)
    }
  }
  elapsed_ordinary <- proc.time()[["elapsed"]]-start
  # The response-curve cache must share the same policy even outside training.
  for (model in names(fit$models)) {
    predict <- castSDM:::.cast_response_predictor(fit$models[[model]], reference, clamp)
    difference <- max(abs(predict(future)-accelerated$total[, model]))
    stopifnot(is.finite(difference), difference < 1e-10)
    results[[length(results)+1L]] <- data.frame(model = model, clamp = clamp,
      combination = "response_api", n = nrow(future),
      n_future_outside_training_range = sum(outside), max_difference = difference)
  }
  timing[[length(timing)+1L]] <- data.frame(clamp = clamp, n = nrow(current),
    combinations = length(scenarios), ordinary_s = elapsed_ordinary, paired_s = elapsed_fast,
    threads = getOption("castSDM.predict_threads"), R = as.character(getRversion()),
    package = as.character(packageVersion("castSDM")),
    interpretation = "Matched engineering timings on saved fits; not full-workflow or ecological validation")
}
write.csv(do.call(rbind, results), file.path(qa, "prediction_contract.csv"), row.names = FALSE)
write.csv(do.call(rbind, timing), file.path(qa, "prediction_timings.csv"), row.names = FALSE)
inputs <- c(file.path(case, c("training/result.rds", "training/stages/01_prepare.rds", "scenario_table.csv")), mapping$path)
write.csv(data.frame(file = inputs, md5 = unname(tools::md5sum(inputs))),
  file.path(qa, "prediction_inputs.csv"), row.names = FALSE)
cat("All actual hold-out locations and declared climate combinations passed both clamp policies.\n")
