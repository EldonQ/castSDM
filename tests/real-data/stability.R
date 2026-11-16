# Conditional response sensitivity on every saved hold-out location.
# Rscript tests/real-data/stability.R <completed_real_case> <separate_QA_dir>
# Uses existing fits and actual future rasters. No refitting or row subsampling.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Supply a completed real case and a separate QA directory.")
library(castSDM)
case <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2L], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(tolower(qa), "/"), paste0(tolower(case), "/")))
  stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
options(castSDM.predict_threads = 2L)
run <- readRDS(file.path(case, "training/result.rds"))
comparison <- readRDS(file.path(case, "comparison/comparison.rds"))
split <- readRDS(file.path(case, "training/stages/01_prepare.rds"))
fits <- list(training = run$fit, refit = run$fit_full)
stopifnot(identical(fits$training$env_vars, fits$refit$env_vars))
vars <- fits$refit$env_vars
drivers <- as.character(unlist(comparison$groups, use.names = FALSE))
stopifnot(length(drivers) > 0L)
case_info <- readRDS(file.path(case, "case.rds"))
scenario <- case_info$report_scenario
mapping <- utils::read.csv(file.path(case, "scenario_table.csv"))
mapping <- mapping[mapping$scenario == scenario & mapping$layer %in% drivers, ]
stopifnot(setequal(mapping$layer, drivers), all(file.exists(mapping$path)))
current <- split$test[, vars, drop = FALSE]
future <- current
r <- terra::rast(mapping$path); names(r) <- mapping$layer
future[mapping$layer] <- terra::extract(r, as.matrix(split$test[, c("lon", "lat")]))
stopifnot(nrow(current) == nrow(split$test), nrow(future) == nrow(current))
weights <- comparison$weights
common <- castSDM:::.cast_finite_rows(current) & castSDM:::.cast_finite_rows(future)
for (fit in fits) {
  bounds <- castSDM:::.cast_support_bounds(fit$scaling$reference, comparison$support_probs)
  common <- common & castSDM:::.cast_attribution_inside(current, bounds) &
    castSDM:::.cast_attribution_inside(future, bounds)
}
common[is.na(common)] <- FALSE
stopifnot(any(common))
effects <- lapply(fits, function(fit) {
  context <- castSDM:::.cast_prediction_context(fit, weights, current, future, comparison$clamp)
  base <- context$predict()
  stats::setNames(lapply(drivers, function(v) {
    value <- context$predict(v)-base
    cbind(value, ensemble = as.numeric(value %*% weights))
  }), drivers)
})
rows <- locations <- list()
for (v in drivers) for (m in colnames(effects$training[[v]])) {
  a <- effects$training[[v]][, m]; b <- effects$refit[[v]][, m]
  use <- common & is.finite(a) & is.finite(b)
  correlation <- if (length(unique(a[use])) > 1L && length(unique(b[use])) > 1L)
    stats::cor(a[use], b[use], method = "spearman") else NA_real_
  rows[[length(rows)+1L]] <- data.frame(driver = v, model = m,
    scenario = scenario, n_input = nrow(current), n_common = sum(use),
    mean_training = mean(a[use]), mean_refit = mean(b[use]),
    mean_abs_training = mean(abs(a[use])), mean_abs_refit = mean(abs(b[use])),
    mean_absolute_difference = mean(abs(a[use]-b[use])),
    sign_agreement = mean(castSDM:::.cast_response_sign(a[use]) == castSDM:::.cast_response_sign(b[use])),
    spearman = correlation)
  locations[[length(locations)+1L]] <- data.frame(split$test[, c("lon", "lat")],
    driver = v, model = m, common_support = use, training_effect = a, refit_effect = b)
}
result <- do.call(rbind, rows)
utils::write.csv(result, file.path(qa, "response_stability.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, locations), file.path(qa, "matched_locations.csv"), row.names = FALSE)
inputs <- c(file.path(case, c("training/result.rds", "training/stages/01_prepare.rds", "scenario_table.csv")), mapping$path)
utils::write.csv(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))),
  file.path(qa, "input_provenance.csv"), row.names = FALSE)
cat("All", nrow(current), "saved hold-out locations retained;", sum(common),
  "share paired marginal support under both fits.\n")
print(result[result$model == "ensemble", c("driver", "mean_training", "mean_refit", "sign_agreement", "spearman")], row.names = FALSE)
cat("Conditional sensitivity to adding the hold-out records to refitting; not independent validation, area estimates or causal effects.\n")
