# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# Rscript benchmarks/paired_prediction.R completed_case_directory output_directory
# Uses at most 4,000 sampled eligible cells and existing fitted models. No refit.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Supply completed case directory and benchmark output directory.")
root <- args[1L]; out <- args[2L]
source_dir <- Sys.getenv("CAST_SOURCE_DIR", "")
if (nzchar(source_dir)) pkgload::load_all(source_dir, quiet = TRUE) else library(castSDM)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
internal <- function(name) getFromNamespace(name, "castSDM")
run <- readRDS(file.path(root, "training/result.rds"))
comparison <- readRDS(file.path(root, "comparison/comparison.rds"))
fit <- run$fit_full
case <- readRDS(file.path(root, "case.rds"))
scenario <- case$report_scenario
domain <- terra::rast(file.path(root, "comparison/.cache/common_domain.tif"))
set.seed(20261107)
draw <- sample.int(terra::ncell(domain), min(40000L, terra::ncell(domain)))
cells <- head(draw[which(terra::extract(domain, draw)[[1L]] == 1)], 4000L)
stopifnot(length(cells) >= 100L)
metadata <- read.csv(file.path(root, "training_pool.csv"))
vars <- fit$env_vars
current_paths <- metadata$path[match(vars, metadata$layer)]
table <- read.csv(file.path(root, "scenario_table.csv"))
table <- table[table$scenario == scenario, ]
future_paths <- current_paths
for (i in seq_len(nrow(table))) future_paths[match(table$layer[i], vars)] <- table$path[i]
current <- terra::rast(current_paths); names(current) <- vars
future <- terra::rast(future_paths); names(future) <- vars
a <- terra::extract(current, cells)[, vars, drop = FALSE]
b <- terra::extract(future, cells)[, vars, drop = FALSE]
old <- options(castSDM.predict_threads = 4L)
clamp <- internal(".cast_clamp")
model_predict <- internal("predict_single_model")
base <- vapply(fit$models[names(comparison$weights)], model_predict, numeric(nrow(a)),
  X_raw = if (isTRUE(run$protocol$clamp)) clamp(a, fit$scaling$reference) else a)
weights <- comparison$weights
groups <- comparison$groups
bounds <- internal(".cast_support_bounds")(fit$scaling$reference, comparison$support_probs)
reference <- function() {
  context <- internal(".cast_prediction_context")(fit, weights, a, b, run$protocol$clamp)
  preds <- list(base, context$predict(names(b)), context$predict(groups[[1L]]), context$predict(groups[[2L]]))
  faces <- vapply(preds, function(p) as.numeric(p %*% weights), numeric(nrow(a)))
  effects <- internal(".cast_factorial_values")(faces)
  per_model <- lapply(seq_along(weights), function(j) internal(".cast_factorial_values")(
    matrix(vapply(preds, function(p) p[, j], numeric(nrow(a))), nrow(a), 4)))
  sd <- vapply(1:8, function(k) internal(".ensemble_rowsd")(
    matrix(vapply(per_model, function(p) p[, k], numeric(nrow(a))), nrow(a))), numeric(nrow(a)))
  list(faces = faces, effects = effects, sd = sd)
}
optimized <- function() internal(".cast_compare_block")(fit, a, b, groups,
  comparison$static, weights, base, bounds, run$protocol$clamp)
ref <- reference(); value <- optimized()
errors <- vapply(c("faces", "effects", "sd"), function(nm)
  max(abs(ref[[nm]]-value[[nm]])), numeric(1))
stopifnot(all(errors < 1e-10))
context <- internal(".cast_prediction_context")(fit, weights, a, b, run$protocol$clamp)
changed <- lapply(groups, function(vars) Reduce(`|`, lapply(vars,
  function(v) context$current[[v]] != context$future[[v]])))
counts <- data.frame(n_sample = nrow(a), first_changed = sum(changed[[1L]]),
  second_changed = sum(changed[[2L]]), both_changed = sum(changed[[1L]] & changed[[2L]]),
  reference_rows_per_nonadditive_engine = 3*nrow(a),
  reused_rows_per_nonadditive_engine = sum(changed[[1L]] | changed[[2L]]) +
    2*sum(changed[[1L]] & changed[[2L]]))
methods <- list(reference = reference, optimized = optimized)
timings <- list()
for (replicate in 1:3) for (method in if (replicate %% 2L) names(methods) else rev(names(methods))) {
  gc()
  t <- system.time(value <- methods[[method]]())
  stopifnot(all(vapply(c("faces", "effects", "sd"), function(nm)
    max(abs(ref[[nm]]-value[[nm]])) < 1e-10, logical(1))))
  timings[[length(timings)+1L]] <- data.frame(replicate, method, elapsed_s = unname(t[["elapsed"]]))
}
options(old)
utils::write.csv(do.call(rbind, timings), file.path(out, "timings.csv"), row.names = FALSE)
utils::write.csv(counts, file.path(out, "prediction_rows.csv"), row.names = FALSE)
utils::write.csv(data.frame(quantity = names(errors), max_absolute_error = errors),
  file.path(out, "numerical_errors.csv"), row.names = FALSE)
saveRDS(list(case = normalizePath(root, winslash = "/"), scenario = scenario,
  sample_cells = cells, seed = 20261107, sys = Sys.info(),
  session = utils::capture.output(sessionInfo()),
  limitations = "Warm-filesystem prediction-only sample on one machine. Same models, trees, variables, weights, clamp and inputs; no whole-workflow or predictive-quality claim."),
  file.path(out, "manifest.rds"))
print(counts); print(do.call(rbind, timings)); print(errors)
