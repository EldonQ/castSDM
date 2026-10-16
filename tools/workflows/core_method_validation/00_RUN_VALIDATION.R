# Run six prespecified scenarios concurrently; each worker uses one RF thread.
ROOT <- "E:/Package/castSDM_workflows/core_method_validation"
Sys.setenv(R_LIBS_USER = "E:/Package/.rlib_cast", LANG = "C", LC_ALL = "C")
.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
source(file.path(ROOT, "01_attribution_validation.R"), encoding = "UTF-8")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "production"
out <- file.path(ROOT, "attribution_output", paste0("complete_", mode, "_",
  format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")))
dir.create(out, recursive = TRUE)
writeLines(out, file.path(ROOT, "latest_attribution_active.txt"))
manifest <- source_manifest()
if (mode == "smoke") {
  shards <- run_validation(mode)
} else {
  cl <- parallel::makeCluster(6L, outfile = file.path(out, "workers.log"))
  parallel::clusterEvalQ(cl, {
    .libPaths(c("E:/Package/.rlib_cast", .libPaths()))
    Sys.setenv(LANG = "C", LC_ALL = "C")
    source("E:/Package/castSDM_workflows/core_method_validation/01_attribution_validation.R", encoding = "UTF-8")
    NULL
  })
  shards <- tryCatch(unlist(parallel::parLapply(cl,
    c("linear", "nonlinear", "joint", "null", "correlated", "extrapolation"),
    function(s) run_validation("production", s))), finally = parallel::stopCluster(cl))
}
receipts <- lapply(file.path(shards, "receipt.rds"), readRDS)
if (!all(vapply(receipts, function(x) x$status == "completed", logical(1)))) stop("Incomplete scenario shard")
if (!identical(manifest, source_manifest())) stop("Sources changed during validation")
utils::write.csv(data.frame(scenario_run = shards), file.path(out, "scenario_runs.csv"), row.names = FALSE)
for (file in c("advantage_metrics.csv", "reliability_metrics.csv", "numerical_checks.csv", "errors.csv")) {
  d <- do.call(rbind, lapply(file.path(shards, file), read.csv))
  utils::write.csv(d, file.path(out, file), row.names = FALSE)
}
utils::write.csv(manifest, file.path(out, "source_manifest.csv"), row.names = FALSE)
utils::capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
saveRDS(list(status = "completed", mode = mode, source_manifest = manifest, shards,
  expected_datasets = sum(vapply(receipts, `[[`, numeric(1), "expected_datasets")),
  completed_datasets = sum(vapply(receipts, `[[`, numeric(1), "completed_datasets")),
  checks_passed = TRUE, error_count = 0L), file.path(out, "receipt.rds"))
writeLines(out, file.path(ROOT, paste0("latest_attribution_", mode, "completed.txt")))
source(file.path(ROOT, "02_attribution_figures.R"), encoding = "UTF-8")
make_attribution_figures(out)
message("COMPLETE BENCHMARK: ", out)
