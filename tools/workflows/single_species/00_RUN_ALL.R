ROOT <- "E:/Package/castSDM_workflows/single_species"
source(file.path(ROOT, "01_config.R"), encoding = "UTF-8")
source(file.path(ROOT, "02_pipeline.R"), encoding = "UTF-8")
source(file.path(ROOT, "03_figures.R"), encoding = "UTF-8")

required_workflow_packages()
options(ranger.num.threads = 1L)
taxon <- Sys.getenv("CAST_SPECIES", "Schizopygopsis")
cfg <- castsdm_species_configs()[[taxon]]
if (is.null(cfg)) stop("Unknown CAST_SPECIES: ", taxon)
# Method validation is probed, never required: see core_validation_status().
# Opt in with cfg$require_core_validation = TRUE when the delivery must assert it.
core_validation <- core_validation_status(cfg)
if (isTRUE(cfg$require_core_validation)) {
  if (!isTRUE(core_validation$available)) {
    stop("cfg$require_core_validation is TRUE but the core validation is unavailable: ",
         core_validation$reason)
  }
  core_validation$receipt <- require_core_validation(cfg)
} else if (!isTRUE(core_validation$available)) {
  message("Core-method validation not available (", core_validation$reason,
          "). This does not block the case run: the run records its own source ",
          "provenance, which is what makes it reproducible.")
}
devtools::load_all(cfg$cast_root, quiet = TRUE)
source(cfg$adapter, encoding = "UTF-8")
# Preflight BEFORE any expensive stage: prove the loaded castSDM provides every
# entry point, argument and internal this workflow reaches for, and that it is
# at least cfg$min_cast_version. Without it an API drift surfaces deep inside a
# cached step -- or, worse, as a silently ignored call whose product no longer
# matches the delivered figures. The record goes into the run provenance.
cast_api <- assert_cast_api(cfg)
cfg$output_root <- file.path(ROOT, "output", paste0("production_",
  format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_", Sys.getpid()))
# recursive = TRUE: output/ itself may be absent (a previous delivery may have
# been archived or deleted), and dir.create() without recursion fails on a
# missing parent instead of creating the whole path. The second clause treats
# "already exists" as success and still fails loudly if the path is unusable.
if (!dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE) &&
    !dir.exists(cfg$output_root)) {
  stop("Cannot create a new isolated run directory: ", cfg$output_root)
}
provenance <- record_run_provenance(cfg, core_validation)
# Bind the resolved API surface into the provenance record so every delivered
# bundle states the exact castSDM it was produced against, including the
# version actually loaded by load_all().
provenance$cast_api <- cast_api
saveRDS(provenance, file.path(cfg$output_root, "run_provenance.rds"))
message("\n===== ", cfg$species_name, " -> ", cfg$output_root, " =====")
run_species(cfg)
message("\nCase computation and exports finished; visual QA and acceptance are still required.")
