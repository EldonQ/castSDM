# Small file manifests keep results tied to the source and input artifacts.
workflow_source_signature <- function(cfg) {
  paths <- c(list.files(file.path(cfg$cast_root, "R"), pattern = "\\.[Rr]$", full.names=TRUE),
             file.path(cfg$cast_root, c("DESCRIPTION", "NAMESPACE")), cfg$adapter,
             file.path(dirname(cfg$adapter), "workflow_provenance.R"),
             file.path(dirname(dirname(cfg$adapter)), "single_species",
                       c("00_RUN_ALL.R", "01_config.R", "02_pipeline.R", "03_figures.R")))
  if (any(!file.exists(paths))) stop("Missing source files in provenance manifest")
  hashes <- tools::md5sum(paths)
  list(signature=digest::digest(hashes), hashes=hashes,
       package=as.character(utils::packageVersion("castSDM")))
}

raster_input_manifest <- function(stacks) {
  if (inherits(stacks,"SpatRaster")) stacks <- list(stacks)
  paths <- sort(unique(unlist(lapply(stacks, function(x) as.character(terra::sources(x))))))
  paths <- paths[nzchar(paths)]
  if (!length(paths)) stop("Production raster inputs must have file-backed sources")
  if (any(!file.exists(paths))) stop("A raster source file is missing")
  info <- file.info(paths)
  data.frame(path=normalizePath(paths,winslash="/"),size=info$size,
             modified=as.character(info$mtime),stringsAsFactors=FALSE)
}

require_current_bundle <- function(bundle, cfg) {
  signature <- workflow_source_signature(cfg)$signature
  if (is.null(bundle$provenance$signature) ||
      !identical(bundle$provenance$signature,signature)) {
    stop("Case bundle is missing or has stale source provenance; rerun 00_RUN_ALL.R: ",
         cfg$species_name)
  }
  invisible(TRUE)
}

# Non-fatal probe: what does the core-method validation currently say?
# The case workflow must not depend on this directory existing. Method
# validation is a paper activity with its own inputs and its own figures; a case
# run is a use of the package. Coupling them meant that deleting validation
# artefacts locked every case run, which is the wrong direction of dependency.
# Returns a list with available = FALSE and a reason instead of erroring, so the
# caller can record the state in the run provenance and carry on.
core_validation_status <- function(cfg) {
  root <- file.path(dirname(dirname(cfg$adapter)), "core_method_validation")
  out <- list(available = FALSE, reason = NA_character_, receipt = NULL,
              root = root)
  pointer <- file.path(root, "latest_attribution_productioncompleted.txt")
  if (!file.exists(pointer)) {
    out$reason <- "no completed attribution-validation pointer"
    return(out)
  }
  target <- readLines(pointer, warn = FALSE)
  if (length(target) != 1L || !file.exists(file.path(target, "receipt.rds"))) {
    out$reason <- "completion pointer does not resolve to a receipt"
    return(out)
  }
  receipt <- tryCatch(readRDS(file.path(target, "receipt.rds")), error = function(e) NULL)
  if (is.null(receipt) || !identical(receipt$status, "completed") ||
      !identical(receipt$mode, "production") || !isTRUE(receipt$checks_passed) ||
      !isTRUE(receipt$error_count == 0L) ||
      !isTRUE(receipt$completed_datasets == receipt$expected_datasets)) {
    out$reason <- "attribution validation is incomplete"
    return(out)
  }
  manifest <- receipt$source_manifest
  if (is.null(manifest) || any(!file.exists(manifest$file)) ||
      !identical(as.character(manifest$sha256), unname(vapply(manifest$file,
        function(f) digest::digest(file = f, algo = "sha256"), character(1))))) {
    out$reason <- "completed validation uses different sources"
    return(out)
  }
  out$available <- TRUE
  out$reason <- "completed source-matched attribution validation"
  out$receipt <- receipt
  out$target <- target
  out
}

# What a run needs in order to be reproducible later: the exact sources it was
# built from, the version, the git commit if the tree is versioned, and the state
# of the core validation at run time. Recording is always possible; blocking is
# a policy decision the configuration makes, not a side effect of a directory
# being present.
record_run_provenance <- function(cfg, core_validation = NULL) {
  sig <- workflow_source_signature(cfg)
  cast_root <- cfg$cast_root
  commit <- NA_character_
  if (file.exists(file.path(cast_root, ".git"))) {
    commit <- tryCatch({
      r <- suppressWarnings(system2("git", c("-C", shQuote(cast_root), "rev-parse", "HEAD"),
                                    stdout = TRUE, stderr = FALSE))
      if (length(r)) trimws(r[[1]]) else NA_character_
    }, error = function(e) NA_character_)
  }
  list(
    recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    species = cfg$species_name,
    package_version = sig$package,
    source_signature = sig$signature,
    source_hashes = sig$hashes,
    git_commit = commit,
    core_validation_available = isTRUE(core_validation$available),
    core_validation_reason = if (is.null(core_validation)) "not probed" else core_validation$reason,
    core_validation_run = if (isTRUE(core_validation$available))
      basename(core_validation$target) else NA_character_
  )
}

# Strict variant, kept for the paper side: when preparing a delivery we want to
# assert that the validation on disk was produced by exactly this source. It is
# opt-in through cfg$require_core_validation and is never called implicitly.
require_core_validation <- function(cfg) {
  status <- core_validation_status(cfg)
  if (!isTRUE(status$available)) stop("Complete current-source attribution validation: ", status$reason)
  invisible(status$receipt)
}
