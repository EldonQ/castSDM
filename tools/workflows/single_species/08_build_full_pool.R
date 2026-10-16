# 08_build_full_pool.R ----------------------------------------------------------
# Build the CURRENT-pool manifest (pool_layers_full.csv) from every PRESENT
# covariate under the N-SDM data root, with TWO quality gates:
#
# Gate 1 (scenario exclusion): any layer whose path contains a /scenario/
#   directory segment is a future-scenario product and must never enter the
#   current pool. This catches the bioclim GCM/ensemble products AND the
#   lulc/glu1km scenario tree (2030/2050/2100 x ssp) -- the 20261001 run
#   failed with "Presences: 2" because 486 glu1km scenario layers plus 81
#   low-coverage layers were admitted and their NA holes blanked 99.9% of
#   occurrence cells during environmental extraction.
#
# Gate 2 (coverage gate): every admitted layer must have finite values on
#   >= 95% of the China land mask cells. The audit numbers come from
#   pool_coverage_audit.csv (written by 10_coverage_audit.R, full-grid
#   terra::global per layer); layers below threshold are listed and excluded.
#   The distribution is cleanly bimodal (1018 layers >= 0.98 vs 81 < 0.95),
#   so the threshold is insensitive in [0.95, 0.98].
#
# Layer naming: bioclim present keeps the bare variable id (bio01 ... tmin)
# because nsdm_future_bioclim_path() resolves future files by that exact id;
# every other layer takes its file stem with the "glo_" prefix dropped.
#
# Geometry gate: metadata-only compareGeom vs mask_china.tif.
#
# Compatibility check: all 68 layers of the legacy pool must reappear.

root <- "E:/CHECO26/N-SDM/data"
glo <- file.path(root, "covariates", "glo")
groups <- c("anthro", "bioclim", "dist", "hydro", "lulc", "soil", "topo", "vege")
out_csv <- "E:/Package/castSDM_workflows/single_species/pool_layers_full.csv"
audit_csv <- "E:/Package/castSDM_workflows/single_species/pool_coverage_audit.csv"
COVER_MIN <- 0.95

rows <- list()
for (g in groups) {
  tifs <- list.files(file.path(glo, g), pattern = "[.]tif$",
                     full.names = TRUE, recursive = TRUE)
  tifs <- tifs[!grepl("/scenario/", tifs, fixed = TRUE)]
  for (p in tifs) {
    stem <- sub("[.]tif$", "", basename(p))
    if (identical(g, "bioclim")) {
      layer <- basename(dirname(p))
      is_climate <- TRUE
    } else {
      layer <- sub("^glo_", "", stem)
      is_climate <- FALSE
    }
    rows[[length(rows) + 1L]] <- data.frame(group = g, layer = layer, path = p,
                                            is_climate = is_climate,
                                            stringsAsFactors = FALSE)
  }
}
pool <- do.call(rbind, rows)
message("Present layers after scenario exclusion: ", nrow(pool))

dup <- pool$layer[duplicated(pool$layer)]
if (length(dup)) stop("Duplicate layer names in candidate pool: ",
                      paste(unique(dup), collapse = ", "))

message("Geometry gate (metadata-only compareGeom vs mask_china.tif) ...")
ref <- terra::rast(file.path(root, "mask_china.tif"))
ok <- logical(nrow(pool))
for (i in seq_len(nrow(pool))) {
  ok[i] <- isTRUE(tryCatch(
    terra::compareGeom(ref, terra::rast(pool$path[i]),
                       stopOnError = FALSE, messages = FALSE),
    error = function(e) FALSE))
  if (i %% 200L == 0L) message("  ", i, "/", nrow(pool))
}
bad_geom <- pool[!ok, , drop = FALSE]
pool <- pool[ok, , drop = FALSE]
message("Geometry-aligned: ", nrow(pool), " | mismatched: ", nrow(bad_geom))

message("Coverage gate (finite share on land mask >= ", COVER_MIN, ") ...")
audit <- utils::read.csv(audit_csv, stringsAsFactors = FALSE)
m <- match(pool$layer, audit$layer)
if (anyNA(m)) stop("Coverage audit is missing layers from the candidate pool: ",
                   paste(pool$layer[is.na(m)], collapse = ", "))
pool$cover <- audit$cover[m]
low_cov <- pool$cover < COVER_MIN
bad_cov <- pool[low_cov, c("group", "layer", "path", "cover"), drop = FALSE]
pool <- pool[!low_cov, , drop = FALSE]
message("Coverage-passed: ", nrow(pool), " | excluded: ", nrow(bad_cov))

legacy <- utils::read.csv("E:/Package/castSDM_workflows/single_species/pool_layers.csv",
                          stringsAsFactors = FALSE)
legacy <- legacy[legacy$status %in% c("kept", NA), , drop = FALSE]
m2 <- match(normalizePath(legacy$path, winslash = "/"),
            normalizePath(pool$path, winslash = "/"))
# The guard above is a v060-era invariant ("every legacy layer must still
# exist"). It is still worth keeping, because an UNEXPLAINED disappearance means
# the pool silently shrank. A deliberate removal is a decision, not a fault, so
# it is declared in pool_removed_layers.csv and subtracted here; anything left
# over still aborts the build.
removed_file <- "E:/Package/castSDM_workflows/single_species/pool_removed_layers.csv"
# fileEncoding = "UTF-8-BOM": a BOM makes read.csv name the first column
# "\ufefffile", so $file comes back NULL and every declared removal looks
# undeclared. Editors that add a BOM (Excel among them) are a normal source of
# these files here, so the reader absorbs it rather than the author avoiding it.
declared <- if (file.exists(removed_file)) {
  d <- utils::read.csv(removed_file, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8-BOM")
  if (!"file" %in% names(d)) character() else trimws(d$file)
} else character()
gap <- legacy$path[is.na(m2)]
gap_unexplained <- gap[!(basename(gap) %in% declared)]
if (length(gap_unexplained)) {
  stop("Legacy pool layers missing (or coverage-excluded) in the rebuilt pool: ",
       paste(basename(gap_unexplained), collapse = ", "),
       "\n  If their removal was deliberate, add the file names to ",
       basename(removed_file), " with a reason.")
}
n_removed <- sum(basename(gap) %in% declared)
message(sprintf("Legacy pool: %d layers, %d reappear, %d declared removals (see %s)",
                length(legacy$path), sum(!is.na(m2)), n_removed, basename(removed_file)))
if (n_removed) {
  message("  declared removals: ", paste(basename(gap[nzchar(basename(gap) %in% declared)]),
                                         collapse = ", "))
}

pool$status <- "kept"
pool <- pool[, c("group", "layer", "path", "is_climate", "status", "cover")]
pool <- pool[order(pool$group, pool$layer), , drop = FALSE]
rownames(pool) <- NULL
utils::write.csv(pool, out_csv, row.names = FALSE)
excl <- bad_cov[, c("group", "layer", "path", "cover"), drop = FALSE]
if (nrow(bad_geom)) {
  excl <- rbind(transform(bad_geom, cover = NA_real_)[, c("group", "layer", "path", "cover")],
                excl)
}
utils::write.csv(excl, sub("[.]csv$", "_excluded.csv", out_csv), row.names = FALSE)

cat("\nrebuilt pool: ", nrow(pool), " layers -> ", out_csv, "\n", sep = "")
print(table(pool$group))
cat("is_climate TRUE: ", sum(pool$is_climate), "\n", sep = "")
cat("n_stack_layers for 01_config.R: ", nrow(pool), "\n", sep = "")
