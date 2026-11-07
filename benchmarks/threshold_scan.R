# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# Rscript benchmarks/threshold_scan.R [existing_comparison.tif] [output_directory]
# Bounded I/O benchmark only. No fitting or ecological projection is performed.
args <- commandArgs(TRUE)
input_path <- if (length(args)) args[1L] else ""
out <- if (length(args) >= 2L) args[2L] else file.path("benchmarks", "output")
source_dir <- Sys.getenv("CAST_SOURCE_DIR", "")
if (nzchar(source_dir)) pkgload::load_all(source_dir, quiet = TRUE) else library(castSDM)
stopifnot(utils::packageVersion("castSDM") >= "0.27.0")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
internal <- function(name) getFromNamespace(name, "castSDM")
cuts <- c(.3, .5, .7)
columns <- c("current", "joint", "climate", "lulc", "support")
if (nzchar(input_path)) {
  r <- terra::rast(input_path)
  stopifnot(all(columns %in% names(r)))
  # One fixed 256-row strip, including the original raster's full width.
  first <- min(1600L, max(1L, terra::nrow(r)-255L))
  nr <- min(256L, terra::nrow(r)-first+1L)
  terra::readStart(r)
  values <- terra::readValues(r, row = first, nrows = nr, mat = TRUE)
  terra::readStop(r)
  template <- terra::rast(nrows = nr, ncols = terra::ncol(r), nlyrs = terra::nlyr(r))
  names(template) <- names(r)
  provenance <- list(source = normalizePath(input_path, winslash = "/"),
    source_md5 = unname(tools::md5sum(input_path)), first_row = first, rows = nr)
} else {
  set.seed(62027)
  template <- terra::rast(nrows = 256, ncols = 2048, nlyrs = 21)
  names(template) <- c("current", "joint", "climate", "lulc",
    internal(".cast_factorial_names")(), paste0("sd_", internal(".cast_factorial_names")()), "support")
  values <- matrix(runif(terra::ncell(template)*21), ncol = 21)
  values[, 21L] <- as.integer(values[, 21L] > .25)
  values[seq.int(1L, nrow(values), by = 7L), ] <- NA_real_
  provenance <- list(source = "synthetic_uniform_seed_62027", first_row = 1, rows = 256)
}
terra::values(template) <- values
rm(values)
paths <- c(pixel = file.path(out, "strip_pixel.tif"), band = file.path(out, "strip_band.tif"))
for (layout in names(paths)) terra::writeRaster(template, paths[[layout]], overwrite = TRUE,
  wopt = list(datatype = "FLT8S", gdal = c("COMPRESS=LZW", "TILED=YES",
    paste0("INTERLEAVE=", toupper(layout)))))
rm(template); gc()

# The 0.26.0 reference: values() per short strip, ~100,000 cells per call.
legacy <- function(path) {
  r <- terra::rast(path)[[columns]]
  counts <- array(0, c(length(cuts), 4L, 2L)); n <- numeric(2L)
  rows <- max(1L, as.integer(100000L/terra::ncol(r)))
  for (first in seq.int(1L, terra::nrow(r), by = rows)) {
    x <- terra::values(r, row = first, nrows = min(rows, terra::nrow(r)-first+1L), mat = TRUE)
    valid <- rowSums(!is.finite(x[, 1:4, drop = FALSE])) == 0L
    supported <- valid & is.finite(x[, 5L]) & x[, 5L] != 0
    n <- n + c(sum(valid), sum(supported))
    for (i in seq_along(cuts)) {
      counts[i, , 1L] <- counts[i, , 1L] + colSums(x[valid, 1:4, drop = FALSE] >= cuts[i])
      counts[i, , 2L] <- counts[i, , 2L] + colSums(x[supported, 1:4, drop = FALSE] >= cuts[i])
    }
  }
  rows <- list()
  for (i in seq_along(cuts)) for (d in 1:2) {
    b <- internal(".cast_factorial_binary")(counts[i, , d], n[d])
    b$threshold <- cuts[i]; b$domain <- c("full", "supported")[d]
    rows[[length(rows)+1L]] <- b
  }
  do.call(rbind, rows)
}
optimized <- function(path) internal(".cast_threshold_scan")(
  list(raster_file = path), cuts, TRUE, FALSE)
methods <- list(legacy_pixel = function() legacy(paths[["pixel"]]),
  persistent_pixel = function() optimized(paths[["pixel"]]),
  persistent_band = function() optimized(paths[["band"]]))
reference <- methods$legacy_pixel()
timings <- list()
for (replicate in 1:3) {
  order <- if (replicate %% 2L) names(methods) else rev(names(methods))
  for (method in order) {
    gc()
    t <- system.time(value <- methods[[method]]())
    stopifnot(isTRUE(all.equal(value, reference, tolerance = 0)))
    timings[[length(timings)+1L]] <- data.frame(replicate, method,
      elapsed_s = unname(t[["elapsed"]]), user_s = unname(t[["user.self"]]),
      system_s = unname(t[["sys.self"]]))
  }
}
allocations <- lapply(names(methods), function(method) {
  trace <- file.path(out, paste0(method, "-alloc.txt"))
  utils::Rprofmem(trace)
  value <- tryCatch(methods[[method]](), finally = utils::Rprofmem(NULL))
  stopifnot(isTRUE(all.equal(value, reference, tolerance = 0)))
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(trace, warn = FALSE))))
  unlink(trace)
  data.frame(method, r_allocated_mb = sum(bytes, na.rm = TRUE)/1e6)
})
timings <- do.call(rbind, timings)
utils::write.csv(timings, file.path(out, "timings.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, allocations), file.path(out, "allocations.csv"), row.names = FALSE)
utils::write.csv(reference, file.path(out, "verified_counts.csv"), row.names = FALSE)
saveRDS(list(provenance = provenance, thresholds = cuts,
  artifacts = data.frame(layout = names(paths), bytes = file.info(paths)$size,
    md5 = unname(tools::md5sum(paths))), sys = Sys.info(),
  gdal = terra::gdal(), session = utils::capture.output(sessionInfo()),
  limitations = "Three alternating warm-filesystem trials on one bounded strip and one machine; R allocations are not peak process RSS. No whole-workflow speed claim."),
  file.path(out, "manifest.rds"))
print(timings)
print(do.call(rbind, allocations))
