# Opt-in checks on saved real predictions. No fitting or manufactured data.
# Rscript tests/real-data/display.R <completed_case> <separate_QA_dir> [national]
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Supply a completed real case and a separate QA directory.")
library(castSDM)
case <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2L], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(tolower(qa), "/"), paste0(tolower(case), "/")))
  stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
path <- file.path(case, "figures", ".data", "figure_data.rds")
receipt <- readRDS(path)
d <- castSDM:::.cast_cache_read(path, receipt$key)
stopifnot(!is.null(d), identical(d$comparison$domain, "common"))
att <- d$attribution
source <- att$raster_file
before <- tools::md5sum(source)
r <- terra::rast(source)
national <- length(args) >= 3L && identical(args[3L], "national")
at <- which(terra::values(r[["support_total"]], mat = FALSE) == 1L)[1L]
stopifnot(is.finite(at))
xy <- terra::xyFromCell(r, at)
width <- 20*terra::res(r)
window <- terra::ext(xy[1]-width[1], xy[1]+width[1], xy[2]-width[2], xy[2]+width[2])
if (!national) {
  file <- file.path(qa, "observed_response_crop.tif")
  terra::writeRaster(terra::crop(r, window), file, overwrite = TRUE, datatype = "FLT8S")
  att$raster_file <- file
  r <- terra::rast(file)
}
out <- file.path(qa, "driver_supported.tif")
castSDM:::.cast_figure_driver_display(att, out, .02)
display <- terra::rast(out)
coverage <- utils::read.csv(paste0(out, ".csv"))
contrasts <- c("total_change", paste0("effect_", att$drivers), paste0("retained_", att$drivers), "remainder")
support <- terra::values(r[["support_total"]], mat = FALSE)
inside <- which(support == 1)
x <- as.matrix(terra::extract(r[[contrasts]], inside))
ok <- rowSums(!is.finite(x)) == 0L
stopifnot(coverage$native_valid_cells[1L] == sum(ok),
  all(coverage$mask_order == "native mask before aggregation"))
if (national) stopifnot(all(att$common_summary$n == sum(ok)),
  coverage$native_centres_represented_fraction[1L] > .99)
# Independent native-mask-first oracle, using the same observed effects.
v <- att$drivers[1L]
support[] <- NA_real_
support[inside[ok]] <- 1
mask <- terra::setValues(terra::rast(r[[1L]]), support)
oracle <- terra::project(terra::mask(r[[paste0("effect_", v)]], mask), display[[1L]], method = "average")
stopifnot(max(abs(terra::values(oracle)-terra::values(display[[paste0("effect_", v)]])),
  na.rm = TRUE) < 1e-7)
for (v in att$drivers) {
  layer <- paste0("context_direction_", v)
  if (!layer %in% names(display)) next
  flags <- terra::values(r[[c(paste0("support_effect_", v), paste0("support_retained_", v))]])
  at <- which(flags[, 1L] == 1 & flags[, 2L] == 1)
  n <- sum(is.finite(terra::extract(r[[layer]], at)[[1L]]))
  stopifnot(coverage$native_valid_cells[coverage$reference_layer == layer] == n)
}
stopifnot(identical(before, tools::md5sum(source)))
cat("Native masks, common-summary denominator, per-driver scopes and unchanged source verified.\n")

# Shared baseline/display coverage must equal independent reprojections.
cache <- new.env(parent = emptyenv())
old <- options(castSDM.figure_display_cache = cache)
results <- utils::head(d$comparison$results, 2L)
for (i in seq_along(results)) {
  input <- terra::rast(results[[i]]$raster_file)
  input <- terra::crop(input, window)
  a <- file.path(qa, paste0("shared_", i, ".tif"))
  b <- file.path(qa, paste0("independent_", i, ".tif"))
  castSDM:::.cast_figure_display(input, a, .02, shared_key = "actual_common_domain")
  options(castSDM.figure_display_cache = NULL)
  castSDM:::.cast_figure_display(input, b, .02)
  options(castSDM.figure_display_cache = cache)
  stopifnot(isTRUE(all.equal(terra::values(terra::rast(a)), terra::values(terra::rast(b)), tolerance = 0)))
  flag <- terra::values(terra::rast(a)[["support"]], mat = FALSE)
  stopifnot(all(flag[is.finite(flag)] %in% c(0, 1, 2)))
  bound <- castSDM:::.cast_figure_bound(a, "joint_change", TRUE)
  stopifnot(identical(bound, castSDM:::.cast_figure_bound(a, "joint_change", TRUE)))
}
stopifnot(length(grep("^template_", ls(cache))) == 1L,
  identical(read.csv(file.path(qa, "shared_1.tif.csv")), read.csv(file.path(qa, "shared_2.tif.csv"))))
options(old)
cat("Shared display baseline/coverage equals independent output; mixed support and template/bound reuse verified.\n")
