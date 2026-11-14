# Verify frozen-model optimisations on an existing actual case; no study run.
# Rscript tests/real-data/maintenance.R <actual_case> <writable_QA_directory>
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Supply an actual case and a separate QA directory.")
suppressPackageStartupMessages(library(castSDM))
case <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(qa, "/"), paste0(case, "/"))) stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
r <- readRDS(file.path(case, "training/result.rds"))
p <- readRDS(file.path(case, "training/stages/01_prepare.rds"))
cc <- readRDS(file.path(case, "comparison/comparison.rds"))
cfg <- readRDS(file.path(case, "case.rds"))
fit <- r$fit_full
options(castSDM.predict_threads = 2L)
pool <- read.csv(file.path(case, "training_pool.csv"))
paths <- pool$path[match(fit$env_vars, pool$layer)]
stopifnot(!anyNA(paths), all(file.exists(paths)))
report <- read.csv(file.path(case, "scenario_table.csv"))
report <- report[report$scenario == cfg$report_scenario, , drop = FALSE]
future_paths <- paths
future_paths[match(report$layer, fit$env_vars)] <- report$path
# A native-grid crop centred on a genuinely observed training presence.
pres <- which(p$train$presence == 1L)
xy <- unlist(p$train[pres[ceiling(length(pres)/2)], c("lon", "lat")])
full <- terra::rast(paths); names(full) <- fit$env_vars
cell <- terra::cellFromXY(full, matrix(xy, nrow = 1L))
rc <- terra::rowColFromCell(full, cell)
resolution <- terra::res(full)
window <- terra::ext(xy[1]-3*resolution[1], xy[1]+3*resolution[1],
  xy[2]-3*resolution[2], xy[2]+3*resolution[2])
current <- terra::crop(full, window)
future <- terra::crop(terra::rast(future_paths), window); names(future) <- fit$env_vars
source <- terra::crop(terra::rast(cc$results[[cfg$report_scenario]]$raster_file)[["current"]], window)
mask <- terra::ifel(is.finite(source), 1, NA)
a <- as.data.frame(terra::values(current)); b <- as.data.frame(terra::values(future))
names(a) <- names(b) <- fit$env_vars
keep <- complete.cases(a) & complete.cases(b) & is.finite(terra::values(mask, mat = FALSE))
a <- a[keep, , drop = FALSE]; b <- b[keep, , drop = FALSE]
stopifnot(nrow(a) > 0L)
drivers <- unlist(cc$groups, use.names = FALSE)
expected <- cast_attribution(fit, a, b, drivers = drivers, weights = cc$weights,
  clamp = r$protocol$clamp, verbose = FALSE, diagnostics = TRUE)
if ("gam" %in% names(fit$models)) {
  generic <- fit
  generic$models$gam$type <- "external"
  generic$models$gam$predict_function <- function(model, newdata)
    as.numeric(stats::predict(model, newdata, type = "response"))
  reference <- cast_attribution(generic, a, b, drivers = drivers, weights = cc$weights,
    clamp = r$protocol$clamp, verbose = FALSE, diagnostics = TRUE)
  stopifnot(isTRUE(all.equal(expected$data, reference$data, tolerance = 1e-12)))
}
# Fault injection exercises cancellation, using actual computed raster batches.
.cast_real_check_writes <- 0L
trace(".cast_cache_write", where = asNamespace("castSDM"), print = FALSE,
  exit = quote({
    if (grepl("^block_[0-9]+_[0-9]+[.]rds$", basename(path))) {
      .GlobalEnv$.cast_real_check_writes <- .GlobalEnv$.cast_real_check_writes + 1L
      if (.GlobalEnv$.cast_real_check_writes == 2L) stop("QA checkpoint cancellation")
    }
  }))
output <- tempfile("actual_attribution_", tmpdir = qa, fileext = ".tif")
job <- function() cast_attribution(fit, current, future, drivers = drivers,
  weights = cc$weights, clamp = r$protocol$clamp, diagnostics = TRUE,
  mask = mask, filename = output, block_size = ncol(current), overwrite = TRUE,
  resources = list(cores = 2L, workers = 1L, memory_mb = 512))
cancelled <- tryCatch(job(), error = identity)
untrace(".cast_cache_write", where = asNamespace("castSDM"))
stopifnot(inherits(cancelled, "error"), grepl("QA checkpoint cancellation", conditionMessage(cancelled)))
actual <- job()
z <- terra::values(terra::rast(actual$raster_file))
stopifnot(isTRUE(all.equal(unname(z[keep, , drop = FALSE]),
  unname(as.matrix(expected$data[, actual$layers])), tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(actual$summary, expected$summary, tolerance = 1e-12)))
# Basic plots must never call a fitted engine, PDP or ALE.
trace("predict_single_model", where = asNamespace("castSDM"), print = FALSE,
  tracer = quote(stop("Unexpected engine prediction during basic plots")))
r$output_dir <- file.path(case, "training")
f <- cast_figures(r, cc, file.path(qa, "basic"), include = c("07", "11", "13"),
  dpi = 80, formats = "png", basemap = "china", verbose = TRUE)
untrace("predict_single_model", where = asNamespace("castSDM"))
stopifnot(all(f$status$status == "written_not_visually_checked"))
# Published ALE runs on actual saved reference rows and preserves centred units.
variable <- fit$env_vars[vapply(fit$scaling$reference, function(x)
  length(unique(x[is.finite(x)])) > 1L, logical(1))][1L]
curve <- cast_response_curves(fit, variables = variable, type = "ale", models = names(fit$models)[1L])
stopifnot(identical(attr(curve, "response_type"), "ale"), all(is.finite(curve$prediction)))
invisible(ggplot2::ggplot_build(plot(curve)))
cat("Actual crop parity, interrupted/resumed attribution, prediction-free figures and ALE passed.\n")
