# Opt-in regression on a completed real case. No generated measurements/labels.
# Rscript tests/real-data/scenario-io.R <completed_case> <separate_QA_directory>
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Supply a real completed case and a separate QA directory.")
library(castSDM)
case <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2L], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(tolower(qa), "/"), paste0(tolower(case), "/")))
  stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
p <- castSDM:::.cast_explore_project(case)
info <- p$run$cv$fold_info
at <- which(info$obs == 1L)[1L]
stopifnot(is.finite(at))
xy <- unlist(info[at, c("lon", "lat")])
scenarios <- utils::head(names(p$comparison$results), 3L)
pair <- castSDM:::.cast_explore_stacks(p, scenarios[1L])
res <- terra::res(pair$current)
window <- terra::ext(xy[1]-3*res[1], xy[1]+3*res[1], xy[2]-3*res[2], xy[2]+3*res[2])
current_file <- file.path(qa, "observed_current_crop.tif")
terra::writeRaster(terra::crop(pair$current, window), current_file, datatype = "FLT8S", overwrite = TRUE)
current <- terra::rast(current_file)
changing <- as.character(unlist(p$comparison$groups, use.names = FALSE))
future <- stats::setNames(lapply(scenarios, function(s) {
  value <- current
  if (length(changing)) {
    x <- castSDM:::.cast_explore_stacks(p, s)$future[[changing]]
    path <- file.path(qa, paste0(s, "_actual_crop.tif"))
    terra::writeRaster(terra::crop(x, window), path, datatype = "FLT8S", overwrite = TRUE)
    value[[changing]] <- terra::rast(path)
  }
  value
}), scenarios)
compare <- function(path, workers) cast_compare(p$fit, current, future,
  groups = p$comparison$groups, weights = p$comparison$weights,
  threshold = p$comparison$threshold, clamp = p$run$protocol$clamp,
  support_probs = p$comparison$support_probs, output_dir = path,
  resources = list(cores = 3L, workers = workers, memory_mb = 512))
a <- compare(file.path(qa, "serial"), 1L)
b <- compare(file.path(qa, "parallel"), 3L)
for (s in scenarios) stopifnot(isTRUE(all.equal(
  terra::values(terra::rast(a$results[[s]]$raster_file)),
  terra::values(terra::rast(b$results[[s]]$raster_file)), tolerance = 1e-12)))
c <- compare(file.path(qa, "serial"), 1L)
stopifnot(all(c$events$cache_hit), isTRUE(all.equal(a$summary, c$summary)))
payload <- file.path(qa, "cache_payload.rds")
key <- castSDM:::.cast_digest(a$summary)
castSDM:::.cast_cache_write(a$summary, payload, key)
stopifnot(identical(a$summary, castSDM:::.cast_cache_read(payload, key)))
cat("Actual-case serial/parallel, resume and cache-integrity checks passed.\n")

# Same saved coefficients and observed rows must be independent of batch size.
reference <- p$fit$scaling$reference
reference <- reference[seq_len(min(2000L, nrow(reference))), , drop = FALSE]
chunks <- split(seq_len(nrow(reference)), ceiling(seq_len(nrow(reference))/1000L))
for (model in names(p$fit$models)) {
  whole <- castSDM:::predict_single_model(p$fit$models[[model]], reference)
  parts <- unlist(lapply(chunks, function(i)
    castSDM:::predict_single_model(p$fit$models[[model]], reference[i, , drop = FALSE])), use.names = FALSE)
  stopifnot(max(abs(whole-parts)) < 1e-12)
}

# A new comparison provides checked inputs and paired predictions for explanation.
scenario <- tail(scenarios, 1L)
prepared <- castSDM:::.cast_compare_reuse(p$fit, current, future[[scenario]], a,
  a$weights, isTRUE(p$run$protocol$clamp))
stopifnot(!is.null(prepared$blocks))
response <- castSDM:::.cast_attribution(p$fit, current, future[[scenario]],
  drivers = changing, weights = a$weights, clamp = isTRUE(p$run$protocol$clamp),
  support_probs = a$support_probs, prepared = prepared, diagnostics = TRUE,
  filename = file.path(qa, "actual_response.tif"), resources = list(cores = 2L, memory_mb = 512))
projection <- terra::values(terra::rast(a$results[[scenario]]$raster_file))
response_values <- terra::values(terra::rast(response$raster_file))
if (identical(a$design, "paired")) stopifnot(max(abs(
  projection[, c("current", "joint", "joint_change")]-
  response_values[, c("current", "future", "total_change")]), na.rm = TRUE) < 1e-12)
cat("Saved-model batch invariance and paired response reuse passed.\n")
