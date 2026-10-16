# Same-task comparison with an explicit iml prediction loop, not an estimator contest.
.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
pkgload::load_all("E:/Package/cast", quiet = TRUE)
root <- "E:/Package/castSDM_workflows/core_method_validation"
output <- file.path(root, "workflow_contract_comparison")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("iml", quietly = TRUE), requireNamespace("terra", quietly = TRUE))
set.seed(1610)
reference <- data.frame(temperature = runif(1000, -2, 2), moisture = runif(1000, -1, 1))
y <- rbinom(1000, 1, plogis(reference$temperature * reference$moisture + reference$temperature))
model <- stats::glm(y ~ temperature * moisture, reference, family = stats::binomial())
callback <- function(model, newdata) as.numeric(stats::predict(model, newdata, type = "response"))
explainer <- iml::Predictor$new(model, data = reference, y = y, predict.function = callback)
fit <- cast_as_fit(list(iml = explainer), reference,
  function(model, newdata) as.numeric(model$predict(newdata)[[1]]))
xy <- expand.grid(lon = seq(-1, 1, length.out = 50), lat = seq(-1, 1, length.out = 50))
current <- transform(xy, temperature = 1.5 * lon, moisture = .8 * lat)
future <- current
future$temperature <- future$temperature + .3 + .2 * future$lat
future$moisture <- future$moisture + .15 * future$lon
vars <- names(reference)
manual <- function(a, b) {
  predict <- function(x) as.numeric(explainer$predict(x[, vars, drop = FALSE])[[1]])
  p0 <- predict(a); p1 <- predict(b)
  result <- data.frame(current = p0, future = p1, total_change = p1 - p0)
  for (v in vars) {
    only <- a; only[[v]] <- b[[v]]
    held <- b; held[[v]] <- a[[v]]
    result[[paste0("effect_", v)]] <- predict(only) - p0
    result[[paste0("retained_", v)]] <- p1 - predict(held)
  }
  result$remainder <- result$total_change - rowSums(result[paste0("effect_", vars)])
  result
}
expected <- manual(current, future)
result <- cast_attribution(fit, current, future, verbose = FALSE)
cols <- names(expected)
difference <- function(x) max(abs(as.matrix(x[, cols]) - as.matrix(expected)))
stopifnot(difference(result$data) < 1e-12)

# Both aligned implementations are correct. iml does not promise spatial pairing.
permutation <- sample.int(nrow(future))
shuffled <- future[permutation, ]
matched <- match(paste(current$lon, current$lat), paste(shuffled$lon, shuffled$lat))
stopifnot(!anyNA(matched))
explicitly_aligned <- manual(current, shuffled[matched, ])
automatically_aligned <- cast_attribution(fit, current, shuffled, verbose = FALSE)
naive_error <- difference(manual(current, shuffled))
stopifnot(difference(explicitly_aligned) < 1e-12,
  difference(automatically_aligned$data) < 1e-12, naive_error > .01)
checks <- data.frame(
  check = c("aligned_iml_loop_equivalence", "explicitly_aligned_iml_loop",
    "cast_coordinate_alignment", "naive_rowwise_loop_after_shuffle"),
  max_abs_difference = c(difference(result$data), difference(explicitly_aligned),
    difference(automatically_aligned$data), naive_error),
  interpretation = c(rep("equivalent_correct_comparison", 3),
    "missing_spatial_alignment_in_naive_example_not_an_iml_bug"))

r <- terra::rast(nrows = 50, ncols = 50, nlyrs = 2,
  xmin = -1, xmax = 1, ymin = -1, ymax = 1, crs = "EPSG:4326")
names(r) <- vars
terra::values(r) <- as.matrix(current[, vars])
s <- r; terra::values(s) <- as.matrix(future[, vars])
raster <- cast_attribution(fit, r, s, block_size = 450, verbose = FALSE)
mapped <- as.data.frame(terra::values(terra::rast(raster$raster_file)))
stopifnot(difference(mapped) < 1e-12)
rejected <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
geometry_rejected <- rejected(cast_attribution(fit, r, terra::shift(s, dx = .04), verbose = FALSE))
unmatched <- future; unmatched$lon[1] <- 99
coordinates_rejected <- rejected(cast_attribution(fit, current, unmatched, verbose = FALSE))
stopifnot(geometry_rejected, coordinates_rejected)
guards <- data.frame(check = c("shifted_raster_geometry", "unmatched_coordinates"),
  rejected = c(geometry_rejected, coordinates_rejected))

# Five repeated executions are timing repeats, not independent scientific datasets.
timings <- do.call(rbind, lapply(seq_len(5), function(i) {
  t1 <- system.time(manual(current, future))[["elapsed"]]
  t2 <- system.time(cast_attribution(fit, current, future, verbose = FALSE))[["elapsed"]]
  data.frame(repetition = i, rows = nrow(current), changing_predictors = length(vars),
    iml_manual_seconds = t1, cast_attribution_seconds = t2)
}))
utils::write.csv(checks, file.path(output, "same_task_equivalence.csv"), row.names = FALSE)
utils::write.csv(guards, file.path(output, "spatial_contract_checks.csv"), row.names = FALSE)
utils::write.csv(timings, file.path(output, "execution_times.csv"), row.names = FALSE)
files <- c(list.files("E:/Package/cast/R", "[.]R$", full.names = TRUE),
  "E:/Package/cast/DESCRIPTION", "E:/Package/cast/NAMESPACE",
  file.path(root, "05_workflow_contract_comparison.R"))
utils::write.csv(data.frame(file = files, sha256 = vapply(files,
  function(f) digest::digest(file = f, algo = "sha256"), character(1))),
  file.path(output, "source_manifest.csv"), row.names = FALSE)
utils::write.csv(data.frame(package = c("castSDM", "iml", "terra"),
  version = vapply(c("castSDM", "iml", "terra"), function(p)
    as.character(utils::packageVersion(p)), character(1))),
  file.path(output, "package_versions.csv"), row.names = FALSE)
utils::capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
writeLines(c("# Same-task workflow comparison", "",
  "A correctly aligned iml prediction loop and castSDM reproduce the same paired contrasts.",
  "castSDM also checks raster geometry and pairs table coordinates automatically.",
  "The naive rowwise example omits alignment; this is not a defect in iml.",
  "The raster result agrees with the table result across block boundaries.",
  "Five timing repetitions quantify execution cost only; no speed, usability or memory advantage is inferred.",
  "This uses one external GLM and two changing predictors, not an ALE/SHAP or package-wide benchmark."),
  file.path(output, "RESULTS.md"))
print(checks); print(guards); print(timings)
