# Reuse is optional and strictly tied to the frozen model and physical inputs.
# This is internal plumbing; public attribution calls need no cache settings.
.cast_compare_reuse <- function(fit, current, future, comparison, weights, clamp) {
  if (!inherits(comparison, "cast_comparison") || !inherits(current, "SpatRaster") ||
      !inherits(future, "SpatRaster") || !identical(comparison$clamp, clamp) ||
      !identical(comparison$weights, weights) ||
      !identical(comparison$provenance$model, .cast_fit_identity(fit))) return(NULL)
  memo <- new.env(parent = emptyenv())
  current_key <- .cast_data_identity(current, memo)
  future_key <- .cast_data_identity(future, memo)
  if (!identical(current_key, comparison$provenance$current)) return(NULL)
  scenario <- names(comparison$provenance$future)[match(future_key, comparison$provenance$future)]
  if (is.na(scenario)) return(NULL)
  blocks <- comparison$current_blocks
  if (!is.null(blocks)) {
    directory <- file.path(comparison$output_dir, ".cache", "current_blocks")
    manifest <- .cast_cache_read(file.path(directory, "manifest.rds"), blocks$key)
    expected <- .cast_digest(list(list(current_key, comparison$provenance$mask),
      blocks$rows, fit$env_vars, "current_blocks/1", .cast_code_identity(.cast_raster_block)))
    if (is.null(manifest) || !identical(blocks$key, expected)) blocks <- NULL else {
      blocks <- manifest; blocks$directory <- directory
    }
  }
  prediction <- comparison$results[[scenario]]$raster_file
  layers <- unlist(lapply(c("current_", "future_"), paste0, names(weights)))
  if (!is.character(prediction) || !file.exists(prediction) ||
      !all(layers %in% names(terra::rast(prediction)))) prediction <- NULL
  if (!is.null(prediction)) {
    receipt <- sub("[.]tif$", ".rds", prediction)
    saved <- tryCatch(readRDS(receipt), error = function(e) NULL)
    if (is.null(saved) || is.null(.cast_cache_read(receipt, saved$key))) prediction <- NULL
  }
  gam <- if (!is.null(blocks) && !is.null(comparison$gam_current))
    .cast_gam_current_cache(fit, weights, comparison$changing_predictors, blocks,
      file.path(comparison$output_dir, ".cache", "gam_current"), clamp)
  if (!isTRUE(attr(gam, "ready"))) gam <- NULL
  list(blocks = blocks, gam_current = gam, prediction = prediction, layers = layers,
    current_key = current_key, future_key = future_key)
}

.cast_reuse_predictions <- function(prepared, first, count, positions, weights) {
  if (is.null(prepared$prediction)) return(NULL)
  value <- .cast_raster_block(terra::rast(prepared$prediction)[[prepared$layers]], first, count)
  value <- as.matrix(value[positions, , drop = FALSE])
  k <- length(weights)
  list(current = value[, seq_len(k), drop = FALSE],
    future = value[, k+seq_len(k), drop = FALSE])
}
