# Reuse is optional and strictly tied to the frozen model and physical inputs.
# This is internal plumbing; public attribution calls need no cache settings.
.cast_scenario_bands <- function(paired, models) {
  layers <- c(if (paired) c("current", "joint") else c("current", "joint", "climate", "lulc"),
    if (paired) "joint_change" else .cast_factorial_names())
  effects <- if (paired) "joint_change" else .cast_factorial_names()
  valid <- c(layers, if (paired)
    unlist(lapply(c("current_", "future_", "change_"), paste0, models)))
  expected <- c(layers, paste0("sd_", effects), "support", "common_support",
    "sign_agreement_joint", "sign_conflict_joint", setdiff(valid, layers))
  list(layers = expected, finite = valid)
}

.cast_scenario_readback <- function(value, path, template, models) {
  bands <- .cast_scenario_bands(identical(value$design, "paired"), models)
  n <- unique(value$summary$n_valid)
  if (length(n) != 1L || !is.finite(n)) cli::cli_abort("Saved scenario lacks valid-cell counts.")
  .cast_raster_validate(path, template, bands$layers, valid_layers = bands$finite, n_valid = n)
}

.cast_compare_input_receipt <- function(receipt, key, template, bands) {
  value <- .cast_cache_read(receipt, key)
  if (is.null(value) || isTRUE(attr(value, "raster_verified"))) return(value)
  paths <- file.path(dirname(receipt), unlist(value, use.names = FALSE))
  ok <- tryCatch({
    for (i in seq_along(paths)) .cast_raster_validate(paths[i], template, bands[[i]])
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  attr(value, "raster_verified") <- TRUE
  .cast_cache_write(value, receipt, key, paths)
  value
}

.cast_legacy_scenarios <- function(path, model, current, future, mask, groups, weights,
    threshold, support_probs, clamp, domain, gam_key) {
  if (!file.exists(path)) return(NULL)
  old <- tryCatch(readRDS(path), error = function(e) NULL)
  proven <- old$provenance
  if (!inherits(old, "cast_comparison") || !identical(proven$schema, "cast_comparison/2") ||
      !any(grepl("castSDM_0[.]44[.]0([^0-9.]|$)", proven$session)) ||
      !identical(old$prediction_policy, "regular_gam/1") ||
      !identical(proven$model, model) || !identical(proven$current, current) ||
      !identical(proven$future, future) || !identical(proven$mask, mask) ||
      !identical(old$groups, groups) || !identical(old$weights, weights) ||
      !identical(old$threshold, threshold) || !identical(old$support_probs, support_probs) ||
      !identical(old$clamp, clamp) || !identical(old$domain, domain) ||
      !identical(old$gam_current$key, gam_key)) return(NULL)
  old$results
}

.cast_scenario_resume <- function(receipt, key, file, template, models, legacy_key = NULL) {
  value <- .cast_cache_read(receipt, key)
  if (is.null(file)) return(list(value = value, repair = FALSE))
  if (!file.exists(receipt)) return(list(value = NULL, repair = FALSE))
  if (is.null(value) && !is.null(legacy_key)) value <- .cast_cache_read(receipt, legacy_key)
  owner <- tryCatch(readRDS(receipt)$key, error = function(e) NULL)
  if (is.null(value) || !identical(value$artifact, basename(file))) {
    # Authorize replacement only when a receipt identifies this exact run.
    return(list(value = NULL, repair = !is.null(owner) && owner %in% c(key, legacy_key)))
  }
  needs_validation <- !isTRUE(value$raster_verified)
  if (needs_validation) {
    error <- tryCatch({ .cast_scenario_readback(value, file, template, models); NULL },
      error = conditionMessage)
    if (!is.null(error)) return(list(value = NULL, repair = TRUE, reason = error))
    value$raster_verified <- TRUE
  }
  if (needs_validation || !identical(owner, key)) .cast_cache_write(value, receipt, key, file)
  list(value = value, repair = FALSE)
}

.cast_compare_reuse <- function(fit, current, future, comparison, weights, clamp,
                                memo = new.env(parent = emptyenv())) {
  if (!inherits(comparison, "cast_comparison") || !inherits(current, "SpatRaster") ||
      !inherits(future, "SpatRaster")) return(NULL)
  current_key <- .cast_data_identity(current, memo)
  future_key <- .cast_data_identity(future, memo)
  if (!identical(current_key, comparison$provenance$current)) return(NULL)
  scenario <- names(comparison$provenance$future)[match(future_key, comparison$provenance$future)]
  if (is.na(scenario)) return(NULL)
  blocks <- comparison$current_blocks
  if (!is.null(blocks)) {
    directory <- file.path(comparison$output_dir, ".cache", "current_blocks")
    manifest <- .cast_cache_read(file.path(directory, "manifest.rds"), blocks$key)
    expected <- .cast_current_blocks_key(list(current_key, comparison$provenance$mask),
      blocks$rows, fit$env_vars)
    if (is.null(manifest) || !identical(blocks$key, expected)) blocks <- NULL else {
      blocks <- manifest; blocks$directory <- directory
    }
  }
  reuse_model <- identical(comparison$clamp, clamp) && identical(comparison$weights, weights) &&
    identical(comparison$provenance$model, .cast_fit_identity(fit))
  prediction <- if (reuse_model) comparison$results[[scenario]]$raster_file
  layers <- unlist(lapply(c("current_", "future_"), paste0, names(weights)))
  if (!is.character(prediction) || !file.exists(prediction) ||
      !all(layers %in% names(terra::rast(prediction)))) prediction <- NULL
  if (!is.null(prediction)) {
    receipt <- sub("[.]tif$", ".rds", prediction)
    expected <- comparison$results[[scenario]]$cache_key
    result <- if (!is.null(expected)) .cast_cache_read(receipt, expected)
    if (is.null(result)) prediction <- NULL else if (!isTRUE(result$raster_verified))
      .cast_scenario_readback(result, prediction, current, names(weights))
  }
  gam <- if (reuse_model && !is.null(blocks) && !is.null(comparison$gam_current))
    .cast_gam_current_cache(fit, weights, comparison$changing_predictors, blocks,
      file.path(comparison$output_dir, ".cache", "gam_current"), clamp)
  if (!isTRUE(attr(gam, "ready"))) gam <- NULL
  list(blocks = blocks, gam_current = gam, prediction = prediction, layers = layers,
    current_key = current_key, future_key = future_key, mask_key = comparison$provenance$mask)
}

.cast_reuse_predictions <- function(prepared, first, count, positions, weights) {
  if (is.null(prepared$prediction)) return(NULL)
  value <- .cast_raster_block(terra::rast(prepared$prediction)[[prepared$layers]], first, count)
  value <- as.matrix(value[positions, , drop = FALSE])
  k <- length(weights)
  list(current = value[, seq_len(k), drop = FALSE],
    future = value[, k+seq_len(k), drop = FALSE])
}
