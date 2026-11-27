# Compact, content-verified current GAM contributions shared by future scenarios.
# Predictions use the same fitted, batch-invariant mgcv bases.
.cast_gam_current_cache <- function(fit, weights, variables, blocks, directory, clamp) {
  models <- names(weights)[vapply(names(weights), function(m)
    !is.null(.cast_gam_terms(fit$models[[m]])), logical(1))]
  if (!length(models) || is.null(blocks)) return(NULL)
  key <- .cast_digest(list(.cast_fit_identity(fit), blocks$key, models,
    variables, clamp, .cast_code_identity(.cast_gam_current_values), "gam_current/1"))
  value <- list(directory = directory, key = key, models = models,
    variables = variables, rows = blocks$rows, starts = blocks$starts)
  attr(value, "ready") <- !is.null(.cast_cache_read(file.path(directory, "manifest.rds"), key))
  value
}

.cast_gam_current_values <- function(fit, current, models, variables) {
  stats::setNames(lapply(models, function(m) {
    labels <- .cast_gam_terms(fit$models[[m]])
    cols <- unname(labels[intersect(variables, names(labels))])
    if (!nrow(current)) return(list(eta = numeric(),
      terms = matrix(numeric(), 0L, length(cols), dimnames = list(NULL, cols))))
    a <- .cast_selected_gam_terms(fit$models[[m]]$model, current)
    list(eta = rowSums(a) + attr(a, "constant"), terms = a[, cols, drop = FALSE])
  }), models)
}

.cast_gam_current_path <- function(cache, first) {
  file.path(cache$directory, paste0("rows_", first, ".rds"))
}

.cast_gam_current_read <- function(cache, first, count, indices, required = TRUE) {
  if (is.null(cache)) return(NULL)
  path <- .cast_gam_current_path(cache, first)
  value <- .cast_cache_read(path, .cast_digest(list(cache$key, first, count)))
  pos <- if (!is.null(value)) match(indices, value$indices) else rep(NA_integer_, length(indices))
  if (is.null(value) || anyNA(pos) || !identical(names(value$models), cache$models)) {
    if (!required) return(NULL)
    cli::cli_abort("Current GAM contribution block failed integrity checks: {.path {path}}. Rerun the comparison to rebuild it.")
  }
  lapply(value$models, function(x) list(eta = x$eta[pos], terms = x$terms[pos, , drop = FALSE]))
}

.cast_gam_current_write <- function(cache, first, count, indices, models) {
  if (!is.null(cache)) .cast_cache_write(list(indices = indices, models = models),
    .cast_gam_current_path(cache, first), .cast_digest(list(cache$key, first, count)))
}

.cast_gam_current_finish <- function(cache) {
  if (is.null(cache)) return(invisible(NULL))
  paths <- vapply(cache$starts, function(first) .cast_gam_current_path(cache, first), character(1))
  .cast_cache_write(cache, file.path(cache$directory, "manifest.rds"), cache$key, paths)
  invisible(NULL)
}
