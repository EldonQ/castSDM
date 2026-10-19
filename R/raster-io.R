# Shared row-block I/O. Model-specific calculations stay in their callers.
.cast_finite_rows <- function(x) rowSums(!is.finite(as.matrix(x))) == 0L

# Hash callback dependencies rather than unrelated objects in its enclosing
# environment (which can change merely because a caller stores the result).
.cast_callback_identity <- function(fun, depth = 0L) {
  if (depth > 5L) cli::cli_abort("Prediction callback dependencies are recursive; caching cannot fingerprint them.")
  if (is.primitive(fun)) return(fun)
  env <- environment(fun)
  if (isNamespace(env) || identical(env, baseenv())) {
    return(list(body = body(fun), formals = formals(fun), environment = environmentName(env)))
  }
  globals <- codetools::findGlobals(fun, merge = TRUE)
  values <- lapply(globals, function(name) {
    if (!exists(name, envir = env, inherits = TRUE)) return(NULL)
    value <- get(name, envir = env, inherits = TRUE)
    if (is.function(value)) .cast_callback_identity(value, depth + 1L) else value
  })
  list(body = body(fun), formals = formals(fun), dependencies = stats::setNames(values, globals))
}

.cast_raster_mask <- function(mask, template) {
  if (is.null(mask)) return(NULL)
  if (!inherits(mask, "SpatRaster")) cli::cli_abort("{.arg mask} must be a SpatRaster.")
  if (terra::nlyr(mask) != 1L) cli::cli_abort("{.arg mask} must have exactly one layer.")
  if (!isTRUE(terra::compareGeom(template, mask, lyrs = FALSE, stopOnError = FALSE))) {
    cli::cli_abort(c("Mask geometry must match the prediction grid.",
      i = "Align the mask explicitly with nearest-neighbour resampling before prediction."))
  }
  mask
}

.cast_raster_block <- function(r, row, nrows) {
  values <- terra::values(r, row = row, nrows = nrows, mat = TRUE)
  if (nrow(values) != nrows * terra::ncol(r)) {
    cli::cli_abort("Raster block read returned an unexpected number of cells.")
  }
  as.data.frame(values, check.names = FALSE)
}

.cast_mask_rows <- function(mask, row, nrows, n) {
  if (is.null(mask)) return(rep(TRUE, n))
  m <- tryCatch(terra::values(mask, row = row, nrows = nrows, mat = FALSE),
    error = function(e) cli::cli_abort("Mask read failed: {conditionMessage(e)}"))
  if (length(m) != n) cli::cli_abort("Mask block read returned an unexpected number of cells.")
  is.finite(m) & m != 0
}

# Publish only a closed sibling file. Never copy into a visible final path.
# On Windows an existing destination must first move aside; restore it if
# publication fails. This is recoverable replacement, not a multi-file transaction.
.cast_publish_raster <- function(temp, path, overwrite = FALSE) {
  if (!file.exists(temp)) cli::cli_abort("Completed temporary raster is missing.")
  backup <- NULL
  if (file.exists(path)) {
    if (!overwrite) cli::cli_abort("Raster output already exists: {.path {path}}.")
    backup <- tempfile(".cast-backup-", tmpdir = dirname(path), fileext = ".tif")
    if (!file.rename(path, backup)) cli::cli_abort("Cannot replace locked raster: {.path {path}}.")
  }
  if (!file.rename(temp, path)) {
    restored <- is.null(backup) || file.rename(backup, path)
    if (!restored) cli::cli_abort("Raster publication failed; original output retained at {.path {backup}}.")
    cli::cli_abort("Cannot publish completed raster: {.path {path}}.")
  }
  if (!is.null(backup)) unlink(backup)
  invisible(path)
}

.cast_raster_writer <- function(template, path, layers, datatype = "FLT4S",
                                 compression = "LZW", overwrite = FALSE) {
  if (file.exists(path) && !overwrite) cli::cli_abort("Raster output already exists: {.path {path}}.")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer <- new.env(parent = emptyenv())
  writer$path <- path
  writer$temp <- tempfile(".cast-raster-", tmpdir = dirname(path), fileext = ".tif")
  writer$raster <- terra::rast(template[[1]], nlyrs = length(layers))
  names(writer$raster) <- layers
  writer$open <- FALSE
  writer$overwrite <- overwrite
  tryCatch({
    terra::writeStart(writer$raster, writer$temp, overwrite = TRUE,
      wopt = list(datatype = datatype, gdal = c(paste0("COMPRESS=", compression), "TILED=YES")))
    writer$open <- TRUE
  }, error = function(e) {
    unlink(writer$temp)
    stop(e)
  })
  writer
}

.cast_writer_block <- function(writer, values, row, nrows) {
  terra::writeValues(writer$raster, values, row, nrows)
  invisible(NULL)
}

.cast_writer_finish <- function(writer) {
  terra::writeStop(writer$raster)
  writer$open <- FALSE
  # Release the temporary dataset before renaming it on Windows.
  writer$raster <- NULL
  .cast_publish_raster(writer$temp, writer$path, writer$overwrite)
  invisible(writer$path)
}

.cast_writer_abort <- function(writer) {
  if (isTRUE(writer$open)) try(terra::writeStop(writer$raster), silent = TRUE)
  writer$raster <- NULL
  writer$open <- FALSE
  if (file.exists(writer$temp)) unlink(writer$temp)
  invisible(NULL)
}
