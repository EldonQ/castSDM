#' Tile-based Spatial Prediction for Large Rasters
#'
#' Memory-bounded variant of [cast_predict()] for prediction grids whose
#' full pixel count would otherwise exceed available RAM. Reads the
#' covariate raster stack tile by tile, applies a fitted model, and writes
#' the per-model habitat-suitability surface (HSS) to disk via
#' streamed row-block writing. Output buffers are bounded by tile height
#' times grid width, not total grid size. Model/GDAL memory is additional.
#'
#' Per-tile checkpoints are written eagerly under a per-fit directory
#' (`<output_dir>/.cast_tiles_<fingerprint>/`) so an interrupted run can be
#' resumed without recomputing finished tiles, and a new fit can never
#' replay tiles predicted by an older model. The input-content fingerprint
#' requires one additional streaming read pass before prediction.
#' `overwrite = TRUE` bypasses and refreshes the checkpoints.
#'
#' Cells with missing or non-finite covariates stay `NA` in the output rasters -
#' only finite predictor rows are sent to the fitted engine, matching
#' the missing-cell semantics of
#' [cast_ensemble_raster()]. No suitability values are fabricated outside
#' the data mask (e.g. ocean pixels).
#'
#' @param fit A [cast_fit] object.
#' @param raster A `terra::SpatRaster` whose layer names cover the
#'   environmental variables in `fit$env_vars`. May be a path passed to
#'   `terra::rast()`.
#' @param output_dir Directory in which per-model TIFFs and tile
#'   checkpoints are written. Created if it does not exist.
#' @param tile_size Integer of length 1 or 2. Tile dimensions in pixels
#'   (height x width). Default `512`.
#' @param models Character vector. Which fitted models to predict. Default
#'   all.
#' @param overwrite Logical. Overwrite existing output rasters. Default
#'   `FALSE` - existing files skip only when their receipt matches this fit,
#'   input content and output file. Stale or unverified outputs cause an error.
#' @param compression Character. GeoTIFF compression. Default `"LZW"`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list of class `cast_predict_tiled` with components:
#' \describe{
#'   \item{rasters}{Named list of file paths, one per model.}
#'   \item{n_tiles}{Number of tiles processed.}
#'   \item{models}{Character vector of model names.}
#' }
#'
#' @seealso [cast_predict()]
#' @export
cast_predict_tiled <- function(fit, raster,
                               output_dir,
                               tile_size = 512L,
                               models    = NULL,
                               overwrite = FALSE,
                               compression = "LZW",
                               verbose   = TRUE) {
  check_suggested("terra", "for tile-based raster prediction")
  .cast_check_fit(fit)
  if (is.character(raster)) raster <- terra::rast(raster)
  if (!inherits(raster, "SpatRaster")) cli::cli_abort("{.arg raster} must be a SpatRaster or readable path.")
  mdl_names <- models %||% names(fit$models)
  if (!length(mdl_names) || anyNA(mdl_names) || anyDuplicated(mdl_names) ||
      !all(mdl_names %in% names(fit$models))) cli::cli_abort("{.arg models} must name unique fitted models.")
  env_vars <- fit$env_vars
  if (!all(env_vars %in% names(raster))) cli::cli_abort("Raster is missing fitted predictors.")
  raster <- raster[[env_vars]]
  if (!is.numeric(tile_size) || !length(tile_size) %in% c(1L, 2L) ||
      any(!is.finite(tile_size)) || any(tile_size < 1) ||
      any(tile_size != floor(tile_size)) || any(tile_size > .Machine$integer.max)) {
    cli::cli_abort("{.arg tile_size} must be one or two positive whole numbers of pixels.")
  }
  ts <- rep(as.integer(tile_size), length.out = 2L)
  nr <- terra::nrow(raster); nc <- terra::ncol(raster)
  starts <- seq.int(1L, nr, by = ts[1]); columns <- seq.int(1L, nc, by = ts[2])
  n_tiles <- length(starts) * length(columns)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(file.path(output_dir, paste0("HSS_", mdl_names, ".tif")), mdl_names)
  # Content hashes distinguish spatial permutations with the same summary stats.
  input_hash <- vapply(starts, function(row) .cast_digest(.cast_raster_block(raster,
    row, min(ts[1], nr - row + 1L))), character(1))
  grid <- list(extent = as.vector(terra::ext(raster)), dimensions = dim(raster),
    crs = terra::crs(raster), variables = env_vars, tile_size = ts, input = input_hash)
  signatures <- stats::setNames(vapply(mdl_names, function(m) {
    info <- fit$models[[m]]
    if (is.function(info$predict_function)) info$predict_function <- .cast_callback_identity(info$predict_function)
    .cast_digest(list(algorithm = "streamed_tiles_v1", grid = grid, model = info,
      scaling = fit$scaling, predictor = body(predict_single_model)))
  }, character(1)), mdl_names)
  receipt_path <- function(m) paste0(paths[[m]], ".receipt.rds")
  todo <- mdl_names
  for (m in mdl_names[file.exists(paths) & !overwrite]) {
    receipt <- tryCatch(readRDS(receipt_path(m)), error = function(e) NULL)
    if (is.null(receipt) || !identical(receipt$signature, signatures[[m]]) ||
        !identical(receipt$md5, unname(tools::md5sum(paths[[m]])))) {
      cli::cli_abort("Existing raster does not match this fit/input: {.path {paths[[m]]}}. Use a new directory or overwrite = TRUE.")
    }
    todo <- setdiff(todo, m)
  }
  writers <- list()
  on.exit(for (writer in writers) .cast_writer_abort(writer), add = TRUE)
  for (m in todo) writers[[m]] <- .cast_raster_writer(raster, paths[[m]], paste0("HSS_", m),
    compression = compression, overwrite = overwrite)
  checkpoints <- stats::setNames(lapply(todo, function(m) {
    path <- file.path(output_dir, paste0(".cast_tiles_", signatures[[m]]))
    dir.create(path, showWarnings = FALSE)
    path
  }), todo)
  if (length(todo)) for (row in starts) {
    count <- min(ts[1], nr - row + 1L)
    X <- .cast_raster_block(raster, row, count)
    if (!identical(.cast_digest(X), input_hash[match(row, starts)])) {
      cli::cli_abort("Raster input changed during prediction; no completed output is published.")
    }
    stripe <- matrix(NA_real_, nrow(X), length(todo), dimnames = list(NULL, todo))
    for (col in columns) {
      width <- min(ts[2], nc - col + 1L)
      index <- as.vector(t(outer(seq_len(count) - 1L, seq.int(col, col + width - 1L),
                                function(r, c) r * nc + c)))
      tile <- X[index, , drop = FALSE]
      ok <- .cast_finite_rows(tile)
      for (m in todo) {
        file <- file.path(checkpoints[[m]], sprintf("%d_%d.rds", row, col))
        checkpoint <- if (!overwrite && file.exists(file)) tryCatch(readRDS(file), error = function(e) NULL) else NULL
        vals <- if (is.list(checkpoint) &&
          identical(checkpoint$digest, .cast_digest(checkpoint$values))) checkpoint$values else NULL
        if (!is.numeric(vals) || length(vals) != nrow(tile) ||
            any(is.infinite(vals)) || any(!is.na(vals[!ok]))) {
          vals <- rep(NA_real_, nrow(tile))
          if (any(ok)) {
            pred <- tryCatch(predict_single_model(fit$models[[m]], tile[ok, , drop = FALSE]),
              error = function(e) cli::cli_abort("Prediction failed for {.val {m}} in tile {row}/{col}: {conditionMessage(e)}"))
            if (!is.numeric(pred) || length(pred) != sum(ok)) cli::cli_abort("Invalid prediction vector for {.val {m}}.")
            pred[!is.finite(pred)] <- NA_real_
            vals[ok] <- pred
          }
          saveRDS(list(values = vals, digest = .cast_digest(vals)), file)
        }
        stripe[index, m] <- vals
      }
    }
    for (m in todo) .cast_writer_block(writers[[m]], stripe[, m], row, count)
    if (verbose) cli::cli_inform("Prediction rows {row}-{row + count - 1L}/{nr}")
  }
  for (m in todo) {
    .cast_writer_finish(writers[[m]])
    saveRDS(list(signature = signatures[[m]], md5 = unname(tools::md5sum(paths[[m]]))), receipt_path(m))
  }
  structure(list(rasters = paths, n_tiles = n_tiles, models = mdl_names), class = "cast_predict_tiled")
}

#' @export
print.cast_predict_tiled <- function(x, ...) {  cat("<cast_predict_tiled>\n")
  cat("  models  :", paste(x$models, collapse = ", "), "\n")
  cat("  tiles   :", x$n_tiles, "\n")
  cat("  rasters :\n")
  for (m in names(x$rasters))
    cat("    -", m, "->", x$rasters[[m]], "\n")
  invisible(x)
}
