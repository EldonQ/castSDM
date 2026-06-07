#' Tile-based Spatial Prediction for Large Rasters
#'
#' Memory-bounded variant of [cast_predict()] for prediction grids whose
#' full pixel count would otherwise exceed available RAM. Reads the
#' covariate raster stack tile by tile, applies a fitted model, and writes
#' the per-model habitat-suitability surface (HSS) to disk via
#' `terra::writeRaster()`.
#'
#' Per-tile checkpoints (`<output_dir>/.cast_tile_<model>_<tileid>.rds`) are
#' written eagerly so that an interrupted run can be resumed without
#' recomputing already-finished tiles.
#'
#' Tiles are dispatched through [future.apply::future_lapply()] when a
#' `cast_worker_budget` is given, using the Windows-friendly
#' `multisession` (PSOCK) backend.
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
#' @param budget A [cast_worker_budget] (optional). When supplied,
#'   `future_lapply()` is used to evaluate tiles in parallel.
#' @param overwrite Logical. Overwrite existing output rasters. Default
#'   `FALSE` - existing files trigger a `cli_inform()` and skip.
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
#' @seealso [cast_predict()], [cast_worker_budget()]
#' @export
cast_predict_tiled <- function(fit, raster,
                               output_dir,
                               tile_size = 512L,
                               models    = NULL,
                               budget    = NULL,
                               overwrite = FALSE,
                               compression = "LZW",
                               verbose   = TRUE) {
  check_suggested("terra", "for tile-based raster prediction")
  if (inherits(raster, "character")) raster <- terra::rast(raster)
  if (!inherits(raster, "SpatRaster"))
    cli::cli_abort("{.arg raster} must be a SpatRaster or readable path.")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  mdl_names <- models %||% names(fit$models)
  mdl_names <- intersect(mdl_names, names(fit$models))
  if (length(mdl_names) == 0)
    cli::cli_abort("No matching fitted models in fit.")

  env_vars <- fit$env_vars
  missing <- setdiff(env_vars, names(raster))
  if (length(missing) > 0)
    cli::cli_abort("Raster missing required layer{?s}: {.val {missing}}.")

  ts <- if (length(tile_size) == 1L) c(tile_size, tile_size) else tile_size
  nrow_r <- terra::nrow(raster); ncol_r <- terra::ncol(raster)
  row_starts <- seq(1L, nrow_r, by = ts[1])
  col_starts <- seq(1L, ncol_r, by = ts[2])
  tile_ids <- expand.grid(r0 = row_starts, c0 = col_starts,
                          KEEP.OUT.ATTRS = FALSE)
  n_tiles <- nrow(tile_ids)
  if (verbose) cli::cli_inform(
    "cast_predict_tiled: {n_tiles} tiles ({ts[1]}x{ts[2]} px) for {length(mdl_names)} model(s)."
  )

  # Output rasters: one per model, init NA, written tile by tile.
  out_paths <- stats::setNames(
    file.path(output_dir, paste0("HSS_", mdl_names, ".tif")),
    mdl_names
  )
  for (mdl in mdl_names) {
    op <- out_paths[[mdl]]
    if (file.exists(op) && !overwrite) {
      cli::cli_inform("  {basename(op)} exists; skipping (overwrite = FALSE).")
      next
    }
    proto <- terra::rast(raster, nlyrs = 1L)
    names(proto) <- paste0("HSS_", mdl)
    terra::writeRaster(proto, op, overwrite = TRUE,
                       gdal = c(paste0("COMPRESS=", compression),
                                "TILED=YES"),
                       wopt = list(datatype = "FLT4S"))
  }

  do_tile <- function(k) {
    r0 <- tile_ids$r0[k]; c0 <- tile_ids$c0[k]
    nr <- min(ts[1], nrow_r - r0 + 1L)
    nc <- min(ts[2], ncol_r - c0 + 1L)
    ext_tile <- terra::ext(
      terra::xFromCol(raster, c0) - terra::xres(raster) / 2,
      terra::xFromCol(raster, c0 + nc - 1L) + terra::xres(raster) / 2,
      terra::yFromRow(raster, r0 + nr - 1L) - terra::yres(raster) / 2,
      terra::yFromRow(raster, r0) + terra::yres(raster) / 2
    )
    tile_r <- terra::crop(raster, ext_tile)
    df <- terra::as.data.frame(tile_r, xy = TRUE, na.rm = FALSE)
    if (!"x" %in% names(df)) names(df)[1:2] <- c("x", "y")
    names(df)[names(df) == "x"] <- "lon"
    names(df)[names(df) == "y"] <- "lat"
    list(tile_id = k, r0 = r0, c0 = c0, nr = nr, nc = nc, df = df)
  }

  predict_tile <- function(tile, mdl_name) {
    df <- tile$df
    nr <- tile$nr; nc <- tile$nc
    if (nrow(df) == 0L) return(matrix(NA_real_, nr, nc))
    chk_path <- file.path(output_dir,
                          sprintf(".cast_tile_%s_%05d.rds",
                                  mdl_name, tile$tile_id))
    if (file.exists(chk_path)) {
      vals <- tryCatch(readRDS(chk_path), error = function(e) NULL)
      if (!is.null(vals)) return(matrix(vals, nr, nc, byrow = TRUE))
    }
    pred_obj <- tryCatch(
      cast_predict(fit, df, models = mdl_name),
      error = function(e) NULL
    )
    if (is.null(pred_obj)) {
      vals <- rep(NA_real_, nr * nc)
    } else {
      vals <- pred_obj$predictions[[paste0("HSS_", mdl_name)]]
      if (length(vals) != nr * nc) {
        # nrow(df) may differ if cropping clipped pixels; pad NAs.
        v2 <- rep(NA_real_, nr * nc); v2[seq_along(vals)] <- vals
        vals <- v2
      }
    }
    saveRDS(vals, chk_path)
    matrix(vals, nr, nc, byrow = TRUE)
  }

  parallel <- !is.null(budget) &&
    inherits(budget, "cast_worker_budget") && budget$total > 1L &&
    requireNamespace("future.apply", quietly = TRUE)

  for (mdl in mdl_names) {
    op <- out_paths[[mdl]]
    if (file.exists(op) && !overwrite) next
    if (verbose) cli::cli_inform("  predicting tiles for {.val {mdl}}...")

    process_one <- function(k) {
      tile <- do_tile(k)
      list(tile = tile, vals = predict_tile(tile, mdl))
    }

    results <- if (parallel) {
      old_plan <- future::plan(future::multisession,
                               workers = budget$total)
      on.exit(future::plan(old_plan), add = TRUE)
      future.apply::future_lapply(seq_len(n_tiles), process_one,
                                  future.seed = TRUE)
    } else {
      lapply(seq_len(n_tiles), process_one)
    }

    # Write tiles back to the persistent raster.
    out_r <- terra::rast(op)
    for (res in results) {
      tile <- res$tile
      vals <- as.numeric(t(res$vals))
      cells <- terra::cellFromRowColCombine(
        out_r,
        seq(tile$r0, tile$r0 + tile$nr - 1L),
        seq(tile$c0, tile$c0 + tile$nc - 1L)
      )
      out_r[cells] <- vals
    }
    terra::writeRaster(
      out_r, op, overwrite = TRUE,
      gdal = c(paste0("COMPRESS=", compression), "TILED=YES"),
      wopt = list(datatype = "FLT4S")
    )
  }

  out <- list(
    rasters = out_paths,
    n_tiles = n_tiles,
    models  = mdl_names
  )
  class(out) <- "cast_predict_tiled"
  out
}


#' @export
print.cast_predict_tiled <- function(x, ...) {
  cat("<cast_predict_tiled>\n")
  cat("  models  :", paste(x$models, collapse = ", "), "\n")
  cat("  tiles   :", x$n_tiles, "\n")
  cat("  rasters :\n")
  for (m in names(x$rasters))
    cat("    -", m, "->", x$rasters[[m]], "\n")
  invisible(x)
}
