#' Tile-based Spatial Prediction for Large Rasters
#'
#' Memory-bounded variant of [cast_predict()] for prediction grids whose
#' full pixel count would otherwise exceed available RAM. Reads the
#' covariate raster stack tile by tile, applies a fitted model, and writes
#' the per-model habitat-suitability surface (HSS) to disk via
#' `terra::writeRaster()`.
#'
#' Per-tile checkpoints are written eagerly under a per-fit directory
#' (`<output_dir>/.cast_tiles_<fingerprint>/`) so an interrupted run can be
#' resumed without recomputing finished tiles, and a new fit can never
#' replay tiles predicted by an older model. `overwrite = TRUE` bypasses
#' and refreshes the checkpoints.
#'
#' Cells with missing (`NA`) covariates stay `NA` in the output rasters -
#' they are recorded before [cast_predict()] imputes them and written back
#' as `NA` afterwards, matching the NA semantics of
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
  if (inherits(raster, "character")) raster <- terra::rast(raster)
  if (!inherits(raster, "SpatRaster"))
    cli::cli_abort("{.arg raster} must be a SpatRaster or readable path.")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  mdl_names <- if (is.null(models)) names(fit$models) else models
  mdl_names <- intersect(mdl_names, names(fit$models))
  if (length(mdl_names) == 0)
    cli::cli_abort("No matching fitted models in fit.")

  env_vars <- fit$env_vars
  missing <- setdiff(env_vars, names(raster))
  if (length(missing) > 0)
    cli::cli_abort("Raster missing required layer{?s}: {.val {missing}}.")

  ts <- if (length(tile_size) == 1L) c(tile_size, tile_size) else tile_size
  # Unvalidated, `tile_size = 0` produced `seq(1, n, by = 0)` (an R-internal
  # error) and a length-3 vector silently dropped its third element.
  if (!is.numeric(ts) || !length(ts) %in% c(1L, 2L) ||
      any(!is.finite(ts)) || any(ts < 1)) {
    cli::cli_abort(
      "{.arg tile_size} must be one or two positive whole numbers of pixels."
    )
  }
  ts <- as.integer(ts)
  nrow_r <- terra::nrow(raster); ncol_r <- terra::ncol(raster)
  row_starts <- seq(1L, nrow_r, by = ts[1])
  col_starts <- seq(1L, ncol_r, by = ts[2])
  tile_ids <- expand.grid(r0 = row_starts, c0 = col_starts,
                          KEEP.OUT.ATTRS = FALSE)
  n_tiles <- nrow(tile_ids)
  if (verbose) cli::cli_inform(
    "cast_predict_tiled: {n_tiles} tiles ({ts[1]}x{ts[2]} px) for {length(mdl_names)} model(s)."
  )

  # Output rasters: one per model, written once at the end of the sweep.
  # `pre_existing` records which outputs existed BEFORE this call. No NA
  # prototype is written up front: an interrupted run used to leave a
  # pure-NA TIFF that the next `overwrite = FALSE` run accepted as complete.
  out_paths <- stats::setNames(
    file.path(output_dir, paste0("HSS_", mdl_names, ".tif")),
    mdl_names
  )
  pre_existing <- stats::setNames(file.exists(out_paths), mdl_names)
  for (mdl in mdl_names) {
    if (pre_existing[[mdl]] && !overwrite) {
      cli::cli_inform("  {basename(out_paths[[mdl]])} exists; skipping (overwrite = FALSE).")
    }
  }

  do_tile <- function(k, rast) {
    r0 <- tile_ids$r0[k]; c0 <- tile_ids$c0[k]
    nr <- min(ts[1], nrow_r - r0 + 1L)
    nc <- min(ts[2], ncol_r - c0 + 1L)
    ext_tile <- terra::ext(
      terra::xFromCol(rast, c0) - terra::xres(rast) / 2,
      terra::xFromCol(rast, c0 + nc - 1L) + terra::xres(rast) / 2,
      terra::yFromRow(rast, r0 + nr - 1L) - terra::yres(rast) / 2,
      terra::yFromRow(rast, r0) + terra::yres(rast) / 2
    )
    tile_r <- terra::crop(rast, ext_tile)
    df <- terra::as.data.frame(tile_r, xy = TRUE, na.rm = FALSE)
    if (!"x" %in% names(df)) names(df)[1:2] <- c("x", "y")
    names(df)[names(df) == "x"] <- "lon"
    names(df)[names(df) == "y"] <- "lat"
    list(tile_id = k, r0 = r0, c0 = c0, nr = nr, nc = nc, df = df)
  }

  # Tile checkpoints live in a per-fit directory so a new fit can never
  # replay predictions of an older model, and `overwrite = TRUE` bypasses
  # (and refreshes) the checkpoints instead of silently reusing them.
  predict_tile <- function(tile, mdl_name, chk_root) {
    df <- tile$df
    nr <- tile$nr; nc <- tile$nc
    if (nrow(df) == 0L) return(matrix(NA_real_, nr, nc))
    # `chk_root` is NA when fingerprinting failed: checkpoints are then
    # disabled rather than shared, so no tile can cross models.
    chk_path <- if (is.na(chk_root)) NULL else
      file.path(chk_root, sprintf("%s_%05d.rds", mdl_name, tile$tile_id))
    if (!is.null(chk_path) && !isTRUE(overwrite) && file.exists(chk_path)) {
      vals <- tryCatch(readRDS(chk_path), error = function(e) NULL)
      if (!is.null(vals) && length(vals) == nr * nc) {
        return(matrix(vals, nr, nc, byrow = TRUE))
      }
    }
    pred_obj <- tryCatch(
      cast_predict(fit, df, models = mdl_name),
      error = function(e) NULL
    )
    if (is.null(pred_obj)) {
      vals <- rep(NA_real_, nr * nc)
    } else {
      vals <- pred_obj$predictions[[paste0("HSS_", mdl_name)]]
      # Write NA back over cells with NA covariates: cast_predict() imputes
      # them with training medians, but those pixels must stay NA instead
      # of carrying fabricated suitabilities (cast_ensemble_raster NA
      # semantics).
      na_rows <- !stats::complete.cases(df[, env_vars, drop = FALSE])
      if (length(vals) == length(na_rows)) vals[na_rows] <- NA_real_
      if (length(vals) != nr * nc) {
        # nrow(df) may differ if cropping clipped pixels; pad NAs.
        v2 <- rep(NA_real_, nr * nc); v2[seq_along(vals)] <- vals
        vals <- v2
      }
    }
    if (!is.null(chk_path)) {
      tryCatch(saveRDS(vals, chk_path), error = function(e) NULL)
    }
    matrix(vals, nr, nc, byrow = TRUE)
  }

  # Tiles are processed sequentially: ranger/gbm predictions dominate the
  # cost and a serial sweep bounds peak memory to one tile at a time.
  #
  # The model loop is INSIDE the tile loop on purpose. The previous order
  # (model outside, tile inside) re-read the whole predictor stack from disk
  # once per model, because `do_tile()` crops `raster` on every call -- with
  # four models that is four redundant passes over a multi-gigabyte stack.
  # Reading each tile once and predicting every model from it cuts that I/O
  # to a single pass. Peak memory still holds one tile at a time: each tile
  # keeps one prediction vector per model, not one full raster.
  todo <- mdl_names[!pre_existing | overwrite]
  if (!length(todo)) {
    cli::cli_inform("All output rasters exist; nothing to do (overwrite = FALSE).")
  }

  # Fit fingerprints (one per model) so a new fit can never replay tiles
  # predicted by an older model.
  mdl_sig <- stats::setNames(
    vapply(mdl_names, function(m) {
      tryCatch(.cast_digest(fit$models[[m]]$model, len = 10),
               error = function(e) NA_character_)
    }, character(1)),
    mdl_names
  )
  if (anyNA(mdl_sig)) {
    # A digest failure must not collapse distinct models onto one shared
    # checkpoint root, which would let tiles cross models.
    cli::cli_warn(c(
      "Could not fingerprint {sum(is.na(mdl_sig))} fitted model{?s};",
      i = "checkpoints are disabled for those models so no tile can be replayed across models."
    ))
  }

  # Raster fingerprint: geometry + layer names + summary stats, so a different
  # environmental stack can never silently replay old tiles. Computed ONCE
  # (it does not depend on the model) instead of once per model.
  rast_sig <- tryCatch(
    .cast_digest(list(
      as.vector(terra::ext(raster)), dim(raster), names(raster),
      terra::global(raster, c("mean", "min", "max"), na.rm = TRUE)
    ), len = 10),
    error = function(e) NA_character_
  )
  if (is.na(rast_sig)) {
    cli::cli_warn(
      "Could not fingerprint the raster stack; checkpoints are disabled for this run."
    )
  }

  chk_root_of <- function(m) {
    if (is.na(mdl_sig[[m]]) || is.na(rast_sig)) return(NA_character_)
    root <- file.path(output_dir, paste0(".cast_tiles_", mdl_sig[[m]], "_", rast_sig))
    dir.create(root, showWarnings = FALSE, recursive = TRUE)
    root
  }

  if (length(todo)) {
    n_cells_out <- as.double(nrow_r) * as.double(ncol_r)
    # One accumulator per model, filled cell by cell.
    acc <- stats::setNames(
      lapply(todo, function(m) rep(NA_real_, n_cells_out)),
      todo
    )
    # Output templates share the INPUT stack geometry (never re-opening the
    # output file, which would lock it on Windows). `setNames()` is required:
    # `todo` is an unnamed character vector, so `lapply()` returns an UNNAMED
    # list and `templates[[mdl]]` would be NULL.
    templates <- stats::setNames(
      lapply(todo, function(m) {
        tpl <- terra::rast(raster, nlyrs = 1L)
        names(tpl) <- paste0("HSS_", m)
        tpl
      }),
      todo
    )

    for (k in seq_len(n_tiles)) {
      # Read the tile ONCE and reuse it for every model.
      tile <- do_tile(k, raster)
      for (mdl in todo) {
        vals <- as.numeric(t(predict_tile(tile, mdl, chk_root_of(mdl))))
        cells <- terra::cellFromRowColCombine(
          templates[[mdl]],
          seq(tile$r0, tile$r0 + tile$nr - 1L),
          seq(tile$c0, tile$c0 + tile$nc - 1L)
        )
        acc[[mdl]][cells] <- vals
        rm(vals, cells)
      }
      rm(tile)
      if (k %% 20 == 0) invisible(gc())
    }

    for (mdl in todo) {
      if (verbose) cli::cli_inform("  writing tiles for {.val {mdl}}...")
      # Atomic write: a truncated GeoTIFF left behind by an interrupted run
      # used to be accepted by the next `overwrite = FALSE` run as complete.
      .cast_write_raster_atomic(
        terra::setValues(templates[[mdl]], acc[[mdl]]),
        out_paths[[mdl]], compression = compression
      )
    }
    rm(acc, templates)
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
print.cast_predict_tiled <- function(x, ...) {  cat("<cast_predict_tiled>\n")
  cat("  models  :", paste(x$models, collapse = ", "), "\n")
  cat("  tiles   :", x$n_tiles, "\n")
  cat("  rasters :\n")
  for (m in names(x$rasters))
    cat("    -", m, "->", x$rasters[[m]], "\n")
  invisible(x)
}
