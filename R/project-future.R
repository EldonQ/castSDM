#' Project Species Distribution Under Future Climate Scenarios
#'
#' Applies the same ensemble model(s) and weights from current-climate
#' predictions to one or more future climate scenarios, then computes
#' range change statistics (gain, loss, stable).
#'
#' @param fit A [cast_fit] object.
#' @param cv A [cast_cv] object for computing ensemble weights and
#'   thresholds.
#' @param current_env A `data.frame` with `lon`, `lat`, and environmental
#'   variables for the current (baseline) period.
#' @param future_envs A named list of `data.frame`s, each with the same
#'   `lon`, `lat` grid and environmental variables for a future scenario
#'   (e.g., `list(ssp245_2070 = df1, ssp585_2070 = df2)`).
#' @param method Character. Ensemble strategy: `"weighted"` (default),
#'   `"best"`, `"equal"`. See [cast_ensemble()].
#' @param models Character vector. Models to use. Default `NULL` (all).
#' @param save_dir Optional character. When provided, saves CSV and (if
#'   `terra` is available) GeoTIFF outputs to this directory:
#'   - `current_prediction.csv` / `.tif`
#'   - `{scenario}_prediction.csv` / `.tif`
#'   - `{scenario}_change.csv` / `.tif`
#'   - `projection_stats.csv`
#'   Default `NULL` (no file output).
#'
#' @return A `cast_project` object with components:
#' \describe{
#'   \item{current}{A `cast_ensemble` object for the current climate.}
#'   \item{future}{Named list of `cast_ensemble` objects for each scenario.}
#'   \item{changes}{Named list of change `data.frame`s with `lon`, `lat`,
#'     `change` columns.}
#'   \item{stats}{A `data.frame` with columns `scenario`, `n_gain`,
#'     `n_loss`, `n_stable_present`, `n_stable_absent`, `pct_change`,
#'     `centroid_shift_km`.}
#' }
#'
#' @details
#' Change categories per grid cell:
#' - **gain**: absent now, present under future scenario.
#' - **loss**: present now, absent under future scenario.
#' - **stable_present**: present in both.
#' - **stable_absent**: absent in both.
#'
#' Every future grid must correspond row-for-row to `current_env`. A grid
#' with a different number of rows aborts; a grid with the same rows in a
#' different order is re-aligned by matching on `lon`/`lat`. A scenario whose
#' prediction fails is skipped with a warning and an `NA` statistics row.
#'
#' Centroid shift is computed as the great-circle distance between the
#' weighted centroid of current presence and the weighted centroid of
#' future presence. GeoTIFF outputs infer their resolution from the
#' `lon`/`lat` spacing and are written as EPSG:4326.
#'
#' @seealso [cast_ensemble()], [cast_fit()], [cast_cv()]
#'
#' @export
cast_project <- function(fit, cv, current_env, future_envs,
                         method = c("weighted", "best", "equal"),
                         models = NULL,
                         save_dir = NULL) {
  method <- match.arg(method)

  if (!is.list(future_envs) || length(future_envs) == 0) {
    cli::cli_abort("{.arg future_envs} must be a non-empty named list of data.frames.")
  }
  if (is.null(names(future_envs)) || any(names(future_envs) == "")) {
    cli::cli_abort("All elements of {.arg future_envs} must be named (scenario names).")
  }

  # Align every future grid to the current grid before predicting (M16):
  # mismatched row counts abort, re-ordered grids are matched on lon/lat.
  for (scen in names(future_envs)) {
    future_envs[[scen]] <- .cast_align_future_env(current_env,
                                                  future_envs[[scen]], scen)
  }

  # ---- Current prediction -------------------------------------------------
  current <- cast_ensemble(fit, cv, current_env,
                           method = method,
                           models = models)

  # ---- Future predictions -------------------------------------------------
  future_list <- list()
  changes_list <- list()
  stats_rows <- list()
  cur_bin <- current$predictions$binary_ensemble
  has_coords <- all(c("lon", "lat") %in% names(current$predictions))

  for (scen in names(future_envs)) {
    sc_res <- tryCatch({
      fut_env <- future_envs[[scen]]

      # Predict using the same ensemble configuration
      fut <- cast_ensemble(fit, cv, fut_env,
                           method = method,
                           models = models)

      # ---- Compute change map -----------------------------------------------
      fut_bin <- fut$predictions$binary_ensemble

      change <- .cast_change_classes(cur_bin, fut_bin)

      change_df <- if (has_coords) {
        data.frame(
          lon = current$predictions$lon,
          lat = current$predictions$lat,
          change = change,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          site = seq_along(change),
          change = change,
          stringsAsFactors = FALSE
        )
      }

      # ---- Summary statistics -----------------------------------------------
      n_gain   <- sum(change == "gain", na.rm = TRUE)
      n_loss   <- sum(change == "loss", na.rm = TRUE)
      n_stable <- sum(change == "stable_present", na.rm = TRUE)
      n_absent <- sum(change == "stable_absent", na.rm = TRUE)
      total_present_now <- sum(cur_bin == 1, na.rm = TRUE)
      pct_change <- if (total_present_now > 0) {
        100 * (n_gain - n_loss) / total_present_now
      } else {
        NA_real_
      }

      # Centroid shift (great-circle distance in km)
      centroid_km <- NA_real_
      if (has_coords) {
        centroid_km <- tryCatch({
          ok_cur <- !is.na(cur_bin) & cur_bin == 1
          ok_fut <- !is.na(fut_bin) & fut_bin == 1
          cur_pres <- current$predictions[ok_cur, ]
          fut_pres <- fut$predictions[ok_fut, ]
          if (nrow(cur_pres) > 0 && nrow(fut_pres) > 0) {
            # Weight by HSS for more stable centroids
            w_cur <- cur_pres$hss_ensemble
            w_fut <- fut_pres$hss_ensemble
            c_lon <- stats::weighted.mean(cur_pres$lon, w_cur)
            c_lat <- stats::weighted.mean(cur_pres$lat, w_cur)
            f_lon <- stats::weighted.mean(fut_pres$lon, w_fut)
            f_lat <- stats::weighted.mean(fut_pres$lat, w_fut)
            .haversine_km(c_lat, c_lon, f_lat, f_lon)
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_)
      }

      list(
        fut = fut,
        change = change_df,
        stats = data.frame(
          scenario        = scen,
          n_gain          = n_gain,
          n_loss          = n_loss,
          n_stable_present = n_stable,
          n_stable_absent  = n_absent,
          pct_change      = round(pct_change, 2),
          centroid_shift_km = round(centroid_km, 1),
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
      cli::cli_warn("Scenario {.val {scen}} failed: {conditionMessage(e)}")
      NULL
    })

    if (is.null(sc_res)) {
      stats_rows[[scen]] <- data.frame(
        scenario = scen, n_gain = NA_integer_, n_loss = NA_integer_,
        n_stable_present = NA_integer_, n_stable_absent = NA_integer_,
        pct_change = NA_real_, centroid_shift_km = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      future_list[[scen]] <- sc_res$fut
      changes_list[[scen]] <- sc_res$change
      stats_rows[[scen]] <- sc_res$stats
    }
  }

  stats_df <- do.call(rbind, stats_rows)
  rownames(stats_df) <- NULL

  # ---- Save outputs (optional) ---------------------------------------------
  if (!is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

    # CSV outputs
    utils::write.csv(
      current$predictions, file.path(save_dir, "current_prediction.csv"),
      row.names = FALSE
    )
    utils::write.csv(stats_df, file.path(save_dir, "projection_stats.csv"),
                     row.names = FALSE)
    for (scen in names(future_list)) {
      utils::write.csv(
        future_list[[scen]]$predictions,
        file.path(save_dir, paste0(scen, "_prediction.csv")),
        row.names = FALSE
      )
      utils::write.csv(
        changes_list[[scen]],
        file.path(save_dir, paste0(scen, "_change.csv")),
        row.names = FALSE
      )
    }

    # GeoTIFF outputs (if terra available and data has coordinates)
    has_coords <- all(c("lon", "lat") %in% names(current$predictions))
    if (has_coords && requireNamespace("terra", quietly = TRUE)) {
      .save_prediction_tif(
        current$predictions, "lon", "lat", "hss_ensemble",
        file.path(save_dir, "current_prediction.tif")
      )
      for (scen in names(future_list)) {
        .save_prediction_tif(
          future_list[[scen]]$predictions, "lon", "lat", "hss_ensemble",
          file.path(save_dir, paste0(scen, "_prediction.tif"))
        )
        # Change map as numeric: gain=1, loss=-1, stable_present=2, stable_absent=0
        ch <- changes_list[[scen]]
        ch$change_num <- ifelse(ch$change == "gain", 1L,
                         ifelse(ch$change == "loss", -1L,
                         ifelse(ch$change == "stable_present", 2L, 0L)))
        .save_prediction_tif(
          ch, "lon", "lat", "change_num",
          file.path(save_dir, paste0(scen, "_change.tif"))
        )
      }
    }

    cli::cli_inform("Projection outputs saved to {.path {save_dir}}.")
  }

  new_cast_project(
    current = current,
    future  = future_list,
    changes = changes_list,
    stats   = stats_df
  )
}


#' Validate that a future grid corresponds row-for-row to the current grid
#'
#' Aborts on different row counts; re-aligns a grid holding the same rows in
#' a different order by matching on `lon`/`lat`. Silently recycling or
#' pairing mismatched rows would corrupt every change map (review M16).
#'
#' @param current_env Current-grid data.frame.
#' @param fut_env Future-grid data.frame for scenario `scen`.
#' @param scen Scenario name (for error messages).
#' @return The (possibly re-ordered) future grid.
#' @keywords internal
#' @noRd
.cast_align_future_env <- function(current_env, fut_env, scen) {
  if (!is.data.frame(fut_env)) {
    cli::cli_abort("Scenario {.val {scen}}: future grid must be a data.frame.")
  }
  if (nrow(fut_env) != nrow(current_env)) {
    cli::cli_abort(c(
      "Scenario {.val {scen}}: future grid must have one row per current-grid row.",
      "x" = "Got {nrow(fut_env)} future rows vs {nrow(current_env)} current rows."
    ))
  }
  has_ll <- all(c("lon", "lat") %in% names(current_env)) &&
    all(c("lon", "lat") %in% names(fut_env))
  if (has_ll && (!identical(fut_env$lon, current_env$lon) ||
                 !identical(fut_env$lat, current_env$lat))) {
    key_cur <- paste(current_env$lon, current_env$lat)
    key_fut <- paste(fut_env$lon, fut_env$lat)
    if (anyDuplicated(key_cur) || anyDuplicated(key_fut)) {
      cli::cli_abort(c(
        "Scenario {.val {scen}}: {.val lon}/{.val lat} pairs are not unique.",
        i = "Coordinate matching is ambiguous on duplicated cells; supply unique grid coordinates."
      ))
    }
    idx <- match(key_cur, key_fut)
    if (any(is.na(idx))) {
      cli::cli_abort(
        "Scenario {.val {scen}}: future {.val lon}/{.val lat} rows cannot be matched to the current grid."
      )
    }
    fut_env <- fut_env[idx, , drop = FALSE]
    rownames(fut_env) <- NULL
  }
  fut_env
}


#' Classify current vs future binary predictions into change categories
#'
#' Maps each cell to `"gain"`, `"loss"`, `"stable_present"`, or
#' `"stable_absent"`; cells with a non-binary value on either side (failed
#' model, masked cell) receive `NA` rather than an empty-string class.
#'
#' @param cur_bin,fut_bin Integer vectors (0/1, possibly with `NA`).
#' @return Character vector of change classes.
#' @keywords internal
#' @noRd
.cast_change_classes <- function(cur_bin, fut_bin) {
  if (length(fut_bin) != length(cur_bin)) {
    cli::cli_abort(
      "Internal: current and future binary vectors differ in length ({length(cur_bin)} vs {length(fut_bin)})."
    )
  }
  change <- character(length(cur_bin))
  change[cur_bin == 0 & fut_bin == 1] <- "gain"
  change[cur_bin == 1 & fut_bin == 0] <- "loss"
  change[cur_bin == 1 & fut_bin == 1] <- "stable_present"
  change[cur_bin == 0 & fut_bin == 0] <- "stable_absent"
  unknown <- !(cur_bin %in% c(0, 1)) | !(fut_bin %in% c(0, 1))
  change[unknown] <- NA_character_
  change
}


#' Haversine distance in km
#' @keywords internal
#' @noRd
.haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371  # Earth radius in km
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  R * 2 * atan2(sqrt(a), sqrt(1 - a))
}


#' Save a prediction data.frame as GeoTIFF via terra
#'
#' Resolution is inferred from the smallest `lon`/`lat` spacing (no
#' hard-coded grid); `crs` defaults to EPSG:4326, matching the lon/lat
#' prediction grids produced by the package. Cells receiving more than one
#' point (duplicate coordinates) are averaged with `fun = "mean"`; the
#' regular one-point-per-cell path writes values directly.
#'
#' @param df Data.frame with coordinate and value columns.
#' @param lon_col,lat_col,val_col Column names in `df`.
#' @param path Output `.tif` path.
#' @param crs CRS string for the output raster. Default `"EPSG:4326"`.
#' @param res Optional `c(xres, yres)`; inferred from the data when `NULL`.
#' @keywords internal
#' @noRd
.save_prediction_tif <- function(df, lon_col, lat_col, val_col, path,
                                 crs = "EPSG:4326", res = NULL) {
  tryCatch({
    lon <- df[[lon_col]]
    lat <- df[[lat_col]]
    if (is.null(res)) {
      ux <- sort(unique(lon))
      uy <- sort(unique(lat))
      rx <- if (length(ux) > 1L) min(diff(ux)) else 1
      ry <- if (length(uy) > 1L) min(diff(uy)) else 1
      res <- c(rx, ry)
    }
    e <- terra::ext(
      min(lon) - res[1] / 2, max(lon) + res[1] / 2,
      min(lat) - res[2] / 2, max(lat) + res[2] / 2
    )
    tmpl <- terra::rast(e, res = res, crs = crs)
    xy <- cbind(lon, lat)
    cell_id <- terra::cellFromXY(tmpl, xy)
    if (anyDuplicated(cell_id)) {
      # Several points in one cell: aggregate (mean) instead of dropping.
      pts <- terra::vect(xy, type = "points",
                         atts = data.frame(value = df[[val_col]]), crs = crs)
      r <- terra::rasterize(pts, tmpl, field = "value", fun = "mean")
    } else {
      r <- terra::rast(tmpl)
      r[cell_id] <- df[[val_col]]
    }
    terra::writeRaster(r, path, overwrite = TRUE)
  }, error = function(e) {
    cli::cli_warn("Failed to save GeoTIFF {.path {path}}: {e$message}")
  })
}


#' Raster-Based Future Climate Projection
#'
#' Projects species distribution under multiple future climate scenarios
#' using raster inputs and the ensemble model. Generates HSS rasters,
#' binary rasters, change-class rasters, and summary statistics for each
#' scenario. This is the raster-native equivalent of [cast_project()].
#'
#' @param fit A [cast_fit] object.
#' @param cv A [cast_cv] object for computing ensemble weights and
#'   thresholds.
#' @param current_raster A `terra::SpatRaster` stack of current
#'   environmental variables (layer names must match `fit$env_vars`).
#' @param future_rasters A named list of `terra::SpatRaster` stacks,
#'   each for a future scenario (e.g., `list(ssp126_2050 = rast1, ...)`).
#' @param output_dir Character. Output directory for all rasters and CSV.
#' @param method Character. Ensemble method. Default `"weighted"`.
#' @param models Character vector or `NULL`. Models to use.
#' @param mask A `terra::SpatRaster` or `NULL`. Prediction mask.
#' @param overwrite Logical. Overwrite existing outputs. Default `FALSE`.
#' @param compression Character. GeoTIFF compression. Default `"LZW"`.
#' @param clamp Reserved for extrapolation control: forwarded to
#'   [cast_ensemble_raster()] when the installed version supports a `clamp`
#'   argument; otherwise ignored with a warning. Default `NULL`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{current}{List with `hss_path`, `binary_path`, etc.}
#'   \item{future}{Named list of per-scenario results.}
#'   \item{stats}{A `data.frame` with per-scenario statistics.}
#'   \item{output_dir}{The output directory path.}
#' }
#'
#' @details
#' For each scenario, a change-class raster is computed:
#' - `1` = gain (absent now, present in future)
#' - `-1` = loss (present now, absent in future)
#' - `2` = stable present
#' - `0` = stable absent
#'
#' A scenario that fails is skipped with a warning and recorded as an `NA`
#' row in the statistics table; remaining scenarios continue.
#'
#' Centroid shift is computed as the HSS-weighted great-circle distance.
#'
#' @seealso [cast_ensemble_raster()], [cast_project()]
#'
#' @export
cast_project_raster <- function(fit, cv,
                                current_raster,
                                future_rasters,
                                output_dir,
                                method = c("weighted", "best", "equal"),
                                models = NULL,
                                mask = NULL,
                                overwrite = FALSE,
                                compression = "LZW",
                                clamp = NULL,
                                verbose = TRUE) {
  check_suggested("terra", "for raster projection")
  method <- match.arg(method)

  if (!is.list(future_rasters) || length(future_rasters) == 0) {
    cli::cli_abort("{.arg future_rasters} must be a non-empty named list of SpatRasters.")
  }
  if (is.null(names(future_rasters)) || any(names(future_rasters) == "")) {
    cli::cli_abort("All elements of {.arg future_rasters} must be named.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  raster_dir <- file.path(output_dir, "rasters")
  table_dir <- file.path(output_dir, "tables")
  dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  # Reserved passthrough: forward `clamp` only when the installed
  # cast_ensemble_raster() actually supports it.
  er_formals <- names(formals(cast_ensemble_raster))
  if (!is.null(clamp) && !("clamp" %in% er_formals)) {
    cli::cli_warn(
      "{.arg clamp} ignored: this {.fn cast_ensemble_raster} has no clamp support."
    )
  }
  er_call <- function(r, prefix) {
    args <- list(
      fit = fit, cv = cv, raster_stack = r, output_dir = raster_dir,
      method = method, models = models, mask = mask, prefix = prefix,
      overwrite = overwrite, compression = compression, verbose = verbose
    )
    if (!is.null(clamp) && "clamp" %in% er_formals) args$clamp <- clamp
    do.call(cast_ensemble_raster, args)
  }

  # ---- Current prediction -----------------------------------------------------
  if (verbose) cli::cli_h2("Current prediction")

  current_result <- er_call(current_raster, "current")

  # Read current binary for change computation
  cur_bin <- terra::rast(current_result$binary_path)
  cur_hss <- terra::rast(current_result$hss_path)

  # Compute current centroid once
  cur_centroid <- .weighted_centroid_raster(cur_hss, cur_bin)

  # ---- Future predictions -----------------------------------------------------
  future_results <- list()
  stats_rows <- list()

  for (scen in names(future_rasters)) {
    if (verbose) cli::cli_h2("Future: {scen}")

    sc_res <- tryCatch({
      fut_raster <- future_rasters[[scen]]
      if (is.character(fut_raster)) fut_raster <- terra::rast(fut_raster)

      fut_result <- er_call(fut_raster, scen)

      # ---- Change map ---------------------------------------------------------
      change_path <- file.path(raster_dir, paste0(scen, "_change_class.tif"))

      if (!overwrite && file.exists(change_path)) {
        if (verbose) cli::cli_inform("Change raster exists; skipping.")
      } else {
        fut_bin <- terra::rast(fut_result$binary_path)
        fut_hss <- terra::rast(fut_result$hss_path)

        # Change class: gain=1, loss=-1, stable_present=2, stable_absent=0
        change_r <- terra::lapp(
          c(cur_bin, fut_bin),
          fun = function(cur, fut) {
            out <- rep(NA_integer_, length(cur))
            valid <- !is.na(cur) & !is.na(fut)
            out[valid & cur == 0 & fut == 1] <- 1L    # gain
            out[valid & cur == 1 & fut == 0] <- -1L   # loss
            out[valid & cur == 1 & fut == 1] <- 2L    # stable_present
            out[valid & cur == 0 & fut == 0] <- 0L    # stable_absent
            out
          }
        )
        names(change_r) <- "change_class"
        terra::writeRaster(change_r, change_path, overwrite = TRUE,
                           gdal = c(paste0("COMPRESS=", compression)),
                           wopt = list(datatype = "INT2S"))
      }

      # ---- Statistics (from disk, also on the skip branch) --------------------
      change_vals <- terra::values(terra::rast(change_path), mat = FALSE)
      change_vals <- change_vals[!is.na(change_vals)]

      n_gain   <- sum(change_vals == 1L)
      n_loss   <- sum(change_vals == -1L)
      n_stable <- sum(change_vals == 2L)
      n_absent <- sum(change_vals == 0L)
      total_present_now <- n_loss + n_stable

      pct_change <- if (total_present_now > 0) {
        100 * (n_gain - n_loss) / total_present_now
      } else {
        NA_real_
      }

      # Future centroid
      fut_centroid <- .weighted_centroid_raster(
        terra::rast(fut_result$hss_path),
        terra::rast(fut_result$binary_path))

      # Centroid shift
      shift_km <- tryCatch(
        .haversine_km(
          cur_centroid$lat, cur_centroid$lon,
          fut_centroid$lat, fut_centroid$lon
        ),
        error = function(e) NA_real_
      )

      if (verbose) {
        cli::cli_inform(c(
          " " = "Gain: {n_gain} | Loss: {n_loss} | Stable: {n_stable}",
          " " = "Change: {round(pct_change, 1)}% | Shift: {round(shift_km, 1)} km"
        ))
      }

      list(
        result = fut_result,
        stats = data.frame(
          scenario          = scen,
          n_gain            = n_gain,
          n_loss            = n_loss,
          n_stable_present  = n_stable,
          n_stable_absent   = n_absent,
          pct_change        = round(pct_change, 2),
          current_centroid_lon = round(cur_centroid$lon, 4),
          current_centroid_lat = round(cur_centroid$lat, 4),
          future_centroid_lon  = round(fut_centroid$lon, 4),
          future_centroid_lat  = round(fut_centroid$lat, 4),
          centroid_shift_km = round(shift_km, 1),
          stringsAsFactors  = FALSE
        )
      )
    }, error = function(e) {
      cli::cli_warn("Scenario {.val {scen}} failed: {conditionMessage(e)}")
      NULL
    })

    if (is.null(sc_res)) {
      stats_rows[[scen]] <- data.frame(
        scenario = scen, n_gain = NA_integer_, n_loss = NA_integer_,
        n_stable_present = NA_integer_, n_stable_absent = NA_integer_,
        pct_change = NA_real_, current_centroid_lon = NA_real_,
        current_centroid_lat = NA_real_, future_centroid_lon = NA_real_,
        future_centroid_lat = NA_real_, centroid_shift_km = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      future_results[[scen]] <- sc_res$result
      stats_rows[[scen]] <- sc_res$stats
    }
  }

  # ---- Save projection statistics ---------------------------------------------
  if (length(stats_rows) > 0) {
    stats_df <- do.call(rbind, stats_rows)
    rownames(stats_df) <- NULL
    utils::write.csv(
      stats_df,
      file.path(table_dir, "projection_stats.csv"),
      row.names = FALSE
    )
    if (verbose) {
      cli::cli_inform(c("v" = "Projection stats saved to {.path {file.path(table_dir, 'projection_stats.csv')}}"))
    }
  } else {
    stats_df <- data.frame()
  }

  invisible(list(
    current    = current_result,
    future     = future_results,
    stats      = stats_df,
    output_dir = output_dir
  ))
}


#' Compute HSS-weighted centroid from rasters
#' @keywords internal
#' @noRd
.weighted_centroid_raster <- function(hss_raster, binary_raster) {
  # Extract values
  hss_vals <- terra::values(hss_raster, mat = FALSE)
  bin_vals <- terra::values(binary_raster, mat = FALSE)

  # Cells where species is present
  present <- which(!is.na(bin_vals) & bin_vals == 1)

  if (length(present) == 0) {
    return(list(lon = NA_real_, lat = NA_real_))
  }

  xy <- terra::xyFromCell(hss_raster, present)
  w <- hss_vals[present]
  w[is.na(w)] <- 0

  if (sum(w) == 0) w <- rep(1, length(w))

  list(
    lon = stats::weighted.mean(xy[, 1], w),
    lat = stats::weighted.mean(xy[, 2], w)
  )
}
