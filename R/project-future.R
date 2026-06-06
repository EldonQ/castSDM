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
#' @param threshold_method Character. Binary threshold method. Default
#'   `"maxTSS"`.
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
#' Centroid shift is computed as the great-circle distance between the
#' weighted centroid of current presence and the weighted centroid of
#' future presence.
#'
#' @seealso [cast_ensemble()], [cast_fit()], [cast_cv()]
#'
#' @export
cast_project <- function(fit, cv, current_env, future_envs,
                         method = c("weighted", "best", "equal"),
                         threshold_method = "maxTSS",
                         models = NULL,
                         save_dir = NULL) {
  method <- match.arg(method)

  if (!is.list(future_envs) || length(future_envs) == 0) {
    cli::cli_abort("{.arg future_envs} must be a non-empty named list of data.frames.")
  }
  if (is.null(names(future_envs)) || any(names(future_envs) == "")) {
    cli::cli_abort("All elements of {.arg future_envs} must be named (scenario names).")
  }

  # ---- Current prediction -------------------------------------------------
  current <- cast_ensemble(fit, cv, current_env,
                           method = method,
                           threshold_method = threshold_method,
                           models = models)

  # ---- Future predictions -------------------------------------------------
  future_list <- list()
  changes_list <- list()
  stats_rows <- list()

  for (scen in names(future_envs)) {
    fut_env <- future_envs[[scen]]

    # Predict using the same ensemble configuration
    fut <- cast_ensemble(fit, cv, fut_env,
                         method = method,
                         threshold_method = threshold_method,
                         models = models)
    future_list[[scen]] <- fut

    # ---- Compute change map -----------------------------------------------
    cur_bin <- current$predictions$binary_ensemble
    fut_bin <- fut$predictions$binary_ensemble

    change <- character(length(cur_bin))
    change[cur_bin == 0 & fut_bin == 1] <- "gain"
    change[cur_bin == 1 & fut_bin == 0] <- "loss"
    change[cur_bin == 1 & fut_bin == 1] <- "stable_present"
    change[cur_bin == 0 & fut_bin == 0] <- "stable_absent"

    has_coords <- all(c("lon", "lat") %in% names(current$predictions))
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
    changes_list[[scen]] <- change_df

    # ---- Summary statistics -----------------------------------------------
    n_gain   <- sum(change == "gain")
    n_loss   <- sum(change == "loss")
    n_stable <- sum(change == "stable_present")
    n_absent <- sum(change == "stable_absent")
    total_present_now <- sum(cur_bin == 1)
    pct_change <- if (total_present_now > 0) {
      100 * (n_gain - n_loss) / total_present_now
    } else {
      NA_real_
    }

    # Centroid shift (great-circle distance in km)
    centroid_km <- NA_real_
    if (has_coords) {
      centroid_km <- tryCatch({
        cur_pres <- current$predictions[cur_bin == 1, ]
        fut_pres <- fut$predictions[fut_bin == 1, ]
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

    stats_rows[[scen]] <- data.frame(
      scenario        = scen,
      n_gain          = n_gain,
      n_loss          = n_loss,
      n_stable_present = n_stable,
      n_stable_absent  = n_absent,
      pct_change      = round(pct_change, 2),
      centroid_shift_km = round(centroid_km, 1),
      stringsAsFactors = FALSE
    )
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
#' @keywords internal
#' @noRd
.save_prediction_tif <- function(df, lon_col, lat_col, val_col, path) {
  tryCatch({
    pts <- terra::vect(
      cbind(df[[lon_col]], df[[lat_col]]),
      type = "points",
      atts = data.frame(value = df[[val_col]]),
      crs = "EPSG:4326"
    )
    r <- terra::rasterize(pts, terra::rast(pts, res = 0.05), field = "value",
                          fun = "mean")
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
#' @param threshold_method Character. Threshold method. Default `"maxTSS"`.
#' @param models Character vector or `NULL`. Models to use.
#' @param mask A `terra::SpatRaster` or `NULL`. Prediction mask.
#' @param overwrite Logical. Overwrite existing outputs. Default `FALSE`.
#' @param compression Character. GeoTIFF compression. Default `"LZW"`.
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
                                threshold_method = "maxTSS",
                                models = NULL,
                                mask = NULL,
                                overwrite = FALSE,
                                compression = "LZW",
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

  # ---- Current prediction -----------------------------------------------------
  if (verbose) cli::cli_h2("Current prediction")

  current_result <- cast_ensemble_raster(
    fit, cv, current_raster,
    output_dir = raster_dir,
    method = method,
    threshold_method = threshold_method,
    models = models,
    mask = mask,
    prefix = "current",
    overwrite = overwrite,
    compression = compression,
    verbose = verbose
  )

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

    fut_raster <- future_rasters[[scen]]
    if (is.character(fut_raster)) fut_raster <- terra::rast(fut_raster)

    fut_result <- cast_ensemble_raster(
      fit, cv, fut_raster,
      output_dir = raster_dir,
      method = method,
      threshold_method = threshold_method,
      models = models,
      mask = mask,
      prefix = scen,
      overwrite = overwrite,
      compression = compression,
      verbose = verbose
    )
    future_results[[scen]] <- fut_result

    # ---- Change map -----------------------------------------------------------
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

      # ---- Statistics ---------------------------------------------------------
      change_vals <- terra::values(change_r, mat = FALSE)
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
      fut_centroid <- .weighted_centroid_raster(fut_hss, fut_bin)

      # Centroid shift
      shift_km <- tryCatch(
        .haversine_km(
          cur_centroid$lat, cur_centroid$lon,
          fut_centroid$lat, fut_centroid$lon
        ),
        error = function(e) NA_real_
      )

      stats_rows[[scen]] <- data.frame(
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

      if (verbose) {
        cli::cli_inform(c(
          " " = "Gain: {n_gain} | Loss: {n_loss} | Stable: {n_stable}",
          " " = "Change: {round(pct_change, 1)}% | Shift: {round(shift_km, 1)} km"
        ))
      }
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
      cli::cli_inform("v" = "Projection stats saved to {.path {file.path(table_dir, 'projection_stats.csv')}}")
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
