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
                         models = NULL) {
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
