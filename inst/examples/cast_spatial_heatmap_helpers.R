# cast_spatial_heatmap_helpers.R
# ------------------------------------------------------------------------------
# Post-process cast_predict / cast_cate as filled heatmaps (replot only;
# no model refit). The defaults follow fig7_spatial_prediction_maps.py:
# - China extent: [73.5, 135, 18, 53.5]
# - Regular lon/lat grid at 0.06 degree
# - Nearest-neighbor interpolation on grid nodes
# - Optional display downsampling
# - China boundary shap mask (when basemap = "china")
#
# Requires: ggplot2, sf, scales, terra
# For nearest interpolation: FNN (preferred) or gstat fallback.
# ------------------------------------------------------------------------------

.cast_spatial_extent_default <- function() {
  c(lon_min = 73.5, lon_max = 135, lat_min = 18, lat_max = 53.5)
}

.cast_spatial_load_basemap <- function(type) {
  if (!requireNamespace("castSDM", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch(
    utils::getFromNamespace("load_basemap", "castSDM")(type),
    error = function(e) NULL
  )
}

.cast_spatial_hss_palette <- function() {
  c("#0d0887", "#3b049a", "#7201a8", "#a52c60", "#d44842", "#ed7953", "#f0f921")
}

.cast_spatial_diff_palette <- function() {
  c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b")
}

.cast_spatial_build_grid <- function(extent, res_deg) {
  extv <- if (is.null(extent)) {
    .cast_spatial_extent_default()
  } else {
    stats::setNames(as.numeric(extent[seq_len(4)]), c("lon_min", "lon_max", "lat_min", "lat_max"))
  }
  if (!is.finite(res_deg) || res_deg <= 0) {
    stop("res_deg must be a positive number.", call. = FALSE)
  }

  lon_seq <- seq(extv[["lon_min"]], extv[["lon_max"]] + res_deg * 0.5, by = res_deg)
  lat_seq <- seq(extv[["lat_min"]], extv[["lat_max"]] + res_deg * 0.5, by = res_deg)
  list(extent = extv, lon = lon_seq, lat = lat_seq)
}

.cast_spatial_interp_nearest <- function(obs_xy, obs_z, query_xy) {
  if (nrow(obs_xy) < 8L) {
    return(rep(NA_real_, nrow(query_xy)))
  }

  if (requireNamespace("FNN", quietly = TRUE)) {
    idx <- integer(nrow(query_xy))
    chunk <- 120000L
    starts <- seq.int(1L, nrow(query_xy), by = chunk)
    for (st in starts) {
      ed <- min(st + chunk - 1L, nrow(query_xy))
      kn <- FNN::get.knnx(data = obs_xy, query = query_xy[st:ed, , drop = FALSE], k = 1L)
      idx[st:ed] <- kn$nn.index[, 1L]
    }
    return(obs_z[idx])
  }

  if (requireNamespace("gstat", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
    obs <- sf::st_as_sf(
      data.frame(lon = obs_xy[, 1L], lat = obs_xy[, 2L], vv = obs_z),
      coords = c("lon", "lat"), crs = "EPSG:4326"
    )
    grd <- sf::st_as_sf(
      data.frame(lon = query_xy[, 1L], lat = query_xy[, 2L]),
      coords = c("lon", "lat"), crs = "EPSG:4326"
    )
    out <- gstat::idw(vv ~ 1, obs, grd, idp = 0, nmax = 1)
    nm <- names(out)
    i <- grep("\\.pred$", nm)[1L]
    if (is.na(i)) i <- which(nm == "var1.pred")[1L]
    if (is.na(i)) i <- length(nm)
    return(as.numeric(out[[nm[i]]]))
  }

  stop(
    "Nearest interpolation needs package 'FNN' (preferred) or gstat+sf fallback.\n",
    "Install: install.packages('FNN')",
    call. = FALSE
  )
}

.cast_spatial_interp_idw <- function(obs_xy, obs_z, query_xy, nmax_idw = 24L) {
  if (!requireNamespace("gstat", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE)) {
    stop("IDW interpolation requires packages 'gstat' and 'sf'.", call. = FALSE)
  }
  obs <- sf::st_as_sf(
    data.frame(lon = obs_xy[, 1L], lat = obs_xy[, 2L], vv = obs_z),
    coords = c("lon", "lat"), crs = "EPSG:4326"
  )
  grd <- sf::st_as_sf(
    data.frame(lon = query_xy[, 1L], lat = query_xy[, 2L]),
    coords = c("lon", "lat"), crs = "EPSG:4326"
  )
  nm <- min(max(6L, as.integer(nmax_idw)), nrow(obs_xy))
  out <- gstat::idw(vv ~ 1, obs, grd, idp = 2, nmax = nm)
  out_nm <- names(out)
  i <- grep("\\.pred$", out_nm)[1L]
  if (is.na(i)) i <- which(out_nm == "var1.pred")[1L]
  if (is.na(i)) i <- length(out_nm)
  as.numeric(out[[out_nm[i]]])
}

.cast_spatial_to_raster <- function(grid_lon, grid_lat, grid_z, extent, res_deg) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for heatmap post-processing.", call. = FALSE)
  }

  xmn <- extent[["lon_min"]] - res_deg / 2
  xmx <- extent[["lon_max"]] + res_deg / 2
  ymn <- extent[["lat_min"]] - res_deg / 2
  ymx <- extent[["lat_max"]] + res_deg / 2

  r <- terra::rast(
    xmin = xmn, xmax = xmx,
    ymin = ymn, ymax = ymx,
    resolution = res_deg,
    crs = "EPSG:4326"
  )

  v <- terra::vect(data.frame(lon = grid_lon, lat = grid_lat), geom = c("lon", "lat"), crs = "EPSG:4326")
  v$z <- as.numeric(grid_z)
  terra::rasterize(v, r, field = "z", fun = "mean", background = NA_real_)
}

.cast_spatial_apply_china_mask <- function(r, basemap = "china") {
  if (!identical(basemap, "china")) {
    return(r)
  }
  if (!requireNamespace("terra", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE)) {
    return(r)
  }

  china <- .cast_spatial_load_basemap("china")
  if (is.null(china) || !inherits(china, "sf")) {
    return(r)
  }

  poly <- tryCatch(sf::st_union(sf::st_make_valid(china)), error = function(e) NULL)
  if (is.null(poly)) {
    return(r)
  }

  v <- tryCatch(terra::vect(sf::st_as_sf(poly)), error = function(e) NULL)
  if (is.null(v)) {
    return(r)
  }
  terra::mask(r, v)
}

.cast_spatial_raster_to_df <- function(r, display_res_deg = NULL, base_res_deg = NULL) {
  if (!is.null(display_res_deg) && is.finite(display_res_deg) && display_res_deg > 0 &&
      is.finite(base_res_deg) && display_res_deg > base_res_deg) {
    fact <- max(1L, as.integer(round(display_res_deg / base_res_deg)))
    if (fact > 1L) {
      r <- terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
    }
  }

  nm <- names(r)[1L]
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(d) <- c("lon", "lat", "z")
  d
}

.cast_spatial_interpolate_grid <- function(lon,
                                           lat,
                                           z,
                                           extent = NULL,
                                           res_deg = 0.06,
                                           interp_method = "nearest",
                                           display_res_deg = NULL,
                                           nmax_idw = 24L,
                                           basemap = "china") {
  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  z <- as.numeric(z)
  ok <- is.finite(lon) & is.finite(lat) & is.finite(z)
  lon <- lon[ok]
  lat <- lat[ok]
  z <- z[ok]
  if (length(z) < 8L) {
    return(NULL)
  }

  g <- .cast_spatial_build_grid(extent = extent, res_deg = res_deg)
  query <- expand.grid(lon = g$lon, lat = g$lat)
  obs_xy <- cbind(lon, lat)
  query_xy <- as.matrix(query[, c("lon", "lat")])

  interp_method <- tolower(trimws(as.character(interp_method)[1L]))
  pred_z <- switch(
    interp_method,
    nearest = .cast_spatial_interp_nearest(obs_xy, z, query_xy),
    idw = .cast_spatial_interp_idw(obs_xy, z, query_xy, nmax_idw = nmax_idw),
    stop("interp_method must be one of: 'nearest', 'idw'.", call. = FALSE)
  )

  r <- .cast_spatial_to_raster(
    grid_lon = query$lon,
    grid_lat = query$lat,
    grid_z = pred_z,
    extent = g$extent,
    res_deg = res_deg
  )
  r <- .cast_spatial_apply_china_mask(r, basemap = basemap)
  .cast_spatial_raster_to_df(r, display_res_deg = display_res_deg, base_res_deg = res_deg)
}

.cast_spatial_mask_cate_by_hss <- function(df, pred, hss_model, hss_threshold) {
  if (is.null(pred) || !inherits(pred, "cast_predict")) {
    return(df)
  }
  pr <- pred$predictions
  hss_col <- paste0("HSS_", hss_model)
  if (!hss_col %in% names(pr)) {
    return(df)
  }
  key_df <- paste0(signif(df$lon, 8), "_", signif(df$lat, 8))
  key_pr <- paste0(signif(pr$lon, 8), "_", signif(pr$lat, 8))
  mi <- match(key_df, key_pr)
  hss_vals <- pr[[hss_col]][mi]
  low <- is.na(hss_vals) | (as.numeric(hss_vals) < hss_threshold)
  df$cate[low] <- NA_real_
  df
}

.cast_spatial_gg_hss <- function(df_raster,
                                 title,
                                 basemap = "china",
                                 extent = NULL) {
  if (is.null(df_raster) || nrow(df_raster) < 2L) {
    return(NULL)
  }
  if (is.null(extent)) {
    extent <- .cast_spatial_extent_default()
  }

  bm <- if (basemap != "none") .cast_spatial_load_basemap(basemap) else NULL
  dash <- if (identical(basemap, "china")) .cast_spatial_load_basemap("dashline") else NULL

  p <- ggplot2::ggplot()
  if (!is.null(bm)) {
    p <- p + ggplot2::geom_sf(data = bm, fill = "white", color = "black", linewidth = 0.25)
  }

  p <- p +
    ggplot2::geom_raster(
      ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$z),
      data = df_raster,
      interpolate = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colours = .cast_spatial_hss_palette(),
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "transparent",
      name = "HSS"
    ) +
    ggplot2::coord_sf(
      crs = sf::st_crs(4326),
      expand = FALSE,
      xlim = c(extent[["lon_min"]], extent[["lon_max"]]),
      ylim = c(extent[["lat_min"]], extent[["lat_max"]])
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12, color = "black"),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(2.5, "cm"),
      legend.key.height = ggplot2::unit(0.35, "cm"),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8)
    )

  if (!is.null(dash)) {
    p <- p + ggplot2::geom_sf(data = dash, fill = NA, color = "black", linewidth = 0.5)
  }

  p
}

.cast_spatial_gg_cate <- function(df_raster,
                                  title,
                                  subtitle = NULL,
                                  basemap = "china",
                                  extent = NULL,
                                  zlim = NULL) {
  if (is.null(df_raster) || nrow(df_raster) < 2L) {
    return(NULL)
  }
  if (is.null(extent)) {
    extent <- .cast_spatial_extent_default()
  }
  if (is.null(zlim)) {
    zz <- df_raster$z[is.finite(df_raster$z)]
    if (length(zz) < 5L) {
      return(NULL)
    }
    q98 <- as.numeric(stats::quantile(abs(zz), 0.98, na.rm = TRUE))
    s <- stats::sd(zz, na.rm = TRUE)
    zq <- suppressWarnings(max(q98, stats::mad(zz, na.rm = TRUE) * 2, s * 1.5, na.rm = TRUE))
    if (!is.finite(zq) || zq < 1e-20) {
      zq <- max(abs(zz), na.rm = TRUE)
    }
    if (!is.finite(zq) || zq < 1e-20) {
      zq <- 0.05
    }
    zlim <- c(-1, 1) * zq * 1.05
  }

  bm <- if (basemap != "none") .cast_spatial_load_basemap(basemap) else NULL
  dash <- if (identical(basemap, "china")) .cast_spatial_load_basemap("dashline") else NULL

  p <- ggplot2::ggplot()
  if (!is.null(bm)) {
    p <- p + ggplot2::geom_sf(data = bm, fill = "white", color = "black", linewidth = 0.25)
  }

  p <- p +
    ggplot2::geom_raster(
      ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$z),
      data = df_raster,
      interpolate = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colours = .cast_spatial_diff_palette(),
      limits = zlim,
      oob = scales::squish,
      na.value = "transparent",
      values = scales::rescale(c(zlim[1], zlim[1] * 0.5, zlim[1] * 0.2, zlim[1] * 0.05, 0, zlim[2] * 0.05, zlim[2] * 0.2, zlim[2] * 0.5, zlim[2])),
      name = "CATE",
      labels = scales::label_scientific(digits = 2)
    ) +
    ggplot2::coord_sf(
      crs = sf::st_crs(4326),
      expand = FALSE,
      xlim = c(extent[["lon_min"]], extent[["lon_max"]]),
      ylim = c(extent[["lat_min"]], extent[["lat_max"]])
    ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12, color = "black"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey30", size = 9),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(2.5, "cm"),
      legend.key.height = ggplot2::unit(0.35, "cm"),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8)
    )

  if (!is.null(dash)) {
    p <- p + ggplot2::geom_sf(data = dash, fill = NA, color = "black", linewidth = 0.5)
  }

  p
}

#' Re-save HSS maps (heatmap style) and CATE maps from existing objects.
#'
#' @param pred `cast_predict` object or NULL.
#' @param cate `cast_cate` object or NULL.
#' @param nmax_idw Local neighbourhood size for `gstat::idw` (only used when
#'   `interp_method = "idw"`).
#' @param species_label used in CATE subtitle.
#' @param ovis_style logical. If TRUE, filenames match run_ovis_ammon.R
#'   (`ovis_predict_<model>_china.png`, `ovis_cate_<var>.png`); else batch
#'   style (`HSS_<model>.png`, `CATE_<var>.png`).
cast_spatial_replot_hss_cate_heatmaps <- function(pred,
                                                  cate,
                                                  fig_dir,
                                                  fig_dpi,
                                                  var_labels = NULL,
                                                  basemap = "china",
                                                  res_deg = 0.06,
                                                  interp_method = "nearest",
                                                  display_res_deg = 0.02,
                                                  nmax_idw = 24L,
                                                  hss_model = "cast",
                                                  hss_threshold = 0.1,
                                                  species_label = NULL,
                                                  ovis_style = FALSE,
                                                  extent = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("terra", quietly = TRUE)) {
    stop("ggplot2, sf, and terra are required.", call. = FALSE)
  }

  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  ext <- if (is.null(extent)) {
    .cast_spatial_extent_default()
  } else {
    stats::setNames(as.numeric(extent[seq_len(4)]), c("lon_min", "lon_max", "lat_min", "lat_max"))
  }

  if (!is.null(pred) && inherits(pred, "cast_predict")) {
    pr <- pred$predictions
    for (md in pred$models) {
      col <- paste0("HSS_", md)
      if (!col %in% names(pr)) {
        next
      }
      df_r <- .cast_spatial_interpolate_grid(
        pr$lon,
        pr$lat,
        as.numeric(pr[[col]]),
        extent = ext,
        res_deg = res_deg,
        interp_method = interp_method,
        display_res_deg = display_res_deg,
        nmax_idw = nmax_idw,
        basemap = basemap
      )
      ttl <- sprintf("Habitat Suitability (%s)", toupper(md))
      ph <- .cast_spatial_gg_hss(df_r, ttl, basemap = basemap, extent = ext)
      if (!is.null(ph)) {
        fn <- if (isTRUE(ovis_style)) {
          sprintf("ovis_predict_%s_china.png", md)
        } else {
          sprintf("HSS_%s.png", md)
        }
        ggplot2::ggsave(
          file.path(fig_dir, fn),
          ph,
          width = 10,
          height = 7,
          dpi = fig_dpi,
          bg = "white",
          limitsize = FALSE
        )
        message("Heatmap saved: ", fn)
      }
    }
  }

  if (!is.null(cate) && inherits(cate, "cast_cate")) {
    sp_disp <- if (!is.null(species_label)) gsub("_", " ", species_label) else NULL

    for (v in cate$variables) {
      df0 <- cate$effects[cate$effects$variable == v, , drop = FALSE]
      if (nrow(df0) == 0L || !all(c("lon", "lat", "cate") %in% names(df0))) {
        next
      }
      df1 <- if (!is.null(pred)) {
        .cast_spatial_mask_cate_by_hss(df0, pred, hss_model, hss_threshold)
      } else {
        df0
      }

      df_r <- .cast_spatial_interpolate_grid(
        df1$lon,
        df1$lat,
        as.numeric(df1$cate),
        extent = ext,
        res_deg = res_deg,
        interp_method = interp_method,
        display_res_deg = display_res_deg,
        nmax_idw = nmax_idw,
        basemap = basemap
      )

      # Re-apply HSS mask on the interpolated raster grid so that
      # nearest-neighbor infill does not spread CATE into non-habitat areas
      if (!is.null(df_r) && !is.null(pred) && inherits(pred, "cast_predict")) {
        hss_col <- paste0("HSS_", hss_model)
        pr <- pred$predictions
        if (hss_col %in% names(pr)) {
          pr_sf <- sf::st_as_sf(pr, coords = c("lon", "lat"), crs = 4326)
          grid_sf <- sf::st_as_sf(df_r, coords = c("lon", "lat"), crs = 4326)
          nn_idx <- sf::st_nearest_feature(grid_sf, pr_sf)
          hss_at_grid <- as.numeric(pr[[hss_col]][nn_idx])
          df_r$z[is.na(hss_at_grid) | hss_at_grid < hss_threshold] <- NA_real_
          df_r <- df_r[!is.na(df_r$z), , drop = FALSE]
        }
      }

      vlab <- if (!is.null(var_labels) && v %in% names(var_labels)) var_labels[[v]] else v
      ttl <- sprintf("Spatial CATE: %s", vlab)
      nvis <- sum(!is.na(df1$cate))
      sub <- sprintf("Causal forest (%d trees) | n = %d cells displayed", cate$n_trees, nvis)
      if (!is.null(pred)) {
        sub <- paste0(sub, sprintf(" | HSS_%s >= %.2f", hss_model, hss_threshold))
      }
      if (!is.null(sp_disp)) {
        sub <- paste0(sp_disp, " | ", sub)
      }

      zc <- df1$cate[is.finite(df1$cate)]
      zlim_pts <- if (length(zc) >= 5L) {
        q98 <- as.numeric(stats::quantile(abs(zc), 0.98, na.rm = TRUE))
        s <- stats::sd(zc, na.rm = TRUE)
        zq <- suppressWarnings(max(q98, stats::mad(zc, na.rm = TRUE) * 2, s * 1.5, na.rm = TRUE))
        if (!is.finite(zq) || zq < 1e-20) zq <- max(abs(zc), na.rm = TRUE)
        if (is.finite(zq) && zq > 0) c(-1, 1) * zq * 1.05 else NULL
      } else {
        NULL
      }

      pc <- .cast_spatial_gg_cate(df_r, ttl, subtitle = sub, basemap = basemap, extent = ext, zlim = zlim_pts)
      if (!is.null(pc)) {
        fn <- if (isTRUE(ovis_style)) sprintf("ovis_cate_%s.png", v) else sprintf("CATE_%s.png", v)
        ggplot2::ggsave(
          file.path(fig_dir, fn),
          pc,
          width = 10,
          height = 7,
          dpi = fig_dpi,
          bg = "white",
          limitsize = FALSE
        )
        message("Heatmap saved: ", fn)
      }
    }
  }

  invisible(TRUE)
}
