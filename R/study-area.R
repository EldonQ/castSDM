#' Define Species-Specific Study Area (Accessible Area M)
#'
#' Creates a study area mask following the BAM framework (Barve et al. 2011),
#' which recommends using the species' accessible area (M) rather than
#' administrative boundaries. The study area constrains both background
#' point sampling and spatial prediction extent.
#'
#' @param occurrences A `data.frame` with at least `lon` and `lat` columns,
#'   or an `sf` points object.
#' @param raster_template A `terra::SpatRaster` used as the reference grid
#'   (e.g., one layer of the environmental stack). Only cells that are
#'   non-NA in this template can be part of the study area.
#' @param method Character. How to define the study area:
#'   - `"buffer"` (default): convex hull around occurrences + buffer.
#'   - `"bbox_buffer"`: bounding box of occurrences + buffer (simpler shape).
#'   - `"full"`: entire non-NA extent of `raster_template` (no restriction).
#' @param buffer_km Numeric. Buffer distance in kilometers around occurrences.
#'   Only used when `method` is `"buffer"` or `"bbox_buffer"`. Default `300`.
#' @param min_cells Integer. Minimum number of valid cells required in the
#'   study area. If the result has fewer cells, the buffer is automatically
#'   expanded. Default `1000L`.
#' @param verbose Logical. Print informational messages. Default `TRUE`.
#'
#' @return A `cast_study_area` object (S3 class) with components:
#' \describe{
#'   \item{mask}{A `terra::SpatRaster` (0/1 or NA) at the same resolution and
#'     CRS as `raster_template`.}
#'   \item{boundary}{An `sf` polygon of the study area boundary.}
#'   \item{n_cells}{Integer. Number of valid (non-NA) cells in the mask.}
#'   \item{method}{Character. The method used.}
#'   \item{buffer_km}{Numeric. The buffer distance used (NA for `"full"`).}
#'   \item{n_occurrences}{Integer. Number of occurrence points provided.}
#' }
#'
#' @references
#' Barve, N., Barve, V., Jimenez-Valverde, A. et al. (2011). The crucial
#' role of the accessible area in ecological niche modeling and species
#' distribution modeling. *Ecological Modelling*, 222(11), 1810-1819.
#'
#' @seealso [cast_background()], [cast_ensemble_raster()]
#'
#' @export
cast_study_area <- function(occurrences,
                            raster_template,
                            method = c("buffer", "bbox_buffer", "full"),
                            buffer_km = 300,
                            min_cells = 1000L,
                            verbose = TRUE) {
  check_suggested("terra", "for raster operations")
  check_suggested("sf", "for spatial operations")
  method <- match.arg(method)

  # ---- Validate inputs -------------------------------------------------------
  if (!inherits(raster_template, "SpatRaster")) {
    cli::cli_abort("{.arg raster_template} must be a {.cls SpatRaster}.")
  }

  # Convert occurrences to sf points
  if (inherits(occurrences, "sf")) {
    occ_sf <- occurrences
  } else if (is.data.frame(occurrences)) {
    if (!all(c("lon", "lat") %in% names(occurrences))) {
      cli::cli_abort("{.arg occurrences} must have {.val lon} and {.val lat} columns.")
    }
    occ_sf <- sf::st_as_sf(
      occurrences,
      coords = c("lon", "lat"),
      crs = 4326
    )
  } else {
    cli::cli_abort("{.arg occurrences} must be a data.frame or sf object.")
  }

  n_occ <- nrow(occ_sf)
  if (n_occ < 3) {
    cli::cli_abort("Need at least 3 occurrence points (got {n_occ}).")
  }

  # Ensure CRS match
  raster_crs <- terra::crs(raster_template, describe = TRUE)$code
  occ_crs <- sf::st_crs(occ_sf)$epsg
  if (!is.na(raster_crs) && !is.na(occ_crs) && raster_crs != occ_crs) {
    occ_sf <- sf::st_transform(occ_sf, terra::crs(raster_template))
  }

  # ---- Build study area polygon -----------------------------------------------
  if (method == "full") {
    # Use entire non-NA extent
    mask_r <- !is.na(raster_template[[1]])
    mask_r[mask_r == 0] <- NA
    boundary_sf <- .raster_extent_to_sf(raster_template)
    used_buffer <- NA_real_

  } else if (method == "buffer") {
    # Convex hull + buffer
    hull <- sf::st_convex_hull(sf::st_union(occ_sf))

    # Buffer in meters (convert km → m)
    buffer_m <- buffer_km * 1000
    if (sf::st_is_longlat(occ_sf)) {
      # For geographic CRS, project to equal-area for accurate buffering
      ea_crs <- .choose_equal_area_crs(occ_sf)
      hull_proj <- sf::st_transform(hull, ea_crs)
      buffered_proj <- sf::st_buffer(hull_proj, dist = buffer_m)
      boundary_sf <- sf::st_transform(buffered_proj, sf::st_crs(occ_sf))
    } else {
      boundary_sf <- sf::st_buffer(hull, dist = buffer_m)
    }

    mask_r <- .mask_raster_by_polygon(raster_template, boundary_sf)
    used_buffer <- buffer_km

  } else if (method == "bbox_buffer") {
    # Bounding box + buffer
    bbox <- sf::st_as_sfc(sf::st_bbox(occ_sf))
    buffer_m <- buffer_km * 1000
    if (sf::st_is_longlat(occ_sf)) {
      ea_crs <- .choose_equal_area_crs(occ_sf)
      bbox_proj <- sf::st_transform(bbox, ea_crs)
      buffered_proj <- sf::st_buffer(bbox_proj, dist = buffer_m)
      boundary_sf <- sf::st_transform(buffered_proj, sf::st_crs(occ_sf))
    } else {
      boundary_sf <- sf::st_buffer(bbox, dist = buffer_m)
    }

    mask_r <- .mask_raster_by_polygon(raster_template, boundary_sf)
    used_buffer <- buffer_km
  }

  # ---- Check minimum cells and expand if needed --------------------------------
  n_cells <- sum(terra::values(mask_r, mat = FALSE) == 1, na.rm = TRUE)

  if (n_cells < min_cells && method != "full") {
    expanded_buffer <- buffer_km
    while (n_cells < min_cells && expanded_buffer < buffer_km * 3) {
      expanded_buffer <- expanded_buffer * 1.5
      if (verbose) {
        cli::cli_inform(
          "Study area has only {n_cells} cells (< {min_cells}); expanding buffer to {round(expanded_buffer)} km."
        )
      }
      buffer_m <- expanded_buffer * 1000
      if (method == "buffer") {
        hull <- sf::st_convex_hull(sf::st_union(occ_sf))
        if (sf::st_is_longlat(occ_sf)) {
          ea_crs <- .choose_equal_area_crs(occ_sf)
          hull_proj <- sf::st_transform(hull, ea_crs)
          buffered_proj <- sf::st_buffer(hull_proj, dist = buffer_m)
          boundary_sf <- sf::st_transform(buffered_proj, sf::st_crs(occ_sf))
        } else {
          boundary_sf <- sf::st_buffer(hull, dist = buffer_m)
        }
      } else {
        bbox <- sf::st_as_sfc(sf::st_bbox(occ_sf))
        if (sf::st_is_longlat(occ_sf)) {
          ea_crs <- .choose_equal_area_crs(occ_sf)
          bbox_proj <- sf::st_transform(bbox, ea_crs)
          buffered_proj <- sf::st_buffer(bbox_proj, dist = buffer_m)
          boundary_sf <- sf::st_transform(buffered_proj, sf::st_crs(occ_sf))
        } else {
          boundary_sf <- sf::st_buffer(bbox, dist = buffer_m)
        }
      }
      mask_r <- .mask_raster_by_polygon(raster_template, boundary_sf)
      n_cells <- sum(terra::values(mask_r, mat = FALSE) == 1, na.rm = TRUE)
    }
    used_buffer <- expanded_buffer
  }

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Study area defined: {.val {method}} method",
      " " = "Buffer: {ifelse(is.na(used_buffer), 'N/A', paste0(round(used_buffer), ' km'))}",
      " " = "Valid cells: {format(n_cells, big.mark = ',')}",
      " " = "Occurrences: {n_occ}"
    ))
  }

  # ---- Return ------------------------------------------------------------------
  structure(
    list(
      mask          = mask_r,
      boundary      = boundary_sf,
      n_cells       = n_cells,
      method        = method,
      buffer_km     = used_buffer,
      n_occurrences = n_occ
    ),
    class = "cast_study_area"
  )
}


#' @export
print.cast_study_area <- function(x, ...) {
  cat("<cast_study_area>\n")
  cat("  method       :", x$method, "\n")
  if (!is.na(x$buffer_km)) {
    cat("  buffer_km    :", round(x$buffer_km), "\n")
  }
  cat("  valid cells  :", format(x$n_cells, big.mark = ","), "\n")
  cat("  occurrences  :", x$n_occurrences, "\n")
  invisible(x)
}


# ---- Internal helpers ---------------------------------------------------------

#' Mask a raster by an sf polygon
#' @keywords internal
#' @noRd
.mask_raster_by_polygon <- function(raster_template, polygon_sf) {
  # Convert sf to terra SpatVector
  poly_v <- terra::vect(polygon_sf)

  # Ensure CRS match
  if (!identical(terra::crs(raster_template), terra::crs(poly_v))) {
    poly_v <- terra::project(poly_v, terra::crs(raster_template))
  }

  # Create mask: 1 where inside polygon AND raster is non-NA, NA elsewhere
  ref_layer <- raster_template[[1]]
  mask_r <- terra::mask(ref_layer, poly_v)
  mask_r <- terra::classify(mask_r, cbind(-Inf, Inf, 1))  # all non-NA → 1
  mask_r
}


#' Convert raster extent to an sf polygon
#' @keywords internal
#' @noRd
.raster_extent_to_sf <- function(r) {
  ext <- terra::ext(r)
  coords <- matrix(c(
    ext[1], ext[3],
    ext[2], ext[3],
    ext[2], ext[4],
    ext[1], ext[4],
    ext[1], ext[3]
  ), ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  sf::st_sfc(poly, crs = terra::crs(r))
}


#' Choose an appropriate equal-area CRS for buffering
#' @keywords internal
#' @noRd
.choose_equal_area_crs <- function(sf_obj) {
  # Use Lambert Azimuthal Equal Area centered on data centroid
  centroid <- sf::st_coordinates(sf::st_centroid(sf::st_union(sf_obj)))
  sprintf(
    "+proj=laea +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +datum=WGS84 +units=m",
    centroid[1, 2], centroid[1, 1]
  )
}
