#' Thin Presence Records by a Minimum Geographic Distance
#'
#' Greedy great-circle thinning: records are visited in a seeded random
#' order and kept only when they lie at least `min_dist_km` (Haversine)
#' away from every record already kept. This reduces spatial duplication
#' and sampling-cluster bias before modeling. The visiting order decides
#' which record of a conflicting cluster survives; the seed makes that
#' choice reproducible. Rows are returned in their original order.
#'
#' @param occurrences A `data.frame` with longitude/latitude columns
#'   (decimal degrees). Additional columns are preserved.
#' @param min_dist_km Numeric > 0. Minimum great-circle distance between
#'   retained records, in kilometres. Default `1`.
#' @param lon_col,lat_col Character. Coordinate column names. Defaults
#'   `"lon"` and `"lat"`.
#' @param seed Integer or `NULL`. Random seed controlling the visiting
#'   order (reproducibility of which records are kept).
#' @param verbose Logical. Report the kept/removed counts. Default `TRUE`.
#'
#' @return The thinned `data.frame` with rows of `occurrences` (original
#'   column values, original row order among kept rows).
#'
#' @details
#' Thinning is a standard pre-processing step for clumped occurrence data:
#' dense clusters over-weight locally abundant areas and inflate apparent
#' model performance under spatial cross-validation. The greedy pass is a
#' single sweep (not an iterative maximum-cardinality search), so the
#' retained set depends on the visiting order; supply a `seed` for
#' reproducibility. Non-finite (missing) coordinates are dropped before
#' thinning and reported when `verbose = TRUE`.
#'
#' @references
#' Aiello-Lammens, M. E. et al. (2015). spThin: an R package for spatial
#' thinning of species occurrence records for use in ecological niche
#' models. *Ecography*, 38(5), 541-545.
#'
#' @examples
#' \dontrun{
#' occ_thin <- cast_thin(occurrences, min_dist_km = 5, seed = 42)
#' }
#' @export
cast_thin <- function(occurrences,
                      min_dist_km = 1,
                      lon_col = "lon",
                      lat_col = "lat",
                      seed = NULL,
                      verbose = TRUE) {
  if (!is.data.frame(occurrences)) {
    cli::cli_abort("{.arg occurrences} must be a data.frame.")
  }
  missing_cols <- setdiff(c(lon_col, lat_col), names(occurrences))
  if (length(missing_cols)) {
    cli::cli_abort(
      "{.arg occurrences} is missing coordinate column{?s}: {.val {missing_cols}}."
    )
  }
  if (!is.numeric(min_dist_km) || length(min_dist_km) != 1L ||
      !is.finite(min_dist_km) || min_dist_km <= 0) {
    cli::cli_abort("{.arg min_dist_km} must be a single positive number (km).")
  }
  lon <- as.numeric(occurrences[[lon_col]])
  lat <- as.numeric(occurrences[[lat_col]])
  # Haversine is defined on decimal degrees; silently treating projected
  # metres as degrees would thin by a factor of ~10^5.
  if (any(lon < -180 | lon > 180, na.rm = TRUE) ||
      any(lat < -90 | lat > 90, na.rm = TRUE)) {
    cli::cli_abort(c(
      "Coordinates must be decimal degrees (lon in [-180, 180], lat in [-90, 90]).",
      i = "Projected coordinates are not supported; reproject to lon/lat first."
    ))
  }

  ok <- is.finite(lon) & is.finite(lat)
  n_na <- sum(!ok)
  if (n_na == nrow(occurrences)) {
    cli::cli_abort("{.arg occurrences} has no finite coordinates.")
  }
  if (n_na > 0 && verbose) {
    cli::cli_warn("{n_na} record{?s} with missing coordinates dropped before thinning.")
  }

  if (!is.null(seed)) set.seed(seed)
  ord <- sample(which(ok))  # random visiting order

  # Haversine great-circle distance, vectorised over already-kept records.
  R <- 6371.0088  # IUGG mean Earth radius, km
  lon_rad <- lon[ord] * pi / 180
  lat_rad <- lat[ord] * pi / 180
  kept <- integer(length(ord))
  kept_n <- 0L
  for (i in seq_along(ord)) {
    if (kept_n == 0L) {
      kept_n <- 1L
      kept[[kept_n]] <- i
      next
    }
    j <- kept[seq_len(kept_n)]
    dlat <- lat_rad[i] - lat_rad[j]
    dlon <- lon_rad[i] - lon_rad[j]
    a <- sin(dlat / 2)^2 +
      cos(lat_rad[j]) * cos(lat_rad[i]) * sin(dlon / 2)^2
    d <- 2 * R * asin(pmin(1, sqrt(a)))
    if (all(d >= min_dist_km)) {
      kept_n <- kept_n + 1L
      kept[[kept_n]] <- i
    }
  }

  keep_rows <- sort(ord[kept[seq_len(kept_n)]])
  out <- occurrences[keep_rows, , drop = FALSE]
  rownames(out) <- NULL
  if (verbose) {
    cli::cli_inform(
      "Thinning: kept {kept_n}/{length(ord)} record{?s} ({length(ord) - kept_n} removed at >= {min_dist_km} km)."
    )
  }
  out
}
