#' Generate Adaptive Background (Pseudo-Absence) Points
#'
#' Creates background points within a study area for species distribution
#' modeling, with adaptive count based on the number of presence records.
#' Extracts environmental variables from a raster stack and merges with
#' presence data to produce a ready-to-model dataset.
#'
#' @param occurrences A `data.frame` with at least `lon` and `lat` columns
#'   (presence locations). May include additional columns that will be
#'   preserved in the output.
#' @param study_area A [cast_study_area] object defining the spatial extent
#'   for background sampling. If `NULL`, backgrounds are sampled from all
#'   non-NA cells of `raster_stack`.
#' @param raster_stack A `terra::SpatRaster` of environmental variables.
#'   Used both for extracting values at occurrence and background points
#'   and for defining valid (non-NA) sampling cells.
#' @param n_bg Integer or `NULL`. Number of background points. If `NULL`
#'   (default), adaptively determined as `clamp(ratio * n_presence,
#'   min_bg, max_bg)`.
#' @param ratio Numeric. Multiplier for adaptive background count.
#'   `n_bg = ratio * n_presence`. Default `2`.
#' @param min_bg Integer. Minimum background points. Default `500L`.
#' @param max_bg Integer. Maximum background points. Default `20000L`.
#' @param strategy Character. Background sampling strategy:
#'   - `"random"` (default): uniform random sampling from valid cells.
#'   - `"environmental"`: stratified random sampling in environmental space
#'     (partitions environmental PCA space into bins and samples uniformly
#'     across bins).
#' @param cell_thin Logical. If `TRUE` (default), ensures only one
#'   occurrence per raster cell (removes spatial duplicates at raster
#'   resolution).
#' @param exclude_presence Logical. If `TRUE` (default), background points
#'   cannot fall in cells occupied by occurrences.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param verbose Logical. Print informational messages. Default `TRUE`.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{lon, lat}{Coordinates.}
#'   \item{presence}{Integer (1 = presence, 0 = background).}
#'   \item{...}{One column per layer in `raster_stack` with extracted values,
#'     followed by any additional columns from `occurrences` (kept for
#'     presence rows, `NA` for background rows).}
#' }
#' Rows with any `NA` in environmental variables are removed.
#'
#' @section Interpretation and sampling bias:
#' Background (pseudo-absence) points define an arbitrary reference prevalence,
#' so downstream model outputs are **relative habitat suitability**, not
#' calibrated probabilities of occurrence: changing `ratio` shifts predicted
#' values up or down without changing their spatial ranking. Presence records
#' are also subject to sampling bias (accessibility, survey effort); unless
#' corrected, the models partly describe where the species was *observed*
#' rather than where it *occurs*. Carry both caveats into any HSS, effect, or
#' projection interpretation.
#'
#' @references
#' Barbet-Massin, M. et al. (2012). Selecting pseudo-absences for species
#' distribution models: how, where and how many?
#' *Methods in Ecology and Evolution*, 3(2), 327-338.
#'
#' Phillips, S. J. et al. (2009). Sample selection bias and presence-only
#' distribution models: implications for background and pseudo-absence data.
#' *Ecological Applications*, 19(1), 181-197.
#'
#' @seealso [cast_study_area()], [cast_prepare()]
#'
#' @export
cast_background <- function(occurrences,
                            study_area = NULL,
                            raster_stack,
                            n_bg = NULL,
                            ratio = 2,
                            min_bg = 500L,
                            max_bg = 20000L,
                            strategy = c("random", "environmental"),
                            cell_thin = TRUE,
                            exclude_presence = TRUE,
                            seed = NULL,
                            verbose = TRUE) {
  check_suggested("terra", "for raster extraction")
  strategy <- match.arg(strategy)

  # ---- Validate inputs -------------------------------------------------------
  if (!is.data.frame(occurrences)) {
    cli::cli_abort("{.arg occurrences} must be a data.frame.")
  }
  if (!all(c("lon", "lat") %in% names(occurrences))) {
    cli::cli_abort("{.arg occurrences} must have {.val lon} and {.val lat} columns.")
  }
  if (!inherits(raster_stack, "SpatRaster")) {
    cli::cli_abort("{.arg raster_stack} must be a {.cls SpatRaster}.")
  }
  if (!is.null(study_area) && !inherits(study_area, "cast_study_area")) {
    cli::cli_abort("{.arg study_area} must be a {.cls cast_study_area} or NULL.")
  }

  # The study-area mask is indexed against raster_stack cell numbers below,
  # so mismatched geometry would silently sample the wrong cells (M19).
  if (!is.null(study_area)) {
    geom_ok <- tryCatch(
      terra::compareGeom(study_area$mask, raster_stack, lyrs = FALSE),
      error = function(e) e
    )
    if (!isTRUE(geom_ok)) {
      cli::cli_abort(c(
        "{.arg study_area} mask and {.arg raster_stack} have incompatible geometry.",
        "x" = if (inherits(geom_ok, "error")) conditionMessage(geom_ok) else
          "compareGeom() mismatch",
        "i" = "Build the study area from the same raster grid you sample from."
      ))
    }
  }

  if (!is.null(seed)) set.seed(seed)

  # ---- Cell-thinning of occurrences -------------------------------------------
  occ_xy <- as.matrix(occurrences[, c("lon", "lat")])
  occ_cells <- terra::cellFromXY(raster_stack, occ_xy)

  if (cell_thin) {
    keep <- !duplicated(occ_cells) & !is.na(occ_cells)
    occurrences <- occurrences[keep, , drop = FALSE]
    occ_xy <- occ_xy[keep, , drop = FALSE]
    occ_cells <- occ_cells[keep]
    if (verbose) {
      cli::cli_inform(
        "Cell-thinned: {sum(!keep)} duplicate cell{?s} removed, {nrow(occurrences)} remain."
      )
    }
  }

  n_pres <- nrow(occurrences)

  # ---- Determine number of background points ----------------------------------
  if (is.null(n_bg)) {
    n_bg <- as.integer(round(ratio * n_pres))
    n_bg <- max(min_bg, min(max_bg, n_bg))
    if (verbose) {
      cli::cli_inform(
        "Adaptive background: {ratio} x {n_pres} = {n_bg} points (clamped to [{min_bg}, {max_bg}])."
      )
    }
  }

  # ---- Identify valid sampling cells -------------------------------------------
  if (!is.null(study_area)) {
    # Mask the raster stack by study area
    ref_mask <- study_area$mask
  } else {
    ref_mask <- !is.na(raster_stack[[1]])
    ref_mask[ref_mask == 0] <- NA
  }

  # Get all valid cell indices: inside the mask AND non-NA in every layer
  # (per-cell anyNA across layers, not just the first).
  n_na_layers <- terra::app(is.na(raster_stack), fun = "sum")
  valid_r <- !is.na(ref_mask) & (n_na_layers == 0)
  valid_cells <- which(as.logical(terra::values(valid_r, mat = FALSE)))

  # Exclude presence cells if requested
  if (exclude_presence) {
    valid_cells <- setdiff(valid_cells, occ_cells)
  }

  if (length(valid_cells) < n_bg) {
    if (verbose) {
      cli::cli_warn(
        "Only {length(valid_cells)} valid cells available; sampling with replacement."
      )
    }
    sample_replace <- TRUE
  } else {
    sample_replace <- FALSE
  }

  # ---- Sample background cells ------------------------------------------------
  if (strategy == "random") {
    bg_cells <- sample(valid_cells, size = n_bg, replace = sample_replace)

  } else if (strategy == "environmental") {
    # Environmental stratification via PCA binning
    bg_cells <- .sample_environmental(
      raster_stack, valid_cells, n_bg, seed, sample_replace
    )
  }

  # ---- Extract coordinates and environmental values ----------------------------
  bg_xy <- terra::xyFromCell(raster_stack, bg_cells)

  # Extract env values for both presence and background
  all_cells <- c(occ_cells, bg_cells)
  all_xy <- rbind(occ_xy, bg_xy)
  env_df <- as.data.frame(raster_stack[all_cells])

  # Build output data.frame
  out_df <- data.frame(
    lon = all_xy[, 1],
    lat = all_xy[, 2],
    presence = c(rep(1L, n_pres), rep(0L, n_bg)),
    stringsAsFactors = FALSE
  )
  out_df <- cbind(out_df, env_df)

  # Remove rows with NA in any environmental variable, then top the
  # background set back up to the requested size. `used_cells` accumulates
  # every background cell ever sampled so no cell is topped up twice.
  used_cells <- bg_cells
  topup_rounds <- 0L
  while (topup_rounds < 5L) {
    complete <- stats::complete.cases(out_df)
    n_bg_ok <- sum(complete & out_df$presence == 0)
    if (n_bg_ok >= n_bg) break
    pool <- setdiff(valid_cells, used_cells)
    if (!length(pool)) break
    topup_rounds <- topup_rounds + 1L
    add <- sample(pool, min(length(pool), n_bg - n_bg_ok))
    used_cells <- c(used_cells, add)
    add_xy <- terra::xyFromCell(raster_stack, add)
    add_env <- as.data.frame(raster_stack[add])
    add_df <- cbind(
      data.frame(lon = add_xy[, 1], lat = add_xy[, 2], presence = 0L,
                 stringsAsFactors = FALSE),
      add_env
    )
    out_df <- rbind(out_df, add_df)
  }
  keep_rows <- stats::complete.cases(out_df)
  out_df <- out_df[keep_rows, , drop = FALSE]
  rownames(out_df) <- NULL

  # Preserve additional occurrence columns: presence rows keep their values,
  # background rows are padded with NA.
  extra_cols <- setdiff(names(occurrences), c("lon", "lat"))
  if (length(extra_cols)) {
    pres_part <- occurrences[keep_rows[seq_len(n_pres)], extra_cols,
                             drop = FALSE]
    bg_part <- as.data.frame(lapply(pres_part, function(v) {
      v[rep(NA_integer_, sum(out_df$presence == 0))]
    }))
    names(bg_part) <- extra_cols
    out_df <- cbind(out_df, rbind(pres_part, bg_part))
    rownames(out_df) <- NULL
  }

  n_removed <- n_bg - sum(out_df$presence == 0)
  n_removed <- max(0L, n_removed)
  n_pres_final <- sum(out_df$presence == 1)
  n_bg_final <- sum(out_df$presence == 0)

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Background sampling complete:",
      " " = "Strategy: {.val {strategy}}",
      " " = "Presences: {n_pres_final} | Backgrounds: {n_bg_final}",
      if (n_removed > 0)
        c("!" = "{n_removed} rows removed due to NA environmental values.")
    ))
  }

  out_df
}


# ---- Internal helpers ---------------------------------------------------------

#' Environmental-space stratified background sampling
#' @keywords internal
#' @noRd
.sample_environmental <- function(raster_stack, valid_cells, n_bg,
                                  seed, replace) {
  # Extract env values at valid cells (sample subset if too many)
  max_extract <- min(length(valid_cells), 50000L)
  if (length(valid_cells) > max_extract) {
    sub_idx <- sample.int(length(valid_cells), max_extract)
    sub_cells <- valid_cells[sub_idx]
  } else {
    sub_cells <- valid_cells
    sub_idx <- seq_along(valid_cells)
  }

  env_vals <- as.data.frame(raster_stack[sub_cells])
  env_complete <- stats::complete.cases(env_vals)
  env_vals <- env_vals[env_complete, , drop = FALSE]
  sub_cells_clean <- sub_cells[env_complete]

  if (nrow(env_vals) < n_bg) {
    # Not enough valid cells after NA removal; fall back to random
    return(sample(sub_cells_clean, size = n_bg, replace = TRUE))
  }

  # PCA on env values
  pca <- tryCatch(
    stats::prcomp(env_vals, center = TRUE, scale. = TRUE, rank. = 2),
    error = function(e) NULL
  )

  if (is.null(pca)) {
    return(sample(sub_cells_clean, size = n_bg, replace = replace))
  }

  # Bin PC1 x PC2 space into grid
  scores <- pca$x[, 1:min(2, ncol(pca$x)), drop = FALSE]
  n_bins <- min(20L, ceiling(sqrt(nrow(scores))))

  # Create bin IDs
  bin_ids <- rep(1L, nrow(scores))
  for (j in seq_len(ncol(scores))) {
    breaks <- seq(
      min(scores[, j]) - 1e-6,
      max(scores[, j]) + 1e-6,
      length.out = n_bins + 1
    )
    bin_ids <- bin_ids + (as.integer(cut(scores[, j], breaks)) - 1L) *
      (n_bins^(j - 1L))
  }

  # Sample evenly across bins
  unique_bins <- unique(bin_ids)
  per_bin <- ceiling(n_bg / length(unique_bins))

  sampled_idx <- integer(0)
  for (b in unique_bins) {
    in_bin <- which(bin_ids == b)
    take <- min(per_bin, length(in_bin))
    sampled_idx <- c(sampled_idx,
                     sample(in_bin, size = take, replace = replace))
  }

  # Trim to exactly n_bg
  if (length(sampled_idx) > n_bg) {
    sampled_idx <- sample(sampled_idx, n_bg)
  } else if (length(sampled_idx) < n_bg) {
    extra <- sample(seq_along(sub_cells_clean),
                    n_bg - length(sampled_idx), replace = TRUE)
    sampled_idx <- c(sampled_idx, extra)
  }

  sub_cells_clean[sampled_idx]
}
