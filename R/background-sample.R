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
#'   - `"sre"`: surface-range-envelope sampling — cells whose environmental
#'     values fall inside the per-variable
#'     `[sre_quantile, 1 - sre_quantile]` quantile range of the presence
#'     environments (Barbet-Massin et al. 2012). Environmentally plausible,
#'     geographically unconstrained pseudo-absences; harder than random.
#'   - `"disk"`: ring sampling — cells at a distance between `disk_min` and
#'     `disk_max` map units from the nearest presence (Barbet-Massin et al.
#'     2012). Geographically close, hard pseudo-absences.
#' @param user_table Optional `data.frame` with `lon` and `lat` columns.
#'   When supplied, the background sample is taken from exactly these points
#'   (Phillips et al. 2009-style user-defined background): `strategy` and the
#'   adaptive count are ignored, `n_bg` defaults to the number of usable
#'   points, and points are only filtered by raster coverage, duplicate cells,
#'   presence-cell exclusion, and `NA` environments. Intended for already
#'   available species records used as background, or expert-chosen controls.
#' @param bias_raster Optional `terra::SpatRaster` (single layer, same grid
#'   as `raster_stack`). When supplied, random sampling draws cells with
#'   probability proportional to the raster values (Phillips et al. 2009
#'   target-group background): e.g. a raster of the sampling effort of the
#'   target taxon group. Non-finite and negative values are treated as zero
#'   weight; requires `strategy = "random"`.
#' @param sre_quantile Numeric in `[0, 0.5)`. Envelope tail probability for
#'   `strategy = "sre"`: each variable's range is trimmed at this quantile on
#'   both sides of the presence distribution. `0` keeps the full min-max
#'   envelope. Default `0.025`.
#' @param disk_min,disk_max Numeric. Distance band (map units; degrees for
#'   lon/lat rasters, consistent with `buffer` in [cast_cv()]) around
#'   presences for `strategy = "disk"`. `disk_min` defaults to `0`;
#'   `disk_max` to half the shortest raster extent side. Must satisfy
#'   `0 <= disk_min < disk_max`.
#' @param bin_method Binning method for `strategy = "environmental"`:
#'   `"pca2d"` (default) stratifies a PC1 x PC2 grid as before; `"kmeans"`
#'   partitions the scaled environment space with k-means (the sdm
#'   convention), which adapts to correlated drivers that collapse the
#'   first two PCs.
#' @param bin_k Integer >= 2. Bin count: the grid side for `"pca2d"`
#'   (capped at `sqrt(n)`), the number of cluster centers for `"kmeans"`.
#'   Default `20`.
#' @param cell_thin Logical. If `TRUE` (default), ensures only one
#'   occurrence per raster cell (removes spatial duplicates at raster
#'   resolution). Occurrences falling outside the raster or on NA-valued
#'   cells are dropped with a warning reporting the lost share.
#' @param exclude_presence Logical. If `TRUE` (default), background points
#'   cannot fall in cells occupied by occurrences.
#' @param n_rep Integer >= 1. Number of independent pseudo-absence
#'   replicate sets (Barbet-Massin et al. 2012 fit separately to each
#'   replicate to propagate PA-selection uncertainty). With `n_rep = 1`
#'   (default) a single `data.frame` is returned as before; with
#'   `n_rep > 1` a named list of such data frames (class
#'   `cast_background_reps`), each sampled under its own sub-seed.
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
                            strategy = c("random", "environmental", "sre", "disk"),
                            user_table = NULL,
                            bias_raster = NULL,
                            sre_quantile = 0.025,
                            disk_min = NULL,
                            disk_max = NULL,
                            bin_method = c("pca2d", "kmeans"),
                            bin_k = 20L,
                            cell_thin = TRUE,
                            exclude_presence = TRUE,
                            n_rep = 1L,
                            seed = NULL,
                            verbose = TRUE) {
  check_suggested("terra", "for raster extraction")
  strategy <- match.arg(strategy)
  bin_method <- match.arg(bin_method)
  if (!is.numeric(bin_k) || length(bin_k) != 1L || !is.finite(bin_k) ||
      bin_k < 2L) {
    cli::cli_abort("{.arg bin_k} must be a single number >= 2.")
  }
  bin_k <- as.integer(bin_k)

  use_user <- !is.null(user_table)
  use_bias <- !is.null(bias_raster)
  if (use_user && use_bias) {
    cli::cli_abort("Pass either {.arg user_table} or {.arg bias_raster}, not both.")
  }
  if (use_bias && strategy != "random") {
    cli::cli_abort(c(
      "{.arg bias_raster} combines with {.code strategy = \"random\"} only.",
      i = "Bias weighting and the environmental/SRE/disk designs are competing sampling designs."
    ))
  }

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
  if (use_user) {
    if (!is.data.frame(user_table) ||
        !all(c("lon", "lat") %in% names(user_table))) {
      cli::cli_abort("{.arg user_table} must be a data.frame with {.val lon} and {.val lat} columns.")
    }
    if (!nrow(user_table)) {
      cli::cli_abort("{.arg user_table} is empty.")
    }
  }
  if (use_bias) {
    if (!inherits(bias_raster, "SpatRaster")) {
      cli::cli_abort("{.arg bias_raster} must be a {.cls SpatRaster}.")
    }
    if (terra::nlyr(bias_raster) > 1L) {
      cli::cli_warn(
        "{.arg bias_raster} has {terra::nlyr(bias_raster)} layers; using the first one."
      )
      bias_raster <- bias_raster[[1L]]
    }
    # Weight lookup indexes cells against raster_stack below, so mismatched
    # geometry would silently weight the wrong cells.
    geom_ok <- tryCatch(
      terra::compareGeom(bias_raster, raster_stack, lyrs = FALSE),
      error = function(e) e
    )
    if (!isTRUE(geom_ok)) {
      cli::cli_abort(c(
        "{.arg bias_raster} and {.arg raster_stack} have incompatible geometry.",
        "x" = if (inherits(geom_ok, "error")) conditionMessage(geom_ok) else
          "compareGeom() mismatch",
        "i" = "Build the bias raster on the same raster grid you sample from."
      ))
    }
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

  if (strategy == "sre") {
    if (!is.numeric(sre_quantile) || length(sre_quantile) != 1L ||
        !is.finite(sre_quantile) || sre_quantile < 0 || sre_quantile >= 0.5) {
      cli::cli_abort("{.arg sre_quantile} must be a single number in [0, 0.5).")
    }
  }
  if (strategy == "disk") {
    if (is.null(disk_min)) disk_min <- 0
    if (is.null(disk_max)) {
      ext <- terra::ext(raster_stack)
      disk_max <- 0.5 * min(ext$xmax - ext$xmin, ext$ymax - ext$ymin)
    }
    if (!is.numeric(disk_min) || length(disk_min) != 1L ||
        !is.finite(disk_min) || disk_min < 0 ||
        !is.numeric(disk_max) || length(disk_max) != 1L ||
        !is.finite(disk_max) || disk_max <= disk_min) {
      cli::cli_abort(c(
        "{.arg disk_min} and {.arg disk_max} must satisfy 0 <= disk_min < disk_max (map units).",
        i = "Both are resolved in the CRS of {.arg raster_stack} (degrees for lon/lat data)."
      ))
    }
  }

  n_rep <- as.integer(n_rep)
  if (is.na(n_rep) || n_rep < 1L) {
    cli::cli_abort("{.arg n_rep} must be an integer >= 1.")
  }

  # ---- Cell-thinning of occurrences -------------------------------------------
  # Deterministic (seed-independent), so it runs once outside the replicate
  # loop: every replicate starts from the same thinned occurrence set.
  occ_xy <- as.matrix(occurrences[, c("lon", "lat")])
  occ_cells <- terra::cellFromXY(raster_stack, occ_xy)

  if (cell_thin) {
    # Accounting note: `duplicated()` treats NAs as matching values, so a
    # duplicate NA cell would be counted both as outside and as duplicate.
    # Split the flags to keep the two counts exact.
    n_outside <- sum(is.na(occ_cells))
    dup_flag <- duplicated(occ_cells) & !is.na(occ_cells)
    n_dup <- sum(dup_flag)
    keep <- !dup_flag & !is.na(occ_cells)
    occurrences <- occurrences[keep, , drop = FALSE]
    occ_xy <- occ_xy[keep, , drop = FALSE]
    occ_cells <- occ_cells[keep]
    if (verbose) {
      cli::cli_inform(
        "Cell-thinned: {n_dup} duplicate cell{?s} removed, {nrow(occurrences)} remain."
      )
    }
    # Records outside the raster or on NA-valued cells used to vanish into
    # the generic drop count; losing coverage is worth a warning with the
    # lost share (S12).
    if (n_outside > 0) {
      pct_out <- round(100 * n_outside / (nrow(occurrences) + n_outside + n_dup), 1)
      cli::cli_warn(c(
        "{n_outside} occurrence{?s} ({pct_out}%) fall outside {.arg raster_stack} or on NA-valued cells and were dropped before sampling.",
        "i" = "Check the raster extent/mask against the occurrence coordinates."
      ))
    }
  }

  n_pres <- nrow(occurrences)

  # One sampling replicate. With n_rep = 1 the plain user seed is kept for
  # backwards compatibility; with n_rep > 1 each replicate draws under its
  # own sub-seed so the sets are independent.
  bg_once <- function(rep_i) {
    if (!is.null(seed)) {
      set.seed(if (n_rep == 1L) seed else seed + rep_i - 1L)
    }

  # ---- Background cell selection -----------------------------------------------
  if (use_user) {
    # User-defined background: exactly the supplied points. Filtered only by
    # raster coverage, duplicate cells, presence-cell exclusion and, later,
    # NA environments; the study-area mask does not apply.
    u_cells <- terra::cellFromXY(raster_stack,
                                 as.matrix(user_table[, c("lon", "lat")]))
    n_outside <- sum(is.na(u_cells))
    if (n_outside > 0 && verbose) {
      cli::cli_warn(
        "{n_outside} user background point{?s} fall outside the raster; dropped."
      )
    }
    u_cells <- unique(stats::na.omit(u_cells))
    if (exclude_presence) u_cells <- setdiff(u_cells, occ_cells)
    if (!length(u_cells)) {
      cli::cli_abort("No usable points in {.arg user_table} after filtering.")
    }
    if (is.null(n_bg)) {
      n_bg <- length(u_cells)
    } else {
      n_bg <- as.integer(n_bg)
      if (n_bg > length(u_cells)) {
        cli::cli_warn(
          "Only {length(u_cells)} usable user background point{?s}; n_bg reduced from {n_bg}."
        )
        n_bg <- length(u_cells)
      }
    }
    valid_cells <- u_cells  # the top-up pool is the remaining user points
    bg_cells <- if (length(u_cells) > n_bg) sample(u_cells, n_bg) else u_cells
  } else {
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

    # ---- Strategy-specific candidate pools ---------------------------------------
    # The pool REPLACES valid_cells so the NA top-up loop below stays inside
    # the same strategy design instead of drifting back to uniform sampling.
    if (strategy == "sre") {
      valid_cells <- .sre_pool(raster_stack, occ_cells, valid_cells,
                               sre_quantile)
      if (verbose) {
        cli::cli_inform(
          "SRE envelope (q = {sre_quantile}): {length(valid_cells)} candidate cells."
        )
      }
    } else if (strategy == "disk") {
      valid_cells <- .disk_pool(raster_stack, occ_xy, occ_cells, valid_cells,
                                disk_min, disk_max)
      if (verbose) {
        cli::cli_inform(
          "Disk band [{disk_min}, {disk_max}] map units: {length(valid_cells)} candidate cells."
        )
      }
    }
    if (!length(valid_cells)) {
      hint <- switch(strategy,
        sre = "Widen {.arg sre_quantile} toward 0.5 or use {.code strategy = \"random\"}.",
        disk = "Lower {.arg disk_min} or raise {.arg disk_max}, or use {.code strategy = \"random\"}.",
        "Check the study-area mask and raster coverage."
      )
      cli::cli_abort(c(
        "No candidate background cells for strategy {.val {strategy}}.",
        "i" = hint
      ))
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
    if (use_bias) {
      # Target-group background (Phillips et al. 2009): cells are drawn with
      # probability proportional to the bias surface.
      bias_vals <- as.numeric(terra::values(bias_raster, mat = FALSE))
      w <- bias_vals[valid_cells]
      w[!is.finite(w) | w < 0] <- 0
      if (sum(w) <= 0) {
        cli::cli_abort(
          "{.arg bias_raster} has no positive weights inside the valid sampling area."
        )
      }
      bg_cells <- sample(valid_cells, size = n_bg, prob = w,
                         replace = sample_replace)
    } else if (strategy == "random" || strategy == "sre" || strategy == "disk") {
      # sre / disk pools are already restricted; the draw is uniform within.
      bg_cells <- sample(valid_cells, size = n_bg, replace = sample_replace)

    } else if (strategy == "environmental") {
      # Environmental stratification via PCA-grid or k-means binning
      bg_cells <- .sample_environmental(
        raster_stack, valid_cells, n_bg, seed, sample_replace,
        bin_method = bin_method, bin_k = bin_k
      )
    }
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
    take_n <- min(length(pool), n_bg - n_bg_ok)
    add <- if (use_bias) {
      # Top-up follows the same bias surface as the initial draw.
      wp <- bias_vals[pool]
      wp[!is.finite(wp) | wp < 0] <- 0
      if (sum(wp) > 0) {
        sample(pool, take_n, prob = wp)
      } else {
        sample(pool, take_n)
      }
    } else {
      sample(pool, take_n)
    }
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
    strategy_label <- if (use_user) "user_table" else
      if (use_bias) "bias-weighted random" else strategy
    cli::cli_inform(c(
      "v" = "Background sampling complete:",
      " " = "Strategy: {.val {strategy_label}}",
      " " = "Presences: {n_pres_final} | Backgrounds: {n_bg_final}",
      if (n_removed > 0)
        c("!" = "{n_removed} rows removed due to NA environmental values.")
    ))
  }

    out_df
  }

  reps <- lapply(seq_len(n_rep), bg_once)
  if (n_rep == 1L) return(reps[[1]])
  reps <- stats::setNames(reps, sprintf("rep%d", seq_len(n_rep)))
  class(reps) <- c("cast_background_reps", "list")
  if (verbose) {
    cli::cli_inform(
      "Generated {n_rep} pseudo-absence replicate sets ({names(reps)})."
    )
  }
  reps
}


# ---- Internal helpers ---------------------------------------------------------

#' Surface Range Envelope candidate pool (SRE)
#'
#' Cells whose environmental values fall inside the per-variable
#' `[q, 1 - q]` quantile range of the presence environments
#' (Barbet-Massin et al. 2012, "sre"): environmentally plausible,
#' geographically unconstrained pseudo-absences.
#' @keywords internal
#' @noRd
.sre_pool <- function(raster_stack, occ_cells, valid_cells, sre_quantile) {
  occ_cells <- occ_cells[!is.na(occ_cells)]
  pres_env <- as.data.frame(raster_stack[occ_cells])
  pres_env <- pres_env[stats::complete.cases(pres_env), , drop = FALSE]
  if (nrow(pres_env) < 2L) {
    cli::cli_abort(
      "SRE needs environmental values at at least 2 presence cells."
    )
  }
  probs <- c(sre_quantile, 1 - sre_quantile)
  bounds <- apply(pres_env, 2, stats::quantile, probs = probs, names = FALSE)
  lo <- bounds[1, ]
  hi <- bounds[2, ]
  pool <- integer(0)
  chunk <- 200000L  # bounded-memory sweep over the candidate cells
  for (s in seq(1L, length(valid_cells), by = chunk)) {
    idx <- s:min(s + chunk - 1L, length(valid_cells))
    cells_chunk <- valid_cells[idx]
    ev <- as.data.frame(raster_stack[cells_chunk])
    ok <- stats::complete.cases(ev)
    ev <- ev[ok, , drop = FALSE]
    inside <- rep(TRUE, nrow(ev))
    for (j in seq_along(ev)) {
      inside <- inside & ev[[j]] >= lo[j] & ev[[j]] <= hi[j]
    }
    pool <- c(pool, cells_chunk[ok][inside])
  }
  pool
}

#' Disk-band candidate pool
#'
#' Cells whose planar Euclidean distance to the nearest presence coordinate
#' lies in `[disk_min, disk_max]` map units (Barbet-Massin et al. 2012,
#' "disk"): geographically close, hard pseudo-absences. Distances are
#' computed between cell centroids and presence coordinates in the same
#' units as `buffer` in [cast_cv()] (degrees for lon/lat data) rather than
#' via [terra::distance()], which switches to geodesic metres on
#' lon/lat geometries.
#' @keywords internal
#' @noRd
.disk_pool <- function(raster_stack, occ_xy, occ_cells, valid_cells,
                       disk_min, disk_max) {
  ok_p <- !is.na(occ_cells)
  if (!any(ok_p)) {
    cli::cli_abort("Disk strategy needs presences inside the raster extent.")
  }
  occ_xy <- occ_xy[ok_p, , drop = FALSE]
  if (!length(valid_cells)) return(integer(0))
  xy_v <- terra::xyFromCell(raster_stack, valid_cells)
  # Running row-min of squared distances, one presence at a time: O(n_pres)
  # vectorised passes, O(n_cells) memory, no terra version-dependent units.
  dmin2 <- rep(Inf, nrow(xy_v))
  for (p in seq_len(nrow(occ_xy))) {
    d2 <- (xy_v[, 1] - occ_xy[p, 1])^2 + (xy_v[, 2] - occ_xy[p, 2])^2
    dmin2 <- pmin(dmin2, d2)
  }
  valid_cells[dmin2 >= disk_min^2 & dmin2 <= disk_max^2]
}

#' Environmental-space stratified background sampling
#' @keywords internal
#' @noRd
.sample_environmental <- function(raster_stack, valid_cells, n_bg,
                                  seed, replace, bin_method = "pca2d",
                                  bin_k = 20L) {
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

  if (identical(bin_method, "kmeans")) {
    # sdm-style k-means partition of the scaled environment space: the
    # cluster id plays the role of the PCA-grid bin id in the stratified
    # draw below. centers are capped at the number of rows; a degenerate
    # partition (zero-variance space, duplicated rows) falls back to random.
    km <- tryCatch(
      stats::kmeans(scale(env_vals),
                    centers = min(as.integer(bin_k), nrow(env_vals)),
                    iter.max = 50, nstart = 3),
      error = function(e) NULL
    )
    if (is.null(km)) {
      return(sample(sub_cells_clean, size = n_bg, replace = replace))
    }
    bin_ids <- as.integer(km$cluster)
  } else {
    # PCA on env values
    pca <- tryCatch(
      stats::prcomp(env_vals, center = TRUE, scale. = TRUE, rank. = 2),
      error = function(e) NULL
    )

    if (is.null(pca)) {
      return(sample(sub_cells_clean, size = n_bg, replace = replace))
    }

    # Bin PC1 x PC2 space into a grid
    scores <- pca$x[, 1:min(2, ncol(pca$x)), drop = FALSE]
    n_bins <- min(as.integer(bin_k), ceiling(sqrt(nrow(scores))))

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
