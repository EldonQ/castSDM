# Factorial causal scenario attribution --------------------------------------
#
# The package's causal layer is deliberately narrow: a 2 x 2 factorial design
# (current, climate-only, LULC-only, joint) defines the counterfactual
# contrasts. No DAG, adjustment set, propensity model or causal-discovery
# algorithm is introduced. The contrasts are model-based causal scenario
# attributions: they are causal only under the stated scenario assumptions and
# do not identify a real-world total effect from observational data.

#' Factorial causal attribution of climate and land-use scenario effects
#'
#' Given four aligned suitability surfaces from a 2 x 2 factorial scenario
#' design, compute the climate main effect, the land-use main effect, the joint
#' change, the interaction, and the two conditional effects at the joint
#' scenario. The design is:
#'
#' \itemize{
#'   \item `current`: current climate and current land use;
#'   \item `climate`: future climate and current land use;
#'   \item `lulc`: current climate and future land use;
#'   \item `joint`: future climate and future land use.
#' }
#'
#' The contrasts are defined on the supplied suitability scale:
#'
#' \deqn{climate main = climate - current}
#' \deqn{land-use main = lulc - current}
#' \deqn{joint change = joint - current}
#' \deqn{interaction = joint - climate - lulc + current}
#' \deqn{climate at joint = joint - lulc}
#' \deqn{land-use at joint = joint - climate}
#'
#' These are model-based causal scenario contrasts, not identified ecological
#' causal effects. The function records the assumptions in `contract`; it does
#' not test them. The factorial design itself fixes the other component, so no
#' user-supplied adjustment set is required.
#'
#' @param current,joint,climate,lulc A single-layer `terra::SpatRaster` or a
#'   path to one. These are the four suitability surfaces of the factorial
#'   design. The layers must share the same geometry.
#' @param current_binary,joint_binary,climate_binary,lulc_binary Optional
#'   binary range rasters (0/1) for the same four scenarios. Supply all four or
#'   none. If none are supplied and `threshold` is given, binary rasters are
#'   derived from the continuous surfaces.
#' @param threshold Optional numeric threshold in the continuous suitability
#'   scale, used only when binary rasters are not supplied.
#' @param mask Optional single-layer mask. Cells that are `NA` or `0` are
#'   excluded from all contrasts and summaries.
#' @param support Optional single-layer support raster (1 = supported). When
#'   supplied, summaries report the supported subset separately. The support
#'   raster is not used to delete the raw contrasts.
#' @param current_sd,joint_sd,climate_sd,lulc_sd Optional cross-model standard
#'   deviation rasters for the four scenarios. When all four are supplied,
#'   an approximate sqrt-sum surface-SD propagation is computed under an independence assumption
#'   and reported. It describes surface-level model disagreement, not the SD of the contrast, and is not a confidence interval.
#' @param scenario,ssp,period Optional scenario metadata stored on the result.
#' @param filename Optional GeoTIFF path for the contrast stack. The stack
#'   contains the four input surfaces, the six contrasts and any SD layers.
#' @param overwrite,compression Passed to [terra::writeRaster()].
#' @param verbose Print a short progress message.
#'
#' @return A `cast_factorial` object with components:
#' \describe{
#'   \item{summary}{Continuous contrast summaries.}
#'   \item{binary}{Binary range counts and percentage changes, or `NULL`.}
#'   \item{raster}{The contrast `SpatRaster`.}
#'   \item{raster_file}{The written file, or `NA`.}
#'   \item{contract}{The claim level, contrast definitions and assumptions.}
#'   \item{scenario,ssp,period}{Scenario metadata.}
#' }
#'
#' @references
#' The factorial contrast is a model-based counterfactual comparison. It does
#' not identify a causal effect without assumptions about the fitted model and
#' the scenario manipulations.
#'
#' @examples
#' \dontrun{
#' cur <- terra::rast(matrix(0.2, 4, 4))
#' clim <- terra::rast(matrix(0.3, 4, 4))
#' lulc <- terra::rast(matrix(0.4, 4, 4))
#' joint <- terra::rast(matrix(0.5, 4, 4))
#' fa <- cast_factorial_attribution(cur, joint, clim, lulc, threshold = 0.5)
#' fa$summary
#' }
#' @export
cast_factorial_attribution <- function(
    current, joint, climate, lulc,
    current_binary = NULL, joint_binary = NULL,
    climate_binary = NULL, lulc_binary = NULL,
    threshold = NULL, mask = NULL, support = NULL,
    current_sd = NULL, joint_sd = NULL,
    climate_sd = NULL, lulc_sd = NULL,
    scenario = NULL, ssp = NULL, period = NULL,
    filename = NULL, overwrite = FALSE, compression = "LZW",
    verbose = TRUE) {
  check_suggested("terra", "for factorial attribution")
  if (!requireNamespace("terra", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg terra} is required for factorial attribution.")
  }

  as_hss <- function(x, arg) {
    if (is.character(x)) {
      if (length(x) != 1L || !file.exists(x)) {
        cli::cli_abort("{.arg {arg}} must be a readable raster path or a SpatRaster.")
      }
      x <- terra::rast(x)
    }
    if (!inherits(x, "SpatRaster")) {
      cli::cli_abort("{.arg {arg}} must be a SpatRaster or a raster path.")
    }
    if (terra::nlyr(x) != 1L) {
      cli::cli_warn("{.arg {arg}} has {terra::nlyr(x)} layers; using the first.")
      x <- x[[1L]]
    }
    x
  }
  check_geom <- function(reference, x, arg) {
    same <- tryCatch(
      terra::compareGeom(reference, x, stopOnError = FALSE, messages = FALSE),
      error = function(e) FALSE)
    if (!isTRUE(same)) {
      cli::cli_abort("{.arg {arg}} is not aligned with {.arg current}.")
    }
    invisible(TRUE)
  }

  cur <- as_hss(current, "current")
  joi <- as_hss(joint, "joint")
  cli <- as_hss(climate, "climate")
  lul <- as_hss(lulc, "lulc")
  check_geom(cur, joi, "joint")
  check_geom(cur, cli, "climate")
  check_geom(cur, lul, "lulc")

  if (!is.null(mask)) {
    mask <- as_hss(mask, "mask")
    check_geom(cur, mask, "mask")
  }
  if (!is.null(support)) {
    support <- as_hss(support, "support")
    check_geom(cur, support, "support")
  }

  binary_args <- list(
    current = current_binary, joint = joint_binary,
    climate = climate_binary, lulc = lulc_binary)
  n_binary <- sum(!vapply(binary_args, is.null, logical(1)))
  if (n_binary %in% c(0L, 4L)) {
    if (n_binary == 0L && !is.null(threshold)) {
      if (!is.numeric(threshold) || length(threshold) != 1L ||
          !is.finite(threshold)) {
        cli::cli_abort("{.arg threshold} must be one finite number.")
      }
      current_binary <- terra::ifel(cur >= threshold, 1L, 0L)
      joint_binary <- terra::ifel(joi >= threshold, 1L, 0L)
      climate_binary <- terra::ifel(cli >= threshold, 1L, 0L)
      lulc_binary <- terra::ifel(lul >= threshold, 1L, 0L)
    }
  } else {
    cli::cli_abort("Supply all four binary rasters, or none of them.")
  }
  if (!is.null(current_binary)) {
    current_binary <- as_hss(current_binary, "current_binary")
    joint_binary <- as_hss(joint_binary, "joint_binary")
    climate_binary <- as_hss(climate_binary, "climate_binary")
    lulc_binary <- as_hss(lulc_binary, "lulc_binary")
    for (b in list(current_binary, joint_binary, climate_binary, lulc_binary)) {
      check_geom(cur, b, "binary")
    }
  }

  sd_args <- list(current_sd, joint_sd, climate_sd, lulc_sd)
  has_sd <- all(!vapply(sd_args, is.null, logical(1)))
  if (any(!vapply(sd_args, is.null, logical(1))) && !has_sd) {
    cli::cli_abort("Supply all four SD rasters, or none of them.")
  }
  if (has_sd) {
    current_sd <- as_hss(current_sd, "current_sd")
    joint_sd <- as_hss(joint_sd, "joint_sd")
    climate_sd <- as_hss(climate_sd, "climate_sd")
    lulc_sd <- as_hss(lulc_sd, "lulc_sd")
    for (s in list(current_sd, joint_sd, climate_sd, lulc_sd)) {
      check_geom(cur, s, "SD")
    }
  }

  finite <- !is.na(cur) & !is.na(joi) & !is.na(cli) & !is.na(lul)
  if (!is.null(mask)) {
    finite <- finite & !is.na(mask) & mask > 0
  }
  apply_finite <- function(r) {
    terra::mask(r, finite, maskvalues = c(0, NA))
  }

  climate_main <- apply_finite(cli - cur)
  lulc_main <- apply_finite(lul - cur)
  joint_change <- apply_finite(joi - cur)
  interaction <- apply_finite(joi - cli - lul + cur)
  climate_conditional <- apply_finite(joi - lul)
  lulc_conditional <- apply_finite(joi - cli)

  contrasts <- c(
    cur, joi, cli, lul,
    climate_main, lulc_main, joint_change, interaction,
    climate_conditional, lulc_conditional)
  names(contrasts) <- c(
    "current", "joint", "climate", "lulc",
    "climate_main", "lulc_main", "joint_change", "interaction",
    "climate_conditional", "lulc_conditional")
  effects <- list(
    climate_main = climate_main, lulc_main = lulc_main,
    joint_change = joint_change, interaction = interaction,
    climate_conditional = climate_conditional,
    lulc_conditional = lulc_conditional)

  contrast_sd <- NULL
  if (has_sd) {
    contrast_sd <- list(
      climate_main = sqrt(current_sd^2 + climate_sd^2),
      lulc_main = sqrt(current_sd^2 + lulc_sd^2),
      joint_change = sqrt(current_sd^2 + joint_sd^2),
      interaction = sqrt(current_sd^2 + climate_sd^2 + lulc_sd^2 + joint_sd^2),
      climate_conditional = sqrt(joint_sd^2 + lulc_sd^2),
      lulc_conditional = sqrt(joint_sd^2 + climate_sd^2))
    sd_stack <- do.call(c, unname(contrast_sd))
    names(sd_stack) <- paste0("sd_", names(contrast_sd))
    contrasts <- c(contrasts, sd_stack)
  }

  summary_one <- function(r, name, sd = NULL) {
    valid <- !is.na(r)
    n_valid <- terra::global(valid, "sum", na.rm = TRUE)[[1]][1]
    mean_v <- terra::global(r, "mean", na.rm = TRUE)[[1]][1]
    mean_abs <- terra::global(abs(r), "mean", na.rm = TRUE)[[1]][1]
    n_supported <- NA_real_
    mean_supported <- NA_real_
    mean_abs_supported <- NA_real_
    support_share <- NA_real_
    if (!is.null(support)) {
      sup <- support > 0
      sup <- terra::ifel(is.na(sup), 0, sup)
      n_supported <- terra::global(valid & sup, "sum", na.rm = TRUE)[[1]][1]
      rs <- terra::mask(r, sup, maskvalues = c(0, NA))
      mean_supported <- terra::global(rs, "mean", na.rm = TRUE)[[1]][1]
      mean_abs_supported <- terra::global(abs(rs), "mean", na.rm = TRUE)[[1]][1]
      support_share <- if (is.finite(n_valid) && n_valid > 0) {
        n_supported / n_valid
      } else {
        NA_real_
      }
    }
    mean_sd <- if (is.null(sd)) NA_real_ else
      terra::global(sd, "mean", na.rm = TRUE)[[1]][1]
    data.frame(
      component = name, scale = "continuous",
      n_valid = n_valid, n_supported = n_supported,
      support_share = support_share,
      mean = mean_v, mean_abs = mean_abs,
      mean_supported = mean_supported,
      mean_abs_supported = mean_abs_supported,
      mean_surface_sd = mean_sd,
      stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, lapply(names(effects), function(nm) {
    summary_one(effects[[nm]], nm,
                if (is.null(contrast_sd)) NULL else contrast_sd[[nm]])
  }))

  binary <- NULL
  if (!is.null(current_binary)) {
    binary_values <- function(r) {
      r <- terra::ifel(is.na(r), 0L, r)
      terra::global(r, "sum", na.rm = TRUE)[[1]][1]
    }
    cb <- binary_values(current_binary)
    jb <- binary_values(joint_binary)
    clb <- binary_values(climate_binary)
    lub <- binary_values(lulc_binary)
    comps <- c(
      climate_main = clb - cb,
      lulc_main = lub - cb,
      joint_change = jb - cb,
      interaction = jb - clb - lub + cb,
      climate_conditional = jb - lub,
      lulc_conditional = jb - clb)
    future_n <- c(
      climate_main = clb, lulc_main = lub, joint_change = jb,
      interaction = NA_real_, climate_conditional = jb,
      lulc_conditional = jb)
    binary <- data.frame(
      component = names(comps),
      n_current_range = cb,
      n_future_range = unname(future_n),
      estimate = unname(comps),
      pct_of_current = if (is.finite(cb) && cb > 0) {
        100 * unname(comps) / cb
      } else {
        rep(NA_real_, length(comps))
      },
      stringsAsFactors = FALSE)
  }

  if (!is.null(filename)) {
    if (!is.character(filename) || length(filename) != 1L ||
        !nzchar(filename)) {
      cli::cli_abort("{.arg filename} must be one non-empty path.")
    }
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(
      contrasts, filename, overwrite = overwrite,
      wopt = list(datatype = "FLT4S",
                  gdal = c(paste0("COMPRESS=", compression), "TILED=YES")))
  }

  contract <- data.frame(
    key = c(
      "claim_level", "design", "current", "joint", "climate_only",
      "lulc_only", "climate_main", "lulc_main", "joint_change",
      "interaction", "climate_conditional", "lulc_conditional",
      "support", "surface_sd", "assumptions", "not_identified"),
    value = c(
      "model_based_causal_scenario_attribution",
      "2x2 factorial: current, climate-only, LULC-only, joint",
      "current climate + current land use",
      "future climate + future land use",
      "future climate + current land use",
      "current climate + future land use",
      "climate - current",
      "lulc - current",
      "joint - current",
      "joint - climate - lulc + current",
      "joint - lulc (climate effect conditional on future land use)",
      "joint - climate (land-use effect conditional on future climate)",
      if (is.null(support)) "not supplied; summaries use finite cells only" else
        "supported subset reported separately; raw contrasts retained",
      if (!has_sd) "not supplied" else
        "approximate sqrt-sum of scenario surface SDs; not a contrast SD or CI",
      paste(
        "Scenario components are well-defined manipulations;",
        "the fitted ensemble is a valid causal response surface for these scenarios;",
        "no omitted scenario component or unmeasured driver changes the contrast;",
        "counterfactual values lie within the training support;",
        "no interference between cells. The function does not test these assumptions."),
      paste(
        "Not identified: real-world total ecological causal effect, occupancy,",
        "absolute occurrence probability, or future realised distribution.",
        "Presence-background suitability is a relative model scale.")),
    stringsAsFactors = FALSE)

  structure(list(
    summary = summary,
    binary = binary,
    raster = contrasts,
    raster_file = if (is.null(filename)) NA_character_ else filename,
    contract = contract,
    scenario = scenario,
    ssp = ssp,
    period = period,
    threshold = threshold
  ), class = "cast_factorial")
}

#' @export
print.cast_factorial <- function(x, ...) {
  cli::cli_h1("castSDM Factorial Causal Attribution")
  cli::cli_ul(c(
    "Scenario: {x$scenario %||% 'unspecified'}",
    "SSP/period: {x$ssp %||% 'unspecified'} / {x$period %||% 'unspecified'}",
    "Claim: {x$contract$value[x$contract$key == 'claim_level']}"
  ))
  if (!is.null(x$summary) && nrow(x$summary)) {
    show <- x$summary[, c("component", "mean", "mean_abs",
                          "support_share", "mean_surface_sd")]
    print(show, row.names = FALSE)
  }
  if (!is.null(x$binary) && nrow(x$binary)) {
    cli::cli_h2("Binary range")
    print(x$binary[, c("component", "estimate", "pct_of_current")],
          row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.cast_factorial <- function(object, ...) {
  print(object)
  invisible(object)
}

#' @export
plot.cast_factorial <- function(x, layers = NULL, ...) {
  check_suggested("terra", "for plotting factorial contrasts")
  r <- x$raster
  if (is.null(layers)) {
    layers <- c("climate_main", "lulc_main", "joint_change", "interaction")
  }
  layers <- intersect(layers, names(r))
  if (!length(layers)) {
    cli::cli_abort("No requested contrast layer is present in the raster stack.")
  }
  terra::plot(r[[layers]], ...)
  invisible(x)
}
