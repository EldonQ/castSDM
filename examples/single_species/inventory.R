# Metadata inventory for the national 1 km reg library (legacy glo supported).
#
# The inventory is the only place that knows the directory layout. It writes
# one current-layer table and one future-layer table; the model, the screen and
# the scenario table consume those tables and never hard-code variable names.
#
# Current layers are classified by data availability, not by ecological
# judgement:
#   climate      bioclim present, future scenarios available
#   lulc         the selected land-use product, future scenarios available
#   static       soil / topography / hydrography / distance: no future product
#   current_only vegetation / anthropogenic / CLCD: no future product
# Only climate, lulc and static layers are eligible for a future projection
# model. current_only layers are inventoried for provenance but are not offered
# to the screen.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}

glo_inventory_files <- function(root) {
  cfg_dir <- file.path(root, "config")
  list(
    current_all = file.path(cfg_dir, "pool_current_all.csv"),
    eligible_glu = file.path(cfg_dir, "pool_current_eligible_glu.csv"),
    eligible_hlu = file.path(cfg_dir, "pool_current_eligible_hlu.csv"),
    future = file.path(cfg_dir, "future_layers.csv"),
    version = file.path(cfg_dir, "inventory_version.csv")
  )
}

case_covariate_scope <- function(data_root) {
  for (scope in c("reg", "glo"))
    if (dir.exists(file.path(data_root, "covariates", scope))) return(scope)
  stop("Environmental library not found below ", file.path(data_root, "covariates"),
    "; expected reg or legacy glo.", call. = FALSE)
}

.glo_scale <- function(layer) {
  if (grepl("_f1000$", layer)) return("f1000")
  if (grepl("_f5000$", layer)) return("f5000")
  if (grepl("_f10000$", layer)) return("f10000")
  if (grepl("_f20000$", layer)) return("f20000")
  "base"
}

.glo_stat <- function(layer) {
  if (grepl("_max$", layer)) return("max")
  if (grepl("_min$", layer)) return("min")
  "base"
}

.glo_strip_scale_stat <- function(x) {
  x <- sub("_f(1000|5000|10000|20000)$", "", x)
  x <- sub("_(max|min)$", "", x)
  x
}

.glo_strip_year <- function(x) {
  gsub("_(19|20|21)[0-9]{2}(?=_|$)", "", x, perl = TRUE)
}

# One current layer. `rel` is the path below covariates/reg/ or legacy glo/.
glo_current_meta_one <- function(path, rel, data_root) {
  seg <- strsplit(rel, "/", fixed = TRUE)[[1]]
  group <- seg[1]
  dataset <- paste(seg[1:2], collapse = "/")
  product <- seg[2]
  stem <- sub("[.]tif$", "", basename(path))
  layer <- if (identical(group, "bioclim")) {
    basename(dirname(path))
  } else {
    gsub("-", "_", sub("^(reg|glo)_", "", stem))
  }
  year <- suppressWarnings(as.integer(seg[grepl("^(19|20|21)[0-9]{2}$", seg)][1]))
  scale <- .glo_scale(layer)
  stat <- .glo_stat(layer)
  base <- .glo_strip_scale_stat(layer)
  construct_id <- if (identical(group, "bioclim")) {
    layer
  } else if (identical(group, "lulc")) {
    base <- sub("^lulc_(glu1km|hlu1km|clcd)_(19|20|21)[0-9]{2}_", "", base)
    paste0("lulc_", .glo_strip_year(base))
  } else if (identical(group, "soil")) {
    paste0("soil_", sub("^soil_csdlv2_", "", base))
  } else if (identical(group, "hydro")) {
    paste0("hydro_", sub("^hydro_", "", base))
  } else if (identical(group, "dist")) {
    paste0("dist_", sub("^dist_", "", base))
  } else if (identical(group, "topo")) {
    paste0("topo_", sub("^topo_chinadem_", "", base))
  } else if (identical(group, "vege")) {
    paste0("vege_", .glo_strip_year(sub("^vege_", "", base)))
  } else if (identical(group, "anthro")) {
    paste0("anthro_", product)
  } else {
    paste0(group, "_", base)
  }
  future_policy <- if (identical(group, "bioclim")) {
    "climate"
  } else if (identical(group, "lulc") && product %in% c("glu1km", "hlu1km")) {
    "lulc"
  } else if (group %in% c("soil", "topo", "hydro", "dist")) {
    "static"
  } else {
    "current_only"
  }
  data.frame(
    group = group, dataset = dataset, product = product,
    layer = layer, path = path, construct_id = construct_id,
    year = year, scale = scale, stat = stat,
    is_climate = identical(group, "bioclim"),
    future_policy = future_policy,
    stringsAsFactors = FALSE
  )
}

glo_future_meta_one <- function(path) {
  rel <- sub("^.*/covariates/(reg|glo)/", "", normalizePath(path, winslash = "/"))
  seg <- strsplit(rel, "/", fixed = TRUE)[[1]]
  group <- seg[1]
  product <- seg[2]
  dataset <- paste(seg[1:2], collapse = "/")
  stem <- sub("[.]tif$", "", basename(path))
  layer_future <- gsub("-", "_", sub("^(reg|glo)_", "", stem))
  year <- suppressWarnings(as.integer(seg[grepl("^(19|20|21)[0-9]{2}$", seg)][1]))
  ssp <- seg[grepl("^ssp[0-9]+$", seg)][1]
  scale <- .glo_scale(layer_future)
  stat <- .glo_stat(layer_future)
  base <- .glo_strip_scale_stat(layer_future)
  if (identical(group, "bioclim")) {
    construct_id <- basename(dirname(path))
    model <- sub("^cnclim1km-", "", product)
    model <- gsub("-", "_", model)
    scale <- "base"
    stat <- "base"
  } else if (identical(group, "lulc")) {
    construct_id <- sub(
      "^lulc_(glu1km|hlu1km|clcd)_scenario_(19|20|21)[0-9]{2}_ssp[0-9]+_",
      "", base)
    construct_id <- paste0("lulc_", .glo_strip_year(construct_id))
    model <- NA_character_
  } else {
    stop("Future inventory only supports bioclim and lulc; got: ", path)
  }
  data.frame(
    group = group, dataset = dataset, product = product, model = model,
    construct_id = construct_id, layer_future = layer_future,
    year = year, ssp = ssp, scale = scale, stat = stat, path = path,
    stringsAsFactors = FALSE
  )
}

glo_scan_current <- function(glo_root) {
  specs <- list(
    list(dir = "bioclim/cnclim1km/present", group = "bioclim"),
    list(dir = "lulc/glu1km/2020", group = "lulc"),
    list(dir = "lulc/hlu1km/2020", group = "lulc"),
    list(dir = "lulc/clcd/2020", group = "lulc"),
    list(dir = "soil/csdlv2", group = "soil"),
    list(dir = "hydro/h90m", group = "hydro"),
    list(dir = "hydro/hydrosheds", group = "hydro"),
    list(dir = "dist", group = "dist"),
    list(dir = "topo/chinadem", group = "topo"),
    list(dir = "vege", group = "vege"),
    list(dir = "anthro/hfp/2020", group = "anthro"),
    list(dir = "anthro/smod/2020", group = "anthro")
  )
  rows <- list()
  for (s in specs) {
    files <- list.files(file.path(glo_root, s$dir), pattern = "[.]tif$",
                        recursive = TRUE, full.names = TRUE)
    files <- files[!grepl("[.]aux[.]xml$", files)]
    if (!length(files)) next
    files <- normalizePath(files, winslash = "/", mustWork = TRUE)
    rel <- sub(paste0("^", normalizePath(glo_root, winslash = "/"), "/"),
               "", files)
    for (i in seq_along(files)) {
      rows[[length(rows) + 1L]] <- glo_current_meta_one(files[i], rel[i],
                                                        dirname(glo_root))
    }
  }
  out <- do.call(rbind, rows)
  if (anyDuplicated(out$layer)) {
    stop("Current inventory has duplicate layer ids: ",
         paste(unique(out$layer[duplicated(out$layer)]), collapse = ", "))
  }
  out
}

glo_scan_future <- function(glo_root) {
  files <- list.files(file.path(glo_root, c("bioclim", "lulc")),
                      pattern = "[.]tif$", recursive = TRUE, full.names = TRUE)
  files <- files[grepl("/scenario/", files, fixed = TRUE)]
  files <- files[!grepl("[.]aux[.]xml$", files)]
  if (!length(files)) stop("No future scenario rasters found below ", glo_root)
  files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  out <- do.call(rbind, lapply(files, glo_future_meta_one))
  key <- paste(out$construct_id, out$product, out$scale, out$stat,
               out$year, out$ssp, sep = "|")
  if (anyDuplicated(key)) {
    dup <- unique(key[duplicated(key)])
    stop("Future inventory has duplicate construct/product/scale/stat/year/ssp keys: ",
         paste(utils::head(dup, 10), collapse = ", "))
  }
  out
}

glo_coverage <- function(paths, mask, verbose = TRUE) {
  n_land <- sum(terra::values(mask, mat = FALSE) == 1, na.rm = TRUE)
  if (!is.finite(n_land) || n_land <= 0) stop("The reference land mask has no valid cells.")
  out <- numeric(length(paths))
  for (i in seq_along(paths)) {
    r <- terra::rast(paths[i])
    n_fin <- terra::global(r, function(x) sum(is.finite(x), na.rm = TRUE))[[1]][1]
    out[i] <- n_fin / n_land
    if (isTRUE(verbose) && i %% 50L == 0L) {
      message("  coverage ", i, "/", length(paths))
    }
  }
  out
}

#' Build the current/future inventory for the reg or legacy glo library.
#'
#' @param root Workflow root containing `config/`.
#' @param data_root Project data root containing `covariates/reg` (or legacy `glo`) and
#'   `mask_china.tif`.
#' @param coverage_min Minimum finite share on the land mask for a current layer.
#' @param coverage_csv Optional CSV with `path`/`cover` columns to reuse; when
#'   missing, coverage is computed from the rasters.
#' @param check_geom Compare every current raster to the land mask.
#' @param verbose Print progress.
#' @return A list of output paths and summary counts.
glo_inventory_build <- function(root, data_root, coverage_min = 0.95,
                                coverage_csv = NULL, check_geom = TRUE,
                                verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to build the environmental inventory.")
  }
  glo_root <- file.path(data_root, "covariates", case_covariate_scope(data_root))
  mask_file <- file.path(data_root, "mask_china.tif")
  if (!file.exists(mask_file)) stop("Land mask not found: ", mask_file)
  files <- glo_inventory_files(root)
  dir.create(dirname(files$current_all), recursive = TRUE, showWarnings = FALSE)
  mask <- terra::rast(mask_file)

  if (isTRUE(verbose)) message("Scanning current layers ...")
  current <- glo_scan_current(glo_root)

  if (isTRUE(check_geom)) {
    if (isTRUE(verbose)) message("Geometry gate against ", basename(mask_file), " ...")
    ok <- vapply(current$path, function(p) {
      isTRUE(tryCatch(
        terra::compareGeom(mask, terra::rast(p), stopOnError = FALSE,
                           messages = FALSE),
        error = function(e) FALSE))
    }, logical(1))
    current$geom_ok <- ok
    if (any(!ok)) {
      message("  ", sum(!ok), " current layer(s) do not match the mask; excluded.")
    }
  } else {
    current$geom_ok <- TRUE
  }

  if (!is.null(coverage_csv) && file.exists(coverage_csv)) {
    if (isTRUE(verbose)) message("Reusing coverage from ", coverage_csv)
    cov <- utils::read.csv(coverage_csv, stringsAsFactors = FALSE)
    if (!all(c("path", "cover") %in% names(cov))) {
      stop("coverage_csv must have path and cover columns: ", coverage_csv)
    }
    current$cover <- cov$cover[match(
      normalizePath(current$path, winslash = "/"),
      normalizePath(cov$path, winslash = "/"))]
  } else {
    if (isTRUE(verbose)) message("Coverage gate (finite share on land mask >= ",
                                 coverage_min, ") ...")
    current$cover <- glo_coverage(current$path, mask, verbose = verbose)
  }
  current$status <- ifelse(current$geom_ok & is.finite(current$cover) &
                             current$cover >= coverage_min, "kept", "excluded")

  current$eligible_glu <- current$status == "kept" &
    current$future_policy %in% c("climate", "lulc", "static") &
    current$product != "hlu1km"
  current$eligible_hlu <- current$status == "kept" &
    current$future_policy %in% c("climate", "lulc", "static") &
    current$product != "glu1km"

  current <- current[, c("group", "dataset", "product", "layer", "path",
                         "construct_id", "year", "scale", "stat", "is_climate",
                         "future_policy", "cover", "status", "geom_ok",
                         "eligible_glu", "eligible_hlu")]
  current <- current[order(current$group, current$product, current$layer), ,
                     drop = FALSE]
  rownames(current) <- NULL
  utils::write.csv(current, files$current_all, row.names = FALSE)
  utils::write.csv(current[current$eligible_glu, , drop = FALSE],
                   files$eligible_glu, row.names = FALSE)
  utils::write.csv(current[current$eligible_hlu, , drop = FALSE],
                   files$eligible_hlu, row.names = FALSE)

  if (isTRUE(verbose)) message("Scanning future layers ...")
  future <- glo_scan_future(glo_root)
  future <- future[order(future$group, future$product, future$construct_id,
                         future$scale, future$year, future$ssp), , drop = FALSE]
  rownames(future) <- NULL
  utils::write.csv(future, files$future, row.names = FALSE)

  version <- data.frame(
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    data_root = normalizePath(data_root, winslash = "/", mustWork = TRUE),
    covariate_scope = basename(glo_root),
    n_current = nrow(current),
    n_current_kept = sum(current$status == "kept"),
    n_eligible_glu = sum(current$eligible_glu),
    n_eligible_hlu = sum(current$eligible_hlu),
    n_future = nrow(future),
    coverage_min = coverage_min,
    coverage_source = if (!is.null(coverage_csv) && file.exists(coverage_csv))
      normalizePath(coverage_csv, winslash = "/") else "computed",
    current_md5 = unname(tools::md5sum(files$current_all)),
    future_md5 = unname(tools::md5sum(files$future)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(version, files$version, row.names = FALSE)

  if (isTRUE(verbose)) {
    cat("Inventory written:\n")
    cat("  current: ", nrow(current), " (kept ", sum(current$status == "kept"),
        "; eligible GLU ", sum(current$eligible_glu),
        "; eligible HLU ", sum(current$eligible_hlu), ")\n", sep = "")
    cat("  future:  ", nrow(future), "\n", sep = "")
  }
  invisible(list(paths = files, current = current, future = future,
                 version = version))
}

glo_inventory_read_future <- function(cfg) {
  p <- cfg$future_layers_csv %||% ""
  if (!length(p) || is.na(p) || !nzchar(p) || !file.exists(p)) {
    stop("Future inventory not found: ", p,
         ". Run: Rscript 00_INVENTORY.R")
  }
  utils::read.csv(p, stringsAsFactors = FALSE)
}
