# CHELSA v2.1 current climate and ensemble-mean future climate inventory.
# Land use, soil, topography and hydrology remain current in every projection.
# Other dynamic products are inventoried but excluded from the projection pool.
# Definitions and units come from the delivered library; no rescaling is done.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}

case_inventory_files <- function(root) {
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
  if (!dir.exists(file.path(data_root, "covariates", "reg")))
    stop("Expected the national environmental library below covariates/reg: ", data_root, call. = FALSE)
  "reg"
}

.case_scale <- function(layer) {
  if (grepl("_f1000$", layer)) return("f1000")
  if (grepl("_f5000$", layer)) return("f5000")
  if (grepl("_f10000$", layer)) return("f10000")
  if (grepl("_f20000$", layer)) return("f20000")
  "base"
}

.case_stat <- function(layer) {
  if (grepl("_max$", layer)) return("max")
  if (grepl("_min$", layer)) return("min")
  "base"
}

.case_strip_scale_stat <- function(x) {
  x <- sub("_f(1000|5000|10000|20000)$", "", x)
  x <- sub("_(max|min)$", "", x)
  x
}

.case_strip_year <- function(x) {
  gsub("_(19|20|21)[0-9]{2}(?=_|$)", "", x, perl = TRUE)
}

# One current layer. `rel` is the path below covariates/reg/ only.
case_current_meta_one <- function(path, rel, data_root) {
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
  scale <- .case_scale(layer)
  stat <- .case_stat(layer)
  base <- .case_strip_scale_stat(layer)
  construct_id <- if (identical(group, "bioclim")) {
    layer
  } else if (identical(group, "lulc")) {
    base <- sub("^lulc_(glu1km|hlu1km|clcd)_(19|20|21)[0-9]{2}_", "", base)
    paste0("lulc_", .case_strip_year(base))
  } else if (identical(group, "soil")) {
    paste0("soil_", sub("^soil_csdlv2_", "", base))
  } else if (identical(group, "hydro")) {
    paste0("hydro_", sub("^hydro_", "", base))
  } else if (identical(group, "dist")) {
    paste0("dist_", sub("^dist_", "", base))
  } else if (identical(group, "topo")) {
    paste0("topo_", sub("^topo_chinadem_", "", base))
  } else if (identical(group, "vege")) {
    paste0("vege_", .case_strip_year(sub("^vege_", "", base)))
  } else if (identical(group, "anthro")) {
    paste0("anthro_", product)
  } else {
    paste0(group, "_", base)
  }
  future_policy <- if (identical(group, "bioclim")) {
    "climate"
  } else if (identical(group, "lulc") && product %in% c("glu1km", "hlu1km")) {
    "static"
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

case_future_meta_one <- function(path) {
  rel <- sub("^.*/covariates/(reg|glo)/", "", normalizePath(path, winslash = "/"))
  seg <- strsplit(rel, "/", fixed = TRUE)[[1]]
  group <- seg[1]
  product <- seg[2]
  dataset <- paste(seg[1:2], collapse = "/")
  stem <- sub("[.]tif$", "", basename(path))
  layer_future <- gsub("-", "_", sub("^(reg|glo)_", "", stem))
  year <- suppressWarnings(as.integer(seg[grepl("^(19|20|21)[0-9]{2}$", seg)][1]))
  ssp <- seg[grepl("^ssp[0-9]+$", seg)][1]
  scale <- .case_scale(layer_future)
  stat <- .case_stat(layer_future)
  base <- .case_strip_scale_stat(layer_future)
  if (identical(group, "bioclim")) {
    construct_id <- basename(dirname(path))
    if (!identical(product, "chelsa21-ensmean")) stop("Only CHELSA ensemble-mean climate scenarios are supported here.")
    model <- "ensmean"
    scale <- "base"
    stat <- "base"
  } else stop("Future predictors must be CHELSA ensemble climate layers: ", path)
  data.frame(
    group = group, dataset = dataset, product = product, model = model,
    construct_id = construct_id, layer_future = layer_future,
    year = year, ssp = ssp, scale = scale, stat = stat, path = path,
    stringsAsFactors = FALSE
  )
}

case_scan_current <- function(cov_root) {
  specs <- list(
    list(dir = "bioclim/chelsa21/present", group = "bioclim"),
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
    files <- list.files(file.path(cov_root, s$dir), pattern = "[.]tif$",
                        recursive = TRUE, full.names = TRUE)
    files <- files[!grepl("[.]aux[.]xml$", files)]
    if (!length(files)) next
    files <- normalizePath(files, winslash = "/", mustWork = TRUE)
    rel <- sub(paste0("^", normalizePath(cov_root, winslash = "/"), "/"),
               "", files)
    for (i in seq_along(files)) {
      rows[[length(rows) + 1L]] <- case_current_meta_one(files[i], rel[i],
                                                        dirname(cov_root))
    }
  }
  out <- do.call(rbind, rows)
  if (anyDuplicated(out$layer)) {
    stop("Current inventory has duplicate layer ids: ",
         paste(unique(out$layer[duplicated(out$layer)]), collapse = ", "))
  }
  out
}

case_scan_future <- function(cov_root) {
  files <- list.files(file.path(cov_root, "bioclim", "chelsa21-ensmean"),
                      pattern = "[.]tif$", recursive = TRUE, full.names = TRUE)
  files <- files[grepl("/scenario/", files, fixed = TRUE)]
  files <- files[!grepl("[.]aux[.]xml$", files)]
  if (!length(files)) stop("No future scenario rasters found below ", cov_root)
  files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  out <- do.call(rbind, lapply(files, case_future_meta_one))
  key <- paste(out$construct_id, out$product, out$scale, out$stat,
               out$year, out$ssp, sep = "|")
  if (anyDuplicated(key)) {
    dup <- unique(key[duplicated(key)])
    stop("Future inventory has duplicate construct/product/scale/stat/year/ssp keys: ",
         paste(utils::head(dup, 10), collapse = ", "))
  }
  out
}

case_coverage <- function(paths, mask, verbose = TRUE) {
  n_land <- terra::global(mask == 1, "sum", na.rm = TRUE)[[1L]][1L]
  if (!is.finite(n_land) || n_land <= 0) stop("The reference land mask has no valid cells.")
  out <- numeric(length(paths))
  for (i in seq_along(paths)) {
    r <- terra::rast(paths[i])
    n_fin <- terra::global(r, "notNA")[[1L]][1L]
    out[i] <- n_fin / n_land
    if (isTRUE(verbose) && i %% 50L == 0L) {
      message("  coverage ", i, "/", length(paths))
    }
  }
  out
}

#' Build the current/future inventory for the CHELSA national library.
#'
#' @param root Workflow root containing `config/`.
#' @param data_root Project data root containing `covariates/reg` and
#'   `mask_china.tif`.
#' @param coverage_min Minimum finite share on the land mask for a current layer.
#' @param coverage_csv Optional CSV with `path`/`cover` columns to reuse; when
#'   missing, coverage is computed from the rasters.
#' @param check_geom Compare every current raster to the land mask.
#' @param verbose Print progress.
#' @return A list of output paths and summary counts.
case_inventory_build <- function(root, data_root, coverage_min = 0.95,
                                coverage_csv = NULL, check_geom = TRUE,
                                verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to build the environmental inventory.")
  }
  cov_root <- file.path(data_root, "covariates", case_covariate_scope(data_root))
  mask_file <- file.path(data_root, "mask_china.tif")
  if (!file.exists(mask_file)) stop("Land mask not found: ", mask_file)
  files <- case_inventory_files(root)
  dir.create(dirname(files$current_all), recursive = TRUE, showWarnings = FALSE)
  mask <- terra::rast(mask_file)

  if (isTRUE(verbose)) message("Scanning current layers ...")
  current <- case_scan_current(cov_root)

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
    dictionary_path <- file.path(dirname(data_root), "00_docs", "variable_dictionary.csv")
    current$cover <- NA_real_
    if (file.exists(dictionary_path)) {
      dictionary <- utils::read.csv(dictionary_path, stringsAsFactors = FALSE, encoding = "UTF-8")
      rel <- sub(paste0(normalizePath(file.path(data_root, "covariates"), winslash = "/"), "/"), "", current$path, fixed = TRUE)
      if (all(c("relpath", "n_valid") %in% names(dictionary))) {
        n_land <- terra::global(mask == 1, "sum", na.rm = TRUE)[[1L]][1L]
        current$cover <- dictionary$n_valid[match(rel, dictionary$relpath)]/n_land
      }
    }
    missing_coverage <- which(!is.finite(current$cover))
    if (length(missing_coverage)) current$cover[missing_coverage] <-
      case_coverage(current$path[missing_coverage], mask, verbose = verbose)
  }
  current$status <- ifelse(current$geom_ok & is.finite(current$cover) &
                             current$cover >= coverage_min, "kept", "excluded")

  current$eligible_glu <- current$status == "kept" &
    current$future_policy %in% c("climate", "static") &
    !(current$group == "lulc" & current$product != "glu1km")
  current$eligible_hlu <- current$status == "kept" &
    current$future_policy %in% c("climate", "static") &
    !(current$group == "lulc" & current$product != "hlu1km")

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
  future <- case_scan_future(cov_root)
  if (isTRUE(check_geom)) {
    aligned <- vapply(future$path, function(p) isTRUE(terra::compareGeom(
      mask, terra::rast(p), stopOnError = FALSE, messages = FALSE)), logical(1))
    if (any(!aligned)) stop("Future climate rasters must match the current grid: ", future$path[which(!aligned)[1L]])
  }
  period_map <- c(`2025` = "2011_2040", `2055` = "2041_2070", `2085` = "2071_2100")
  future$period <- unname(period_map[as.character(future$year)])
  if (anyNA(future$period) || any(!future$ssp %in% c("ssp126", "ssp370", "ssp585"))) stop("Unexpected CHELSA scenario period or SSP; inspect the source inventory.")
  future$baseline_period <- "1981_2010"
  future <- future[order(future$group, future$product, future$construct_id,
                         future$scale, future$year, future$ssp), , drop = FALSE]
  rownames(future) <- NULL
  utils::write.csv(future, files$future, row.names = FALSE)
  expected <- as.vector(outer(c("ssp126", "ssp370", "ssp585"), c(2025, 2055, 2085), paste))
  for (v in current$construct_id[current$future_policy == "climate"]) {
    present <- future[future$construct_id == v, , drop = FALSE]
    if (!setequal(paste(present$ssp, present$year), expected))
      stop("Incomplete future-climate coverage for ", v, "; expected all nine real SSP/period combinations.")
  }
  case_inventory_metadata(root, data_root, current)

  version <- data.frame(
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    data_root = normalizePath(data_root, winslash = "/", mustWork = TRUE),
    covariate_scope = basename(cov_root),
    n_current = nrow(current),
    n_current_kept = sum(current$status == "kept"),
    n_eligible_glu = sum(current$eligible_glu),
    n_eligible_hlu = sum(current$eligible_hlu),
    n_future = nrow(future),
    coverage_min = coverage_min,
    coverage_source = if (!is.null(coverage_csv) && file.exists(coverage_csv))
      normalizePath(coverage_csv, winslash = "/") else "delivered variable dictionary; raster fallback for missing entries",
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

case_inventory_read_future <- function(cfg) {
  p <- cfg$future_layers_csv %||% ""
  if (!length(p) || is.na(p) || !nzchar(p) || !file.exists(p)) {
    stop("Future inventory not found: ", p,
         ". Run: Rscript 00_INVENTORY.R")
  }
  utils::read.csv(p, stringsAsFactors = FALSE)
}

# Source definitions are copied into each case, never inferred from BIO ids.
case_inventory_metadata <- function(root, data_root, current) {
  provenance <- file.path(dirname(data_root), "00_docs", "provenance_table.csv")
  if (!file.exists(provenance)) stop("Source provenance table is required: ", provenance)
  cols <- c("layer_path", "variable_definition", "units", "source_citation", "notes")
  info <- if (requireNamespace("data.table", quietly = TRUE))
    data.table::fread(provenance, select = cols, data.table = FALSE, encoding = "UTF-8") else
    utils::read.csv(provenance, stringsAsFactors = FALSE, encoding = "UTF-8")[, cols]
  rel <- substring(current$path, nchar(normalizePath(file.path(data_root, "covariates"), winslash = "/")) + 2L)
  if (anyDuplicated(info$layer_path)) stop("Duplicate provenance layer paths.")
  at <- match(rel, info$layer_path)
  if (anyNA(at)) stop("Missing provenance for: ", paste(current$layer[is.na(at)], collapse = ", "))
  metadata <- data.frame(layer = current$layer, product = current$product,
    definition = info$variable_definition[at], units = info$units[at],
    source = info$source_citation[at], notes = info$notes[at], label = current$layer)
  climate_labels <- c("Annual mean temperature", "Mean diurnal range", "Isothermality",
    "Temperature seasonality", "Warmest-month maximum temperature", "Coldest-month minimum temperature",
    "Annual temperature range", "Wettest-quarter mean temperature", "Driest-quarter mean temperature",
    "Warmest-quarter mean temperature", "Coldest-quarter mean temperature", "Annual precipitation",
    "Wettest-month precipitation", "Driest-month precipitation", "Precipitation seasonality",
    "Wettest-quarter precipitation", "Driest-quarter precipitation", "Warmest-quarter precipitation",
    "Coldest-quarter precipitation", "Growing degree days above 0 C", "Growing degree days above 5 C", "Growing degree days above 10 C")
  labels <- stats::setNames(climate_labels, c(sprintf("bio%02d", 1:19), "gdd0", "gdd5", "gdd10"))
  hit <- metadata$product == "chelsa21" & metadata$layer %in% names(labels)
  metadata$label[hit] <- paste0(metadata$layer[hit], " | ", labels[metadata$layer[hit]])
  utils::write.csv(metadata, file.path(root, "config", "predictor_metadata.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

case_inventory_signature <- function(data_root) {
  cov <- file.path(data_root, "covariates", case_covariate_scope(data_root))
  current <- case_scan_current(cov)$path
  future <- case_scan_future(cov)$path
  paths <- c(current, future, file.path(data_root, "mask_china.tif"),
    file.path(dirname(data_root), "00_docs", c("variable_dictionary.csv", "provenance_table.csv")))
  info <- file.info(paths)
  if (anyNA(info$size)) stop("Missing input when indexing the environmental library.")
  list(schema = "chelsa_inventory/1", root = normalizePath(data_root, winslash = "/"),
    paths = paths, size = info$size, modified = as.numeric(info$mtime))
}

case_inventory_ensure <- function(root, data_root) {
  signature <- case_inventory_signature(data_root)
  receipt <- file.path(root, "config", "inventory_signature.rds")
  files <- unlist(case_inventory_files(root), use.names = FALSE)
  files <- c(files, file.path(root, "config", "predictor_metadata.csv"))
  signature$script <- unname(tools::md5sum(file.path(root, "inventory.R")))
  previous <- if (file.exists(receipt)) readRDS(receipt) else NULL
  if (all(file.exists(files)) && identical(previous$inputs, signature) &&
      identical(previous$outputs, unname(tools::md5sum(files)))) return(invisible(FALSE))
  case_inventory_build(root, data_root)
  saveRDS(list(inputs = signature, outputs = unname(tools::md5sum(files))), receipt)
  invisible(TRUE)
}
