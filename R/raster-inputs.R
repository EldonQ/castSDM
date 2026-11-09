# AI-assisted maintenance (2026-11-09); scope and verification: inst/AI_ASSISTANCE.md.
.cast_assert_raster_sources <- function(x, label) {
  if (is.null(x) || !inherits(x, "SpatRaster")) return(invisible(TRUE))
  paths <- unique(as.character(terra::sources(x)))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  # GDAL virtual/subdataset connections are not ordinary filesystem paths.
  virtual <- grepl("^/vsi|^[[:alpha:]][[:alnum:]+.-]*:", paths) &
    !grepl("^[[:alpha:]]:/", chartr("\\", "/", paths))
  missing <- paths[!virtual & !file.exists(paths)]
  if (length(missing)) cli::cli_abort(c(
    "Missing file-backed raster source for {label}.",
    i = "{length(missing)} source file{?s} unavailable; first: {.path {missing[1L]}}.",
    i = "Update input paths after moving data. Keep source files available until comparison finishes."))
  invisible(TRUE)
}

.cast_assert_compare_sources <- function(plan) {
  if (!plan$raster) return(invisible(TRUE))
  .cast_assert_raster_sources(plan$current, "current")
  for (scenario in names(plan$future))
    .cast_assert_raster_sources(plan$future[[scenario]], paste("future scenario", scenario))
  .cast_assert_raster_sources(plan$mask, "mask")
  invisible(TRUE)
}
