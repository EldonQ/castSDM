# Complete publication figure collection using the package's established maps.
# Numerical diagnostic products are cached; render_only never predicts/refits.
case_figure_stacks <- function(out_dir, result, scenario) {
  vars <- (result$fit_full %||% result$fit)$env_vars
  pool <- utils::read.csv(file.path(out_dir, "training_pool.csv"), stringsAsFactors = FALSE)
  paths <- pool$path[match(vars, pool$layer)]
  if (anyNA(paths) || any(!file.exists(paths))) stop("Current figure predictor sources are missing; use render-only with a complete figure-data cache.")
  current <- terra::rast(paths); names(current) <- vars
  table <- utils::read.csv(file.path(out_dir, "scenario_table.csv"), stringsAsFactors = FALSE)
  table <- table[table$scenario == scenario, , drop = FALSE]
  if (!nrow(table)) stop("No saved predictor mapping for reporting scenario: ", scenario)
  future_paths <- paths
  for (i in seq_len(nrow(table))) {
    hit <- match(table$layer[i], vars)
    if (!is.na(hit)) future_paths[hit] <- table$path[i]
  }
  if (any(!file.exists(future_paths))) stop("Future figure predictor sources are missing; use render-only with a complete figure-data cache.")
  future <- terra::rast(future_paths); names(future) <- vars
  list(current = current, future = future)
}

case_figures <- function(comparison, out_dir, scenario, result = NULL, cfg = NULL,
                         current = NULL, future = NULL, attribution = NULL,
                         render_only = FALSE, include = "all") {
  if (is.null(result)) result <- readRDS(file.path(out_dir, "training", "result.rds"))
  result$output_dir <- file.path(out_dir, "training")
  if (is.null(cfg)) {
    case_file <- file.path(out_dir, "case.rds")
    cfg <- if (file.exists(case_file)) readRDS(case_file)$cfg else list()
  }
  if (!render_only && (is.null(current) || is.null(future))) {
    stacks <- case_figure_stacks(out_dir, result, scenario)
    current <- stacks$current; future <- stacks$future
  }
  figures <- castSDM::cast_figures(result, comparison,
    output_dir = file.path(out_dir, "figures"), scenario = scenario,
    current = current, future = future, attribution = attribution,
    basemap = "china", include = include, render_only = render_only,
    response_grid = cfg$figure_response_grid %||% 50L,
    display_resolution = cfg$figure_display_resolution %||% .02,
    width_mm = 183, dpi = 600,
    resources = cfg$resources %||% list(cores = 4, memory_mb = 2048))
  # Remove only the obsolete generic previews, after successful replacement.
  if ("02_habitat_change" %in% figures$status$figure &&
      all(figures$status$status[figures$status$figure == "02_habitat_change"] == "written_not_visually_checked")) {
    for (name in c("scenario_comparison.png", "scenario_comparison.pdf",
                   "spatial_validation.png", "spatial_validation.pdf")) {
      path <- file.path(out_dir, "figures", name)
      if (file.exists(path)) unlink(path)
    }
  }
  invisible(figures)
}
