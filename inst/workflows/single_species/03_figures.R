# Complete publication figure collection using the package's established maps.
# Numerical diagnostic products are cached; render_only never predicts/refits.
case_figure_stacks <- function(out_dir, result, comparison, scenario) {
  fit <- result$fit_full %||% result$fit
  pool <- utils::read.csv(file.path(out_dir, "training_pool.csv"), stringsAsFactors = FALSE)
  pool <- pool[match(fit$env_vars, pool$layer), , drop = FALSE]
  if (anyNA(pool$path)) stop("The case inventory lacks fitted predictors.")
  scenarios <- utils::read.csv(file.path(out_dir, "scenario_table.csv"), stringsAsFactors = FALSE)
  castSDM:::.cast_explore_stacks(list(fit = fit, pool = pool, scenarios = scenarios,
    comparison = comparison), scenario)
}

case_figures <- function(comparison, out_dir, scenario, result = NULL, cfg = NULL,
                         current = NULL, future = NULL, attribution = NULL,
                         render_only = FALSE, include = "core", prepare_only = FALSE) {
  if (is.null(result)) result <- readRDS(file.path(out_dir, "training", "result.rds"))
  result$output_dir <- file.path(out_dir, "training")
  if (is.null(cfg)) {
    case_file <- file.path(out_dir, "case.rds")
    cfg <- if (file.exists(case_file)) readRDS(case_file)$cfg else list()
  }
  products <- utils::getFromNamespace(".cast_figure_products", "castSDM")(include)
  if (!render_only && any(products %in% c("attribution", "mess")) &&
      (is.null(current) || is.null(future))) {
    stacks <- case_figure_stacks(out_dir, result, comparison, scenario)
    current <- stacks$current; future <- stacks$future
  }
  figures <- castSDM::cast_figures(result, comparison,
    output_dir = file.path(out_dir, "figures"), scenario = scenario,
    current = current, future = future, attribution = attribution,
    basemap = "china", include = include, render_only = render_only,
    prepare_only = prepare_only,
    response_method = cfg$figure_response_method %||% "ale",
    response_grid = cfg$figure_response_grid %||% 50L,
    display_resolution = cfg$figure_display_resolution %||% .02,
    width_mm = 183, dpi = 600,
    resources = cfg$resources %||% list(cores = 4, memory_mb = 2048))
  invisible(figures)
}
