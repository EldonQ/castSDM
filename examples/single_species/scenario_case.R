# Dataset-specific mapping only; cast_compare constructs the four corners.
case_pool_metadata <- function(cfg, vars) {
  pool <- nsdm_pool_layers(cfg$nsdm_pool_csv)
  hit <- match(vars, pool$layer)
  if (anyNA(hit)) stop("Prepared predictors absent from pool: ", paste(vars[is.na(hit)], collapse = ", "))
  pool[hit, , drop = FALSE]
}

case_scenarios <- function(cfg, selected, current, out_dir) {
  meta <- case_pool_metadata(cfg, selected)
  groups <- list(climate = selected[meta$future_policy == "climate"],
                 lulc = selected[meta$future_policy == "lulc"])
  if (any(!lengths(groups))) stop("The selected model must contain both climate and LULC groups.")
  if (any(!meta$future_policy %in% c("climate", "lulc", "static")))
    stop("Selected pool includes a predictor without a declared future policy.")
  future <- glo_inventory_read_future(cfg)
  required <- c("construct_id", "product", "model", "year", "ssp", "scale", "stat", "path")
  if (!all(required %in% names(future))) stop("Future inventory has missing metadata fields.")
  rows <- list()
  for (ssp in cfg$scenario_ssps) for (period in cfg$scenario_periods) {
    climate_year <- cfg$scenario_climate_years[[period]]
    lulc_year <- cfg$scenario_lulc_years[[period]]
    if (is.null(climate_year) || is.null(lulc_year) || anyNA(c(climate_year, lulc_year)))
      stop("Missing year mapping for ", period)
    scenario <- paste("joint", cfg$scenario_climate_model, ssp, period, sep = "_")
    for (i in which(meta$future_policy != "static")) {
      m <- meta[i, ]; climate <- m$future_policy == "climate"
      hit <- future$construct_id == m$construct_id & future$ssp == ssp &
        future$scale == m$scale & future$stat == m$stat &
        future$year == if (climate) climate_year else lulc_year
      hit <- hit & if (climate) future$model == cfg$scenario_climate_model else
        future$product == paste0(cfg$scenario_lulc_product, "1km")
      match <- future[which(hit %in% TRUE), , drop = FALSE]
      if (nrow(match) != 1L) stop("Expected one future file for ", m$layer, " / ", scenario,
                                  "; found ", nrow(match), ". Rebuild the inventory if necessary.")
      rows[[length(rows)+1L]] <- data.frame(scenario = scenario, family = "joint",
        climate_model = cfg$scenario_climate_model, ssp = ssp, period = period,
        climate_year = climate_year, lulc_year = lulc_year, layer = m$layer,
        path = match$path, source = if (climate) "future_climate" else "future_lulc")
    }
  }
  table <- do.call(rbind, rows)
  path <- file.path(out_dir, "scenario_table.csv")
  utils::write.csv(table, path, row.names = FALSE)
  stacks <- nsdm_build_scenario_stacks(current, selected, path)
  report <- paste("joint", cfg$scenario_climate_model, cfg$scenario_report_ssp,
                  cfg$scenario_report_period, sep = "_")
  if (!report %in% names(stacks)) stop("Report scenario was not requested: ", report)
  list(future = stacks, groups = groups, report = report, table = table)
}
