# CHELSA current versus ensemble-mean future climate; other predictors static.
case_pool_metadata <- function(cfg, vars) {
  pool <- nsdm_pool_layers(cfg$nsdm_pool_csv)
  hit <- match(vars, pool$layer)
  if (anyNA(hit)) stop("Prepared predictors absent from pool: ", paste(vars[is.na(hit)], collapse = ", "))
  pool[hit, , drop = FALSE]
}
case_scenarios <- function(cfg, selected, current, out_dir) {
  meta <- case_pool_metadata(cfg, selected)
  groups <- list(climate = selected[meta$future_policy == "climate"])
  if (!length(groups$climate)) stop("At least one climate predictor is required for climate projections.")
  if (any(!meta$future_policy %in% c("climate", "static")))
    stop("Every nonclimate predictor must have a declared static future policy.")
  future <- case_inventory_read_future(cfg)
  rows <- list()
  for (ssp in cfg$scenario_ssps) for (period in cfg$scenario_periods) {
    scenario <- paste("climate", "ensmean", ssp, period, sep = "_")
    for (i in which(meta$future_policy == "climate")) {
      m <- meta[i, ]
      hit <- future$construct_id == m$construct_id & future$ssp == ssp &
        future$period == period & future$product == "chelsa21-ensmean"
      matched <- future[which(hit %in% TRUE), , drop = FALSE]
      if (nrow(matched) != 1L) stop("Expected one CHELSA ensemble layer: ", m$layer, " / ", scenario)
      rows[[length(rows) + 1L]] <- data.frame(scenario = scenario, family = "climate",
        climate_model = "ensmean", ssp = ssp, period = period, baseline_period = "1981_2010",
        climate_year = matched$year, layer = m$layer, path = matched$path, source = "future_climate")
    }
  }
  table <- do.call(rbind, rows)
  path <- file.path(out_dir, "scenario_table.csv")
  utils::write.csv(table, path, row.names = FALSE)
  utils::write.csv(data.frame(layer = selected, policy = meta$future_policy,
    current_path = meta$path, baseline_period = "1981_2010",
    future_assumption = ifelse(meta$future_policy == "climate", "CHELSA ensemble-mean climate", "held at current value")),
    file.path(out_dir, "projection_assumptions.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    members = "GFDL-ESM4; IPSL-CM6A-LR; MPI-ESM1-2-HR; MRI-ESM2-0; UKESM1-0-LL",
    interpretation = "Prediction at ensemble-mean climate is not the mean of GCM predictions.",
    missingness = "Delivered ensemble uses available-member means; per-cell member counts are not provided.",
    uncertainty = "Algorithm disagreement excludes between-GCM uncertainty; no GCM interval is estimated."),
    file.path(out_dir, "climate_ensemble.csv"), row.names = FALSE)
  stacks <- nsdm_build_scenario_stacks(current, selected, path)
  report <- paste("climate", "ensmean", cfg$scenario_report_ssp, cfg$scenario_report_period, sep = "_")
  if (!report %in% names(stacks)) stop("Report scenario was not requested: ", report)
  list(future = stacks, groups = groups, report = report, table = table)
}
