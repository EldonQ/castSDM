# Known-response environmental scenario benchmark. Independent datasets are units.
ROOT <- "E:/Package/castSDM_workflows/core_method_validation"
PACKAGE <- "E:/Package/cast"
source_files <- function() c(list.files(file.path(PACKAGE, "R"), "[.]R$", full.names = TRUE),
  file.path(PACKAGE, c("DESCRIPTION", "NAMESPACE")),
  file.path(ROOT, c("00_RUN_VALIDATION.R", "01_attribution_validation.R",
    "02_attribution_figures.R", "RESEARCH_PLAN_20261015.md")))
source_manifest <- function() data.frame(file = source_files(),
  sha256 = vapply(source_files(), function(f) digest::digest(file = f, algo = "sha256"), character(1)))

truth <- function(X, scenario) {
  eta <- switch(scenario,
    linear = -.3 + 1.2 * X$temperature + .6 * X$moisture + .4 * X$terrain,
    nonlinear = 1.4 - .8 * X$temperature^2 + .6 * X$moisture + .4 * X$terrain,
    joint = -.2 + 1.1 * X$temperature * X$moisture + .4 * X$terrain,
    null = -.2 + .8 * X$moisture + .4 * X$terrain,
    correlated = -.3 + 1.2 * X$temperature + .3 * X$moisture + .4 * X$terrain,
    extrapolation = 1.4 - .8 * X$temperature^2 + .6 * X$moisture + .4 * X$terrain)
  plogis(eta)
}
make_data <- function(seed, scenario, n = 1500L) {
  set.seed(seed)
  d <- data.frame(lon = runif(n, -1, 1), lat = runif(n, -1, 1),
    temperature = runif(n, -2, 2), moisture = runif(n, -2, 2), terrain = runif(n, -1, 1))
  if (scenario == "correlated") d$moisture <- .9 * d$temperature + rnorm(n, sd = .22)
  d$presence <- rbinom(n, 1, truth(d, scenario))
  d
}
make_grid <- function(scenario) {
  r <- terra::rast(nrows = 40, ncols = 40, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  a <- data.frame(lon = xy[, 1], lat = xy[, 2], temperature = 1.6 * xy[, 1],
    moisture = 1.6 * xy[, 2], terrain = sin(pi * xy[, 1]) * cos(pi * xy[, 2]))
  if (scenario == "correlated") a$moisture <- .9 * a$temperature + .18 * sin(pi * a$lat)
  b <- a; b$temperature <- a$temperature + .35 + .7 * a$lat
  if (scenario == "extrapolation") b$temperature <- a$temperature + 2.2 + .8 * a$lon
  b$moisture <- a$moisture + .25 * a$lon
  list(current = a, future = b, raster = r)
}
pred <- function(fit, X, model = names(fit$models)[1]) {
  castSDM:::predict_single_model(fit$models[[model]], X[, fit$env_vars, drop = FALSE])
}
global_effect <- function(fit, train, a, b, driver) {
  # Marginal response averaged over all training rows, not a median profile.
  values <- seq(min(c(a[[driver]], b[[driver]])), max(c(a[[driver]], b[[driver]])), length.out = 65)
  response <- vapply(values, function(v) {
    d <- train; d[[driver]] <- v
    mean(pred(fit, d))
  }, numeric(1))
  f <- stats::approxfun(values, response, rule = 2)
  f(b[[driver]]) - f(a[[driver]])
}
metric <- function(estimate, target, keep = rep(TRUE, length(target))) {
  k <- keep & is.finite(estimate) & is.finite(target)
  sign_k <- k & abs(target) > .01
  h <- which(k)
  n_hot <- if (length(h)) max(1L, ceiling(length(h) * .2)) else 0L
  truth_hot <- h[order(abs(target[h]), decreasing = TRUE)][seq_len(n_hot)]
  est_hot <- h[order(abs(estimate[h]), decreasing = TRUE)][seq_len(n_hot)]
  c(n = sum(k), rmse = if (any(k)) sqrt(mean((estimate[k] - target[k])^2)) else NA_real_,
    bias = if (any(k)) mean(estimate[k] - target[k]) else NA_real_,
    sign_agreement = if (any(sign_k)) mean(sign(estimate[sign_k]) == sign(target[sign_k])) else NA_real_,
    hotspot_recall = if (n_hot) length(intersect(truth_hot, est_hot)) / n_hot else NA_real_)
}
similarity <- function(a, b) {
  k <- is.finite(a) & is.finite(b)
  c(map_correlation = if (sum(k) > 2 && sd(a[k]) > 0 && sd(b[k]) > 0) cor(a[k], b[k]) else NA_real_,
    map_rmse = if (any(k)) sqrt(mean((a[k] - b[k])^2)) else NA_real_,
    sign_agreement = if (any(k & abs(a) > .01)) mean(sign(a[k & abs(a) > .01]) == sign(b[k & abs(a) > .01])) else NA_real_)
}
run_validation <- function(mode = "production", scenario = NULL) {
  pkgload::load_all(PACKAGE, quiet = TRUE)
  manifest <- source_manifest()
  scenarios <- c("linear", "nonlinear", "joint", "null", "correlated", "extrapolation")
  if (!is.null(scenario)) scenarios <- intersect(scenarios, scenario)
  if (!length(scenarios)) stop("Unknown benchmark scenario")
  reps <- if (mode == "smoke") 1L else 20L
  output <- file.path(ROOT, "attribution_output", paste0(mode, "_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
    if (!is.null(scenario)) paste0("_", scenario) else ""))
  dir.create(output, recursive = TRUE)
  dir.create(file.path(output, "datasets")); dir.create(file.path(output, "maps"))
  if (is.null(scenario)) writeLines(output, file.path(ROOT, paste0("latest_attribution_", mode, "_running.txt")))
  utils::write.csv(manifest, file.path(output, "source_manifest.csv"), row.names = FALSE)
  metrics <- reliability <- checks <- errors <- list()
  add_check <- function(scenario, rep, name, error, tolerance = 1e-10) {
    checks[[length(checks) + 1L]] <<- data.frame(scenario, replicate = rep, check = name,
      max_abs_error = error, tolerance, passed = is.finite(error) && error <= tolerance)
  }
  export_progress <- function() {
    bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
    utils::write.csv(bind(metrics), file.path(output, "advantage_metrics.csv"), row.names = FALSE)
    utils::write.csv(bind(reliability), file.path(output, "reliability_metrics.csv"), row.names = FALSE)
    utils::write.csv(bind(checks), file.path(output, "numerical_checks.csv"), row.names = FALSE)
    utils::write.csv(if (length(errors)) bind(errors) else data.frame(scenario = character(), replicate = integer(), error = character()),
      file.path(output, "errors.csv"), row.names = FALSE)
  }
  for (s in scenarios) for (rep in seq_len(reps)) {
    id <- paste(s, sprintf("%02d", rep), sep = "_")
    message("VALIDATION ", id)
    tryCatch({
      scenario_id <- match(s, c("linear", "nonlinear", "joint", "null", "correlated", "extrapolation"))
      train <- make_data(12000 + scenario_id * 100 + rep, s)
      train2 <- make_data(22000 + scenario_id * 100 + rep, s)
      grid <- make_grid(s); a <- grid$current; b <- grid$future
      saveRDS(list(train = train, resample = train2, current = a, future = b,
        protocol = list(scenario = s, replicate = rep, n = 1500L, grid = c(40L, 40L))),
        file.path(output, "datasets", paste0(id, ".rds")))
      f <- cast_fit(train, models = c("rf", "gam"), rf_ntree = 500, seed = 100 + rep, verbose = FALSE)
      fr <- cast_fit(train2, models = "rf", rf_ntree = 500, seed = 100 + rep, verbose = FALSE)
      fs <- cast_fit(train, models = "rf", rf_ntree = 500, seed = 1100 + rep, verbose = FALSE)
      only_truth <- a; only_truth$temperature <- b$temperature
      held_truth <- b; held_truth$temperature <- a$temperature
      target <- truth(only_truth, s) - truth(a, s)
      retained_truth <- truth(b, s) - truth(held_truth, s)
      total_truth <- truth(b, s) - truth(a, s)
      maps <- list()
      for (engine in c("rf", "gam")) {
        z <- cast_attribution(f, a, b, models = engine, verbose = FALSE)
        native <- z$data$effect_temperature
        fg <- f; fg$models <- f$models[engine]
        global <- global_effect(fg, train, a, b, "temperature")
        mean_future <- a; mean_future$temperature <- a$temperature + mean(b$temperature - a$temperature)
        uniform <- pred(f, mean_future, engine) - pred(f, a, engine)
        supported <- z$data$support_effect_temperature == 1
        for (method in c("local_scenario", "global_response", "mean_shift")) {
          estimate <- switch(method, local_scenario = native, global_response = global, mean_shift = uniform)
          for (domain in c("all", "within_range", "outside_range")) {
            keep <- switch(domain, all = rep(TRUE, length(target)), within_range = supported, outside_range = !supported)
            m <- metric(estimate, target, keep)
            metrics[[length(metrics) + 1L]] <- data.frame(scenario = s, replicate = rep,
              engine, method, domain, t(m), support = mean(supported),
              predictive_rmse = sqrt(mean((pred(f, a, engine) - truth(a, s))^2)), row.names = NULL)
          }
        }
        for (kind in c("retained", "total")) {
          estimate <- if (kind == "retained") z$data$retained_temperature else z$data$total_change
          m <- metric(estimate, if (kind == "retained") retained_truth else total_truth)
          metrics[[length(metrics) + 1L]] <- data.frame(scenario = s, replicate = rep,
            engine, method = kind, domain = "all", t(m), support = mean(supported),
            predictive_rmse = sqrt(mean((pred(f, a, engine) - truth(a, s))^2)), row.names = NULL)
        }
        maps[[engine]] <- z
        if (rep == 1L) {
          display <- data.frame(a[, c("lon", "lat")], truth = target, local_scenario = native,
            global_response = global, mean_shift = uniform, error = native - target,
            retained_truth, total_truth, z$data[, setdiff(names(z$data), c("lon", "lat"))])
          utils::write.csv(display, file.path(output, "maps", paste0(s, "_", engine, ".csv")), row.names = FALSE)
        }
      }
      z <- maps$rf
      zr <- cast_attribution(fr, a, b, verbose = FALSE)
      zs <- cast_attribution(fs, a, b, verbose = FALSE)
      for (type in c("resampling", "forest_seed", "model_family")) {
        other <- switch(type, resampling = zr, forest_seed = zs, model_family = maps$gam)
        sm <- similarity(z$data$effect_temperature, other$data$effect_temperature)
        rank1 <- z$summary$mean_abs[z$summary$contrast == "only"]
        rank2 <- other$summary$mean_abs[other$summary$contrast == "only"]
        reliability[[length(reliability) + 1L]] <- data.frame(scenario = s, replicate = rep, type,
          t(sm), top_driver_agreement = which.max(rank1) == which.max(rank2), row.names = NULL)
      }
      # Independent numerical prediction oracle for both scenario comparisons.
      xo <- a; xo$temperature <- b$temperature
      xh <- b; xh$temperature <- a$temperature
      add_check(s, rep, "only_prediction_oracle", max(abs(z$data$effect_temperature - (pred(f, xo, "rf") - pred(f, a, "rf")))))
      add_check(s, rep, "retained_prediction_oracle", max(abs(z$data$retained_temperature - (pred(f, b, "rf") - pred(f, xh, "rf")))))
      add_check(s, rep, "static_terrain_zero", max(abs(z$data$effect_terrain)))
      add_check(s, rep, "table_summary", abs(z$summary$mean_signed[1] - mean(z$data$total_change)))
      if (rep == 1L) {
        ar <- terra::rast(grid$raster, nlyrs = 3); names(ar) <- f$env_vars
        br <- ar
        terra::values(ar) <- as.matrix(a[, f$env_vars]); terra::values(br) <- as.matrix(b[, f$env_vars])
        mr <- cast_attribution(f, ar, br, models = "rf", block_size = 80L,
          filename = file.path(output, "maps", paste0(s, "_native.tif")), verbose = FALSE)
        mm <- cast_attribution(f, ar, br, models = "rf", block_size = 2000L, verbose = FALSE)
        values <- terra::values(terra::rast(mr$raster_file))
        compare <- function(x, y) {
          x <- as.numeric(x); y <- as.numeric(y)
          if (!identical(is.na(x), is.na(y))) return(Inf)
          max(abs(x[is.finite(x)] - y[is.finite(y)]), 0)
        }
        add_check(s, rep, "table_raster", compare(values, as.matrix(z$data[, mr$layers])))
        add_check(s, rep, "block_size", compare(values, terra::values(terra::rast(mm$raster_file))))
      }
    }, error = function(e) {
      message("ERROR ", id, ": ", conditionMessage(e))
      errors[[length(errors) + 1L]] <<- data.frame(scenario = s, replicate = rep, error = conditionMessage(e))
    })
    export_progress()
  }
  export_progress()
  stable <- identical(manifest, source_manifest())
  passed <- !length(errors) && all(vapply(checks, function(x) x$passed, logical(1))) && stable
  utils::capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  files <- list.files(output, recursive = TRUE, full.names = TRUE)
  artifacts <- data.frame(file = files, sha256 = vapply(files, function(f) digest::digest(file = f, algo = "sha256"), character(1)))
  saveRDS(list(status = if (passed) "completed" else "incomplete", mode = mode, source_manifest = manifest,
    source_stable = stable, expected_datasets = length(scenarios) * reps,
    completed_datasets = length(list.files(file.path(output, "datasets"))),
    error_count = length(errors), checks_passed = passed, artifacts,
    limitations = c("Known fitted-response contrasts, not ecological causal identification.",
      "Independent datasets are replicates; grid cells are dependent prediction targets.",
      "Range flags cannot detect all unsupported joint combinations.",
      "Cross-model spread is not a confidence interval.",
      "Global response and mean-shift are defined baselines, not all competing attribution methods.")),
    file.path(output, "receipt.rds"))
  if (passed && is.null(scenario)) writeLines(output, file.path(ROOT, paste0("latest_attribution_", mode, "completed.txt")))
  if (!passed) stop("Validation incomplete: ", output)
  message("VALIDATION COMPLETED: ", output)
  invisible(output)
}
args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L) run_validation(if (length(args)) args[1] else "production")
