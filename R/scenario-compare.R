.cast_compare_inputs <- function(fit, current, future, groups, mask) {
  vars <- fit$env_vars
  if (!is.list(groups) || length(groups) != 2L || is.null(names(groups)) ||
      anyNA(names(groups)) || any(!nzchar(names(groups))) || anyDuplicated(names(groups)))
    cli::cli_abort("groups must contain exactly two uniquely named predictor groups.")
  changed <- unlist(groups, use.names = FALSE)
  if (!is.character(changed) || anyNA(changed) || any(!lengths(groups)) ||
      anyDuplicated(changed) || !all(changed %in% vars))
    cli::cli_abort("groups must name disjoint, nonempty subsets of fitted predictors.")
  safe <- vapply(names(groups), .cast_safe_name, character(1), arg = "groups")
  if (anyDuplicated(safe)) cli::cli_abort("Group names collide after filename normalization.")
  names(groups) <- safe
  raster <- inherits(current, "SpatRaster") || is.character(current)
  if (raster) {
    check_suggested("terra", "for raster scenario comparisons")
    if (is.character(current)) current <- terra::rast(current)
    if (!all(vars %in% names(current))) cli::cli_abort("current lacks fitted predictors.")
    current <- current[[vars]]
    mask <- .cast_raster_mask(mask, current)
  } else {
    if (!is.data.frame(current) || !nrow(current) || !all(vars %in% names(current)))
      cli::cli_abort("current must be a nonempty predictor frame or raster.")
    .cast_check_numeric_predictors(current[, vars, drop = FALSE])
    if (!is.null(mask) && ((!is.numeric(mask) && !is.logical(mask)) || length(mask) != nrow(current)))
      cli::cli_abort("A table mask must have one numeric/logical value per current row.")
  }
  if (is.data.frame(future) || inherits(future, "SpatRaster") || is.character(future)) future <- list(future = future)
  if (!is.list(future) || !length(future) || is.null(names(future)) || anyNA(names(future)) ||
      any(!nzchar(names(future))) || anyDuplicated(names(future)))
    cli::cli_abort("future must be an environment or a uniquely named list of environments.")
  safe <- vapply(names(future), .cast_safe_name, character(1), arg = "future")
  if (anyDuplicated(safe)) cli::cli_abort("Scenario names collide after filename normalization.")
  names(future) <- safe
  static <- setdiff(vars, changed)
  future <- lapply(future, function(x) {
    if (raster) {
      if (is.character(x)) x <- terra::rast(x)
      if (!inherits(x, "SpatRaster") || !all(changed %in% names(x)))
        cli::cli_abort("Every future raster must contain both changing groups.")
      if (!isTRUE(terra::compareGeom(current, x, lyrs = FALSE, stopOnError = FALSE)))
        cli::cli_abort("Future geometry must match current; align inputs explicitly.")
      missing <- setdiff(vars, names(x))
      if (length(missing)) x <- c(x, current[[missing]])
      x[[vars]]
    } else {
      if (!is.data.frame(x) || !all(changed %in% names(x)))
        cli::cli_abort("Every future frame must contain both changing groups.")
      if (all(c("lon", "lat") %in% names(current)) && all(c("lon", "lat") %in% names(x)))
        x <- .cast_align_future_env(current, x, "scenario comparison")
      if (nrow(x) != nrow(current)) cli::cli_abort("Future and current table row counts differ.")
      for (v in setdiff(vars, names(x))) x[[v]] <- current[[v]]
      .cast_check_numeric_predictors(x[, vars, drop = FALSE])
      x[, vars, drop = FALSE]
    }
  })
  list(current = current, future = future, groups = groups, static = static,
       mask = mask, raster = raster, vars = vars)
}

.cast_compare_common_mask <- function(plan, block_size, path = NULL) {
  if (!plan$raster) {
    ok <- .cast_finite_rows(plan$current[, plan$vars, drop = FALSE])
    for (x in plan$future) ok <- ok & .cast_finite_rows(x)
    if (!is.null(plan$mask)) ok <- ok & is.finite(plan$mask) & plan$mask != 0
    return(ok)
  }
  r <- plan$current
  writer <- .cast_raster_writer(r, path, "common_domain", datatype = "INT1U", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  rows <- max(1L, as.integer(block_size/terra::ncol(r)))
  for (row in seq.int(1L, terra::nrow(r), by = rows)) {
    count <- min(rows, terra::nrow(r)-row+1L)
    ok <- .cast_mask_rows(plan$mask, row, count, count*terra::ncol(r))
    if (any(ok)) {
      ok <- ok & .cast_finite_rows(.cast_raster_block(r, row, count))
      for (x in plan$future) ok <- ok & .cast_finite_rows(.cast_raster_block(x, row, count))
    }
    .cast_writer_block(writer, as.integer(ok), row, count)
  }
  .cast_writer_finish(writer)
  path
}

.cast_compare_baseline <- function(fit, plan, weights, clamp, block_size, path = NULL) {
  current <- plan$current
  old <- options(castSDM.predict_threads = getOption("castSDM.compare_threads", 1L))
  on.exit(options(old), add = TRUE)
  if (!plan$raster) {
    a <- current[, plan$vars, drop = FALSE]
    ok <- .cast_finite_rows(a)
    if (!is.null(plan$mask)) ok <- ok & is.finite(plan$mask) & plan$mask != 0
    p <- matrix(NA_real_, nrow(a), length(weights), dimnames = list(NULL, names(weights)))
    if (any(ok)) {
      x <- a[ok, , drop = FALSE]
      if (clamp && !is.null(fit$scaling$reference)) x <- .cast_clamp(x, fit$scaling$reference)
      for (m in names(weights)) p[ok, m] <- predict_single_model(fit$models[[m]], x)
    }
    return(p)
  }
  writer <- .cast_raster_writer(current, path, names(weights), datatype = "FLT8S", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  rows <- max(1L, as.integer(block_size/terra::ncol(current)))
  for (row in seq.int(1L, terra::nrow(current), by = rows)) {
    count <- min(rows, terra::nrow(current)-row+1L)
    ok <- .cast_mask_rows(plan$mask, row, count, count*terra::ncol(current))
    p <- matrix(NA_real_, length(ok), length(weights))
    if (any(ok)) {
      a <- .cast_raster_block(current, row, count)
      ok <- ok & .cast_finite_rows(a)
      if (any(ok)) {
        x <- a[ok, , drop = FALSE]
        if (clamp && !is.null(fit$scaling$reference)) x <- .cast_clamp(x, fit$scaling$reference)
        for (j in seq_along(weights)) p[ok, j] <- predict_single_model(fit$models[[names(weights)[j]]], x)
      }
    }
    .cast_writer_block(writer, p, row, count)
  }
  .cast_writer_finish(writer)
  path
}

.cast_compare_block <- function(fit, a, b, groups, static, weights, base, bounds, clamp) {
  if (length(static)) for (v in static) {
    if (any(a[[v]] != b[[v]])) cli::cli_abort("Static predictor {.val {v}} changes in future; include it in a group or fix its input.")
  }
  context <- .cast_prediction_context(fit, weights, a, b, clamp)
  predictions <- list(base, context$predict(names(b)), context$predict(groups[[1L]]), context$predict(groups[[2L]]))
  ok <- Reduce(`&`, lapply(predictions, function(x) rowSums(!is.finite(x)) == 0L))
  if (any(!ok)) cli::cli_abort(c(
    "An active model returned non-finite predictions on the eligible comparison domain.",
    i = "Check the model callback or declare an explicit valid mask; paired weights and common domains cannot be silently changed."))
  faces <- vapply(predictions, function(x) as.numeric(x %*% weights), numeric(nrow(a)))
  faces <- matrix(faces, nrow(a), 4L)
  faces[!ok, ] <- NA_real_
  effects <- .cast_factorial_values(faces)
  per_model <- lapply(seq_along(weights), function(j) {
    .cast_factorial_values(matrix(vapply(predictions, function(x) x[, j], numeric(nrow(a))), nrow(a), 4L))
  })
  sd <- matrix(NA_real_, nrow(a), 8L)
  if (length(weights) > 1L) for (k in seq_len(8L)) {
    sd[, k] <- .ensemble_rowsd(matrix(vapply(per_model, function(x) x[, k], numeric(nrow(a))), nrow(a)))
  }
  sd[!ok, ] <- NA_real_
  supported <- if (is.null(bounds)) rep(NA, nrow(a)) else
    .cast_attribution_inside(a, bounds) & .cast_attribution_inside(b, bounds)
  list(faces = faces, effects = effects, sd = sd, support = supported)
}

.cast_compare_one <- function(fit, plan, scenario, weights, threshold, bounds, clamp,
                              baseline, domain, block_size, filename = NULL) {
  current <- plan$current; future <- plan$future[[scenario]]
  layers <- c("current", "joint", "climate", "lulc", .cast_factorial_names(),
              paste0("sd_", .cast_factorial_names()), "support")
  acc <- .cast_factorial_accumulator()
  if (plan$raster) {
    base_r <- terra::rast(baseline)
    domain_r <- if (is.character(domain)) terra::rast(domain) else plan$mask
    writer <- .cast_raster_writer(current, filename, layers, datatype = "FLT8S", overwrite = TRUE)
    on.exit(.cast_writer_abort(writer), add = TRUE)
    rows <- max(1L, as.integer(block_size/terra::ncol(current)))
    starts <- seq.int(1L, terra::nrow(current), by = rows)
  } else {
    starts <- seq.int(1L, nrow(current), by = block_size)
    data <- matrix(NA_real_, nrow(current), length(layers), dimnames = list(NULL, layers))
  }
  for (first in starts) {
    if (plan$raster) {
      count <- min(rows, terra::nrow(current)-first+1L)
      keep <- .cast_mask_rows(domain_r, first, count, count*terra::ncol(current))
      n <- length(keep)
    } else {
      idx <- seq.int(first, min(nrow(current), first+block_size-1L))
      n <- length(idx)
      keep <- if (is.logical(domain)) domain[idx] else rep(TRUE, n)
      if (!is.null(plan$mask)) keep <- keep & is.finite(plan$mask[idx]) & plan$mask[idx] != 0
    }
    z <- matrix(NA_real_, n, length(layers), dimnames = list(NULL, layers))
    if (any(keep)) {
      a <- if (plan$raster) .cast_raster_block(current, first, count) else current[idx, plan$vars, drop = FALSE]
      b <- if (plan$raster) .cast_raster_block(future, first, count) else future[idx, , drop = FALSE]
      base <- if (plan$raster) as.matrix(.cast_raster_block(base_r, first, count)) else baseline[idx, , drop = FALSE]
      keep <- keep & .cast_finite_rows(a) & .cast_finite_rows(b)
      if (any(keep)) {
        p <- .cast_compare_block(fit, a[keep, , drop = FALSE], b[keep, , drop = FALSE],
          plan$groups, plan$static, weights, base[keep, , drop = FALSE], bounds, clamp)
        .cast_factorial_accumulate(acc, p$faces, p$effects, p$support,
          if (!is.null(threshold)) (p$faces >= threshold)*1, model_sd = p$sd)
        z[keep, ] <- cbind(p$faces, p$effects, p$sd, as.numeric(p$support))
      }
    }
    if (plan$raster) .cast_writer_block(writer, z, first, count) else data[idx, ] <- z
  }
  if (plan$raster) .cast_writer_finish(writer)
  list(summary = .cast_factorial_summary(acc, !is.null(bounds)),
    binary = if (!is.null(threshold)) .cast_factorial_binary(acc$binary, acc$n_binary),
    binary_supported = if (!is.null(threshold) && !is.null(bounds))
      .cast_factorial_binary(acc$binary_supported, acc$n_binary_supported),
    data = if (!plan$raster) as.data.frame(data), artifact = if (plan$raster) basename(filename))
}

.cast_compare_task <- function(job) {
  start <- Sys.time()
  old <- options(castSDM.predict_threads = job$threads)
  on.exit(options(old), add = TRUE)
  plan <- job$plan
  if (plan$raster) {
    old_temp <- terra::terraOptions(print = FALSE)$tempdir
    terra::terraOptions(tempdir = job$scratch_dir)
    on.exit(terra::terraOptions(tempdir = old_temp), add = TRUE)
  }
  if (plan$raster && inherits(plan$current, "PackedSpatRaster")) {
    plan$current <- terra::unwrap(plan$current)
    plan$future <- lapply(plan$future, terra::unwrap)
    if (!is.null(plan$mask)) plan$mask <- terra::unwrap(plan$mask)
  }
  value <- .cast_compare_one(job$fit, plan, job$scenario, job$weights,
    job$threshold, job$bounds, job$clamp, job$baseline, job$domain,
    job$block_size, job$file)
  .cast_cache_write(value, job$receipt, job$key, if (plan$raster) job$file else character())
  list(value = value, elapsed = as.numeric(difftime(Sys.time(), start, units = "secs")))
}

.cast_compare_tasks <- function(jobs, workers) {
  if (workers < 2L || length(jobs) < 2L) return(lapply(jobs, .cast_compare_task))
  workers <- min(workers, length(jobs))
  cl <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  initialize <- function(libs, version, code) {
    .libPaths(libs)
    loadNamespace("castSDM")
    if (!identical(as.character(utils::packageVersion("castSDM")), version) ||
        !identical(utils::getFromNamespace(".cast_code_identity", "castSDM")(
          utils::getFromNamespace(".cast_compare_one", "castSDM")), code))
      stop("Worker castSDM differs from the calling session. Install the current package before parallel use.")
    TRUE
  }
  environment(initialize) <- baseenv()
  parallel::clusterCall(cl, initialize, .libPaths(), as.character(utils::packageVersion("castSDM")),
    .cast_code_identity(.cast_compare_one))
  for (i in seq_along(jobs)) {
    if (jobs[[i]]$plan$raster) {
      p <- jobs[[i]]$plan
      p$current <- terra::wrap(p$current, proxy = TRUE)
      p$future <- lapply(p$future, terra::wrap, proxy = TRUE)
      if (!is.null(p$mask)) p$mask <- terra::wrap(p$mask, proxy = TRUE)
      jobs[[i]]$plan <- p
    }
    jobs[[i]]$fit$models <- lapply(jobs[[i]]$fit$models, function(m) {
      if (is.function(m$predict_function)) m$predict_function <- .cast_seal_callback(m$predict_function)
      m
    })
  }
  task <- function(job) utils::getFromNamespace(".cast_compare_task", "castSDM")(job)
  environment(task) <- baseenv()
  parallel::parLapplyLB(cl, jobs, task)
}

# No synthetic threshold or per-pixel weight renormalisation in comparisons.
.cast_compare_threshold <- function(cv, weights) {
  oof <- cv$oof
  cols <- paste0("HSS_", names(weights))
  if (!is.data.frame(oof) || !all(c("obs", cols) %in% names(oof))) return(NULL)
  coverage <- .cast_oof_fold_coverage(cv, names(weights))
  if (!all(coverage$ok)) return(NULL)
  values <- as.matrix(oof[, cols, drop = FALSE])
  ok <- .cast_finite_rows(values) & is.finite(oof$obs) & oof$obs %in% c(0, 1)
  if (sum(ok) < 10L || length(unique(oof$obs[ok])) != 2L) return(NULL)
  cast_threshold(as.numeric(values[ok, , drop = FALSE] %*% weights), oof$obs[ok],
    method = cv$threshold_method %||% "max_tss")
}

#' Compare Grouped Environmental Scenarios with Frozen SDMs
#'
#' Constructs current, first-group-only, second-group-only and joint scenarios
#' without copying complete environmental raster stacks. Every contrast uses
#' the same positive-weight models and common finite domain. Model SD is
#' calculated after within-model paired contrasts, never from surface SDs.
#' Every active model must return finite predictions on eligible inputs;
#' otherwise the comparison fails instead of changing weights or domains.
#' @param fit A [cast_fit] or the `cast_result` returned by [cast_run()].
#' @param current Current predictor frame or `SpatRaster`/path.
#' @param future Future environment or uniquely named list of environments.
#'   All grouped predictors are required. Omitted static layers inherit current;
#'   supplied static layers must remain equal to current.
#' @param groups Named list of exactly two disjoint predictor groups. Other
#'   fitted predictors are explicitly recorded as static.
#' @param cv Optional [cast_cv]. A `cast_result` supplies this automatically.
#' @param weights Optional finite named non-negative model weights. Without CV,
#'   the declared default is equal weights; no validation scores are invented.
#' @param threshold Optional finite threshold. Otherwise a usable pooled OOF
#'   threshold is inherited when available; no binary result is invented.
#' @param mask Optional aligned raster mask, or numeric/logical table-row mask.
#' @param support_probs Training marginal quantiles, default `c(.01, .99)`.
#'   `NULL` disables support. Missing training references give unassessed support.
#' @param clamp Apply the same training-range clamp to every model. Support is
#'   assessed on original inputs. Default `FALSE` unless inherited from a result.
#' @param domain `"common"` uses the intersection across all requested futures;
#'   `"scenario"` uses a separate four-corner finite domain per future.
#' @param output_dir Optional result directory; required for persistent rasters.
#'   A temporary directory is used when omitted for raster inputs.
#' @param resources Optional named list: total `cores` (default 1), `workers`
#'   (default up to 4), total `memory_mb` (default 1024), and existing `scratch_dir`.
#'   The budget determines blocks and learner threads. Scientific definitions
#'   and cache keys do not depend on worker counts or output paths.
#' @param resume Reuse matching receipts only after input/model and output hashes
#'   are verified. Defaults to `TRUE`.
#' @param overwrite Replace stale public artifacts; default `FALSE`. Matching
#'   artifacts can be resumed without overwriting.
#' @param verbose Print the execution plan and progress.
#' @return A `cast_comparison` with combined `summary`, optional binary tables,
#'   named scenario `results`, group/static metadata, frozen `weights`, the
#'   threshold source, artifact paths, stage events and provenance. Legacy
#'   `climate_*` and `lulc_*` columns refer to the first and second named groups;
#'   `*_main` are baseline simple effects; `*_average` are averaged effects.
#' @examples
#' reference <- data.frame(temperature = seq(-1, 1, length.out = 20), land = rep(0:1, 10))
#' model <- list(intercept = -.3)
#' fit <- cast_as_fit(list(example = model), reference, function(model, newdata)
#'   stats::plogis(model$intercept + newdata$temperature - newdata$land))
#' future <- reference; future$temperature <- future$temperature + .2
#' comparison <- cast_compare(fit, reference, future,
#'   groups = list(climate = "temperature", land_use = "land"), verbose = FALSE)
#' comparison$summary
#' @export
cast_compare <- function(fit, current, future, groups, cv = NULL, weights = NULL,
    threshold = NULL, mask = NULL, support_probs = c(.01, .99), clamp = FALSE,
    domain = c("common", "scenario"), output_dir = NULL, resources = list(),
    resume = TRUE, overwrite = FALSE, verbose = TRUE) {
  domain <- match.arg(domain)
  protocol <- NULL
  if (inherits(fit, "cast_result")) {
    run <- fit; fit <- run$fit_full %||% run$fit; cv <- cv %||% run$cv
    protocol <- run$protocol
    if (is.null(weights) && !protocol$ensemble_method %in% c("weighted", "best", "equal"))
      cli::cli_abort("Paired comparisons require explicit linear weights for a nonlinear ensemble result.")
    if (missing(clamp)) clamp <- protocol$clamp %||% FALSE
    weights <- weights %||% run$ensemble$weights
  }
  .cast_check_fit(fit)
  budget <- .cast_resources(resources)
  if (budget$workers > 1L && !file.exists(file.path(
      getNamespaceInfo(asNamespace("castSDM"), "path"), "Meta", "package.rds"))) {
    cli::cli_warn("Development namespace: using one process. Install this package to use scenario workers.")
    budget$workers <- 1L
    budget$threads <- budget$cores
    budget$worker_memory_mb <- budget$memory_mb
  }
  old <- options(castSDM.predict_threads = budget$threads, castSDM.compare_threads = budget$cores)
  on.exit(options(old), add = TRUE)
  plan <- .cast_compare_inputs(fit, current, future, groups, mask)
  if (plan$raster) {
    old_temp <- terra::terraOptions(print = FALSE)$tempdir
    terra::terraOptions(tempdir = budget$scratch_dir)
    on.exit(terra::terraOptions(tempdir = old_temp), add = TRUE)
  }
  method <- if (is.null(cv)) "equal" else protocol$ensemble_method %||% "weighted"
  if (is.null(weights) && !method %in% c("weighted", "best", "equal"))
    cli::cli_abort("Paired comparisons require linear frozen weights; supply explicit weights for this fit.")
  weights <- .cast_attribution_weights(fit, cv, NULL, weights, method,
    protocol$ensemble_min_score %||% 0, protocol$ensemble_decay %||% 1,
    protocol$ensemble_fallback %||% "best", protocol$ensemble_min_metric_folds %||% 2L)
  for (m in names(weights)) {
    pkg <- .cast_engine_package(fit$models[[m]]$name)
    if (!identical(fit$models[[m]]$type, "external") && !is.null(pkg)) check_suggested(pkg, paste("for", m))
  }
  threshold_source <- if (!is.null(threshold)) "supplied" else "not_supplied"
  if (is.null(threshold) && !is.null(cv)) {
    threshold <- .cast_compare_threshold(cv, weights)
    if (!is.null(threshold)) threshold_source <- "out_of_fold_cv"
  }
  if (!is.null(threshold)) threshold <- .cast_check_scalar(threshold, "threshold")
  reference <- fit$scaling$reference
  bounds <- if (!is.null(support_probs) && !is.null(reference)) .cast_support_bounds(reference, support_probs) else NULL
  if (!is.null(bounds) && any(!vapply(bounds, function(b) all(is.finite(b)), logical(1))))
    cli::cli_abort("Training reference has no finite support bounds.")
  if (plan$raster && is.null(output_dir)) output_dir <- tempfile("cast_compare_")
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- if (!is.null(output_dir)) file.path(output_dir, ".cache") else NULL
  raster_dir <- if (!is.null(output_dir)) file.path(output_dir, "rasters") else NULL
  if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(raster_dir)) dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(output_dir)) {
    lock <- .cast_lock(file.path(output_dir, ".comparison.lock"))
    on.exit(filelock::unlock(lock), add = TRUE)
  }
  block_size <- max(1L, as.integer(min(1000000, budget$worker_memory_mb*1e6/
    (8 * (4*length(plan$vars) + 16*length(weights) + 32) * 4))))
  if (verbose) cli::cli_inform("Compare {length(plan$future)} future{?s}; groups {.val {names(groups)}}; {length(weights)} active models; {length(plan$static)} static predictors; {domain} domain.")
  memo <- new.env(parent = emptyenv())
  fit_key <- .cast_fit_identity(fit)
  current_key <- .cast_data_identity(plan$current, memo)
  future_keys <- vapply(plan$future, .cast_data_identity, character(1), memo = memo)
  mask_key <- if (!is.null(plan$mask)) .cast_data_identity(plan$mask, memo) else NULL
  events <- list()
  event <- function(stage, start, reused, scenario = NA_character_) {
    events[[length(events)+1L]] <<- data.frame(stage = stage, scenario = scenario, cache_hit = reused,
      elapsed_s = as.numeric(difftime(Sys.time(), start, units = "secs")),
      finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  }
  start <- Sys.time()
  base_key <- .cast_digest(list(fit_key, current_key, mask_key, names(weights), clamp, .cast_code_identity(.cast_compare_baseline)))
  base_receipt <- if (!is.null(cache_dir)) file.path(cache_dir, "current.rds")
  baseline <- if (resume) .cast_cache_read(base_receipt, base_key) else NULL
  reused <- !is.null(baseline)
  if (reused && plan$raster) baseline <- file.path(cache_dir, baseline)
  if (!reused) {
    baseline <- .cast_compare_baseline(fit, plan, weights, clamp, block_size,
      if (plan$raster) file.path(cache_dir, "current_models.tif"))
    .cast_cache_write(if (plan$raster) basename(baseline) else baseline, base_receipt, base_key,
      if (plan$raster) baseline else character())
  }
  event("current_prediction", start, reused)
  start <- Sys.time()
  common_key <- .cast_digest(list(current_key, future_keys, mask_key, .cast_code_identity(.cast_compare_common_mask)))
  common_receipt <- if (!is.null(cache_dir)) file.path(cache_dir, "domain.rds")
  common <- if (domain == "common" && resume) .cast_cache_read(common_receipt, common_key) else NULL
  reused <- !is.null(common)
  if (reused && plan$raster) common <- file.path(cache_dir, common)
  if (domain == "common" && !reused) {
    common <- .cast_compare_common_mask(plan, block_size,
      if (plan$raster) file.path(cache_dir, "common_domain.tif"))
    .cast_cache_write(if (plan$raster) basename(common) else common, common_receipt, common_key,
      if (plan$raster) common else character())
  }
  event("common_domain", start, reused)
  results <- list()
  jobs <- list()
  for (scenario in names(plan$future)) {
    start <- Sys.time()
    key <- .cast_digest(list(fit_key, current_key, future_keys[[scenario]], mask_key,
      if (domain == "common") common_key, plan$groups, weights, threshold, support_probs, clamp,
      .cast_code_identity(.cast_compare_one), .cast_code_identity(.cast_compare_block), .cast_code_identity(.cast_factorial_accumulate)))
    receipt <- if (!is.null(raster_dir)) file.path(raster_dir, paste0(scenario, ".rds"))
    x <- if (resume) .cast_cache_read(receipt, key) else NULL
    reused <- !is.null(x)
    file <- if (plan$raster) file.path(raster_dir, paste0(scenario, ".tif"))
    if (!reused) {
      if (!is.null(file) && file.exists(file) && !overwrite)
        cli::cli_abort("Stale/unverified scenario {.val {scenario}} exists. Use a new output_dir or overwrite = TRUE.")
      job_plan <- plan
      job_plan$future <- plan$future[scenario]
      jobs[[scenario]] <- list(fit = fit, plan = job_plan, scenario = scenario,
        weights = weights, threshold = threshold, bounds = bounds, clamp = clamp,
        baseline = baseline, domain = common, block_size = block_size, file = file,
        key = key, receipt = receipt, threads = budget$threads, scratch_dir = budget$scratch_dir)
      next
    }
    x$raster_file <- if (plan$raster) file else NULL
    x$scenario <- scenario
    results[[scenario]] <- x
    event("scenario", start, reused, scenario)
    if (verbose) cli::cli_inform("{scenario}: {if (reused) 'verified cache' else 'computed'}")
  }
  computed <- .cast_compare_tasks(jobs, budget$workers)
  for (scenario in names(computed)) {
    x <- computed[[scenario]]$value
    x$raster_file <- if (plan$raster) jobs[[scenario]]$file else NULL
    x$scenario <- scenario
    results[[scenario]] <- x
    events[[length(events)+1L]] <- data.frame(stage = "scenario", scenario = scenario,
      cache_hit = FALSE, elapsed_s = computed[[scenario]]$elapsed,
      finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
    if (verbose) cli::cli_inform("{scenario}: computed")
  }
  results <- results[names(plan$future)]
  combine <- function(field) {
    rows <- lapply(names(results), function(s) {
      x <- results[[s]][[field]]
      if (is.null(x)) return(NULL)
      x$scenario <- s; x
    })
    if (!any(lengths(rows))) return(NULL)
    do.call(rbind, rows)
  }
  result <- structure(list(summary = combine("summary"), binary = combine("binary"),
    binary_supported = combine("binary_supported"), results = results, groups = plan$groups,
    static = plan$static, weights = weights, threshold = threshold, threshold_source = threshold_source,
    domain = domain, support_probs = support_probs, support_assessed = !is.null(bounds),
    resources = budget, events = do.call(rbind, events),
    contract = .cast_factorial_contract(), output_dir = output_dir,
    provenance = list(schema = "cast_comparison/1", model = fit_key, current = current_key,
      future = future_keys, mask = mask_key, session = utils::capture.output(utils::sessionInfo()))),
    class = "cast_comparison")
  if (!is.null(output_dir)) {
    tables <- file.path(output_dir, "tables"); dir.create(tables, showWarnings = FALSE)
    for (field in c("summary", "binary", "binary_supported", "events", "contract")) {
      path <- file.path(tables, paste0(field, ".csv"))
      if (!is.null(result[[field]])) utils::write.csv(result[[field]], path, row.names = FALSE) else
        if (file.exists(path)) unlink(path)
    }
    .cast_atomic_rds(result, file.path(output_dir, "comparison.rds"))
  }
  result
}

#' @export
print.cast_comparison <- function(x, ...) {
  cli::cli_h1("castSDM paired scenario comparison")
  cli::cli_inform("First group: {names(x$groups)[1]}; second: {names(x$groups)[2]}. *_main denotes baseline simple effects.")
  print(x$summary[, c("scenario", "component", "mean", "mean_supported", "mean_model_sd")], row.names = FALSE)
  invisible(x)
}

.cast_compare_thresholds <- function(object, thresholds) {
  if (!is.numeric(thresholds) || !length(thresholds) || any(!is.finite(thresholds)) || anyDuplicated(thresholds))
    cli::cli_abort("thresholds must be unique finite numbers.")
  rows <- list()
  for (scenario in names(object$results)) {
    result <- object$results[[scenario]]
    raster <- !is.null(result$raster_file)
    input <- if (raster) terra::rast(result$raster_file)[[c("current", "joint", "climate", "lulc", "support")]] else
      result$data[, c("current", "joint", "climate", "lulc", "support"), drop = FALSE]
    counts <- array(0, c(length(thresholds), 4L, 2L))
    n <- numeric(2L)
    if (raster) {
      block_rows <- max(1L, as.integer(100000L/terra::ncol(input)))
      starts <- seq.int(1L, terra::nrow(input), by = block_rows)
    } else starts <- seq.int(1L, nrow(input), by = 100000L)
    for (first in starts) {
      x <- if (raster) terra::values(input, row = first,
        nrows = min(block_rows, terra::nrow(input)-first+1L), mat = TRUE) else
        as.matrix(input[seq.int(first, min(nrow(input), first+99999L)), , drop = FALSE])
      valid <- rowSums(!is.finite(x[, 1:4, drop = FALSE])) == 0L
      supported <- valid & is.finite(x[, 5L]) & x[, 5L] != 0
      n <- n + c(sum(valid), sum(supported))
      for (i in seq_along(thresholds)) {
        counts[i, , 1L] <- counts[i, , 1L] + colSums(x[valid, 1:4, drop = FALSE] >= thresholds[i])
        counts[i, , 2L] <- counts[i, , 2L] + colSums(x[supported, 1:4, drop = FALSE] >= thresholds[i])
      }
    }
    for (i in seq_along(thresholds)) for (d in seq_len(if (object$support_assessed) 2L else 1L)) {
      b <- .cast_factorial_binary(counts[i, , d], n[d])
      b$scenario <- scenario; b$threshold <- thresholds[i]
      b$domain <- c("full", "supported")[d]
      rows[[length(rows)+1L]] <- b
    }
  }
  do.call(rbind, rows)
}

#' @rdname cast_compare
#' @param object A comparison object.
#' @param thresholds Optional vector for `summary()`. Computes threshold
#'   sensitivity from saved four-corner predictions without fitting/prediction.
#' @param ... Additional plotting arguments, or ignored by summary/print.
#' @export
summary.cast_comparison <- function(object, thresholds = NULL, ...) {
  if (is.null(thresholds)) object$summary else .cast_compare_thresholds(object, thresholds)
}

#' @export
plot.cast_comparison <- function(x, scenario = names(x$results)[1L],
                                 layers = c("climate_main", "lulc_main", "joint_change", "interaction"), ...) {
  if (!scenario %in% names(x$results)) cli::cli_abort("Unknown scenario.")
  y <- x$results[[scenario]]
  if (!is.null(y$raster_file)) {
    r <- terra::rast(y$raster_file)
    if (!all(layers %in% names(r))) cli::cli_abort("Unknown comparison layers.")
    terra::plot(r[[layers]], ...)
  } else {
    d <- y$summary[y$summary$component %in% layers, ]
    if (!nrow(d)) cli::cli_abort("Unknown comparison layers.")
    graphics::barplot(stats::setNames(d$mean, d$component), ylab = "Mean suitability contrast", ...)
  }
  invisible(x)
}
