# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
.cast_compare_inputs <- function(fit, current, future, groups, mask) {
  vars <- fit$env_vars
  if (!is.list(groups) || !length(groups) %in% 1:2 || is.null(names(groups)) ||
      anyNA(names(groups)) || any(!nzchar(names(groups))) || anyDuplicated(names(groups)))
    cli::cli_abort("groups must contain one or two uniquely named predictor groups.")
  changed <- as.character(unlist(groups, use.names = FALSE))
  if (!is.character(changed) || anyNA(changed) || (length(groups) > 1L && any(!lengths(groups))) ||
      anyDuplicated(changed) || !all(changed %in% vars))
    cli::cli_abort("groups must name disjoint fitted predictor subsets; only a single group may be empty (no fitted predictor changes).")
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
        cli::cli_abort("Every future raster must contain all changing predictors.")
      if (!isTRUE(terra::compareGeom(current, x, lyrs = FALSE, stopOnError = FALSE)))
        cli::cli_abort("Future geometry must match current; align inputs explicitly.")
      missing <- setdiff(vars, names(x))
      if (length(missing)) x <- c(x, current[[missing]])
      x[[vars]]
    } else {
      if (!is.data.frame(x) || !all(changed %in% names(x)))
        cli::cli_abort("Every future frame must contain all changing predictors.")
      if (all(c("lon", "lat") %in% names(current)) && all(c("lon", "lat") %in% names(x)))
        x <- .cast_align_future_env(current, x, "scenario comparison")
      if (nrow(x) != nrow(current)) cli::cli_abort("Future and current table row counts differ.")
      for (v in setdiff(vars, names(x))) x[[v]] <- current[[v]]
      .cast_check_numeric_predictors(x[, vars, drop = FALSE])
      x[, vars, drop = FALSE]
    }
  })
  list(current = current, future = future, groups = groups, static = static,
       design = if (length(groups) == 1L) "paired" else "factorial",
       mask = mask, raster = raster, vars = vars)
}

.cast_shared_raster_layers <- function(current, future, variables) {
  if (!length(variables)) return(TRUE)
  a <- terra::sources(current[[variables]], bands = TRUE)
  b <- terra::sources(future[[variables]], bands = TRUE)
  is.data.frame(a) && is.data.frame(b) && all(nzchar(a$source)) &&
    identical(a, b)
}

.cast_compare_common_mask <- function(plan, block_size, path = NULL,
                                      bounds = NULL, support_path = NULL) {
  if (!plan$raster) {
    ok <- .cast_finite_rows(plan$current[, plan$vars, drop = FALSE])
    supported <- if (!is.null(bounds)) .cast_attribution_inside(plan$current[, plan$vars, drop = FALSE], bounds) else NULL
    for (x in plan$future) {
      ok <- ok & .cast_finite_rows(x)
      if (!is.null(bounds)) supported <- supported & .cast_attribution_inside(x, bounds)
    }
    if (!is.null(plan$mask)) ok <- ok & is.finite(plan$mask) & plan$mask != 0
    return(list(domain = ok, support = if (!is.null(bounds)) ok & supported))
  }
  r <- plan$current
  # Shared file/band references prove that static measurements are unchanged.
  # General inputs without that proof retain the complete verification path.
  future <- lapply(plan$future, function(x) {
    vars <- if (.cast_shared_raster_layers(r, x, plan$static))
      setdiff(plan$vars, plan$static) else plan$vars
    if (length(vars)) x[[vars]] else NULL
  })
  future <- Filter(Negate(is.null), future)
  writer <- .cast_raster_writer(r, path, "common_domain", datatype = "INT1U", overwrite = TRUE)
  on.exit(.cast_writer_abort(writer), add = TRUE)
  support_writer <- if (!is.null(bounds)) .cast_raster_writer(r, support_path,
    "common_support", datatype = "INT1U", overwrite = TRUE)
  if (!is.null(support_writer)) on.exit(.cast_writer_abort(support_writer), add = TRUE)
  rows <- max(1L, as.integer(block_size/terra::ncol(r)))
  for (row in seq.int(1L, terra::nrow(r), by = rows)) {
    count <- min(rows, terra::nrow(r)-row+1L)
    ok <- .cast_mask_rows(plan$mask, row, count, count*terra::ncol(r))
    supported <- rep(FALSE, length(ok))
    if (any(ok)) {
      a <- .cast_raster_block(r, row, count)
      ok <- ok & .cast_finite_rows(a)
      if (!is.null(bounds)) supported <- .cast_attribution_inside(a, bounds)
      for (x in future) {
        b <- .cast_raster_block(x, row, count)
        ok <- ok & .cast_finite_rows(b)
        if (!is.null(bounds)) supported <- supported & .cast_attribution_inside(b, bounds[names(b)])
      }
    }
    .cast_writer_block(writer, as.integer(ok), row, count)
    if (!is.null(support_writer)) .cast_writer_block(support_writer,
      as.integer(ok & !is.na(supported) & supported), row, count)
  }
  .cast_writer_finish(writer)
  if (!is.null(support_writer)) .cast_writer_finish(support_writer)
  list(domain = path, support = if (!is.null(bounds)) support_path)
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
  context <- if (length(groups) == 1L && !length(groups[[1L]])) NULL else .cast_prediction_context(fit, weights, a, b, clamp)
  if (length(groups) == 1L) {
    future <- if (length(groups[[1L]])) context$predict(groups[[1L]]) else base
    if (any(!is.finite(base)) || any(!is.finite(future)))
      cli::cli_abort("An active model returned non-finite paired predictions.")
    change <- future-base
    faces <- cbind(current = as.numeric(base %*% weights), joint = as.numeric(future %*% weights))
    effects <- matrix(as.numeric(change %*% weights), ncol = 1L, dimnames = list(NULL, "joint_change"))
    sd <- matrix(.ensemble_rowsd(change), ncol = 1L, dimnames = list(NULL, "joint_change"))
    supported <- if (is.null(bounds)) rep(NA, nrow(a)) else
      .cast_attribution_inside(a, bounds) & .cast_attribution_inside(b, bounds)
    members <- cbind(base, future, change)
    colnames(members) <- unlist(lapply(c("current_", "future_", "change_"), paste0, names(weights)))
    return(list(faces = faces, effects = effects, sd = sd, support = supported,
      sign = .cast_sign_values(lapply(seq_along(weights), function(j) matrix(change[, j], ncol = 1L))),
      model_layers = members))
  }
  changed <- lapply(groups, function(vars) Reduce(`|`, lapply(vars,
    function(v) context$current[[v]] != context$future[[v]])))
  both <- changed[[1L]] & changed[[2L]]
  joint <- context$predict(names(b), base = base,
    changed_rows = changed[[1L]] | changed[[2L]])
  first <- second <- base
  first[changed[[1L]] & !changed[[2L]], ] <- joint[changed[[1L]] & !changed[[2L]], , drop = FALSE]
  second[!changed[[1L]] & changed[[2L]], ] <- joint[!changed[[1L]] & changed[[2L]], , drop = FALSE]
  predictions <- list(base, joint,
    context$predict(groups[[1L]], base = first, changed_rows = both),
    context$predict(groups[[2L]], base = second, changed_rows = both))
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
  list(faces = faces, effects = effects, sd = sd, support = supported,
       sign = .cast_sign_values(per_model))
}

.cast_compare_one <- function(fit, plan, scenario, weights, threshold, bounds, clamp,
                              baseline, domain, block_size, filename = NULL,
                              common_support = NULL, thresholds = numeric()) {
  current <- plan$current; future <- plan$future[[scenario]]
  shared_static <- plan$raster && .cast_shared_raster_layers(current, future, plan$static)
  changing <- setdiff(plan$vars, plan$static)
  paired <- identical(plan$design, "paired")
  components <- if (paired) "joint_change" else .cast_factorial_names()
  faces <- if (paired) c("current", "joint") else c("current", "joint", "climate", "lulc")
  layers <- c(faces, components, paste0("sd_", components), "support", "common_support",
              "sign_agreement_joint", "sign_conflict_joint",
              if (paired) unlist(lapply(c("current_", "future_", "change_"), paste0, names(weights))))
  acc <- .cast_factorial_accumulator(components, length(faces))
  common_acc <- .cast_factorial_accumulator(components, length(faces))
  cuts <- .cast_threshold_accumulator(unique(c(threshold, thresholds)), length(faces))
  sign_acc <- .cast_sign_accumulator(length(weights), components)
  seconds <- c(read_current = 0, read_future = 0, read_baseline = 0,
    predict = 0, write = 0)
  timed <- function(phase, expr) {
    start <- proc.time()[["elapsed"]]
    value <- force(expr)
    seconds[phase] <<- seconds[phase] + proc.time()[["elapsed"]]-start
    value
  }
  if (plan$raster) {
    base_r <- terra::rast(baseline)
    domain_r <- if (is.character(domain)) terra::rast(domain) else plan$mask
    support_r <- if (is.character(common_support)) terra::rast(common_support) else NULL
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
      a <- timed("read_current", if (plan$raster) .cast_raster_block(current, first, count) else current[idx, plan$vars, drop = FALSE])
      b <- timed("read_future", if (shared_static) {
        value <- a
        if (length(changing)) value[changing] <- .cast_raster_block(future[[changing]], first, count)
        value
      } else if (plan$raster) .cast_raster_block(future, first, count) else future[idx, , drop = FALSE])
      base <- timed("read_baseline", if (plan$raster) as.matrix(.cast_raster_block(base_r, first, count)) else baseline[idx, , drop = FALSE])
      keep <- keep & .cast_finite_rows(a) & .cast_finite_rows(b)
      if (any(keep)) {
        p <- timed("predict", .cast_compare_block(fit, a[keep, , drop = FALSE], b[keep, , drop = FALSE],
          plan$groups, plan$static, weights, base[keep, , drop = FALSE], bounds, clamp))
        .cast_factorial_accumulate(acc, p$faces, p$effects, p$support,
          if (!is.null(threshold)) (p$faces >= threshold)*1, model_sd = p$sd)
        cs <- if (is.null(common_support)) NULL else if (plan$raster)
          .cast_mask_rows(support_r, first, count, n)[keep] else common_support[idx][keep]
        if (!is.null(cs)) .cast_factorial_accumulate(common_acc,
          p$faces[cs, , drop = FALSE], p$effects[cs, , drop = FALSE], rep(TRUE, sum(cs)),
          if (!is.null(threshold)) (p$faces[cs, , drop = FALSE] >= threshold)*1,
          model_sd = p$sd[cs, , drop = FALSE])
        .cast_threshold_accumulate(cuts, p$faces, p$support, cs)
        .cast_sign_accumulate(sign_acc, p$sign, p$support, cs)
        z[keep, ] <- cbind(p$faces, p$effects, p$sd, as.numeric(p$support),
          if (is.null(cs)) rep(NA_real_, sum(keep)) else as.numeric(cs),
          p$sign$agreement[, if (paired) 1L else 3L], p$sign$conflict[, if (paired) 1L else 3L], p$model_layers)
      }
    }
    if (plan$raster) timed("write", .cast_writer_block(writer, z, first, count)) else data[idx, ] <- z
  }
  if (plan$raster) timed("write", .cast_writer_finish(writer))
  list(design = plan$design, summary = .cast_factorial_summary(acc, !is.null(bounds)),
    binary = if (!is.null(threshold)) .cast_factorial_binary(acc$binary, acc$n_binary),
    binary_supported = if (!is.null(threshold) && !is.null(bounds))
      .cast_factorial_binary(acc$binary_supported, acc$n_binary_supported),
    summary_common_support = if (!is.null(common_support)) {
      s <- .cast_factorial_summary(common_acc, FALSE)
      s$n_domain <- acc$moments[, "n"]
      s$support_share <- ifelse(s$n_domain > 0, s$n_valid/s$n_domain, NA_real_)
      s
    },
    binary_common_support = if (!is.null(threshold) && !is.null(common_support))
      .cast_factorial_binary(common_acc$binary, common_acc$n_binary),
    agreement = .cast_sign_summary(sign_acc, !is.null(bounds), !is.null(common_support)),
    threshold_sensitivity = .cast_threshold_summary(cuts, !is.null(bounds), !is.null(common_support)),
    timings = data.frame(phase = names(seconds), elapsed_s = as.numeric(seconds)),
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
  .cast_assert_compare_sources(plan)
  value <- .cast_compare_one(job$fit, plan, job$scenario, job$weights,
    job$threshold, job$bounds, job$clamp, job$baseline, job$domain,
    job$block_size, job$file, job$common_support, job$thresholds)
  .cast_cache_write(value, job$receipt, job$key, if (plan$raster) job$file else character())
  list(value = value, elapsed = as.numeric(difftime(Sys.time(), start, units = "secs")),
    started = format(start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
}

.cast_compare_tasks <- function(jobs, workers, cluster = NULL) {
  if (workers < 2L || length(jobs) < 2L) return(lapply(jobs, .cast_compare_task))
  workers <- min(workers, length(jobs))
  cl <- cluster
  if (is.null(cl)) {
    cl <- .cast_start_compare_cluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
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
  task <- function(job) tryCatch(
    utils::getFromNamespace(".cast_compare_task", "castSDM")(job),
    error = function(e) stop("Scenario ", job$scenario, ": ", conditionMessage(e), call. = FALSE))
  task <- utils::removeSource(task)
  environment(task) <- baseenv()
  parallel::parLapplyLB(cl, jobs, task, chunk.size = 1L)
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
#' One changing group gives a direct current/future comparison. Two groups
#' additionally construct first-group-only and second-group-only scenarios
#' without copying complete environmental raster stacks. Every contrast uses
#' the same positive-weight models and common finite domain. Model SD is
#' calculated after within-model paired contrasts, never from surface SDs.
#' Every active model must return finite predictions on eligible inputs;
#' otherwise the comparison fails instead of changing weights or domains.
#' @param fit A [cast_fit] or the `cast_result` returned by [cast_run()].
#' @param current Current predictor frame or `SpatRaster`/path.
#'   Local file-backed sources must remain available throughout comparison;
#'   missing sources fail before hashing, prediction or output creation.
#' @param future Future environment or uniquely named list of environments.
#'   All grouped predictors are required. Omitted static layers inherit current;
#'   supplied static layers must remain equal to current.
#' @param groups Named list of one or two disjoint predictor groups. Other
#'   fitted predictors are explicitly recorded as static.
#' @param cv Optional [cast_cv]. A `cast_result` supplies this automatically.
#' @param weights Optional finite named non-negative model weights. Without CV,
#'   the declared default is equal weights; no validation scores are invented.
#' @param threshold Optional finite threshold. Otherwise a usable pooled OOF
#'   threshold is inherited when available; no binary result is invented.
#' @param thresholds Optional additional diagnostic thresholds, counted during
#'   prediction. Changing these reuses verified predictions; it never refits.
#'   `summary(x, thresholds = ...)` also supports other thresholds by reading
#'   saved predictions. Counts use `>=` and report their domain denominators.
#' @param mask Optional aligned raster mask, or numeric/logical table-row mask.
#' @param support_probs Training marginal quantiles, default `c(.01, .99)`.
#'   `NULL` disables support. Missing training references give unassessed support.
#' @param clamp Apply the same training-range clamp to every model. Support is
#'   assessed on original inputs. Default `FALSE` unless inherited from a result.
#' @param domain `"common"` uses the intersection across all requested futures;
#'   `"scenario"` uses a separate finite domain per future.
#' @param output_dir Optional result directory; required for persistent rasters.
#'   A temporary directory is used when omitted for raster inputs.
#' @param resources Optional named list: total `cores` (default 1), `workers`
#'   (default up to 4), total `memory_mb` (default 1024), and existing `scratch_dir`.
#'   The budget determines blocks and learner threads. Scientific definitions
#'   and cache keys do not depend on worker counts or output paths.
#'   Local workers start without user/site startup profiles, use the caller's
#'   exact installed castSDM library, and verify their loaded namespace before
#'   input hashing or raster prediction. Version/path/code mismatches are
#'   reported separately; the safety check is not bypassed.
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
#' @details A single empty predictor group explicitly records model non-response:
#'   no fitted variable changes. Identical fitted-input scenarios share verified
#'   raster artifacts; aliases are recorded without duplicate worker timings.
#' @export
cast_compare <- function(fit, current, future, groups, cv = NULL, weights = NULL,
    threshold = NULL, mask = NULL, support_probs = c(.01, .99), clamp = FALSE,
    domain = c("common", "scenario"), output_dir = NULL, resources = list(),
    resume = TRUE, overwrite = FALSE, verbose = TRUE, thresholds = NULL) {
  domain <- match.arg(domain)
  wall_start <- Sys.time()
  thresholds <- .cast_check_thresholds(thresholds)
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
  .cast_assert_compare_sources(plan)
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
  if (!length(unlist(plan$groups))) cli::cli_inform("No fitted predictor changes in these scenarios: paired predictions will be identical. This is model non-response, not evidence of no ecological effect.")
  if (verbose) cli::cli_inform("Compare {length(plan$future)} future{?s}; groups {.val {names(groups)}}; {length(weights)} active models; {length(plan$static)} static predictors; {domain} domain.")
  cluster <- NULL
  worker_start <- Sys.time()
  if (budget$workers > 1L && length(plan$future) > 1L) {
    cluster <- .cast_start_compare_cluster(min(budget$workers, length(plan$future)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    if (verbose) cli::cli_inform("Verified {length(cluster)} scenario workers against the calling installation.")
  }
  worker_elapsed <- as.numeric(difftime(Sys.time(), worker_start, units = "secs"))
  memo <- new.env(parent = emptyenv())
  identity_start <- Sys.time()
  fit_key <- .cast_fit_identity(fit)
  current_key <- .cast_data_identity(plan$current, memo)
  future_keys <- vapply(plan$future, .cast_data_identity, character(1), memo = memo)
  mask_key <- if (!is.null(plan$mask)) .cast_data_identity(plan$mask, memo) else NULL
  identity_elapsed <- as.numeric(difftime(Sys.time(), identity_start, units = "secs"))
  events <- list()
  event <- function(stage, start, reused, scenario = NA_character_) {
    events[[length(events)+1L]] <<- data.frame(stage = stage, scenario = scenario, cache_hit = reused,
      elapsed_s = as.numeric(difftime(Sys.time(), start, units = "secs")),
      started = format(start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
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
  common_key <- .cast_digest(list(current_key, future_keys, mask_key, bounds,
    .cast_code_identity(.cast_compare_common_mask), .cast_code_identity(.cast_shared_raster_layers)))
  common_receipt <- if (!is.null(cache_dir)) file.path(cache_dir, "domain.rds")
  needs_common <- domain == "common" || !is.null(bounds)
  common <- if (needs_common && resume) .cast_cache_read(common_receipt, common_key) else NULL
  reused <- !is.null(common)
  if (reused && plan$raster) common <- lapply(common, function(x) if (!is.null(x)) file.path(cache_dir, x))
  if (needs_common && !reused) {
    common <- .cast_compare_common_mask(plan, block_size,
      if (plan$raster) file.path(cache_dir, "common_domain.tif"), bounds,
      if (plan$raster) file.path(cache_dir, "common_support.tif"))
    .cast_cache_write(if (plan$raster) lapply(common, function(x) if (!is.null(x)) basename(x)) else common,
      common_receipt, common_key, if (plan$raster) unlist(common, use.names = FALSE) else character())
  }
  event("common_domain", start, reused)
  results <- list()
  jobs <- list(); aliases <- list()
  for (scenario in names(plan$future)) {
    first <- names(plan$future)[match(future_keys[[scenario]], future_keys)]
    if (first != scenario) { aliases[[scenario]] <- first; next }
    start <- Sys.time()
    key <- .cast_digest(list(fit_key, current_key, future_keys[[scenario]], mask_key,
      if (needs_common) common_key, plan$groups, weights, threshold, support_probs, clamp,
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
        baseline = baseline, domain = if (domain == "common") common$domain,
        common_support = common$support, thresholds = thresholds,
        block_size = block_size, file = file,
        key = key, receipt = receipt, threads = budget$threads, scratch_dir = budget$scratch_dir)
      next
    }
    x$raster_file <- if (plan$raster) file else NULL
    missing_cuts <- setdiff(thresholds, x$threshold_sensitivity$threshold)
    if (length(missing_cuts)) {
      extra <- .cast_threshold_scan(x, missing_cuts, !is.null(bounds), !is.null(common$support),
        memory_mb = min(256, budget$worker_memory_mb))
      x$threshold_sensitivity <- rbind(x$threshold_sensitivity, extra)
      saved <- x; saved$raster_file <- NULL
      .cast_cache_write(saved, receipt, key, if (plan$raster) file else character())
    }
    x$scenario <- scenario
    results[[scenario]] <- x
    event("scenario", start, reused, scenario)
    if (verbose) cli::cli_inform("{scenario}: {if (reused) 'verified cache' else 'computed'}")
  }
  dispatch_start <- Sys.time()
  computed <- .cast_compare_tasks(jobs, budget$workers, cluster)
  dispatch_elapsed <- as.numeric(difftime(Sys.time(), dispatch_start, units = "secs"))
  for (scenario in names(computed)) {
    x <- computed[[scenario]]$value
    x$raster_file <- if (plan$raster) jobs[[scenario]]$file else NULL
    x$scenario <- scenario
    results[[scenario]] <- x
    events[[length(events)+1L]] <- data.frame(stage = "scenario", scenario = scenario,
      cache_hit = FALSE, elapsed_s = computed[[scenario]]$elapsed,
      started = computed[[scenario]]$started, finished = computed[[scenario]]$finished)
    if (verbose) cli::cli_inform("{scenario}: computed")
  }
  for (scenario in names(aliases)) {
    first <- aliases[[scenario]]
    start <- Sys.time()
    x <- results[[first]]; x$scenario <- scenario; x$alias_of <- first
    x$timings <- NULL
    results[[scenario]] <- x
    event("scenario_alias", start, TRUE, scenario)
    if (verbose) cli::cli_inform("{scenario}: identical fitted inputs; reusing {first}.")
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
    binary_supported = combine("binary_supported"),
    summary_common_support = combine("summary_common_support"),
    binary_common_support = combine("binary_common_support"), agreement = combine("agreement"),
    threshold_sensitivity = combine("threshold_sensitivity"),
    timings = combine("timings"),
    runtime = data.frame(phase = c("worker_startup", "input_identity", "scenario_dispatch_wall", "comparison_wall_before_export"),
      elapsed_s = c(worker_elapsed, identity_elapsed, dispatch_elapsed,
        as.numeric(difftime(Sys.time(), wall_start, units = "secs")))),
    common_support_file = if (plan$raster) common$support,
    results = results, groups = plan$groups, design = plan$design,
    changing_predictors = unlist(plan$groups, use.names = FALSE),
    scenario_aliases = aliases,
    static = plan$static, weights = weights, threshold = threshold, threshold_source = threshold_source,
    domain = domain, support_probs = support_probs, support_assessed = !is.null(bounds),
    resources = budget, events = do.call(rbind, events),
    contract = if (identical(plan$design, "paired")) .cast_paired_contract() else .cast_factorial_contract(), output_dir = output_dir,
    provenance = list(schema = "cast_comparison/2", model = fit_key, current = current_key,
      future = future_keys, mask = mask_key, session = utils::capture.output(utils::sessionInfo()))),
    class = "cast_comparison")
  if (!is.null(output_dir)) {
    tables <- file.path(output_dir, "tables"); dir.create(tables, showWarnings = FALSE)
    for (field in c("summary", "binary", "binary_supported", "summary_common_support",
        "binary_common_support", "agreement", "threshold_sensitivity", "timings", "runtime", "events", "contract")) {
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
  if (identical(x$design, "paired")) cli::cli_inform("Changing group: {names(x$groups)[1]}; all other predictors are static.") else
    cli::cli_inform("First group: {names(x$groups)[1]}; second: {names(x$groups)[2]}. *_main denotes baseline simple effects.")
  print(x$summary[, c("scenario", "component", "mean", "mean_supported", "mean_model_sd")], row.names = FALSE)
  invisible(x)
}

.cast_compare_thresholds <- function(object, thresholds) {
  thresholds <- .cast_check_thresholds(thresholds)
  rows <- list()
  for (scenario in names(object$results)) {
    result <- object$results[[scenario]]
    cached <- result$threshold_sensitivity
    missing <- setdiff(thresholds, cached$threshold)
    if (length(missing)) cached <- rbind(cached, .cast_threshold_scan(result,
      missing, object$support_assessed, !is.null(object$summary_common_support)))
    b <- cached[cached$threshold %in% thresholds, , drop = FALSE]
    b <- b[order(match(b$threshold, thresholds), match(b$domain,
      c("full", "supported", "common_supported"))), , drop = FALSE]
    b$scenario <- scenario
    rows[[length(rows)+1L]] <- b
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname cast_compare
#' @param object A comparison object.
#' @details Marginal support is assessed before clamping. `summary_common_support`
#'   and `binary_common_support` use the intersection of finite inputs and
#'   training quantile bounds over current and all requested futures; their
#'   denominators are identical across scenarios. The legacy `*_supported`
#'   tables use each scenario's own subset. Neither is multivariate AOA.
#'   `agreement` reports the fraction of cells with both positive and negative
#'   model contrasts and mean majority-sign agreement (unweighted across
#'   active models, zeros neutral within 1e-12). One model is unassessed.
#'   These summaries and SDs describe model disagreement, not confidence.
#' @param ... Additional plotting arguments, or ignored by summary/print.
#' @export
summary.cast_comparison <- function(object, thresholds = NULL, ...) {
  if (is.null(thresholds)) object$summary else .cast_compare_thresholds(object, thresholds)
}

#' @rdname cast_compare
#' @param x A saved comparison object to print or plot.
#' @param scenario Saved scenario name for plotting.
#' @param layers Saved prediction/contrast layers to plot.
#' @param basemap Geographic background for plotting: `"world"`, `"china"`
#'   or `"none"`; China retains the established South China Sea inset.
#' @export
plot.cast_comparison <- function(x, scenario = names(x$results)[1L],
                                 layers = NULL,
                                 basemap = c("world", "china", "none"), ...) {
  basemap <- match.arg(basemap)
  layers <- layers %||% if (identical(x$design, "paired")) c("current", "joint", "joint_change") else
    c("climate_main", "lulc_main", "joint_change", "interaction")
  if (!scenario %in% names(x$results)) cli::cli_abort("Unknown scenario.")
  y <- x$results[[scenario]]
  if (!is.null(y$raster_file)) {
    r <- terra::rast(y$raster_file)
    if (!all(layers %in% names(r))) cli::cli_abort("Unknown comparison layers.")
    check_suggested("ggplot2", "for cartographic comparison plots")
    check_suggested("patchwork", "for multi-panel comparison plots")
    effect_layers <- intersect(layers, .cast_factorial_names())
    limits <- if (length(effect_layers)) terra::minmax(r[[effect_layers]], compute = FALSE) else numeric()
    finite <- abs(limits[is.finite(limits)])
    effect_limit <- if (length(finite) && max(finite) > 0) max(finite) else 1
    plots <- lapply(layers, function(layer) {
      surface <- layer %in% c("current", "joint", "climate", "lulc")
      effect <- layer %in% .cast_factorial_names()
      discrete <- layer %in% c("support", "common_support", "sign_conflict_joint")
      cast_plot_effect_map(y$raster_file, layer, gsub("_", " ", layer),
        if (surface) "Suitability" else gsub("_", " ", layer), signed = effect,
        limits = if (surface || discrete || layer == "sign_agreement_joint") c(0, 1) else
          if (effect) c(-effect_limit, effect_limit) else NULL,
        method = if (discrete) "near" else "bilinear", basemap = basemap)
    })
    return(patchwork::wrap_plots(plots, ncol = min(2L, length(plots))))
  } else {
    d <- y$summary[y$summary$component %in% layers, ]
    if (!nrow(d)) cli::cli_abort("Unknown comparison layers.")
    graphics::barplot(stats::setNames(d$mean, d$component), ylab = "Mean suitability contrast", ...)
  }
  invisible(x)
}
