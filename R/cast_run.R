# cast_run(): the single production pipeline ---------------------------------
#
# One orchestration driven by a cast_protocol(): prepare -> select -> fit ->
# nested spatial CV -> one-shot hold-out evaluation -> full-data refit ->
# prediction/ensemble -> report. When output_dir is supplied, every stage
# lands on disk with content digests (receipts) in manifest.rds.
#
# Protocol guarantees (the four evaluation rules):
#   1. selection, cross-validation, thresholds and ensemble weights are
#      estimated on the TRAINING split only;
#   2. the hold-out split is evaluated exactly once, with the pooled
#      out-of-fold thresholds -- never for choosing anything;
#   3. spatial predictions and ensemble maps use the full-data refit;
#   4. threshold and projection surface share one source: the ensemble
#      threshold is optimised on the out-of-fold ensemble surface and
#      applied unchanged to the refit surface.

#' Run the castSDM Pipeline
#'
#' The single production workflow: spatial hold-out split, variable
#' selection, model fitting, nested spatial cross-validation, a one-shot
#' hold-out evaluation, full-data refit, spatial prediction/ensemble (and
#' optional raster projection), and an ODMAP-aligned report -- all driven by
#' one [cast_protocol()]. Every scientific setting is declared once in the
#' protocol; the effective values are echoed into the returned object and,
#' when `output_dir` is supplied, recorded stage by stage with content
#' digests in `manifest.rds`.
#'
#' @section Evaluation protocol:
#' - Variable selection is re-run inside every outer spatial CV fold; the
#'   hold-out split never influences selection, tuning, thresholds or
#'   ensemble weights.
#' - Per-engine binary thresholds and the ensemble threshold are computed
#'   from pooled out-of-fold predictions ([cast_cv()]); the hold-out split
#'   is then evaluated once with those thresholds ([cast_evaluate()]).
#' - The published predictions and ensemble maps use the full-data refit,
#'   carrying the training split's tuned configuration
#'   (`castSDM:::.cast_tuned_args()`).
#' - All stage seeds are applied locally; the caller's RNG state is
#'   restored afterwards.
#'
#' @param species_data A `data.frame` with `lon`, `lat`, `presence` (0/1)
#'   and the environmental predictors.
#' @param env_data Optional `data.frame` of environmental values (same
#'   predictors) for per-site prediction and the ensemble table.
#' @param env_raster Optional `terra::SpatRaster` of current conditions.
#'   When supplied (with `output_dir`), current and optional future
#'   projections are written as GeoTIFFs via [cast_project_raster()].
#' @param future_rasters Optional named list of future `SpatRaster`s; names
#'   become scenario labels. Requires `env_raster`.
#' @param protocol A [cast_protocol()] object.
#' @param select_var_id,select_category Data-aligned selection metadata
#'   passed to [cast_select()] (e.g. variable-level identities and
#'   categories); must have one entry per predictor seen by the screen.
#' @param select_ensure Data-aligned coverage requirements (named list of
#'   logical flag vectors) passed to [cast_select()].
#' @param output_dir Optional directory for stage artifacts, tables, the
#'   ODMAP report and the run manifest. Created (recursively) when missing.
#' @param report_meta Optional named list of report metadata passed to
#'   [cast_report_odmap()] (`taxon`, `authors`, `title`, `purpose`, ...).
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return A `cast_result` object carrying `screen`, `fit` (training split),
#'   `cv`, `eval` (one-shot hold-out), `fit_full` (full-data refit),
#'   `predict`/`ensemble` (when `env_data` is supplied), `projection` (when
#'   `env_raster` is supplied), the effective `protocol` and the receipt
#'   log. Stage files live under `output_dir` when supplied.
#' @seealso [cast_protocol()], [cast_select()], [cast_cv()],
#'   [cast_evaluate()], [cast_report_odmap()]
#' @export
cast_run <- function(species_data,
                     env_data = NULL,
                     env_raster = NULL,
                     future_rasters = NULL,
                     protocol = cast_protocol(),
                     select_var_id = NULL,
                     select_category = NULL,
                     select_ensure = NULL,
                     output_dir = NULL,
                     report_meta = list(),
                     verbose = TRUE) {
  protocol <- .cast_validate_protocol(protocol)
  cl <- match.call()
  if (!is.data.frame(species_data)) {
    cli::cli_abort(
      "{.arg species_data} must be a data.frame, not {.obj_type_friendly {species_data}}."
    )
  }
  if (is.null(env_data) && is.null(env_raster)) {
    cli::cli_inform("{.arg env_data}/{.arg env_raster} not supplied: no spatial prediction stage will run.")
  }
  if (!is.null(future_rasters) && is.null(env_raster)) {
    cli::cli_abort("{.arg future_rasters} requires {.arg env_raster}.")
  }
  if (!is.null(env_raster) && is.null(output_dir)) {
    cli::cli_abort("{.arg env_raster} projection writes GeoTIFFs; supply {.arg output_dir}.")
  }
  out <- .cast_run_dirs(output_dir)
  receipts <- list()

  add_receipt <- function(stage, inputs, value, tables = list(),
                          extra_files = character(0)) {
    elapsed <- Sys.time() - stage_started
    outputs <- data.frame(file = character(0), md5 = character(0),
                          stringsAsFactors = FALSE)
    if (!is.null(out)) {
      files <- character(0)
      stage_file <- file.path(out$stages, paste0(stage, ".rds"))
      saveRDS(value, stage_file)
      files <- c(files, stage_file)
      for (nm in names(tables)) {
        f <- file.path(out$tables, paste0(nm, ".csv"))
        utils::write.csv(tables[[nm]], f, row.names = FALSE)
        files <- c(files, f)
      }
      extra_files <- extra_files[file.exists(extra_files)]
      files <- c(files, extra_files)
      outputs <- data.frame(file = files,
                            md5 = unname(tools::md5sum(files)),
                            stringsAsFactors = FALSE)
    }
    receipts[[length(receipts) + 1L]] <<- list(
      stage = stage,
      inputs = inputs,
      outputs = outputs,
      elapsed_secs = round(as.numeric(elapsed, units = "secs"), 2),
      finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
    invisible(value)
  }
  stage_started <- Sys.time()

  if (verbose) cli::cli_h1("castSDM run ({length(protocol$models)} engine{?s}: {.val {protocol$models}})")

  # ---- 1. Prepare: spatial hold-out split --------------------------------
  if (verbose) cli::cli_h2("Prep: spatial hold-out split")
  stage_started <- Sys.time()
  split <- cast_prepare(
    species_data, train_fraction = protocol$train_fraction,
    seed = protocol$seed, split = "spatial",
    block_method = protocol$block_method, verbose = verbose
  )
  if (verbose) {
    cli::cli_inform("Train: {nrow(split$train)} | Hold-out: {nrow(split$test)} | Predictors: {length(split$env_vars)}")
  }
  add_receipt("01_prepare",
    .cast_digest(list(species_data, protocol$train_fraction,
                      protocol$block_method, protocol$seed)),
    split)

  # ---- 2. Variable selection (training split only) ------------------------
  if (verbose) cli::cli_h2("Select: {protocol$select_method} on the training split")
  stage_started <- Sys.time()
  screen <- cast_select(
    split$train,
    method = protocol$select_method,
    keep = protocol$select_keep,
    corcut = protocol$select_corcut,
    var_id = select_var_id,
    category = select_category,
    ncov = protocol$select_ncov,
    maxncov = protocol$select_maxncov,
    ensure = select_ensure,
    seed = protocol$seed, verbose = verbose
  )
  add_receipt("02_select",
    .cast_digest(list(split$train, protocol$select_method,
                      protocol$select_keep, protocol$select_corcut,
                      protocol$select_ncov, protocol$select_maxncov,
                      protocol$seed)),
    screen, tables = list(selection_scores = screen$scores))

  # ---- 3. Fit on the training split ---------------------------------------
  if (verbose) cli::cli_h2("Fit: engines on the training split")
  stage_started <- Sys.time()
  fit <- cast_fit(
    split$train, screen = screen, models = protocol$models,
    rf_ntree = protocol$rf_ntree, rf_mtry = protocol$rf_mtry,
    maxent_classes = protocol$maxent_classes,
    maxent_regmult = protocol$maxent_regmult,
    tune = protocol$tune, tune_folds = protocol$tune_folds,
    prevalence_target = protocol$prevalence_target,
    num_threads = protocol$num_threads, seed = protocol$seed, verbose = verbose
  )
  add_receipt("03_fit",
    .cast_digest(list(split$train, screen$selected, protocol$models,
                      protocol$rf_ntree, protocol$rf_mtry,
                      protocol$maxent_classes,
                      protocol$maxent_regmult, protocol$tune,
                      protocol$tune_folds, protocol$prevalence_target,
                      protocol$seed)),
    fit)

  # ---- 4. Nested spatial cross-validation (training split only) -----------
  if (verbose) cli::cli_h2("CV: nested spatial cross-validation (k = {protocol$cv_k})")
  stage_started <- Sys.time()
  select_args <- switch(protocol$select_method,
    screen = list(
      corcut = protocol$select_corcut,
      ncov = protocol$select_ncov,
      maxncov = protocol$select_maxncov,
      keep = protocol$select_keep,
      ensure = select_ensure,
      var_id = select_var_id,
      category = select_category
    ),
    full = list(keep = protocol$select_keep)
  )
  cv <- cast_cv(
    split$train, screen = screen,
    select_method = protocol$select_method, select_args = select_args,
    k = protocol$cv_k, models = protocol$models,
    block_method = protocol$block_method, block_size = protocol$block_size,
    buffer = protocol$buffer, buffer_unit = protocol$buffer_unit,
    n_repeat = protocol$n_repeat, tune = protocol$tune,
    tune_folds = protocol$tune_folds,
    threshold_method = protocol$threshold_method,
    aoa = protocol$aoa, prevalence_target = protocol$prevalence_target,
    rf_ntree = protocol$rf_ntree, rf_mtry = protocol$rf_mtry,
    maxent_classes = protocol$maxent_classes,
    maxent_regmult = protocol$maxent_regmult,
    num_threads = protocol$num_threads, seed = protocol$seed, verbose = verbose
  )
  cv_tables <- list(
    cv_metrics = cv$metrics,
    cv_fold_metrics = cv$fold_metrics,
    cv_selection_frequency = cv$selection_freq,
    cv_oof_predictions = cv$oof,
    cv_fold_status = data.frame(
      fold = seq_along(cv$fold_status), status = cv$fold_status,
      selected = vapply(cv$selections %||% vector("list", length(cv$fold_status)),
                        paste, character(1), collapse = ";"),
      stringsAsFactors = FALSE)
  )
  add_receipt("04_cv",
    .cast_digest(list(split$train, screen$selected, select_args,
                      protocol$cv_k, protocol$n_repeat, protocol$block_method,
                      protocol$block_size, protocol$buffer, protocol$buffer_unit,
                      protocol$tune, protocol$tune_folds,
                      protocol$threshold_method, protocol$aoa,
                      protocol$prevalence_target, protocol$seed)),
    cv, tables = cv_tables)

  # ---- 5. One-shot hold-out evaluation (out-of-fold thresholds) -----------
  if (verbose) cli::cli_h2("Evaluate: hold-out split, out-of-fold thresholds")
  stage_started <- Sys.time()
  eval_result <- cast_evaluate(
    fit, split$test,
    threshold_method = protocol$threshold_method,
    thresholds = cv$thresholds
  )
  if (verbose) print(eval_result)
  add_receipt("05_evaluate",
    .cast_digest(list(fit, split$test, cv$thresholds,
                      protocol$threshold_method)),
    eval_result, tables = list(holdout_metrics = eval_result$metrics))

  # ---- 6. Full-data refit (tuned configuration frozen) --------------------
  if (verbose) cli::cli_h2("Refit: full-data models for publication")
  stage_started <- Sys.time()
  tuned <- .cast_tuned_args(fit)
  fit_full <- cast_fit(
    species_data, screen = screen, models = protocol$models,
    rf_ntree = protocol$rf_ntree,
    rf_mtry = tuned$rf_mtry %||% protocol$rf_mtry,
    maxent_classes = tuned$maxent_classes %||% protocol$maxent_classes,
    maxent_regmult = tuned$maxent_regmult %||% protocol$maxent_regmult,
    prevalence_target = protocol$prevalence_target,
    num_threads = protocol$num_threads, seed = protocol$seed, verbose = verbose
  )
  add_receipt("06_refit",
    .cast_digest(list(species_data, screen$selected, tuned, protocol$models,
                      protocol$prevalence_target, protocol$seed)),
    fit_full)

  # ---- 7. Prediction and ensemble (df) / projection (raster) --------------
  predict_result <- NULL
  ensemble_result <- NULL
  projection <- NULL
  if (!is.null(env_data)) {
    if (verbose) cli::cli_h2("Predict: per-site suitability and ensemble")
    stage_started <- Sys.time()
    predict_result <- cast_predict(
      fit_full, env_data, clamp = protocol$clamp,
      cv = cv, aoa = protocol$aoa
    )
    ensemble_result <- cast_ensemble(
      fit_full, cv, env_data,
      method = protocol$ensemble_method,
      min_score = protocol$ensemble_min_score,
      decay = protocol$ensemble_decay,
      fallback = protocol$ensemble_fallback,
      min_metric_folds = protocol$ensemble_min_metric_folds
    )
    ens_summary <- data.frame(
      method = ensemble_result$method,
      threshold = ensemble_result$threshold,
      weights = paste(names(ensemble_result$weights),
                      format(ensemble_result$weights, digits = 4),
                      sep = "=", collapse = "; "),
      weights_fallback = ensemble_result$weights_fallback,
      n = nrow(ensemble_result$predictions),
      mean_hss = mean(ensemble_result$predictions$hss_ensemble, na.rm = TRUE),
      binary_share = mean(ensemble_result$predictions$binary_ensemble, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    add_receipt("07_predict",
      .cast_digest(list(fit_full, env_data, protocol$clamp, protocol$aoa)),
      predict_result)
    add_receipt("08_ensemble",
      .cast_digest(list(fit_full, cv$metrics, env_data,
                        protocol$ensemble_method,
                        protocol$ensemble_min_score, protocol$ensemble_decay,
                        protocol$ensemble_fallback,
                        protocol$ensemble_min_metric_folds)),
      ensemble_result, tables = list(ensemble_summary = ens_summary))
  }
  if (!is.null(env_raster)) {
    if (verbose) cli::cli_h2("Project: current and future raster products")
    stage_started <- Sys.time()
    raster_dir <- file.path(output_dir, "rasters")
    projection <- if (is.null(future_rasters)) {
      current_er <- cast_ensemble_raster(
        fit_full, cv, env_raster, output_dir = raster_dir,
        prefix = "current",
        method = protocol$ensemble_method,
        min_score = protocol$ensemble_min_score,
        decay = protocol$ensemble_decay,
        fallback = protocol$ensemble_fallback,
        min_metric_folds = protocol$ensemble_min_metric_folds,
        clamp = protocol$clamp,
        overwrite = TRUE, verbose = verbose
      )
      list(current = current_er, future = list(), stats = data.frame(),
           output_dir = output_dir)
    } else {
      cast_project_raster(
        fit_full, cv,
        current_raster = env_raster, future_rasters = future_rasters,
        output_dir = file.path(output_dir, "rasters"),
        mask = NULL,
        method = protocol$ensemble_method,
        min_score = protocol$ensemble_min_score,
        decay = protocol$ensemble_decay,
        fallback = protocol$ensemble_fallback,
        min_metric_folds = protocol$ensemble_min_metric_folds,
        clamp = protocol$clamp,
        overwrite = TRUE, verbose = verbose
      )
    }
    raster_outputs <- unlist(lapply(
      c(list(projection$current), projection$future),
      function(x) unlist(x[grep("_path$", names(x))], use.names = FALSE)
    ))
    add_receipt("07_project",
      .cast_digest(list(fit_full, cv$metrics, names(future_rasters),
                        protocol$ensemble_method,
                        protocol$ensemble_min_score, protocol$ensemble_decay,
                        protocol$ensemble_fallback,
                        protocol$ensemble_min_metric_folds, protocol$clamp)),
      list(stats = projection$stats, output_dir = projection$output_dir),
      tables = if (nrow(projection$stats)) {
        list(projection_stats = projection$stats)
      } else list(),
      extra_files = raster_outputs)
  }

  # ---- Result, report and manifest ----------------------------------------
  result <- new_cast_result(
    screen = screen, fit = fit, eval = eval_result, cv = cv,
    predict = predict_result, ensemble = ensemble_result,
    fit_full = fit_full, call = cl, protocol = protocol,
    projection = projection, output_dir = output_dir, receipts = receipts
  )
  if (!is.null(out)) {
    stage_started <- Sys.time()
    report_path <- file.path(output_dir, "odmap_report.md")
    cast_report_odmap(result, path = report_path, meta = report_meta,
                      verbose = FALSE)
    result$report_path <- report_path
    manifest <- list(
      cast_version = as.character(utils::packageVersion("castSDM")),
      protocol = protocol,
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      receipts = receipts,
      finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      session = utils::capture.output(utils::sessionInfo())
    )
    saveRDS(manifest, file.path(output_dir, "manifest.rds"))
    receipt_rows <- do.call(rbind, lapply(receipts, function(r) {
      data.frame(stage = r$stage, inputs = r$inputs,
                 outputs = nrow(r$outputs), elapsed_secs = r$elapsed_secs,
                 finished = r$finished, stringsAsFactors = FALSE)
    }))
    receipt_rows <- rbind(receipt_rows, data.frame(
      stage = "09_report", inputs = .cast_digest(report_meta), outputs = 1L,
      elapsed_secs = round(as.numeric(Sys.time() - stage_started,
                                      units = "secs"), 2),
      finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      stringsAsFactors = FALSE))
    utils::write.csv(receipt_rows,
                     file.path(out$tables, "stage_receipts.csv"),
                     row.names = FALSE)
    if (verbose) cli::cli_inform("Artifacts and manifest under {.path {output_dir}}.")
  }
  if (verbose) cli::cli_h1("castSDM run complete")
  result
}

#' Initialise the cast_run Output Directories
#'
#' @param output_dir Optional directory.
#' @return `NULL`, or a list of subdirectory paths.
#' @keywords internal
#' @noRd
.cast_run_dirs <- function(output_dir) {
  if (is.null(output_dir)) return(NULL)
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    cli::cli_abort("{.arg output_dir} must be a single directory path.")
  }
  dirs <- list(
    root = output_dir,
    stages = file.path(output_dir, "stages"),
    tables = file.path(output_dir, "tables")
  )
  for (d in dirs) {
    if (!dir.exists(d) && !dir.create(d, recursive = TRUE)) {
      cli::cli_abort("Cannot create output directory {.path {d}}.")
    }
  }
  dirs
}
