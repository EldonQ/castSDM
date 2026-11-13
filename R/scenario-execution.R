# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# Shared execution contracts; no scientific defaults depend on output paths.
.cast_code_identity <- function(fun) {
  # AST serialization can change when R JIT compiles embedded closures.
  # Source text is stable across compilation and carries no source-file env.
  paste(deparse(fun, width.cutoff = 500L), collapse = "\n")
}

.cast_lock <- function(path, timeout = 0) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(path, timeout = timeout)
  if (is.null(lock)) cli::cli_abort("An active process holds {.path {path}}. Wait for that run before writing here.")
  lock
}
.cast_resources <- function(resources = list()) {
  allowed <- c("cores", "workers", "memory_mb", "scratch_dir")
  if (!is.list(resources) || (length(resources) && (is.null(names(resources)) ||
      anyNA(names(resources)) || any(!nzchar(names(resources))) || anyDuplicated(names(resources)))) ||
      length(setdiff(names(resources), allowed)))
    cli::cli_abort("resources must be a named list of cores, workers, memory_mb and scratch_dir.")
  cores <- .cast_check_scalar(resources$cores %||% 1L, "resources$cores", lower = 1, integer = TRUE)
  workers <- .cast_check_scalar(resources$workers %||% min(cores, 4L), "resources$workers", lower = 1, integer = TRUE)
  if (workers > cores) cli::cli_abort("resources$workers cannot exceed resources$cores.")
  memory <- .cast_check_scalar(resources$memory_mb %||% 1024, "resources$memory_mb", lower = 32)
  scratch <- resources$scratch_dir %||% tempdir()
  if (!is.character(scratch) || length(scratch) != 1L || is.na(scratch) || !dir.exists(scratch))
    cli::cli_abort("resources$scratch_dir must be an existing writable directory.")
  if (file.access(scratch, 2L) != 0L) cli::cli_abort("The scratch directory is not writable.")
  list(cores = cores, workers = workers, threads = max(1L, as.integer(cores/workers)),
       memory_mb = memory, worker_memory_mb = memory/workers, scratch_dir = scratch)
}

.cast_seal_callback <- function(fun, depth = 0L) {
  if (depth > 5L) cli::cli_abort("Recursive external callback dependencies cannot be sent to workers.")
  env <- environment(fun)
  if (is.primitive(fun) || isNamespace(env) || identical(env, baseenv())) return(fun)
  target <- new.env(parent = baseenv())
  for (name in codetools::findGlobals(fun, merge = TRUE)) {
    if (!exists(name, env, inherits = TRUE)) next
    value <- get(name, env, inherits = TRUE)
    if (is.function(value)) value <- .cast_seal_callback(value, depth + 1L)
    assign(name, value, target)
  }
  environment(fun) <- target
  fun
}

.cast_atomic_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(".cast-stage-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  saveRDS(value, temp, version = 3)
  .cast_publish_raster(temp, path, overwrite = TRUE)
  invisible(path)
}

.cast_cache_read <- function(path, key) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(x) || !identical(x$schema, "cast_cache/2") || !identical(x$key, key) ||
      !identical(x$digest, .cast_digest(x$value))) return(NULL)
  if (length(x$files)) {
    files <- file.path(dirname(path), names(x$files))
    if (any(!file.exists(files)) || !identical(unname(tools::md5sum(files)), unname(x$files))) return(NULL)
  }
  x$value
}

.cast_cache_write <- function(value, path, key, files = character()) {
  if (is.null(path)) return(invisible(value))
  # Artifacts are siblings of the receipt. Relative paths allow relocation.
  if (length(files) && any(dirname(normalizePath(files, winslash = "/")) !=
      normalizePath(dirname(path), winslash = "/", mustWork = FALSE)))
    cli::cli_abort("Cached artifacts must share their receipt directory.")
  hashes <- if (length(files)) stats::setNames(unname(tools::md5sum(files)), basename(files)) else character()
  .cast_atomic_rds(list(schema = "cast_cache/2", key = key, value = value,
    digest = .cast_digest(value), files = hashes), path)
  invisible(value)
}

.cast_fit_identity <- function(fit) {
  package <- as.character(utils::packageVersion("castSDM"))
  # Reviewed maintenance/interface releases keep 0.28.0 prediction/cache
  # semantics. Later releases invalidate conservatively unless reviewed here.
  if (package %in% c("0.28.1", "0.28.2", "0.28.3", "0.29.0", "0.29.1")) package <- "0.28.0"
  models <- lapply(fit$models, function(m) {
    if (is.function(m$predict_function)) m$predict_function <- .cast_callback_identity(m$predict_function)
    m
  })
  engines <- unique(vapply(names(models), function(m) .cast_engine_package(m) %||% "", character(1)))
  versions <- vapply(engines[nzchar(engines)], function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "missing"
  }, character(1))
  .cast_digest(list(package = package,
    models = models, env_vars = fit$env_vars, scaling = fit$scaling,
    calibration = fit$calibration, engines = versions, predictor = .cast_code_identity(predict_single_model),
    contrast = .cast_code_identity(.cast_factorial_values), context = .cast_code_identity(.cast_prediction_context),
    gam = .cast_code_identity(.cast_gam_terms), offset = .cast_code_identity(.cast_apply_offset),
    clamp = .cast_code_identity(.cast_clamp), schema = "paired_scenarios/1"))
}

.cast_data_identity <- function(x, memo = new.env(parent = emptyenv())) {
  if (!inherits(x, "SpatRaster")) return(.cast_digest(x))
  geometry <- list(dim = dim(x), extent = as.vector(terra::ext(x)), crs = terra::crs(x), names = names(x))
  src <- terra::sources(x, bands = TRUE)
  if (is.data.frame(src) && "source" %in% names(src) && all(nzchar(src$source)) && all(file.exists(src$source))) {
    hashes <- vapply(src$source, function(p) {
      p <- normalizePath(p, winslash = "/")
      if (!exists(p, memo, inherits = FALSE)) assign(p, unname(tools::md5sum(p)), memo)
      get(p, memo, inherits = FALSE)
    }, character(1))
    return(.cast_digest(list(geometry, content = unname(hashes), bands = src$bands)))
  }
  # In-memory data are hashed in bounded stripes, including spatial ordering.
  rows <- max(1L, as.integer(100000L/terra::ncol(x)))
  parts <- vapply(seq.int(1L, terra::nrow(x), by = rows), function(r) {
    .cast_digest(terra::values(x, row = r, nrows = min(rows, terra::nrow(x)-r+1L), mat = TRUE))
  }, character(1))
  .cast_digest(list(geometry, parts))
}

.cast_gam_terms <- function(info) {
  if (!identical(info$name, "gam") || identical(info$type, "external")) return(NULL)
  if (!requireNamespace("mgcv", quietly = TRUE)) return(NULL)
  model <- info$model
  if (!inherits(model, "gam") || !identical(model$family$link, "logit") ||
      length(attr(model$pterms, "term.labels")) || length(attr(model$terms, "offset")) ||
      !is.null(model$call$offset) || !length(model$smooth)) return(NULL)
  good <- vapply(model$smooth, function(s) length(s$term) == 1L &&
    identical(s$by, "NA") && !inherits(s, "tensor.smooth"), logical(1))
  if (!all(good)) return(NULL)
  vars <- vapply(model$smooth, function(s) s$term, character(1))
  if (anyDuplicated(vars)) return(NULL)
  stats::setNames(vapply(model$smooth, function(s) s$label, character(1)), vars)
}

.cast_prediction_context <- function(fit, weights, current, future, clamp = FALSE) {
  reference <- fit$scaling$reference
  raw_current <- current; raw_future <- future
  if (clamp && !is.null(reference)) {
    current <- .cast_clamp(current, reference); future <- .cast_clamp(future, reference)
  }
  active <- names(weights)
  terms <- stats::setNames(lapply(active, function(m) .cast_gam_terms(fit$models[[m]])), active)
  prepared <- lapply(active, function(m) {
    if (is.null(terms[[m]])) return(NULL)
    info <- fit$models[[m]]
    a <- stats::predict(info$model, current, type = "terms", newdata.guaranteed = TRUE)
    b <- stats::predict(info$model, future, type = "terms", newdata.guaranteed = TRUE)
    list(eta = rowSums(a) + attr(a, "constant"), delta = b - a)
  })
  names(prepared) <- active
  predict <- function(switched = character(), base = NULL, changed_rows = NULL) {
    X <- current
    X[switched] <- future[switched]
    out <- matrix(NA_real_, nrow(current), length(active), dimnames = list(NULL, active))
    idx <- changed_rows %||% rep(TRUE, nrow(X))
    if (!is.null(base)) out[,] <- base
    for (m in active) {
      info <- fit$models[[m]]
      if (!any(idx)) next
      if (!is.null(prepared[[m]])) {
        p <- prepared[[m]]; cols <- unname(terms[[m]][intersect(switched, names(terms[[m]]))])
        eta <- p$eta[idx] + if (length(cols)) rowSums(p$delta[idx, cols, drop = FALSE]) else 0
        out[idx, m] <- .cast_apply_offset(info$model$family$linkinv(eta), info$offset)
      } else {
        if (any(idx)) {
          p <- tryCatch(predict_single_model(info, X[idx, , drop = FALSE]),
            error = function(e) cli::cli_abort("Prediction failed for {.val {m}}: {conditionMessage(e)}"))
          if (!is.numeric(p) || length(p) != sum(idx)) cli::cli_abort("Model {.val {m}} returned invalid prediction dimensions.")
          out[idx, m] <- p
        }
      }
    }
    out
  }
  list(predict = predict, current = current, future = future,
       raw_current = raw_current, raw_future = raw_future)
}
