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
  if (!is.list(x) || !identical(x$schema, "cast_cache/3") || !identical(x$key, key) ||
      !is.raw(x$payload) || !identical(x$digest,
        digest::digest(x$payload, algo = "sha256", serialize = FALSE))) return(NULL)
  if (length(x$files)) {
    files <- file.path(dirname(path), names(x$files))
    if (any(!file.exists(files)) || !identical(unname(tools::md5sum(files)), unname(x$files))) return(NULL)
  }
  tryCatch(unserialize(x$payload), error = function(e) NULL)
}

.cast_cache_write <- function(value, path, key, files = character()) {
  if (is.null(path)) return(invisible(value))
  # Artifacts are siblings of the receipt. Relative paths allow relocation.
  if (length(files) && any(dirname(normalizePath(files, winslash = "/")) !=
      normalizePath(dirname(path), winslash = "/", mustWork = FALSE)))
    cli::cli_abort("Cached artifacts must share their receipt directory.")
  hashes <- if (length(files)) stats::setNames(unname(tools::md5sum(files)), basename(files)) else character()
  # Verify stored bytes before reconstructing the value. Re-hashing an
  # unserialized object made integrity checks depend on the reader's locale.
  payload <- serialize(value, NULL, version = 3L)
  .cast_atomic_rds(list(schema = "cast_cache/3", key = key, payload = payload,
    digest = digest::digest(payload, algo = "sha256", serialize = FALSE), files = hashes), path)
  invisible(value)
}

# Package-generated GAM formulae contain only data columns and mgcv::s.
# Their fitting-frame environments capture mutable bookkeeping (including
# the fitted object itself), not prediction inputs. Canonicalise only this
# restricted grammar; external/custom formula dependencies remain hashed.
.cast_gam_identity <- function(x, variables) {
  if (inherits(x, "formula") || inherits(x, "terms")) {
    allowed <- c(variables, "presence", "~", "+", "s")
    if (all(all.names(x, functions = TRUE) %in% allowed))
      attr(x, ".Environment") <- asNamespace("mgcv")
  } else if (is.list(x)) {
    for (i in seq_along(x)) x[i] <- list(.cast_gam_identity(x[[i]], variables))
  }
  x
}

.cast_fit_identity <- function(fit) {
  memo <- getOption("castSDM.fit_identity_cache")
  if (is.environment(memo) && identical(memo$fit, fit)) return(memo$key)
  models <- lapply(fit$models, function(m) {
    m$fit_elapsed_s <- NULL
    m$warnings <- NULL
    if (is.function(m$predict_function)) m$predict_function <- .cast_callback_identity(m$predict_function)
    if (identical(m$type, "traditional") && identical(m$name, "gam"))
      m$model <- .cast_gam_identity(m$model, fit$env_vars)
    m
  })
  engines <- unique(vapply(names(models), function(m) .cast_engine_package(m) %||% "", character(1)))
  versions <- vapply(engines[nzchar(engines)], function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "missing"
  }, character(1))
  key <- .cast_digest(list(
    models = models, env_vars = fit$env_vars, scaling = fit$scaling,
    calibration = fit$calibration, engines = versions,
    formula_identity = .cast_code_identity(.cast_gam_identity),
    code = .cast_code_dependencies(c("predict_single_model", ".cast_prediction_context", ".cast_factorial_values")),
    schema = "paired_scenarios/3"))
  # Scoped by the public operation. Mutation is checked with identical(),
  # never with object size, a user-supplied ID or a stale saved digest.
  if (is.environment(memo)) { memo$fit <- fit; memo$key <- key }
  key
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
  # GAM and BAM both use regular, batch-invariant prediction throughout.
  if (!inherits(model, "gam") || !identical(model$family$link, "logit") ||
      length(attr(model$terms, "offset")) ||
      !is.null(model$call$offset) || !length(model$smooth)) return(NULL)
  good <- vapply(model$smooth, function(s) length(s$term) == 1L &&
    identical(s$by, "NA") && !inherits(s, "tensor.smooth"), logical(1))
  if (!all(good)) return(NULL)
  vars <- vapply(model$smooth, function(s) s$term, character(1))
  labels <- vapply(model$smooth, function(s) s$label, character(1))
  # Ordinary numeric/binary linear terms are additive too. Reject factors,
  # transformations and interactions rather than applying an invalid shortcut.
  linear <- attr(model$pterms, "term.labels")
  if (length(linear)) {
    parsed <- lapply(linear, function(x) tryCatch(str2lang(x), error = function(e) NULL))
    if (any(!vapply(parsed, is.symbol, logical(1)))) return(NULL)
    linear_vars <- vapply(parsed, as.character, character(1))
    if (any(!vapply(model$model[linear_vars], is.numeric, logical(1)))) return(NULL)
    vars <- c(vars, linear_vars); labels <- c(labels, linear)
  }
  if (anyDuplicated(vars)) return(NULL)
  stats::setNames(labels, vars)
}

.cast_prediction_context <- function(fit, weights, current, future, clamp = FALSE, gam_current = NULL,
    clamp_ranges = NULL, measure = function(phase, expr) force(expr)) {
  reference <- fit$scaling$reference
  raw_current <- current; raw_future <- future
  changed_between <- function(a, b, vars) vars[!vapply(vars, function(v) identical(a[[v]], b[[v]]), logical(1))]
  changed <- changed_between(current, future, names(current))
  if (clamp && !is.null(reference)) {
    ranges <- clamp_ranges %||% lapply(reference, range, na.rm = TRUE)
    current <- .cast_clamp(current, reference, ranges)
    future <- current
    future[changed] <- .cast_clamp(raw_future[changed], reference, ranges[intersect(changed, names(ranges))])
    changed <- changed_between(current, future, changed)
  }
  active <- names(weights)
  terms <- stats::setNames(lapply(active, function(m) .cast_gam_terms(fit$models[[m]])), active)
  prepared <- lapply(active, function(m) {
    if (is.null(terms[[m]])) return(NULL)
    info <- fit$models[[m]]
    cols <- unname(terms[[m]][intersect(changed, names(terms[[m]]))])
    saved <- gam_current[[m]]
    if (is.null(saved) || !all(cols %in% colnames(saved$terms))) {
      a <- measure("predict_gam", .cast_predict_gam(info$model, current, type = "terms"))
      saved <- list(eta = rowSums(a) + attr(a, "constant"), terms = a)
    }
    delta <- if (length(cols)) {
      b <- measure("predict_gam", .cast_predict_gam(info$model, future, type = "terms", terms = cols))
      b - saved$terms[, colnames(b), drop = FALSE]
    } else matrix(numeric(), nrow(current), 0L)
    list(eta = saved$eta, delta = delta)
  })
  names(prepared) <- active
  # A cached partial prediction cannot stand in for rows never computed.
  # Include both the environmental combination and requested row subset.
  predictions <- new.env(parent = emptyenv())
  predict <- function(switched = character(), base = NULL, changed_rows = NULL) {
    switched <- intersect(names(current), intersect(switched, changed))
    idx <- changed_rows %||% rep(TRUE, nrow(current))
    key <- .cast_digest(list(switched, idx))
    if (exists(key, predictions, inherits = FALSE)) {
      cached <- get(key, predictions, inherits = FALSE)
      if (is.null(base)) cached[!idx, ] <- NA_real_ else cached[!idx, ] <- base[!idx, , drop = FALSE]
      return(cached)
    }
    X <- current
    X[switched] <- future[switched]
    out <- matrix(NA_real_, nrow(current), length(active), dimnames = list(NULL, active))
    if (!is.null(base)) out[,] <- base
    for (m in active) {
      info <- fit$models[[m]]
      if (!any(idx)) next
      if (!is.null(prepared[[m]])) {
        p <- prepared[[m]]; cols <- intersect(unname(terms[[m]][intersect(switched, names(terms[[m]]))]), colnames(p$delta))
        eta <- p$eta[idx] + if (length(cols)) rowSums(p$delta[idx, cols, drop = FALSE]) else 0
        out[idx, m] <- .cast_apply_offset(info$model$family$linkinv(eta), info$offset)
      } else {
        p <- measure(paste0("predict_", m), tryCatch(predict_single_model(info, X[idx, , drop = FALSE]),
          error = function(e) cli::cli_abort("Prediction failed for {.val {m}}: {conditionMessage(e)}")))
        if (!is.numeric(p) || length(p) != sum(idx)) cli::cli_abort("Model {.val {m}} returned invalid prediction dimensions.")
        out[idx, m] <- p
      }
    }
    assign(key, out, predictions)
    out
  }
  list(predict = predict, current = current, future = future,
       raw_current = raw_current, raw_future = raw_future)
}
