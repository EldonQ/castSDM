# AI-assisted maintenance (2026-11-08); scope and verification: inst/AI_ASSISTANCE.md.
# Worker bootstrap transports values; compatibility is checked in the caller.
.cast_compare_worker_identity <- function() {
  ns <- asNamespace("castSDM")
  list(version = as.character(getNamespaceVersion(ns)),
    path = normalizePath(getNamespaceInfo(ns, "path"), winslash = "/", mustWork = TRUE),
    R = as.character(getRversion()),
    code = .cast_code_identity(get(".cast_compare_one", envir = ns, inherits = FALSE)))
}

.cast_check_compare_worker <- function(expected, actual, worker) {
  fields <- c("version", "path", "R", "code")
  changed <- fields[!vapply(fields, function(field)
    identical(expected[[field]], actual[[field]]), logical(1))]
  if (!length(changed)) return(invisible(TRUE))
  labels <- c(version = "package version", path = "package path", R = "R version", code = "core scenario code")
  label <- paste(unname(labels[changed]), collapse = ", ")
  cli::cli_abort(c(
    "Scenario worker {worker} differs in {label}.",
    i = "Caller: castSDM {expected$version}, R {expected$R}, {.path {expected$path}}.",
    i = "Worker: castSDM {actual$version}, R {actual$R}, {.path {actual$path}}.",
    i = "Reload the calling session from the intended installed package; in-memory modifications cannot be sent as an installed namespace."))
}

.cast_start_compare_cluster <- function(workers) {
  # Materialize the loaded caller's identity before crossing the socket.
  # Reading an installed DESCRIPTION alone cannot identify a loaded namespace.
  expected <- .cast_compare_worker_identity()
  libraries <- .libPaths()
  cl <- parallel::makePSOCKcluster(workers, rscript_args = "--vanilla")
  ready <- FALSE
  on.exit(if (!ready) parallel::stopCluster(cl), add = TRUE)
  initialize <- function(libs, package_path) {
    preloaded <- "castSDM" %in% loadedNamespaces()
    library <- dirname(package_path)
    .libPaths(unique(c(library, libs)))
    ns <- loadNamespace("castSDM", lib.loc = library)
    functions <- c(".cast_code_identity", ".cast_compare_one")
    available <- all(vapply(functions, exists, logical(1), envir = ns, inherits = FALSE))
    code <- if (available) get(functions[1L], envir = ns)(get(functions[2L], envir = ns)) else NULL
    list(version = as.character(getNamespaceVersion(ns)),
      path = normalizePath(getNamespaceInfo(ns, "path"), winslash = "/", mustWork = TRUE),
      R = as.character(getRversion()), code = code, preloaded = preloaded)
  }
  # Source/bytecode metadata can itself refer to a package namespace. Strip it
  # before serialization so the worker cannot load a default-library copy
  # while deserializing this closure, before .libPaths/lib.loc take effect.
  initialize <- utils::removeSource(initialize)
  environment(initialize) <- baseenv()
  actual <- parallel::clusterCall(cl, initialize, libraries, expected$path)
  for (i in seq_along(actual)) .cast_check_compare_worker(expected, actual[[i]], i)
  ready <- TRUE
  cl
}
