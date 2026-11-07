# Runtime setup shared by the example entry points. No fitting or raster work.
case_runtime_setup <- function(root, packages, minimum_cast = NULL) {
  override <- Sys.getenv("CAST_R_LIBRARY", "")
  if (nzchar(override) && !dir.exists(override))
    stop("CAST_R_LIBRARY does not exist: ", override, call. = FALSE)
  candidates <- c(override, file.path(root, ".library"))
  ancestor <- root
  repeat {
    candidates <- c(candidates, file.path(ancestor, ".rlib_cast"))
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) break
    ancestor <- parent
  }
  libraries <- unique(candidates[nzchar(candidates) & dir.exists(candidates)])
  if (length(libraries)) .libPaths(c(libraries, .libPaths()))

  packages <- unique(packages)
  minimum <- c(castSDM = minimum_cast, filelock = "1.0.2", digest = "0.6.30",
               terra = "1.7.0", sf = "1.0.0")
  installed <- lapply(packages, function(package) {
    path <- find.package(package, lib.loc = .libPaths(), quiet = TRUE)
    if (!length(path)) return(NULL)
    description <- read.dcf(file.path(path, "DESCRIPTION"), fields = "Version")
    list(path = normalizePath(path, winslash = "/"), version = description[[1L]])
  })
  names(installed) <- packages
  unusable <- vapply(packages, function(package) {
    info <- installed[[package]]
    is.null(info) || (package %in% names(minimum) &&
      utils::compareVersion(info$version, minimum[[package]]) < 0)
  }, logical(1))
  if (any(unusable)) {
    missing <- packages[unusable]
    library <- .libPaths()[1L]
    cran <- setdiff(missing, "castSDM")
    instructions <- if (length(cran)) paste0("install.packages(c(",
      paste(sprintf('"%s"', cran), collapse = ", "), '), lib = "', library,
      '", repos = "https://cloud.r-project.org")') else character()
    if ("castSDM" %in% missing) instructions <- c(instructions,
      'remotes::install_github("EldonQ/castSDM", lib = .libPaths()[1L])')
    stop("Missing/outdated packages: ", paste(missing, collapse = ", "),
      "\nR libraries: ", paste(.libPaths(), collapse = "; "),
      "\nInstall in this R session, then source the entry again:\n",
      paste(instructions, collapse = "\n"), call. = FALSE)
  }

  # Do not mix a previously loaded development/old castSDM with worker installs.
  conflicting <- vapply(intersect(packages, loadedNamespaces()), function(package) {
    selected <- installed[[package]]
    different_version <- !identical(as.character(getNamespaceVersion(package)), selected$version)
    different_cast <- package == "castSDM" && !identical(
      normalizePath(getNamespaceInfo(package, "path"), winslash = "/"), selected$path)
    different_version || different_cast
  }, logical(1))
  if (any(conflicting)) stop("Restart R, then source the entry again. This session has a different loaded copy of: ",
    paste(names(conflicting)[conflicting], collapse = ", "), call. = FALSE)
  failures <- vapply(packages, function(package) {
    tryCatch({ loadNamespace(package); "" }, error = function(e) conditionMessage(e))
  }, character(1))
  if (any(nzchar(failures))) stop("Cannot load required packages:\n",
    paste(paste(names(failures)[nzchar(failures)], failures[nzchar(failures)], sep = ": "),
      collapse = "\n"), call. = FALSE)
  versions <- data.frame(package = packages,
    version = vapply(installed, `[[`, character(1), "version"),
    path = vapply(installed, `[[`, character(1), "path"), row.names = NULL)
  if ("castSDM" %in% packages)
    message("castSDM ", installed$castSDM$version, " | R library: ", dirname(installed$castSDM$path))
  invisible(list(libraries = .libPaths(), dependencies = versions))
}
