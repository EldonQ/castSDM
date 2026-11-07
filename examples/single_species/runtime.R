# Runtime setup shared by the example entry points. No fitting or raster work.
case_runtime_setup <- function(root, packages, minimum_cast = NULL) {
  override <- Sys.getenv("CAST_R_LIBRARY", "")
  if (nzchar(override) && !dir.exists(override))
    stop("CAST_R_LIBRARY does not exist: ", override, call. = FALSE)
  candidates <- c(override, file.path(root, ".library"))
  ancestor <- root
  repeat {
    candidates <- c(candidates, file.path(ancestor, ".rlib_cast"))
    versions <- file.path(ancestor, ".rlib_cast_versions")
    if (dir.exists(versions)) candidates <- c(candidates,
      list.dirs(versions, recursive = FALSE, full.names = TRUE))
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) break
    ancestor <- parent
  }
  libraries <- unique(candidates[nzchar(candidates) & dir.exists(candidates)])
  # Do not overwrite a package library used by a running R process. An
  # explicitly selected library retains precedence over side-by-side versions.
  if (!nzchar(override) && !dir.exists(file.path(root, ".library")) && length(libraries)) {
    with_cast <- libraries[file.exists(file.path(libraries, "castSDM", "DESCRIPTION"))]
    if (length(with_cast)) {
      newest <- Reduce(function(a, b) {
        av <- read.dcf(file.path(a, "castSDM", "DESCRIPTION"), fields = "Version")[[1L]]
        bv <- read.dcf(file.path(b, "castSDM", "DESCRIPTION"), fields = "Version")[[1L]]
        if (utils::compareVersion(bv, av) > 0) b else a
      }, with_cast)
      libraries <- c(newest, setdiff(libraries, newest))
    }
  }
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

  # castSDM is a pure-R package. Switch just this namespace when an old/dev
  # copy is loaded; do not require an interactive R restart for a path change.
  reattach_cast <- FALSE
  if ("castSDM" %in% packages && "castSDM" %in% loadedNamespaces()) {
    selected <- installed$castSDM
    active_path <- normalizePath(getNamespaceInfo("castSDM", "path"), winslash = "/")
    active_version <- as.character(getNamespaceVersion("castSDM"))
    code_files <- file.path("R", paste0("castSDM.", c("rdb", "rdx")))
    active_files <- file.path(active_path, code_files)
    selected_files <- file.path(selected$path, code_files)
    same_code <- identical(active_version, selected$version) &&
      (identical(active_path, selected$path) ||
        (all(file.exists(c(active_files, selected_files))) &&
          identical(unname(tools::md5sum(active_files)), unname(tools::md5sum(selected_files)))))
    if (!same_code) {
      users <- getNamespaceUsers("castSDM")
      if (length(users)) stop("Cannot switch castSDM while it is imported by: ",
        paste(users, collapse = ", "),
        ". Run 00_RUN_ALL.R with Rscript in a separate R process.", call. = FALSE)
      has_native <- any(vapply(getLoadedDLLs(), function(dll)
        startsWith(normalizePath(dll[["path"]], winslash = "/", mustWork = FALSE),
          paste0(active_path, "/")), logical(1)))
      if (has_native)
        stop("The loaded castSDM has native code; run the entry with Rscript in a separate R process.", call. = FALSE)
      message("Switching castSDM ", active_version, " (", active_path,
        ") -> ", selected$version, " (", selected$path, ")")
      reattach_cast <- "package:castSDM" %in% search()
      if ("pkgload" %in% loadedNamespaces() && pkgload::is_dev_package("castSDM")) {
        pkgload::unload("castSDM", quiet = TRUE)
      } else {
        if (reattach_cast) detach("package:castSDM", unload = FALSE, character.only = TRUE)
        unloadNamespace("castSDM")
      }
    }
  }
  # Leave already loaded native dependencies in place. A genuine version
  # mismatch matters for workers, whereas a same-version path difference does not.
  conflicting <- vapply(intersect(packages, loadedNamespaces()), function(package) {
    selected <- installed[[package]]
    !identical(as.character(getNamespaceVersion(package)), selected$version)
  }, logical(1))
  if (any(conflicting)) stop("Loaded dependency versions differ from the project library: ",
    paste(names(conflicting)[conflicting], collapse = ", "),
    ". Run the entry with Rscript in a separate R process.", call. = FALSE)
  failures <- vapply(packages, function(package) {
    tryCatch({ loadNamespace(package); "" }, error = function(e) conditionMessage(e))
  }, character(1))
  if (any(nzchar(failures))) stop("Cannot load required packages:\n",
    paste(paste(names(failures)[nzchar(failures)], failures[nzchar(failures)], sep = ": "),
      collapse = "\n"), call. = FALSE)
  if (reattach_cast) library("castSDM", lib.loc = dirname(installed$castSDM$path), character.only = TRUE)
  versions <- data.frame(package = packages,
    version = vapply(packages, function(package) as.character(getNamespaceVersion(package)), character(1)),
    path = vapply(packages, function(package)
      normalizePath(getNamespaceInfo(package, "path"), winslash = "/"), character(1)), row.names = NULL)
  if ("castSDM" %in% packages)
    message("castSDM ", getNamespaceVersion("castSDM"), " | R library: ",
      dirname(getNamespaceInfo("castSDM", "path")))
  invisible(list(libraries = .libPaths(), dependencies = versions))
}
