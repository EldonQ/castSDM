# Explicitly run when the environmental library changes; not part of each fit.
local({
  sources <- vapply(sys.frames(), function(frame)
    if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else
      if (is.character(frame$file) && length(frame$file) == 1L &&
          basename(frame$file) == "00_INVENTORY.R") frame$file else "", character(1))
  script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
  override <- Sys.getenv("CAST_EXAMPLE_ROOT", "")
  roots <- if (nzchar(override)) override else c(
    dirname(tail(sources[nzchar(sources)], 1L)), dirname(script), getwd())
  roots <- roots[file.exists(file.path(roots, "runtime.R"))]
  if (!length(roots)) stop("Cannot locate runtime.R; source 00_INVENTORY.R using its full path.")
  root <- normalizePath(roots[1L], winslash = "/", mustWork = TRUE)
  sys.source(file.path(root, "runtime.R"), envir = environment())
  case_runtime_setup(root, "terra")
  sys.source(file.path(root, "inventory.R"), envir = environment())
  case_inventory_build(root, Sys.getenv("CAST_DATA_ROOT", "E:/CHECO26/N-SDM/data"))
})
