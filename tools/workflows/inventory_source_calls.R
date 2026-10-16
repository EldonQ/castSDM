.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
root <- "E:/Package/cast"
files <- list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)
functions <- list()
for (f in files) {
  expressions <- parse(f, keep.source = TRUE)
  refs <- attr(expressions, "srcref")
  for (i in seq_along(expressions)) {
    e <- expressions[[i]]
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        is.symbol(e[[2]]) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
      name <- as.character(e[[2]])
      functions[[name]] <- list(file = basename(f), line = as.integer(refs[[i]][1]),
        arguments = length(e[[3]][[2]]), body = e[[3]][[3]])
    }
  }
}
calls <- function(expr) {
  if (!is.call(expr) && !is.expression(expr) && !is.pairlist(expr)) return(character())
  result <- if (is.call(expr) && is.symbol(expr[[1]])) as.character(expr[[1]]) else character()
  items <- as.list(expr)
  for (i in seq_along(items)) {
    if (identical(items[[i]], quote(expr = ))) next
    result <- c(result, calls(items[[i]]))
  }
  unique(result)
}
edges <- do.call(rbind, lapply(names(functions), function(n) {
  targets <- intersect(calls(functions[[n]]$body), names(functions))
  if (length(targets)) data.frame(caller = n, callee = targets) else NULL
}))
namespace <- readLines(file.path(root, "NAMESPACE"))
public <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))
methods <- grep("^S3method\\(", namespace, value = TRUE)
methods <- sub("^S3method\\(([^,]+),([^\\)]+)\\).*$", "\\1.\\2", methods)
workflow_files <- c(list.files("E:/Package/castSDM_workflows/single_species", "[.]R$", full.names = TRUE),
  list.files("E:/Package/castSDM_workflows/R", "[.]R$", full.names = TRUE))
workflow_text <- paste(unlist(lapply(workflow_files, readLines, warn = FALSE)), collapse = "\n")
catalog <- do.call(rbind, lapply(names(functions), function(n) {
  entry <- functions[[n]]
  data.frame(function_name = n, file = entry$file, line = entry$line,
    arguments = entry$arguments, exported = n %in% public, S3_method = n %in% methods,
    internal_callers = sum(edges$callee == n), internal_callees = sum(edges$caller == n),
    mentioned_in_case_workflow = grepl(n, workflow_text, fixed = TRUE))
}))
output <- file.path(root, "tools/workflows")
utils::write.csv(catalog, file.path(output, "function_call_inventory.csv"), row.names = FALSE)
utils::write.csv(edges, file.path(output, "internal_call_edges.csv"), row.names = FALSE)
writeLines(c("Static direct-call inventory; reflective get(), dispatch and serialized callbacks are not resolved.",
  "Zero static callers does not prove unused code. Exported functions and S3 methods remain public contracts.",
  "Inspect tests, examples and downstream workflows before removing any candidate."),
  file.path(output, "CALL_INVENTORY_LIMITS.txt"))
print(catalog[!catalog$exported & !catalog$S3_method & catalog$internal_callers == 0, ])
