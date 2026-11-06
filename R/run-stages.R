# Training dependencies are tracked independently of projection and plotting.
.cast_run_stage_identity <- function(stage, models = character()) {
  roots <- switch(stage,
    `01_prepare` = c("cast_prepare", ".cast_verify_prepared"), `02_select` = "cast_select",
    `03_fit` = "cast_fit", `04_cv` = "cast_cv",
    `05_evaluate` = "cast_evaluate", `06_refit` = "cast_fit",
    cli::cli_abort("Unknown training stage {.val {stage}}."))
  ns <- asNamespace("castSDM")
  seen <- character(); code <- list()
  visit <- function(name) {
    if (name %in% seen || !exists(name, ns, inherits = FALSE)) return(invisible(NULL))
    fun <- get(name, ns, inherits = FALSE)
    if (!is.function(fun) || !identical(environment(fun), ns)) return(invisible(NULL))
    seen <<- c(seen, name)
    code[[name]] <<- .cast_code_identity(fun)
    for (dep in codetools::findGlobals(fun, merge = TRUE)) visit(dep)
    invisible(NULL)
  }
  for (root in roots) visit(root)
  packages <- switch(stage,
    `01_prepare` = c("blockCV", "sf"), `02_select` = "glmnet",
    `04_cv` = c("blockCV", "sf", "FNN", "glmnet", "pROC"),
    `05_evaluate` = "pROC", character())
  if (stage %in% c("03_fit", "04_cv", "06_refit"))
    packages <- c(packages, vapply(models, .cast_engine_package, character(1)))
  versions <- vapply(unique(packages), function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not_installed"
  }, character(1))
  list(code = code[sort(names(code))], versions = versions, R = as.character(getRversion()))
}

.cast_verify_prepared <- function(data, prepared) {
  required <- c("train", "test", "train_index", "test_index", "env_vars", "split")
  if (!is.list(prepared) || !all(required %in% names(prepared)))
    cli::cli_abort("prepared must be a current cast_prepare() split with row indices.")
  a <- prepared$train_index; b <- prepared$test_index
  indices <- c(a, b)
  if (!is.numeric(indices) || anyNA(indices) || anyDuplicated(indices) ||
      !setequal(indices, seq_len(nrow(data))) || !length(a) || !length(b))
    cli::cli_abort("prepared row indices must partition species_data exactly once.")
  for (part in c("train", "test")) {
    frame <- prepared[[part]]
    idx <- prepared[[paste0(part, "_index")]]
    if (!is.data.frame(frame) || !all(names(frame) %in% names(data)) ||
        !identical(frame, data[idx, names(frame), drop = FALSE]))
      cli::cli_abort("prepared {part} values do not match species_data at the recorded indices.")
  }
  prepared
}
