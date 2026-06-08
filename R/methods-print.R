# Print Methods ---------------------------------------------------------------

#' @export
print.cast_dag <- function(x, ...) {
  n_edges <- if (is.null(x$edges)) 0L else nrow(x$edges)
  n_nodes <- length(x$nodes)
  sm <- x$structure_method %||% "pc"
  resp <- x$response_node
  cli::cli_h1("castSDM DAG")
  cli::cli_ul(c(
    "Method: {sm}",
    "Nodes: {n_nodes}",
    "Edges: {n_edges} (strength >= {x$strength_threshold})",
    if (isTRUE(x$metadata$dense_mb)) {
      "Dense MB fallback: response-focused screening graph"
    } else NULL,
    if (!is.null(resp)) "Response node: {resp}" else NULL,
    if (isTRUE(!is.na(x$boot_R))) "Bootstrap replicates: {x$boot_R}" else
      "Bootstrap replicates: n/a (not used for this method)",
    "Score / criterion label: {x$score}"
  ))
  invisible(x)
}

#' @export
print.cast_select <- function(x, ...) {
  n_selected <- length(x$selected)
  cli::cli_h1("castSDM Variable Selection")
  cli::cli_ul(c(
    "Selected variables: {n_selected}"
  ))
  if (!is.null(x$roles) && nrow(x$roles) > 0) {
    role_tbl <- table(x$roles$role)
    cli::cli_text("Roles: {paste(names(role_tbl), role_tbl, sep='=', collapse=', ')}")
  }
  cli::cli_text("Variables: {.val {x$selected}}")
  invisible(x)
}

#' @export
print.cast_refute <- function(x, ...) {
  cli::cli_h1("castSDM Screen Refutation")
  if (!is.null(x$settings)) {
    cli::cli_ul(c(
      "Repetitions: {x$settings$reps %||% NA_integer_}",
      "Permuted-response mean overlap: {round(x$settings$permuted_response_mean_overlap %||% NA_real_, 3)}",
      "Random common-cause hit rate: {round(x$settings$random_common_cause_hit_rate %||% NA_real_, 3)}"
    ))
  }
  if (!is.null(x$summary) && nrow(x$summary) > 0) {
    print(x$summary)
  }
  invisible(x)
}

#' @export
print.cast_fit <- function(x, ...) {
  model_names <- names(x$models)
  cli::cli_h1("castSDM Model Fit")
  cli::cli_ul(c(
    "Models: {.val {model_names}}",
    "Variables: {length(x$cast_vars)}"
  ))
  invisible(x)
}

#' @export
print.cast_eval <- function(x, ...) {
  src <- if (isTRUE(x$cv_source)) "Spatial CV" else "Hold-out test set"
  cli::cli_h1("castSDM Model Evaluation ({src})")
  print(x$metrics)
  invisible(x)
}

#' @export
print.cast_cv <- function(x, ...) {
  cli::cli_h1("castSDM Spatial Cross-Validation")
  cli::cli_ul(c(
    "Folds (k): {x$k}",
    "Block method: {x$block_method}",
    "Models: {.val {x$metrics$model}}"
  ))
  cli::cli_h2("Aggregated metrics (mean +/- SD)")
  m <- x$metrics
  for (i in seq_len(nrow(m))) {
    cli::cli_li(paste0(
      m$model[i],
      " | AUC=", round(m$auc_mean[i], 3),
      " (", round(m$auc_sd[i], 3), ")",
      " TSS=", round(m$tss_mean[i], 3),
      " CBI=", round(m$cbi_mean[i], 3)
    ))
  }
  cli::cli_h2("Optimal thresholds (max TSS)")
  for (nm in names(x$thresholds)) {
    cli::cli_li("{nm}: {round(x$thresholds[nm], 3)}")
  }
  invisible(x)
}

#' @export
print.cast_predict <- function(x, ...) {
  n_sites <- nrow(x$predictions)
  cli::cli_h1("castSDM Spatial Predictions")
  cli::cli_ul(c(
    "Sites: {n_sites}",
    "Models: {.val {x$models}}"
  ))
  invisible(x)
}

#' @export
print.cast_ensemble <- function(x, ...) {
  cli::cli_h1("castSDM Ensemble Prediction")
  cli::cli_ul(c(
    "Method: {x$method}",
    "Threshold: {round(x$threshold, 3)}",
    "Sites: {nrow(x$predictions)}"
  ))
  if (!is.null(x$weights) && length(x$weights) > 0) {
    cli::cli_text("Weights: {paste(names(x$weights), round(x$weights, 3), sep='=', collapse=', ')}")
  }
  invisible(x)
}

#' @export
print.cast_project <- function(x, ...) {
  n_scenarios <- length(x$future)
  cli::cli_h1("castSDM Future Projection")
  cli::cli_ul(c(
    "Scenarios: {n_scenarios} ({.val {names(x$future)}})",
    "Current range cells: {sum(x$current$predictions$binary_ensemble == 1, na.rm = TRUE)}"
  ))
  if (!is.null(x$stats) && nrow(x$stats) > 0) {
    cli::cli_h2("Range change summary")
    for (i in seq_len(nrow(x$stats))) {
      s <- x$stats[i, ]
      cli::cli_li("{s$scenario}: gain={s$gain} loss={s$loss} stable={s$stable} shift={round(s$centroid_shift_km, 1)}km")
    }
  }
  invisible(x)
}

#' @export
print.cast_cate <- function(x, ...) {
  n_sites <- length(unique(paste(x$effects$lon, x$effects$lat)))
  cli::cli_h1("castSDM Spatial CATE")
  cli::cli_ul(c(
    "Sites: {n_sites}",
    "Variables: {.val {x$variables}}",
    "Causal forest trees: {x$n_trees}"
  ))
  invisible(x)
}

#' @export
print.cast_result <- function(x, ...) {
  cli::cli_h1("castSDM Pipeline Result")
  cli::cli_ul(c(
    "DAG: {nrow(x$dag$edges)} edges",
    "Selected: {length(x$screen$selected)} variables",
    "Models: {.val {names(x$fit$models)}}",
    "Predictions: {if (!is.null(x$predict)) 'Yes' else 'No'}",
    "Ensemble: {if (!is.null(x$ensemble)) 'Yes' else 'No'}",
    "CATE: {if (!is.null(x$cate)) 'Yes' else 'No'}"
  ))
  invisible(x)
}

#' @export
print.cast_batch <- function(x, ...) {
  cli::cli_h1("castSDM Batch Results")
  cli::cli_ul(c(
    "Species: {.val {x$species}}",
    "Models: {.val {x$models}}"
  ))
  if (!is.null(x$output_dir)) {
    cli::cli_text("  Output: {x$output_dir}")
  }
  if (!is.null(x$species_metrics) && nrow(x$species_metrics) > 0) {
    cli::cli_h2("Summary metrics")
    print(x$species_metrics[, intersect(
      c("species", "model", "auc_mean", "tss_mean", "cbi_mean"),
      names(x$species_metrics)
    ), drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}
