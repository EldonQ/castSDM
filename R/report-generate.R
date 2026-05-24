#' Generate HTML Summary Report
#'
#' Creates an HTML report summarizing the castSDM pipeline results,
#' including DAG visualization, variable selection, model evaluation,
#' and optional SHAP explanations. Uses the rmarkdown package to render
#' a self-contained HTML document.
#'
#' @param result A `cast_result` object from [cast()].
#' @param output_file Character. Output file path. Default
#'   `"cast_report.html"`.
#' @param species Character. Species name for the report title.
#'   Default `"Species"`.
#' @param var_labels Optional named character vector mapping variable names
#'   to display labels.
#' @param open Logical. Open the report in a browser after generation.
#'   Default `TRUE`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return The output file path (invisibly).
#'
#' @seealso [cast()]
#'
#' @export
cast_report <- function(result,
                        output_file = "cast_report.html",
                        species = "Species",
                        var_labels = NULL,
                        open = TRUE,
                        verbose = TRUE) {
  if (!inherits(result, "cast_result")) {
    cli::cli_abort("{.arg result} must be a {.cls cast_result} object.")
  }
  check_suggested("rmarkdown", "for HTML report generation")
  check_suggested("knitr", "for HTML report generation")
  check_suggested("ggplot2", "for report figures")

  # Create a temporary Rmd from template
  rmd_template <- system.file(
    "report", "cast_report_template.Rmd", package = "castSDM"
  )

  if (rmd_template == "") {
    # Use inline template if the bundled one isn't found
    rmd_content <- .cast_report_rmd_content(species)
    rmd_tmp <- tempfile(fileext = ".Rmd")
    writeLines(rmd_content, rmd_tmp)
  } else {
    rmd_tmp <- rmd_template
  }

  # Resolve output path
  output_file <- normalizePath(output_file, mustWork = FALSE)
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (verbose) cli::cli_inform("Rendering castSDM report to {.file {output_file}}...")

  rmarkdown::render(
    input = rmd_tmp,
    output_file = basename(output_file),
    output_dir = output_dir,
    params = list(
      result = result,
      species = species,
      var_labels = var_labels
    ),
    envir = new.env(parent = globalenv()),
    quiet = !verbose
  )

  if (verbose) cli::cli_inform("Report saved to {.file {output_file}}")

  if (open && interactive()) utils::browseURL(output_file)
  invisible(output_file)
}


#' Generate inline Rmd template
#' @keywords internal
#' @noRd
.cast_report_rmd_content <- function(species = "Species") {
  paste0('---
title: "castSDM Pipeline Report: ', species, '"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    self_contained: true
params:
  result: NULL
  species: "Species"
  var_labels: NULL
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
                      fig.width = 10, fig.height = 6)
library(castSDM)
result <- params$result
var_labels <- params$var_labels
```
## 1. DAG Structure

```{r dag-plot, fig.height=8}
if (!is.null(result$dag)) {
  plot(result$dag, roles = result$roles, screen = result$screen,
       var_labels = var_labels, species = params$species)
}
```

## 2. Variable Selection

```{r screen-plot}
if (!is.null(result$screen)) {
  plot(result$screen, var_labels = var_labels)
}
```

```{r roles-table}
if (!is.null(result$screen) && !is.null(result$screen$roles)) {
  knitr::kable(result$screen$roles, caption = "Causal Roles (Markov Blanket)")
}
```

## 3. Model Evaluation

```{r eval-plot}
if (!is.null(result$eval)) {
  plot(result$eval)
}
```

```{r cv-plot}
if (!is.null(result$cv)) {
  plot(result$cv)
}
```

## 4. Prediction Maps

```{r pred-plot, fig.height=8}
if (!is.null(result$predict)) {
  for (mdl in result$predict$models) {
    print(plot(result$predict, model = mdl))
  }
}
```

```{r ensemble-plot, fig.height=8}
if (!is.null(result$ensemble)) {
  plot(result$ensemble)
}
```

## 5. SHAP Explanations

```{r shap-section, results="asis"}
if (!is.null(result$shap)) {
  for (nm in names(result$shap)) {
    if (!is.null(result$shap[[nm]])) {
      cat(sprintf("\\n### SHAP: %s\\n\\n", nm))
      print(plot(result$shap[[nm]]))
      cat("\\n")
    }
  }
} else {
  cat("*SHAP not computed in this run.*\\n")
}
```

## Session Info

```{r session}
sessionInfo()
```
')
}
