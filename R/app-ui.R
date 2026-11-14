# AI-assisted development (2026-11-13); see inst/AI_ASSISTANCE.md.

.cast_app_scenario_label <- function(x) {
  pattern <- "^.*_(ssp[0-9]+)_([0-9]{4})_([0-9]{4})$"
  matched <- grepl(pattern, x)
  out <- x
  ssp <- sub(pattern, "\\1", x[matched])
  label <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp370 = "SSP3-7.0", ssp585 = "SSP5-8.5")[ssp]
  label[is.na(label)] <- ssp[is.na(label)]
  out[matched] <- paste(label, paste0(sub(pattern, "\\2", x[matched]), "\u2013", sub(pattern, "\\3", x[matched])), sep = " \u00b7 ")
  out
}

.cast_app_driver_label <- function(project, variables) {
  out <- variables
  if (is.null(project) || !"label" %in% names(project$metadata)) return(out)
  labels <- project$metadata$label[match(variables, project$metadata$layer)]
  known <- !is.na(labels) & nzchar(labels)
  out[known] <- labels[known]
  out
}

.cast_app_metric <- function(label, value, suffix = NULL, title = NULL) {
  shiny::tags$div(class = "cast-metric", title = title, shiny::tags$span(class = "cast-metric-label", label),
    shiny::tags$div(class = "cast-metric-value", value, shiny::tags$small(suffix)))
}

.cast_app_point_contrast <- function(point, layer, focus) {
  v <- point$values
  result <- list(before = v$current, after = v$future, change = v$total_change,
    supported = isTRUE(v$support_total == 1), kind = "joint", selected = NULL)
  if (layer %in% names(v)) result$selected <- v[[layer]]
  if (!is.null(focus) && !is.null(point$driver_response)) {
    d <- point$driver_response; s <- point$driver_support
    if (layer %in% names(d)) result$selected <- d[[layer]]
    only <- d[[paste0("effect_", focus)]]; retained <- d[[paste0("retained_", focus)]]
    if (identical(layer, paste0("effect_", focus))) {
      result$change <- only; result$after <- v$current + only; result$kind <- "only"
      result$supported <- isTRUE(s[[1L]] == 1)
    } else if (identical(layer, paste0("retained_", focus))) {
      result$change <- retained; result$before <- v$future - retained; result$kind <- "retained"
      result$supported <- isTRUE(s[[2L]] == 1)
    } else if (identical(layer, paste0("context_gap_", focus))) {
      result$before <- only; result$after <- retained; result$change <- retained - only
      result$kind <- "gap"; result$supported <- all(unlist(s) == 1)
    }
  }
  result
}

.cast_app_asset <- function(name) {
  path <- system.file("app", name, package = "castSDM")
  if (!nzchar(path) || !file.exists(path))
    cli::cli_abort("Reinstall castSDM: missing application asset {.val {name}}.")
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

.cast_app_driver_table <- function(project, data) {
  value <- function(x) if (is.finite(x)) sprintf("%+.4f", x) else "\u2014"
  percent <- function(x) if (is.finite(x)) sprintf("%.1f%%", 100 * x) else "\u2014"
  cell <- function(x) shiny::tags$td(class = if (is.finite(x) && x > 0) "cast-value-positive" else
    if (is.finite(x) && x < 0) "cast-value-negative" else NULL, value(x))
  headers <- c("\u73af\u5883\u56e0\u5b50", "\u53d8\u5316\u6837\u70b9", "\u6bd4\u8f83\u524d\u54cd\u5e94", "\u6bd4\u8f83\u540e\u54cd\u5e94", "\u80cc\u666f\u53cd\u53f7", "\u6a21\u578b\u540c\u5411\u00b7\u524d", "\u6a21\u578b\u540c\u5411\u00b7\u540e")
  shiny::tags$table(class = "cast-table",
    shiny::tags$thead(shiny::tags$tr(lapply(headers, function(x) shiny::tags$th(scope = "col", x)))),
    shiny::tags$tbody(lapply(seq_len(nrow(data)), function(i) {
      d <- data[i, ]
      shiny::tags$tr(shiny::tags$td(shiny::tags$button(type = "button", class = "cast-driver-button", `data-driver` = d$driver,
        title = paste(d$driver, "| Model agreement: changed supported locations only"), shiny::tags$span(class = "cast-rank", sprintf("%02d", d$rank)),
        shiny::tags$span(class = "cast-driver-name", .cast_app_driver_label(project, d$driver)))),
        shiny::tags$td(d$n_changed_common), cell(d$mean_only), cell(d$mean_retained), shiny::tags$td(percent(d$context_reversal_fraction)),
        shiny::tags$td(percent(d$all_models_agree_changed_only)), shiny::tags$td(percent(d$all_models_agree_changed_retained)))
    })))
}

.cast_app_ui <- function(case_dir, workflow = NULL) {
  css <- .cast_app_asset("app.css")
  script <- .cast_app_asset("app.js")
  shiny::fluidPage(
    theme = bslib::bs_theme(version = 5, bg = "#FFFFFF", fg = "#202B32", primary = "#206B61"),
    shiny::tags$head(shiny::tags$style(shiny::HTML(paste(css, collapse = "\n"))),
      shiny::tags$script(shiny::HTML(paste(script, collapse = "\n")))),
    shiny::div(class = "cast-shell",
      shiny::tags$header(class = "cast-header",
        shiny::div(class = "cast-brand", shiny::tags$span(class = "cast-wordmark", "castSDM"),
          shiny::tags$span(class = "cast-brand-description", "\u73af\u5883\u54cd\u5e94\u5206\u6790")),
        shiny::uiOutput("case_label", inline = TRUE)),
      shiny::div(class = "cast-layout",
        shiny::tags$aside(class = "cast-sidebar",
          shiny::tags$p(class = "cast-eyebrow", "\u5206\u6790\u65b9\u6848"),
          shiny::selectInput("scenario", "\u672a\u6765\u60c5\u666f", choices = character(), width = "100%"),
          shiny::actionButton("view", "\u67e5\u770b\u5df2\u6709\u5730\u56fe", class = "btn-outline-secondary cast-wide", `data-cast-action` = ""),
          shiny::actionButton("analyse", "\u5206\u6790\u73af\u5883\u54cd\u5e94", class = "btn-primary cast-wide", `data-cast-action` = ""),
          shiny::div(class = "cast-status", role = "status", `aria-live` = "polite", shiny::uiOutput("status")),
          shiny::actionButton("cancel_job", "\u505c\u6b62\u8ba1\u7b97", class = "btn-link cast-cancel", hidden = "hidden"),
          shiny::tags$details(class = "cast-settings", shiny::tags$summary("\u66f4\u591a\u8bbe\u7f6e"),
            shiny::tags$p(class = "cast-section-label", "\u73af\u5883\u6539\u53d8"),
            shiny::radioButtons("mode", NULL, choices = stats::setNames(c("scenario", "shift"), c("\u4f7f\u7528\u5b9e\u9645\u60c5\u666f", "\u5047\u8bbe\u5355\u56e0\u5b50\u6539\u53d8"))),
            shiny::selectInput("variables", "\u53c2\u4e0e\u6bd4\u8f83\u7684\u56e0\u5b50", choices = character(), width = "100%"),
            shiny::uiOutput("variable_definition"),
            shiny::conditionalPanel("input.mode == 'shift'",
              shiny::numericInput("shift", "\u6539\u53d8\u5e45\u5ea6 \u00b7 \u539f\u59cb\u5355\u4f4d", value = NA_real_),
              shiny::tags$p(class = "cast-note", "\u6d3e\u751f\u6307\u6807\u6539\u53d8\u5c5e\u4e8e\u6a21\u578b\u654f\u611f\u6027\u5206\u6790\u3002\u571f\u5730\u5229\u7528\u6539\u53d8\u4f18\u5148\u4f7f\u7528\u4e00\u81f4\u7684\u60c5\u666f\u56fe\u5c42\u3002")),
            shiny::conditionalPanel("input.mode == 'shift'", shiny::selectInput("context", "\u5176\u4f59\u56e0\u5b50\u7684\u73af\u5883", choices = stats::setNames(c("current", "future"), c("\u4fdd\u6301\u5f53\u524d\u503c", "\u4f7f\u7528\u672a\u6765\u60c5\u666f\u503c")))),
            shiny::tags$p(class = "cast-section-label", "\u8ba1\u7b97"),
            shiny::numericInput("sample_n", "\u6837\u70b9\u4e0a\u9650", 5000L, min = 4L, step = 100L),
            shiny::checkboxInput("grouping", "\u53ef\u9009\uff1a\u63cf\u8ff0\u6027\u5206\u533a", FALSE),
            shiny::conditionalPanel("input.grouping", shiny::numericInput("maxdepth", "\u5206\u7ec4\u6df1\u5ea6\u4e0a\u9650", 3L, min = 1L, max = 8L)),
            shiny::numericInput("cores", "\u8ba1\u7b97\u7ebf\u7a0b", 4L, min = 1L),
            shiny::checkboxInput("maps", "\u540c\u65f6\u8ba1\u7b97\u56e0\u5b50\u539f\u751f\u5730\u56fe", FALSE),
            shiny::tags$p(class = "cast-note", "\u9ed8\u8ba4\u5148\u8ba1\u7b97\u6837\u70b9\u6982\u89c8\uff0c\u518d\u6309\u9700\u751f\u6210\u56e0\u5b50\u5730\u56fe\u3002"),
            shiny::tags$p(class = "cast-section-label", "\u6848\u4f8b"),
            shiny::textInput("case_dir", "\u7ed3\u679c\u76ee\u5f55", value = case_dir),
            shiny::actionButton("load_case", "\u6253\u5f00\u6848\u4f8b", class = "btn-outline-secondary cast-wide", `data-cast-action` = ""))),
        shiny::tags$main(class = "cast-main", id = "analysis",
          shiny::tabsetPanel(id = "result_tab", type = "tabs", selected = if (file.exists(file.path(case_dir, "completion.rds"))) "map" else "setup",
            .cast_app_setup_ui(.cast_app_setup_defaults(case_dir, workflow)),
            shiny::tabPanel("\u5730\u56fe", value = "map",
              shiny::div(class = "cast-view-heading", shiny::uiOutput("map_caption")),
              shiny::div(class = "cast-map-toolbar",
                shiny::selectInput("layer", "\u663e\u793a\u5185\u5bb9", choices = stats::setNames(c("total_change", "current", "future", "support_total", "total_sd"),
                  c("\u9002\u5b9c\u5ea6\u53d8\u5316", "\u6539\u53d8\u524d", "\u6539\u53d8\u540e", "\u8303\u56f4\u652f\u6301", "\u6a21\u578b\u5206\u6b67")), width = "100%"),
                shiny::div(class = "cast-supported-filter", shiny::checkboxInput("supported_only", "\u4ec5\u663e\u793a\u8303\u56f4\u652f\u6301\u4f4d\u7f6e", TRUE))),
              shiny::div(class = "cast-map-frame", leaflet::leafletOutput("map", height = 520), shiny::uiOutput("map_hint")),
              shiny::tags$p(class = "cast-note cast-map-note", "\u9884\u89c8\u7528\u4e8e\u6d4f\u89c8\uff1b\u70b9\u51fb\u8bfb\u53d6\u539f\u751f\u7f51\u683c\u3002"),
              shiny::tags$details(id = "point_details", class = "cast-detail", shiny::tags$summary("\u4f4d\u7f6e\u67e5\u8be2"),
                shiny::uiOutput("point_summary"),
                shiny::selectInput("profile_driver", "\u5c40\u5730\u54cd\u5e94\u5256\u9762", choices = character()),
                shiny::plotOutput("profile", height = 230),
                shiny::tags$p(class = "cast-note", "\u5256\u9762\u56fa\u5b9a\u6bd4\u8f83\u524d\u7684\u5176\u4ed6\u56e0\u5b50\u3002"),
                shiny::tags$details(shiny::tags$summary("\u6240\u9009\u56e0\u5b50\u7684\u4e24\u79cd\u80cc\u666f\u54cd\u5e94"), shiny::tableOutput("point_driver")),
                shiny::tags$details(shiny::tags$summary("\u539f\u59cb\u6570\u503c\u4e0e\u73af\u5883"), shiny::tableOutput("point"), shiny::tableOutput("environment"))),
              shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u5982\u4f55\u7406\u89e3\u5730\u56fe"),
                shiny::tags$p(class = "cast-note", "\u7ed3\u679c\u662f\u62df\u5408\u6a21\u578b\u7684\u76f8\u5bf9\u9002\u5b9c\u5ea6\u54cd\u5e94\u3002\u8303\u56f4\u652f\u6301\u4e0d\u4ee3\u8868\u8054\u5408\u652f\u6301\u6216\u56e0\u679c\u8bc6\u522b\uff1b\u6a21\u578b\u5206\u6b67\u4e0d\u662f\u7f6e\u4fe1\u533a\u95f4\u3002"))),
            shiny::tabPanel("\u73af\u5883\u54cd\u5e94", value = "drivers",
              shiny::div(class = "cast-view-heading", shiny::tags$h2("\u54ea\u4e9b\u53d8\u5316\u5f71\u54cd\u66f4\u660e\u663e"),
                shiny::tags$p(class = "cast-subtitle", "\u5728\u540c\u4e00\u6279\u652f\u6301\u6837\u70b9\u4e0a\u6bd4\u8f83\u73af\u5883\u56e0\u5b50\uff0c\u67e5\u770b\u80cc\u666f\u4f9d\u8d56\u4e0e\u6a21\u578b\u4e00\u81f4\u6027\u3002")),
              shiny::uiOutput("response_summary"), shiny::uiOutput("evidence_note"),
              shiny::div(class = "cast-table-scroll", shiny::uiOutput("drivers")),
              shiny::div(class = "cast-focus-toolbar", shiny::selectInput("focus_driver", "\u67e5\u770b\u56e0\u5b50\u5730\u56fe", choices = character(), width = "100%"),
                shiny::actionButton("map_driver", "\u751f\u6210\u56e0\u5b50\u5730\u56fe", class = "btn-primary", `data-cast-action` = "")),
              shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u6837\u70b9\u4e2d\u672a\u53d8\u5316\u6216\u672a\u6392\u5e8f\u7684\u56e0\u5b50"), shiny::uiOutput("drivers_inactive")),
              shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u54cd\u5e94\u5bf9\u7167\u56fe"), shiny::plotOutput("driver_plot", height = 300)),
              shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u81ea\u52a8\u5206\u7ec4\u4e0e\u8bca\u65ad"), shiny::uiOutput("group_note"),
                shiny::tableOutput("groups"), shiny::plotOutput("group_plot", height = 230), shiny::tableOutput("tuning")),
              shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u6307\u6807\u5b9a\u4e49\u4e0e\u6765\u6e90"),
                shiny::div(class = "cast-table-scroll", shiny::tableOutput("metadata"))),
              shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u89e3\u91ca\u8fb9\u754c"),
                shiny::tags$p(class = "cast-note", "\u6392\u5e8f\u9488\u5bf9\u5b9e\u9645\u60c5\u666f\u6539\u53d8\uff0c\u4e0d\u4ee3\u8868\u56fa\u6709\u91cd\u8981\u6027\u6216\u6cbb\u7406\u6210\u672c\u3002\u672a\u6539\u53d8\u4e0d\u4ee3\u8868\u4e0d\u91cd\u8981\uff1b\u5355\u56e0\u5b50\u54cd\u5e94\u4e0d\u5fc5\u52a0\u603b\u4e3a\u603b\u53d8\u5316\u3002\u80cc\u666f\u548c\u6a21\u578b\u540c\u5411\u6027\u662f\u8bca\u65ad\uff0c\u4e0d\u662f\u56e0\u679c\u6548\u5e94\u6216\u7f6e\u4fe1\u5ea6\u3002"))),
            shiny::tabPanel("\u7ed3\u679c", value = "figures",
              shiny::div(class = "cast-view-heading", shiny::tags$h2("\u4fdd\u5b58\u4e0e\u590d\u73b0")),
              shiny::div(class = "cast-export-actions",
                shiny::downloadButton("download_case", "\u5b8c\u6574\u7ed3\u679c", class = "btn-outline-secondary"),
                shiny::downloadButton("download_bundle", "\u4e0b\u8f7d\u5206\u6790\u5305", class = "btn-primary"),
                shiny::downloadButton("download_points", "\u6837\u70b9\u6570\u636e", class = "btn-outline-secondary"),
                shiny::downloadButton("download_figure", "\u5f53\u524d\u56fe\u4ef6", class = "btn-outline-secondary")),
              shiny::tags$p(class = "cast-note", "\u5206\u6790\u5305\u5305\u542b\u7ed3\u679c\u3001\u8bbe\u7f6e\u4e0e\u53ef\u6267\u884c\u7684\u754c\u9762\u590d\u73b0\u811a\u672c\u3002\u9700\u8981\u539f\u59cb\u6848\u4f8b\u6587\u4ef6\u3002"),
              shiny::tags$details(class = "cast-detail", open = "open", shiny::tags$summary("\u539f\u6709\u6b63\u5f0f\u56fe\u4ef6"),
                shiny::selectInput("figure", "\u9009\u62e9\u56fe\u4ef6", choices = character()), shiny::imageOutput("publication", height = "auto")),
              shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u6a21\u578b\u8bc4\u4ef7"),
                shiny::div(class = "cast-table-scroll", shiny::tableOutput("evaluation")))))))))

}
