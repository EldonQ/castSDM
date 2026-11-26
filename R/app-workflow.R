.cast_app_setup_defaults <- function(case_dir, workflow) {
  root <- if (!is.null(workflow)) dirname(workflow) else file.path(getwd(), "castSDM_analysis")
  selection <- .cast_workflow_selection(root, allow_unconfigured = TRUE)
  settings <- selection$settings
  settings <- settings %||% list()
  utils::modifyList(list(path = root, data_dir = Sys.getenv("CAST_DATA_ROOT", ""),
    occurrence_file = "", species = "", domain = "national",
    input_type = "gbif", lulc_product = "glu", models = c("rf", "maxent", "gam"), cv_k = 5L,
    resources = list(cores = 4L, workers = 1L, memory_mb = 4096)), settings)
}

.cast_app_setup_ui <- function(defaults) {
  shiny::tabPanel("\u5f00\u59cb", value = "setup",
    shiny::div(class = "cast-view-heading", shiny::tags$h2("\u4ece\u6570\u636e\u5230\u73af\u5883\u54cd\u5e94"),
      shiny::tags$p(class = "cast-subtitle", "\u4e2d\u56fd\u73af\u5883\u5e93\u4e0e CHELSA \u6848\u4f8b\uff1a\u68c0\u67e5\u6570\u636e \u2192 \u7a7a\u95f4\u9a8c\u8bc1 \u2192 \u60c5\u666f\u9884\u6d4b \u2192 \u73af\u5883\u54cd\u5e94")),
    shiny::div(class = "cast-setup-form",
      shiny::textInput("setup_data", "\u73af\u5883\u6570\u636e\u76ee\u5f55", defaults$data_dir),
      shiny::textInput("setup_occurrence", "\u7269\u79cd\u8bb0\u5f55\u6587\u4ef6", defaults$occurrence_file),
      shiny::div(class = "cast-setup-grid",
        shiny::textInput("setup_species", "\u5206\u6790\u7c7b\u7fa4\uff08\u82f1\u6587\u6807\u8bc6\uff09", defaults$species),
        shiny::selectInput("setup_domain", "\u5206\u6790\u57df", stats::setNames(c("national", "custom_mask"), c("\u5168\u56fd\u9646\u57df", "\u81ea\u5b9a\u4e49\u7814\u7a76\u57df")), defaults$domain)),
      shiny::conditionalPanel("input.setup_domain == 'custom_mask'",
        shiny::textInput("setup_domain_mask", "\u7814\u7a76\u57df\u6805\u683c\u8def\u5f84\uff08\u6b63\u503c\u7f51\u683c\u6709\u6548\uff09", defaults$domain_mask %||% "")),
      shiny::textInput("setup_path", "\u5de5\u4f5c\u4e0e\u8f93\u51fa\u76ee\u5f55", defaults$path),
      shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u5206\u6790\u8bbe\u7f6e"),
        shiny::checkboxGroupInput("setup_models", "\u6a21\u578b", choices = stats::setNames(c("rf", "maxent", "gam"), c("Random forest", "MaxEnt", "GAM")), selected = defaults$models, inline = TRUE),
        shiny::div(class = "cast-setup-grid",
          shiny::selectInput("setup_validation_target", "\u9a8c\u8bc1\u76ee\u6807",
            stats::setNames(c("within_domain", "spatial_transfer"), c("\u57df\u5185\u9884\u6d4b", "\u8de8\u533a\u57df\u5916\u63a8")),
            defaults$validation_target %||% "spatial_transfer"),
          shiny::textInput("setup_predictor_file", "\u5019\u9009\u53d8\u91cf CSV\uff08\u53ef\u9009\uff09", defaults$predictor_file %||% ""),
          shiny::textInput("setup_target_group", "\u540c\u7c7b\u8c03\u67e5\u80cc\u666f\u70b9\u6587\u4ef6\uff08\u53ef\u9009\uff09", defaults$target_group_file %||% ""),
          shiny::numericInput("setup_background", "\u80cc\u666f\u70b9\u6570\u91cf\uff08\u7559\u7a7a\u81ea\u52a8\uff09", defaults$background_points %||% NA_real_, min = 1L, step = 1000L),
          shiny::numericInput("setup_thin", "\u8ddd\u79bb\u7a00\u758f\u5316\uff08km\uff1b0 \u4e3a\u4e0d\u7a00\u758f\uff09", defaults$thin_km %||% 0, min = 0),
          shiny::numericInput("setup_coverage", "\u5019\u9009\u53d8\u91cf\u6700\u4f4e\u57df\u5185\u8986\u76d6\u7387\uff080\u20131\uff09", defaults$predictor_min_coverage %||% .95, min = 0, max = 1),
          shiny::numericInput("setup_vif", "VIF \u9608\u503c", defaults$vif_threshold %||% 10, min = 1, step = 1),
          shiny::selectInput("setup_missing", "\u7f3a\u5931\u73af\u5883\u6307\u6807", stats::setNames(c("complete_cases", "drop_predictors"), c("\u4fdd\u7559\u6307\u6807\uff0c\u53ea\u7528\u5b8c\u6574\u6837\u672c", "\u4e25\u683c\u5254\u9664\u7f3a\u5931\u6307\u6807\uff0c\u4fdd\u7559\u4f4d\u7f6e")), defaults$predictor_missing %||% "complete_cases"),
          shiny::selectInput("setup_block_method", "\u7a7a\u95f4\u5206\u5757",
            stats::setNames(c("auto", "blockcv", "grid"), c("\u81ea\u52a8", "blockCV", "\u5730\u7406\u7f51\u683c")),
            defaults$block_method %||% "auto"),
          shiny::conditionalPanel("input.setup_block_method == 'blockcv'",
            shiny::numericInput("setup_block_size", "\u5185\u5c42\u5757\u8fb9\u957f\uff08\u7c73\uff1b\u7559\u7a7a\u81ea\u52a8\uff09", defaults$block_size %||% NA_real_, min = 1),
            shiny::numericInput("setup_holdout_block_size", "\u5916\u5c42\u5757\u8fb9\u957f\uff08\u7c73\uff1b\u53ef\u7559\u7a7a\uff09", defaults$holdout_block_size %||% NA_real_, min = 1)),
          shiny::numericInput("setup_cv", "\u7a7a\u95f4\u4ea4\u53c9\u9a8c\u8bc1\u6298\u6570", defaults$cv_k, min = 3L, step = 1L),
          shiny::numericInput("setup_cores", "\u603b\u8ba1\u7b97\u7ebf\u7a0b", defaults$resources$cores, min = 1L, step = 1L),
          shiny::numericInput("setup_projection_cores", "\u60c5\u666f\u9884\u6d4b\u7ebf\u7a0b\uff08\u53ef\u7559\u7a7a\uff09", defaults$projection_cores %||% NA_real_, min = 1L, step = 1L),
          shiny::numericInput("setup_workers", "\u5e76\u884c\u60c5\u666f\u6570", defaults$resources$workers, min = 1L, step = 1L),
          shiny::numericInput("setup_memory", "\u603b\u5185\u5b58\u9884\u7b97 MB", defaults$resources$memory_mb, min = 256L, step = 512L),
          shiny::selectInput("setup_input", "\u8bb0\u5f55\u683c\u5f0f", c("GBIF / literature" = "gbif", "N-SDM occurrence" = "nsdm_occ"), defaults$input_type),
          shiny::selectInput("setup_lulc", "\u5f53\u524d\u571f\u5730\u5229\u7528", c("GLU 2020" = "glu", "HLU 2020" = "hlu"), defaults$lulc_product))),
      shiny::tags$p(class = "cast-note", "CHELSA 1981\u20132010 \u2192 SSP1-2.6 / SSP3-7.0 / SSP5-8.5\uff1b2011\u20132040\u30012041\u20132070\u30012071\u20132100\u3002\u672a\u6765\u4ec5\u66ff\u6362\u96c6\u6210\u6c14\u5019\uff0c\u5176\u4ed6\u56e0\u5b50\u4fdd\u6301\u5f53\u524d\u503c\u3002"),
      shiny::div(class = "cast-export-actions",
        shiny::actionButton("check_inputs", "\u68c0\u67e5\u6570\u636e\u4e0e\u7a7a\u95f4\u6298", class = "btn-outline-secondary", `data-cast-action` = ""),
        shiny::actionButton("start_analysis", "\u8fd0\u884c\u5b8c\u6574\u5206\u6790", class = "btn-primary", `data-cast-action` = "")),
      shiny::uiOutput("setup_summary"),
      shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u8fd0\u884c\u65e5\u5fd7"), shiny::uiOutput("log"))))
}
