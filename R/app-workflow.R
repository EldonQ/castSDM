.cast_app_setup_defaults <- function(case_dir, workflow) {
  root <- if (!is.null(workflow)) dirname(workflow) else file.path(getwd(), "castSDM_analysis")
  settings_file <- file.path(root, "workflow_settings.rds")
  settings <- if (file.exists(settings_file)) readRDS(settings_file) else NULL
  if (is.null(settings) && !is.null(workflow)) {
    env <- new.env(parent = environment())
    for (f in c("inventory.R", "01_config.R")) sys.source(file.path(root, f), env)
    cfg <- tryCatch(env$castsdm_species_configs(root)[[Sys.getenv("CAST_SPECIES", "Schizothoracinae")]], error = function(e) NULL)
    if (!is.null(cfg)) settings <- list(data_dir = cfg$nsdm_data_root,
      occurrence_file = cfg$occurrence_file, species = cfg$species_name,
      domain = cfg$domain_mode, input_type = cfg$input_type, lulc_product = cfg$current_lulc_product,
      models = cfg$protocol$models, cv_k = cfg$protocol$cv_k, background_points = cfg$bg_n, resources = cfg$resources)
  }
  settings <- settings %||% list()
  utils::modifyList(list(path = root, data_dir = Sys.getenv("CAST_DATA_ROOT", ""),
    occurrence_file = "", species = "Schizothoracinae", domain = "freshwater_river",
    input_type = "gbif", lulc_product = "glu", models = c("rf", "maxent", "gam"), cv_k = 5L,
    resources = list(cores = 4L, workers = 1L, memory_mb = 4096)), settings)
}

.cast_app_setup_ui <- function(defaults) {
  shiny::tabPanel("\u5f00\u59cb", value = "setup",
    shiny::div(class = "cast-view-heading", shiny::tags$h2("\u4ece\u6570\u636e\u5230\u73af\u5883\u54cd\u5e94"),
      shiny::tags$p(class = "cast-subtitle", "\u68c0\u67e5\u6570\u636e \u2192 \u8bad\u7ec3\u4e0e\u7a7a\u95f4\u9a8c\u8bc1 \u2192 \u672a\u6765\u9884\u6d4b \u2192 \u56fe\u4ef6\u4e0e\u73af\u5883\u54cd\u5e94")),
    shiny::div(class = "cast-setup-form",
      shiny::textInput("setup_data", "\u73af\u5883\u6570\u636e\u76ee\u5f55", defaults$data_dir),
      shiny::textInput("setup_occurrence", "\u7269\u79cd\u8bb0\u5f55\u6587\u4ef6", defaults$occurrence_file),
      shiny::div(class = "cast-setup-grid",
        shiny::textInput("setup_species", "\u5206\u6790\u7c7b\u7fa4\uff08\u82f1\u6587\u6807\u8bc6\uff09", defaults$species),
        shiny::selectInput("setup_domain", "\u5206\u6790\u57df", stats::setNames(c("freshwater_river", "national"), c("\u6cb3\u7f51\u4e0e\u8bb0\u5f55\u9644\u8fd1", "\u5168\u56fd\u9646\u57df")), defaults$domain)),
      shiny::textInput("setup_path", "\u5de5\u4f5c\u4e0e\u8f93\u51fa\u76ee\u5f55", defaults$path),
      shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u5206\u6790\u8bbe\u7f6e"),
        shiny::checkboxGroupInput("setup_models", "\u6a21\u578b", choices = stats::setNames(c("rf", "maxent", "gam"), c("Random forest", "MaxEnt", "GAM")), selected = defaults$models, inline = TRUE),
        shiny::div(class = "cast-setup-grid",
          shiny::numericInput("setup_background", "\u80cc\u666f\u70b9\u6570\u91cf", defaults$background_points %||% 10000L, min = 1L, step = 1000L),
          shiny::numericInput("setup_cv", "\u7a7a\u95f4\u4ea4\u53c9\u9a8c\u8bc1\u6298\u6570", defaults$cv_k, min = 3L, step = 1L),
          shiny::numericInput("setup_cores", "\u603b\u8ba1\u7b97\u7ebf\u7a0b", defaults$resources$cores, min = 1L, step = 1L),
          shiny::numericInput("setup_workers", "\u5e76\u884c\u60c5\u666f\u6570", defaults$resources$workers, min = 1L, step = 1L),
          shiny::numericInput("setup_memory", "\u603b\u5185\u5b58\u9884\u7b97 MB", defaults$resources$memory_mb, min = 256L, step = 512L),
          shiny::selectInput("setup_input", "\u8bb0\u5f55\u683c\u5f0f", c("GBIF / literature" = "gbif", "N-SDM occurrence" = "nsdm_occ"), defaults$input_type),
          shiny::selectInput("setup_lulc", "\u5f53\u524d\u571f\u5730\u5229\u7528", c("GLU 2020" = "glu", "HLU 2020" = "hlu"), defaults$lulc_product))),
      shiny::tags$p(class = "cast-note", "CHELSA 1981\u20132010 \u2192 SSP1-2.6 / SSP3-7.0 / SSP5-8.5\uff1b2011\u20132040\u30012041\u20132070\u30012071\u20132100\u3002\u672a\u6765\u4ec5\u66ff\u6362\u96c6\u6210\u6c14\u5019\uff0c\u5176\u4ed6\u56e0\u5b50\u4fdd\u6301\u5f53\u524d\u503c\u3002"),
      shiny::div(class = "cast-export-actions",
        shiny::actionButton("check_inputs", "\u68c0\u67e5\u6570\u636e", class = "btn-outline-secondary", `data-cast-action` = ""),
        shiny::actionButton("start_analysis", "\u8fd0\u884c\u5b8c\u6574\u5206\u6790", class = "btn-primary", `data-cast-action` = "")),
      shiny::uiOutput("setup_summary"),
      shiny::tags$details(class = "cast-detail", shiny::tags$summary("\u8fd0\u884c\u65e5\u5fd7"), shiny::uiOutput("log"))))
}
