# AI-assisted development (2026-11-11); see inst/AI_ASSISTANCE.md.

.cast_app_identity <- function() {
  functions <- c(".cast_app_task", "cast_explore", "cast_attribution", "cast_response_groups")
  list(version = as.character(utils::packageVersion("castSDM")),
    path = normalizePath(getNamespaceInfo(asNamespace("castSDM"), "path"), winslash = "/"),
    code = vapply(functions, function(n) .cast_code_identity(get(n, asNamespace("castSDM"))), character(1)))
}

.cast_app_start <- function(request, logfile) {
  worker <- function(request, libraries, expected) {
    .libPaths(libraries)
    ns <- loadNamespace("castSDM")
    actual <- get(".cast_app_identity", ns)()
    if (!identical(actual, expected))
      stop("Background castSDM differs from the calling installation. Reinstall and reopen the app.")
    get(".cast_app_task", ns)(request)
  }
  # No package namespace/source reference may be deserialised before libPaths.
  environment(worker) <- baseenv()
  callr::r_bg(worker, args = list(request = request, libraries = .libPaths(),
    expected = .cast_app_identity()), libpath = .libPaths(),
    stdout = logfile, stderr = "2>&1", supervise = TRUE,
    system_profile = FALSE, user_profile = FALSE)
}

.cast_app_layer <- function(result, layer) {
  r <- terra::rast(result$raster_file)
  name <- if (identical(result$schema, "cast_app_view/1"))
    switch(layer, future = "joint", total_change = "joint_change",
      support_total = "support", total_sd = "sd_joint_change", layer) else layer
  r[[name]]
}

.cast_app_point <- function(project, result, longitude, latitude, driver) {
  settings <- result$settings
  pair <- .cast_explore_pair(project, settings$scenario, settings$variables,
    settings$shifts, settings$context)
  xy <- matrix(c(longitude, latitude), ncol = 2L)
  xy <- terra::project(xy, "EPSG:4326", terra::crs(pair$current))
  before <- terra::extract(pair$current, xy)
  after <- terra::extract(pair$future, xy)
  for (v in names(settings$shifts)) after[[v]] <- after[[v]] + settings$shifts[[v]]
  vars <- project$fit$env_vars
  if (!all(is.finite(unlist(before[vars]))) || !all(is.finite(unlist(after[vars]))))
    cli::cli_abort("\u8be5\u4f4d\u7f6e\u7f3a\u5c11\u5b8c\u6574\u73af\u5883\u503c\uff0c\u8bf7\u5728\u5206\u6790\u57df\u5185\u70b9\u51fb\u3002")
  focus <- settings$focus
  paired <- cast_attribution(project$fit, before, after, drivers = focus %||% character(),
    weights = settings$weights, support_probs = project$comparison$support_probs,
    clamp = settings$clamp %||% FALSE, verbose = FALSE, diagnostics = TRUE)
  values <- paired$data[1L, c("current", "future", "total_change", "support_total", "total_sd")]
  # Read the native analysis mask/results, not a smoothed display pixel.
  domain <- if (!is.null(result$raster_file))
    terra::extract(.cast_app_layer(result, "current"), xy)[1L, 1L] else
      terra::extract(terra::rast(file.path(project$case$preparation, "analysis_domain_mask.tif")), xy)[1L, 1L]
  if (!is.finite(domain) || (is.null(result$raster_file) && domain == 0))
    cli::cli_abort("\u8be5\u4f4d\u7f6e\u4e0d\u5728\u8fd9\u9879\u5206\u6790\u7684\u5171\u540c\u6709\u6548\u57df\u5185\u3002")
  if (!is.null(result$raster_file)) {
    fields <- names(values)
    values <- as.data.frame(stats::setNames(lapply(fields, function(n)
      terra::extract(.cast_app_layer(result, n), xy)[1L, 1L]), fields))
  }
  group <- if (!is.null(result$groups) && isTRUE(values$support_total == 1))
    stats::predict(result$groups, before)[1L] else NA_integer_
  reference <- .cast_reference(project$fit, vars)[[driver]]
  grid <- unique(as.numeric(stats::quantile(reference[is.finite(reference)],
    probs = seq(0, 1, length.out = 40L), type = 1L, names = FALSE)))
  profile <- before[rep(1L, length(grid)), , drop = FALSE]
  profile[[driver]] <- grid
  context <- .cast_prediction_context(project$fit, settings$weights,
    profile, profile, settings$clamp %||% FALSE)
  predicted <- context$predict()
  curve <- as.numeric(predicted %*% settings$weights)
  curve[rowSums(!is.finite(predicted)) > 0L] <- NA_real_
  list(values = values, group = group,
    driver_response = if (!is.null(focus)) paired$data[, paste0(c("effect_", "retained_",
      "effect_agreement_", "retained_agreement_", "context_direction_", "context_gap_"), focus), drop = FALSE],
    environment = data.frame(variable = vars, before = unlist(before[vars], use.names = FALSE),
      after = unlist(after[vars], use.names = FALSE), row.names = NULL),
    curve = data.frame(environment = grid, suitability = curve), driver = driver,
    longitude = longitude, latitude = latitude)
}

.cast_app_ui <- function(case_dir) {
  shiny::fluidPage(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#176B87"),
    shiny::tags$head(shiny::tags$style(shiny::HTML(paste(
      ".cast-panel{background:#f5f8fa;padding:18px;border-radius:10px;margin-bottom:16px}",
      ".cast-note{color:#52616b;font-size:13px}.cast-heading{margin:20px 0}",
      ".cast-map-caption{padding:8px 0;font-weight:600}table{font-size:13px}",
      ".shiny-output-error-validation{color:#52616b} .cast-panel .btn{margin:3px 0}")))),
    shiny::tags$div(class = "cast-heading", shiny::h2("castSDM \u00b7 \u73af\u5883\u54cd\u5e94\u5730\u56fe"),
      shiny::tags$p(class = "cast-note", "\u6bd4\u8f83\u660e\u786e\u7684\u73af\u5883\u6539\u53d8\uff0c\u67e5\u770b\u9002\u5b9c\u5ea6\u54cd\u5e94\u53ca\u73af\u5883\u5206\u7ec4\u3002\u7ed3\u679c\u89e3\u91ca\u62df\u5408\u6a21\u578b\uff0c\u4e0d\u4ee3\u8868\u5df2\u8bc6\u522b\u7684\u751f\u6001\u56e0\u679c\u6548\u5e94\u3002")),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width = 3, class = "cast-panel",
        shiny::textInput("case_dir", "\u6848\u4f8b\u7ed3\u679c\u76ee\u5f55", value = case_dir),
        shiny::actionButton("load_case", "\u6253\u5f00 / \u5237\u65b0\u6848\u4f8b", class = "btn-outline-primary"),
        shiny::tags$hr(),
        shiny::selectInput("scenario", "\u672a\u6765\u60c5\u666f", choices = character()),
        shiny::actionButton("view", "\u6d4f\u89c8\u5df2\u6709\u60c5\u666f", class = "btn-outline-primary"),
        shiny::tags$hr(),
        shiny::radioButtons("mode", "\u73af\u5883\u6539\u53d8", choices = c(stats::setNames("scenario", "\u4f7f\u7528\u5df2\u6709\u60c5\u666f"),
          stats::setNames("shift", "\u5047\u8bbe\u5355\u4e2a\u6307\u6807\u6539\u53d8"))),
        shiny::selectInput("variables", "\u6539\u53d8\u54ea\u4e9b\u56e0\u5b50", choices = character()),
        shiny::uiOutput("variable_definition"),
        shiny::conditionalPanel("input.mode == 'shift'",
          shiny::numericInput("shift", "\u6539\u53d8\u5e45\u5ea6\uff08\u6307\u6807\u539f\u59cb\u6570\u503c\u5355\u4f4d\uff09", value = NA_real_),
          shiny::tags$p(class = "cast-note", "\u5355\u4e2a\u6d3e\u751f\u6307\u6807\u7684\u6539\u53d8\u4ec5\u4f5c\u654f\u611f\u6027\u5206\u6790\uff1b\u571f\u5730\u5229\u7528\u6bd4\u4f8b\u6216\u683c\u5c40\u53d8\u5316\u5e94\u4f18\u5148\u4f7f\u7528\u4e00\u81f4\u7684\u60c5\u666f\u56fe\u5c42\u3002")),
        shiny::selectInput("context", "\u5176\u4ed6\u56e0\u5b50\u7684\u80cc\u666f", choices = c(
          stats::setNames("current", "\u4fdd\u6301\u5f53\u524d\u503c"), stats::setNames("future", "\u4f7f\u7528\u9009\u5b9a\u672a\u6765\u60c5\u666f\u7684\u503c"))),
        shiny::tags$details(shiny::tags$summary("\u8ba1\u7b97\u8bbe\u7f6e"),
          shiny::numericInput("sample_n", "\u771f\u5b9e CV \u6837\u70b9\u4e0a\u9650", 5000L, min = 4L, step = 100L),
          shiny::numericInput("maxdepth", "\u5206\u7ec4\u6811\u6700\u5927\u6df1\u5ea6", 3L, min = 1L, max = 8L),
          shiny::numericInput("cores", "\u8ba1\u7b97\u7ebf\u7a0b", 4L, min = 1L),
          shiny::checkboxInput("maps", "\u8ba1\u7b97\u539f\u751f\u7f51\u683c\u53ca\u5730\u56fe\u9884\u89c8", FALSE)),
        shiny::actionButton("analyse", "\u8ba1\u7b97\u54cd\u5e94\u4e0e\u81ea\u52a8\u5206\u7ec4", class = "btn-primary"),
        shiny::tags$hr(), shiny::uiOutput("status"),
        shiny::tags$p(class = "cast-note", "\u957f\u8ba1\u7b97\u5728\u540e\u53f0\u6267\u884c\u3002\u8c03\u6574\u5de6\u4fa7\u9009\u9879\u540e\uff0c\u70b9\u51fb\u8ba1\u7b97\u624d\u66f4\u65b0\u7ed3\u679c\uff1b\u5730\u56fe\u6807\u9898\u59cb\u7ec8\u5bf9\u5e94\u5df2\u5b8c\u6210\u7684\u65b9\u6848\u3002")),
      shiny::mainPanel(width = 9,
        shiny::tabsetPanel(id = "result_tab",
          shiny::tabPanel("\u6548\u5e94\u5730\u56fe", value = "map",
            shiny::div(class = "cast-map-caption", shiny::uiOutput("map_caption")),
            shiny::selectInput("layer", "\u663e\u793a\u56fe\u5c42", choices = c(
              stats::setNames("total_change", "\u9002\u5b9c\u5ea6\u53d8\u5316"), stats::setNames("current", "\u6539\u53d8\u524d"), stats::setNames("future", "\u6539\u53d8\u540e"),
              stats::setNames("support_total", "\u8303\u56f4\u652f\u6301\uff081 / 0\uff09"), stats::setNames("total_sd", "\u6a21\u578b\u95f4\u5206\u6b67"),
              stats::setNames("response_group", "\u73af\u5883\u54cd\u5e94\u5206\u7ec4"))),
            shiny::checkboxInput("supported_only", "\u4ec5\u663e\u793a\u8303\u56f4\u68c0\u67e5\u901a\u8fc7\u7684\u4f4d\u7f6e", TRUE),
            leaflet::leafletOutput("map", height = 510),
            shiny::tags$p(class = "cast-note", "\u9884\u89c8\u7ecf\u8fc7\u964d\u91c7\u6837\uff1b\u70b9\u51fb\u8bfb\u53d6\u539f\u751f\u7f51\u683c\u3002\u8303\u56f4\u652f\u6301\u4e0d\u7b49\u4e8e\u8054\u5408\u652f\u6301\u6216\u56e0\u679c\u8bc6\u522b\uff0c\u6a21\u578b\u95f4\u5206\u6b67\u4e0d\u662f\u7f6e\u4fe1\u533a\u95f4\u3002"),
            shiny::fluidRow(shiny::column(5, shiny::h4("\u70b9\u51fb\u4f4d\u7f6e\u7684\u7ed3\u679c"), shiny::tableOutput("point")),
              shiny::column(7, shiny::selectInput("profile_driver", "\u67e5\u770b\u5c40\u5730\u9884\u6d4b\u5256\u9762", choices = character()),
                shiny::plotOutput("profile", height = 240))),
            shiny::tags$details(shiny::tags$summary("\u6240\u9009\u56e0\u5b50\u7684\u4e24\u79cd\u80cc\u666f\u54cd\u5e94"), shiny::tableOutput("point_driver")),
            shiny::tags$details(shiny::tags$summary("\u8be5\u4f4d\u7f6e\u7684\u73af\u5883\u503c"), shiny::tableOutput("environment"))),
          shiny::tabPanel("\u56e0\u5b50\u54cd\u5e94", value = "drivers", shiny::h4("\u540c\u4e00\u6279\u4f4d\u7f6e\u4e0a\u7684\u56e0\u5b50\u54cd\u5e94"),
            shiny::uiOutput("evidence_note"),
            shiny::div(style = "overflow-x:auto", shiny::tableOutput("drivers")),
            shiny::plotOutput("driver_plot", height = 360),
            shiny::tags$details(shiny::tags$summary("\u6307\u6807\u5b9a\u4e49\u4e0e\u6765\u6e90"),
              shiny::div(style = "overflow-x:auto", shiny::tableOutput("metadata"))),
            shiny::selectInput("focus_driver", "\u67e5\u770b\u56e0\u5b50\u5730\u56fe\uff08\u9ed8\u8ba4\u9009\u62e9\u6700\u5927\u54cd\u5e94\uff09", choices = character()),
            shiny::actionButton("map_driver", "\u8ba1\u7b97\u8be5\u56e0\u5b50\u7684\u4e24\u79cd\u80cc\u666f\u5730\u56fe", class = "btn-outline-primary"),
            shiny::tags$p(class = "cast-note", "\u9010\u6307\u6807\u66ff\u6362\u63cf\u8ff0\u6a21\u578b\u654f\u611f\u6027\uff1b\u4e0d\u540c\u56e0\u5b50\u7684\u5b9e\u9645\u6539\u53d8\u91cf\u548c\u5355\u4f4d\u4e0d\u540c\u3002\u5e45\u5ea6\u6392\u5e8f\u4e0d\u662f\u6cbb\u7406\u6210\u672c\u6216\u56e0\u679c\u6548\u5e94\u6392\u5e8f\uff0c\u5355\u56e0\u5b50\u54cd\u5e94\u4e4b\u548c\u4e5f\u4e0d\u5fc5\u7b49\u4e8e\u603b\u53d8\u5316\u3002\u672a\u5165\u6a21\u56e0\u5b50\u4e0d\u4f5c\u4f30\u8ba1\u3002")),
          shiny::tabPanel("\u54cd\u5e94\u5206\u7ec4", value = "groups", shiny::h4("\u81ea\u52a8\u8bc6\u522b\u7684\u73af\u5883\u5b50\u7fa4"),
            shiny::uiOutput("group_note"), shiny::tableOutput("groups"),
            shiny::plotOutput("group_plot", height = 270),
            shiny::tags$p(class = "cast-note", "\u5206\u7ec4\u6c47\u603b\u51bb\u7ed3\u6a21\u578b\u7684\u54cd\u5e94\uff0c\u53ef\u80fd\u53ea\u6709\u4e00\u4e2a\u7ec4\u3002\u7a7a\u95f4\u6298\u7528\u4e8e\u9009\u62e9\u6811\u590d\u6742\u5ea6\uff1b\u5206\u7ec4\u4e0d\u662f CATE\uff0c\u4e5f\u4e0d\u6784\u6210\u72ec\u7acb\u751f\u6001\u9a8c\u8bc1\u3002\u5b50\u7fa4\u6620\u5c04\u540e\u4e0d\u4e00\u5b9a\u5730\u7406\u8fde\u7eed\u3002"),
            shiny::tags$details(shiny::tags$summary("\u7a7a\u95f4\u6298\u8c03\u53c2\u4e0e\u8bef\u5dee"), shiny::tableOutput("tuning"))),
          shiny::tabPanel("\u6a21\u578b\u4e0e\u6b63\u5f0f\u56fe\u4ef6", value = "figures", shiny::h4("\u4fdd\u5b58\u7684\u7559\u51fa\u8bc4\u4ef7"),
            shiny::tableOutput("evaluation"), shiny::h4("\u539f\u6709\u8bba\u6587\u56fe\u4ef6"),
            shiny::selectInput("figure", "\u56fe\u4ef6", choices = character()),
            shiny::imageOutput("publication", height = "auto")),
          shiny::tabPanel("\u8fd0\u884c\u4e0e\u5bfc\u51fa", value = "export", shiny::h4("\u5df2\u914d\u7f6e\u7684\u5b8c\u6574\u5de5\u4f5c\u6d41"),
            shiny::tags$p("\u4f7f\u7528\u5355\u7269\u79cd\u914d\u7f6e\u6587\u4ef6\u4e2d\u7684\u79d1\u5b66\u8bbe\u7f6e\uff0c\u5728\u72ec\u7acb\u8fdb\u7a0b\u8fd0\u884c\u3002\u5df2\u6709\u6b63\u786e\u7ed3\u679c\u6309\u7ba1\u7ebf\u89c4\u5219\u7eed\u8dd1\uff1b\u5de5\u4f5c\u6d41\u4e0d\u4f1a\u56e0\u6253\u5f00\u5e94\u7528\u800c\u81ea\u52a8\u8fd0\u884c\u3002"),
            shiny::uiOutput("workflow_path"),
            shiny::actionButton("run_workflow", "\u8fd0\u884c\u5df2\u914d\u7f6e\u5de5\u4f5c\u6d41", class = "btn-outline-primary"),
            shiny::h4("\u5f53\u524d\u5206\u6790\u5bfc\u51fa"),
            shiny::downloadButton("download_points", "\u771f\u5b9e\u6837\u70b9\u54cd\u5e94\u8868"),
            shiny::downloadButton("download_bundle", "\u5206\u6790\u5305\uff08\u542b\u6805\u683c\u4e0e\u590d\u73b0\u811a\u672c\uff09"),
            shiny::downloadButton("download_figure", "\u6240\u9009\u6b63\u5f0f\u56fe\u4ef6"),
            shiny::tags$details(shiny::tags$summary("\u540e\u53f0\u8ba1\u7b97\u65e5\u5fd7"), shiny::uiOutput("log")))
        )
      )
    )
  )
}

.cast_app_server <- function(case_dir, workflow, state = NULL) {
  force(case_dir); force(workflow); force(state)
  function(input, output, session) {
    project <- shiny::reactiveVal(NULL)
    restore_state <- state
    active <- shiny::reactiveVal(NULL)
    active_request <- shiny::reactiveVal(NULL)
    point <- shiny::reactiveVal(NULL)
    status <- shiny::reactiveVal("\u6253\u5f00\u5df2\u5b8c\u6210\u6848\u4f8b\uff0c\u6216\u8fd0\u884c\u5df2\u914d\u7f6e\u7684\u5de5\u4f5c\u6d41\u3002")
    log_text <- shiny::reactiveVal("")
    job <- NULL; logfile <- NULL; pending <- NULL
    notify_error <- function(e) {
      status(conditionMessage(e))
      shiny::showNotification(conditionMessage(e), type = "error", duration = 12)
    }
    start <- function(request) {
      if (!is.null(job) && job$is_alive()) {
        shiny::showNotification("\u5df2\u6709\u540e\u53f0\u4efb\u52a1\uff0c\u8bf7\u7b49\u5f85\u5b83\u5b8c\u6210\u3002", type = "message")
        return(invisible(FALSE))
      }
      logfile <<- tempfile("cast-app-", fileext = ".log")
      pending <<- request
      job <<- .cast_app_start(request, logfile)
      status("\u540e\u53f0\u8ba1\u7b97\u4e2d\uff1b\u53ef\u4ee5\u7ee7\u7eed\u6d4f\u89c8\u5df2\u5b8c\u6210\u7ed3\u679c\u3002")
      invisible(TRUE)
    }
    load_project <- function(view = TRUE) {
      p <- .cast_explore_project(input$case_dir %||% case_dir)
      project(p); active(NULL); point(NULL)
      scenarios <- names(p$comparison$results)
      shiny::updateSelectInput(session, "scenario", choices = scenarios, selected = p$case$report_scenario)
      groups <- names(p$comparison$groups)
      choices <- c(stats::setNames(".all", "\u5168\u90e8\u60c5\u666f\u56e0\u5b50"),
        stats::setNames(paste0(".group:", groups), paste0("\u7ec4\uff1a", groups)),
        stats::setNames(p$fit$env_vars, p$fit$env_vars))
      shiny::updateSelectInput(session, "variables", choices = choices)
      shiny::updateSelectInput(session, "profile_driver", choices = p$fit$env_vars)
      files <- list.files(file.path(p$directory, "figures"), pattern = "\\.png$", full.names = FALSE)
      shiny::updateSelectInput(session, "figure", choices = files)
      status(paste("\u5df2\u6253\u5f00", p$case$species, "\uff1b\u8bfb\u53d6\u5df2\u6709\u6a21\u578b\u4e0e\u7ed3\u679c\u3002"))
      if (view) start(list(operation = "view", case_dir = p$directory, scenario = p$case$report_scenario))
      invisible(p)
    }
    session$onFlushed(function() {
      if (file.exists(file.path(case_dir, "case.rds"))) tryCatch(shiny::isolate({
        load_project(view = is.null(restore_state))
        if (!is.null(restore_state)) {
          request <- restore_state$request; request$case_dir <- case_dir
          shiny::updateSelectInput(session, "scenario", selected = request$scenario)
          shiny::updateSelectInput(session, "context", selected = request$context %||% "current")
          vars <- request$variables %||% project()$fit$env_vars
          choice <- if (setequal(vars, project()$fit$env_vars)) ".all" else if (length(vars) == 1L) vars else {
            groups <- project()$comparison$groups
            matched <- names(groups)[vapply(groups, function(x) setequal(x, vars), logical(1))]
            if (!length(matched)) stop("Saved predictor selection no longer matches the case's groups.")
            paste0(".group:", matched[1L])
          }
          shiny::updateSelectInput(session, "variables", selected = choice)
          shiny::updateRadioButtons(session, "mode", selected = if (is.null(request$shifts)) "scenario" else "shift")
          if (!is.null(request$shifts)) shiny::updateNumericInput(session, "shift", value = unname(request$shifts[1L]))
          shiny::updateSelectInput(session, "profile_driver", selected = restore_state$profile_driver)
          shiny::updateCheckboxInput(session, "supported_only", value = isTRUE(restore_state$supported_only))
          if (!is.null(request$sample_n)) shiny::updateNumericInput(session, "sample_n", value = request$sample_n)
          if (!is.null(request$maxdepth)) shiny::updateNumericInput(session, "maxdepth", value = request$maxdepth)
          if (!is.null(request$maps)) shiny::updateCheckboxInput(session, "maps", value = request$maps)
          start(request)
        }
      }), error = notify_error)
    }, once = TRUE)
    shiny::observeEvent(input$load_case, {
      if (!is.null(job) && job$is_alive()) {
        shiny::showNotification("\u8bf7\u7b49\u5f85\u5f53\u524d\u4efb\u52a1\u5b8c\u6210\u540e\u518d\u5207\u6362\u6848\u4f8b\u3002", type = "message")
      } else tryCatch(load_project(), error = notify_error)
    })
    shiny::observeEvent(input$view, {
      shiny::req(project(), input$scenario)
      tryCatch(start(list(operation = "view", case_dir = project()$directory,
        scenario = input$scenario)), error = notify_error)
    })
    shiny::observeEvent(input$analyse, {
      shiny::req(project(), input$scenario, input$variables)
      tryCatch({
        p <- project(); choice <- input$variables
        variables <- if (choice == ".all") NULL else if (startsWith(choice, ".group:"))
          p$comparison$groups[[substring(choice, 8L)]] else choice
        shifts <- NULL
        if (identical(input$mode, "shift")) {
          if (!choice %in% p$fit$env_vars) stop("\u5047\u8bbe\u6539\u53d8\u6a21\u5f0f\u8bf7\u9009\u62e9\u4e00\u4e2a\u5177\u4f53\u6307\u6807\u3002")
          if (is.null(input$shift) || !is.finite(input$shift) || input$shift == 0)
            stop("\u8bf7\u8f93\u5165\u6709\u9650\u3001\u975e\u96f6\u7684\u539f\u59cb\u5355\u4f4d\u6539\u53d8\u5e45\u5ea6\u3002")
          shifts <- stats::setNames(input$shift, choice)
        }
        start(list(operation = "response", case_dir = p$directory, scenario = input$scenario,
          variables = variables, shifts = shifts, context = input$context,
          sample_n = input$sample_n, maxdepth = input$maxdepth, maps = input$maps,
          resources = list(cores = input$cores, workers = 1L, memory_mb = 2048)))
      }, error = notify_error)
    })
    shiny::observeEvent(input$map_driver, {
      shiny::req(active()$evidence, input$focus_driver)
      s <- active()$settings
      tryCatch(start(list(operation = "response", case_dir = active()$case_dir,
        scenario = s$scenario, variables = s$variables, shifts = s$shifts,
        context = s$context, sample_n = s$sample_n, maxdepth = s$maxdepth,
        maps = TRUE, focus = input$focus_driver,
        resources = list(cores = input$cores, workers = 1L, memory_mb = 2048))), error = notify_error)
    })
    shiny::observeEvent(input$run_workflow, {
      tryCatch({
        if (is.null(workflow) || !file.exists(workflow)) stop("\u8bf7\u901a\u8fc7 05_SHINY.R \u6253\u5f00\uff0c\u6216\u5411 cast_app() \u63d0\u4f9b workflow \u8def\u5f84\u3002")
        target <- normalizePath(input$case_dir %||% case_dir, winslash = "/", mustWork = FALSE)
        p <- project()
        environment <- c(CAST_SPECIES = if (!is.null(p)) p$case$species else basename(target),
          CAST_RUN_ID = basename(dirname(target)), CAST_OUTPUT_ROOT = dirname(dirname(target)),
          CAST_EXAMPLE_ROOT = dirname(workflow))
        if (!is.null(p)) environment <- c(environment,
          CAST_DATA_ROOT = p$case$cfg$nsdm_data_root,
          CAST_OCCURRENCE_ROOT = dirname(p$case$cfg$occurrence_file))
        start(list(operation = "workflow", workflow = workflow,
          case_dir = target, environment = environment))
      }, error = notify_error)
    })
    shiny::observe({
      shiny::invalidateLater(750, session)
      if (is.null(job)) return()
      if (!is.null(logfile) && file.exists(logfile))
        log_text(paste(utils::tail(readLines(logfile, warn = FALSE, encoding = "UTF-8"), 60L), collapse = "\n"))
      if (job$is_alive()) return()
      finished <- job; job <<- NULL
      tryCatch({
        value <- finished$get_result()
        if (identical(pending$operation, "workflow")) {
          load_project()
        } else {
          active(value); active_request(pending); point(NULL)
          if (!is.null(value$evidence)) shiny::updateSelectInput(session, "focus_driver",
            choices = value$evidence$summary$driver, selected = value$settings$focus)
          shiny::updateTabsetPanel(session, "result_tab", selected =
            if (!is.null(value$evidence) && !isTRUE(pending$maps)) "drivers" else "map")
          layers <- c("total_change", "current", "future", "support_total", "total_sd")
          labels <- c("\u603b\u9002\u5b9c\u5ea6\u53d8\u5316", "\u6539\u53d8\u524d", "\u6539\u53d8\u540e",
            "\u603b\u53d8\u5316\u8303\u56f4\u652f\u6301", "\u603b\u53d8\u5316\u6a21\u578b\u5206\u6b67")
          focus <- value$settings$focus
          if (!is.null(focus) && !is.null(value$raster_file)) {
            layers <- c(layers, paste0(c("effect_", "retained_", "context_direction_",
              "context_gap_", "effect_agreement_", "retained_agreement_"), focus))
            labels <- c(labels, paste(focus, c("| \u6bd4\u8f83\u524d\u80cc\u666f\u54cd\u5e94", "| \u6bd4\u8f83\u540e\u80cc\u666f\u54cd\u5e94",
              "| \u80cc\u666f\u65b9\u5411\u4e00\u81f4\u6027", "| \u4e24\u80cc\u666f\u54cd\u5e94\u5dee",
              "| \u6bd4\u8f83\u524d\u80cc\u666f\u6a21\u578b\u540c\u5411\u6bd4\u4f8b", "| \u6bd4\u8f83\u540e\u80cc\u666f\u6a21\u578b\u540c\u5411\u6bd4\u4f8b")))
          }
          if (!is.null(value$group_file)) {
            layers <- c(layers, "response_group"); labels <- c(labels, "\u73af\u5883\u54cd\u5e94\u5206\u7ec4")
          }
          selected <- if (!is.null(focus) && !is.null(value$raster_file)) paste0("effect_", focus) else "total_change"
          shiny::updateSelectInput(session, "layer", choices = stats::setNames(layers, labels), selected = selected)
          if (!is.null(restore_state)) {
            if (restore_state$layer %in% layers) shiny::updateSelectInput(session, "layer", selected = restore_state$layer)
            if (!is.null(restore_state$center) && !is.null(restore_state$zoom))
              leaflet::setView(leaflet::leafletProxy("map", session),
                lng = restore_state$center$lng, lat = restore_state$center$lat, zoom = restore_state$zoom)
            if (!is.null(restore_state$click)) point(.cast_app_point(project(), value,
              restore_state$click$lng, restore_state$click$lat, restore_state$profile_driver))
            restore_state <<- NULL
          }
          status("\u8ba1\u7b97\u5b8c\u6210\uff1b\u5f53\u524d\u663e\u793a\u7ed3\u679c\u5bf9\u5e94\u4e0b\u65b9\u5730\u56fe\u6807\u9898\u3002")
        }
      }, error = notify_error)
    })
    session$onSessionEnded(function() {
      if (!is.null(job) && job$is_alive()) job$kill()
    })
    output$status <- shiny::renderUI(shiny::tags$span(status()))
    output$log <- shiny::renderUI(shiny::tags$pre(log_text()))
    output$workflow_path <- shiny::renderUI(shiny::tags$pre(workflow %||% "No configured workflow."))
    output$map_caption <- shiny::renderUI({
      x <- active(); if (is.null(x)) return(shiny::tags$span("\u5c1a\u65e0\u5df2\u663e\u793a\u7684\u65b9\u6848\u3002"))
      vars <- x$settings$variables
      shiny::tags$span(paste(x$settings$scenario, "|", if (is.null(x$settings$shifts)) "\u60c5\u666f\u6bd4\u8f83" else "\u5047\u8bbe\u6307\u6807\u6539\u53d8",
        "|", paste(vars, collapse = ", "), "|", x$settings$context))
    })
    output$map <- leaflet::renderLeaflet({
      m <- leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE))
      m <- leaflet::addTiles(m, options = leaflet::tileOptions(opacity = .3))
      leaflet::setView(m, lng = 104, lat = 35, zoom = 4)
    })
    shiny::observe({
      x <- active(); layer <- input$layer
      proxy <- leaflet::leafletProxy("map", session)
      proxy <- leaflet::clearImages(proxy); proxy <- leaflet::clearControls(proxy)
      if (is.null(x) || is.null(layer)) return()
      if (is.null(x$preview) || !file.exists(x$preview)) return()
      r <- terra::rast(x$preview)
      if (!layer %in% names(r)) return()
      support <- if (startsWith(layer, "effect_")) paste0("support_effect_", x$settings$focus) else
        if (startsWith(layer, "retained_")) paste0("support_retained_", x$settings$focus) else "support_total"
      if ((isTRUE(input$supported_only) || startsWith(layer, "context_")) &&
          !layer %in% c("support_total", "response_group"))
        r[[layer]] <- terra::mask(r[[layer]], r[[support]], maskvalues = 0)
      r <- r[[layer]]
      direction <- startsWith(layer, "context_direction_")
      signed <- layer == "total_change" || startsWith(layer, "context_gap_") ||
        layer %in% paste0(c("effect_", "retained_"), x$settings$focus)
      if (signed) {
        # Shared fixed score scale permits comparisons between scenarios.
        limits <- if (startsWith(layer, "context_gap_")) c(-2, 2) else c(-1, 1)
        palette <- leaflet::colorNumeric(c("#B35806", "#F7F7F7", "#2166AC"), limits, na.color = "#00000000")
      } else if (direction) {
        ids <- c(-1, 0, 1, 2, 3)
        palette <- leaflet::colorFactor(c("#B35806", "#D9D9D9", "#2166AC", "#8C2981", "#F2A541"),
          domain = ids, levels = ids, na.color = "#00000000")
      } else if (layer == "response_group") {
        ids <- x$groups$summary$group
        palette <- leaflet::colorFactor(grDevices::hcl.colors(max(3L, length(ids)), "Dark 3"), domain = ids,
          na.color = "#00000000")
      } else palette <- leaflet::colorNumeric(c("#F5F7FA", "#468FAF", "#014F86"), c(0, 1), na.color = "#00000000")
      proxy <- leaflet::addRasterImage(proxy, r, colors = palette, project = FALSE,
        opacity = .85, layerId = "response", maxBytes = 8 * 1024 * 1024)
      if (direction) leaflet::addLegend(proxy,
        colors = palette(ids), labels = c("Both negative", "Both zero", "Both positive", "Opposite signs", "Zero in one background"),
        title = "Background directions", opacity = .85) else
        leaflet::addLegend(proxy, pal = palette, values = if (layer == "response_group") ids else
          if (signed) limits else c(0, 1), title = if (signed) "Suitability response" else layer, opacity = .85)
    })
    shiny::observeEvent(list(input$map_click, input$profile_driver, active()), {
      shiny::req(project(), active(), input$profile_driver)
      click <- input$map_click %||% if (!is.null(point())) list(lng = point()$longitude, lat = point()$latitude)
      shiny::req(click)
      tryCatch(point(.cast_app_point(project(), active(), click$lng,
        click$lat, input$profile_driver)), error = function(e) {
          point(NULL); shiny::showNotification(conditionMessage(e), type = "message", duration = 6)
        })
    }, ignoreInit = TRUE)
    shiny::observe({
      p <- point()
      proxy <- leaflet::clearGroup(leaflet::leafletProxy("map", session), "query")
      if (!is.null(p)) leaflet::addCircleMarkers(proxy, lng = p$longitude,
        lat = p$latitude, radius = 5, color = "#17252A", fillColor = "#FFFFFF",
        fillOpacity = 1, weight = 2, group = "query")
    })
    output$point <- shiny::renderTable({
      shiny::validate(shiny::need(point(), "\u8bf7\u5728\u5206\u6790\u57df\u5185\u70b9\u51fb\u5730\u56fe\u3002"))
      p <- point(); z <- unlist(p$values)
      data.frame(item = c("Longitude", "Latitude", "Suitability before", "Suitability after", "Change",
        "Range support (1/0)", "Model disagreement", "Response group"),
        value = c(p$longitude, p$latitude, z, p$group))
    }, digits = 5)
    output$environment <- shiny::renderTable({ shiny::req(point()); point()$environment }, digits = 5)
    output$metadata <- shiny::renderTable({ shiny::req(project()); project()$metadata })
    output$variable_definition <- shiny::renderUI({
      shiny::req(project(), input$variables)
      m <- project()$metadata
      m <- m[m$layer == input$variables, , drop = FALSE]
      if (!nrow(m) || is.na(m$definition[1L])) return(NULL)
      shiny::tags$p(class = "cast-note", paste(m$definition[1L], "|", m$units[1L]))
    })
    output$profile <- shiny::renderPlot({
      shiny::req(point())
      p <- point()
      graphics::plot(p$curve$environment, p$curve$suitability, type = "l", col = "#176B87", lwd = 2,
        xlab = p$driver, ylab = "Relative suitability", main = "Fixed-other-predictor profile")
      value <- p$environment$before[p$environment$variable == p$driver]
      graphics::abline(v = value, lty = 2, col = "#B35806")
    })
    output$group_note <- shiny::renderUI({
      shiny::req(active()); shiny::tags$span(active()$group_note %||% paste(
        "\u652f\u6301\u6837\u70b9\uff1a", active()$groups$settings$n_supported,
        "\uff1b\u7a7a\u95f4\u6298\uff1a", active()$groups$settings$folds, "\uff1b\u91c7\u7528\u4e00\u6807\u51c6\u8bef\u5dee\u526a\u679d\u3002"))
    })
    driver_summary <- shiny::reactive({
      shiny::req(active()$evidence)
      active()$evidence$summary
    })
    output$evidence_note <- shiny::renderUI({
      shiny::req(active()$evidence)
      p <- active()$evidence$population
      shiny::tags$p(paste("\u5171\u540c\u652f\u6301\u4f4d\u7f6e\uff1a", p$n_common, "/", p$n_input,
        "\u3002\u6240\u6709\u56e0\u5b50\u4f7f\u7528\u540c\u4e00\u5206\u6bcd\u3002\u672a\u53d8\u5316\u6307\u6807\u4e0d\u53c2\u4e0e\u6392\u540d\uff1b\u4e0d\u8868\u793a\u8be5\u6307\u6807\u4e0d\u91cd\u8981\u3002",
        "\u6a21\u578b\u540c\u5411\u5217\u8868\u793a\u5168\u90e8\u6a21\u578b\u540c\u5411\u7684\u4f4d\u7f6e\u6bd4\u4f8b\uff1b\u4e0d\u662f\u7f6e\u4fe1\u5ea6\u3002"))
    })
    output$drivers <- shiny::renderTable({
      d <- driver_summary()[, c("rank", "driver", "status", "n_changed_common", "mean_only", "mean_retained",
        "context_reversal_fraction", "all_models_agree_only", "all_models_agree_retained")]
      labels <- c(changed = "Changed", unchanged_scenario = "Unchanged",
        unchanged_on_common_domain = "Unchanged here", no_common_support = "No support", zero_model_response = "Zero response",
        missing_predictor_pairs = "Missing pairs", unchanged_known_pairs = "Unchanged where known")
      d$status <- unname(labels[d$status])
      names(d) <- c("Rank", "Driver", "Status", "Changed n", "Before mean", "After mean",
        "Flip fraction", "All agree: before", "All agree: after")
      d
    }, digits = 5)
    output$point_driver <- shiny::renderTable({
      shiny::req(point()$driver_response)
      data.frame(item = c("Before-background response", "After-background response", "Model agreement: before",
        "Model agreement: after", "Background direction class", "Background response gap"),
        value = as.numeric(point()$driver_response[1L, ]))
    }, digits = 5)
    output$driver_plot <- shiny::renderPlot({
      d <- driver_summary()
      d <- d[!is.na(d$rank), , drop = FALSE]
      shiny::req(nrow(d), any(is.finite(d$mean_only)))
      old <- graphics::par(mar = c(4, 17, 2, 1)); on.exit(graphics::par(old))
      graphics::barplot(rbind(rev(d$mean_only), rev(d$mean_retained)), names.arg = rev(d$driver), beside = TRUE,
        horiz = TRUE, las = 1, cex.names = .7,
        col = c("#176B87", "#D17A22"), legend.text = c("Before comparison background", "After comparison background"),
        args.legend = list(x = "bottomright", bty = "n", cex = .8),
        xlab = "Mean fitted response on the same supported saved locations")
      graphics::abline(v = 0, col = "grey60")
    })
    output$groups <- shiny::renderTable({ shiny::req(active()$groups); active()$groups$summary }, digits = 5)
    output$tuning <- shiny::renderTable({ shiny::req(active()$groups); active()$groups$cv }, digits = 6)
    output$group_plot <- shiny::renderPlot({
      shiny::req(active()$groups)
      d <- active()$groups$summary
      graphics::barplot(d$mean_response, names.arg = d$group, col = ifelse(d$mean_response >= 0,
        "#2166AC", "#B35806"), xlab = "Environmental response group", ylab = "Mean fitted response")
      graphics::abline(h = 0, col = "grey60")
    })
    output$evaluation <- shiny::renderTable({
      shiny::req(project())
      tables <- Filter(is.data.frame, project()$run$eval)
      shiny::req(length(tables)); tables[[1L]]
    }, digits = 4)
    publication <- shiny::reactive({
      shiny::req(project(), input$figure)
      # Dropdown values must remain a basename in the case's figure directory.
      shiny::req(identical(basename(input$figure), input$figure))
      path <- file.path(project()$directory, "figures", input$figure)
      shiny::req(file.exists(path)); path
    })
    output$publication <- shiny::renderImage({
      list(src = publication(), contentType = "image/png", width = "100%", alt = input$figure)
    }, deleteFile = FALSE)
    output$download_points <- shiny::downloadHandler(filename = function() "castSDM_responses.csv",
      content = function(file) {
        shiny::req(active()$points_file); file.copy(active()$points_file, file, overwrite = TRUE)
      })
    output$download_figure <- shiny::downloadHandler(filename = function() input$figure,
      content = function(file) file.copy(publication(), file, overwrite = TRUE))
    output$download_bundle <- shiny::downloadHandler(filename = function() "castSDM_exploration.zip",
      content = function(file) {
        shiny::req(active()); directory <- active()$directory
        files <- list.files(directory, full.names = FALSE)
        files <- files[!grepl("\\.lock$|^\\.", files)]
        check_suggested("zip", "for the downloadable analysis bundle")
        bundle <- tempfile("cast-view-export-"); dir.create(bundle)
        file.copy(file.path(directory, files), bundle, recursive = TRUE)
        state <- list(schema = "cast_app_state/1", request = active_request(), layer = input$layer,
          supported_only = input$supported_only, center = input$map_center, zoom = input$map_zoom,
          click = if (!is.null(point())) list(lng = point()$longitude, lat = point()$latitude),
          profile_driver = input$profile_driver)
        saveRDS(state, file.path(bundle, "view_state.rds"))
        writeLines(c("# Run from the extracted bundle; original real-case inputs must remain available.",
          "library(castSDM)", "state <- readRDS('view_state.rds')",
          "cast_app(state$request$case_dir, state = state)"), file.path(bundle, "reproduce_view.R"))
        zip::zipr(file, files = list.files(bundle, full.names = TRUE), root = bundle)
      })
  }
}

#' Open the Environmental Response Shiny Application
#'
#' Browse completed cases, query native-grid responses, explore explicit
#' environmental changes and descriptive response groups, and export results.
#' A configured workflow can be submitted manually to an isolated background
#' R process. Opening the app never starts model fitting.
#'
#' @param case_dir Completed single-species workflow directory. It may be a
#'   not-yet-created output directory when `workflow` is supplied.
#' @param workflow Optional path to a configured `00_RUN_ALL.R` entry.
#' @param run Run the application, default `TRUE`. `FALSE` returns a Shiny app
#'   object without starting a server, computations or model fitting.
#' @param launch.browser Open a browser when running, default `interactive()`.
#' @param host,port Local Shiny server settings. Host defaults to `127.0.0.1`.
#' @param state Optional view-state list from an exported `view_state.rds`.
#'   Restores a saved view/response request, layer, map position and point query.
#'   Original case files must remain available. It never submits a workflow.
#' @return Invisibly the result of `shiny::runApp()`, or the app when `run = FALSE`.
#' @details Shiny, bslib, leaflet and callr are optional interface dependencies.
#'   Background jobs receive paths/settings and pin the calling installation.
#'   Map previews are downsampled; point queries use native-grid results.
#'   Model responses, descriptive groups and model disagreement do not identify
#'   ecological causal effects, CATE, probability changes or sampling intervals.
#'   Original publication figures remain available separately.
#' @export
cast_app <- function(case_dir, workflow = NULL, run = TRUE,
                      launch.browser = interactive(), host = "127.0.0.1", port = NULL,
                      state = NULL) {
  for (pkg in c("shiny", "bslib", "leaflet", "callr", "terra", "rpart"))
    check_suggested(pkg, "for the environmental response application")
  if (!is.character(case_dir) || length(case_dir) != 1L || is.na(case_dir) || !nzchar(case_dir))
    cli::cli_abort("Supply a case directory.")
  case_dir <- normalizePath(case_dir, winslash = "/", mustWork = FALSE)
  if (!is.null(workflow)) workflow <- normalizePath(workflow, winslash = "/", mustWork = TRUE)
  if (!is.null(state) && (!is.list(state) || !identical(state$schema, "cast_app_state/1") ||
      !is.list(state$request) || length(state$request$operation) != 1L ||
      !state$request$operation %in% c("view", "response") || length(state$layer) != 1L ||
      length(state$profile_driver) != 1L))
    cli::cli_abort("Supply a castSDM view/response state exported by the app.")
  app <- shiny::shinyApp(.cast_app_ui(case_dir), .cast_app_server(case_dir, workflow, state))
  if (!run) return(app)
  shiny::runApp(app, launch.browser = launch.browser, host = host, port = port)
}
