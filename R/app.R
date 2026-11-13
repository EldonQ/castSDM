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
  name <- layer
  if (identical(result$schema, "cast_app_view/1")) {
    aliases <- c(future = "joint", total_change = "joint_change", support_total = "support", total_sd = "sd_joint_change")
    at <- match(layer, names(aliases)); found <- !is.na(at)
    name[found] <- unname(aliases[at[found]])
  }
  out <- r[[name]]; names(out) <- layer
  out
}

.cast_app_point <- function(project, result, longitude, latitude, driver) {
  settings <- result$settings
  source_paths <- c(project$pool$path, project$scenarios$path[project$scenarios$scenario == settings$scenario])
  if (any(!file.exists(source_paths)) && !is.null(result$raster_file)) {
    # Native results remain inspectable when historical source inputs are archived.
    basic <- c("current", "future", "total_change", "support_total", "total_sd")
    details <- if (!is.null(settings$focus)) paste0(c("effect_", "retained_", "effect_agreement_",
      "retained_agreement_", "context_direction_", "context_gap_", "support_effect_", "support_retained_"), settings$focus)
    has_details <- length(details) == 8L && all(details %in% names(terra::rast(result$raster_file)))
    raster <- .cast_app_layer(result, c(basic, if (has_details) details))
    xy <- terra::project(matrix(c(longitude, latitude), ncol = 2L), "EPSG:4326", terra::crs(raster))
    native <- terra::extract(raster, xy); values <- native[, basic, drop = FALSE]
    if (!is.finite(values$current)) cli::cli_abort("Click inside the saved analysis domain.")
    return(list(values = values, group = NA_integer_, environment = NULL, curve = NULL,
      driver_response = if (has_details) native[, details[1:6], drop = FALSE],
      driver_support = if (has_details) native[, details[7:8], drop = FALSE],
      driver = driver, longitude = longitude, latitude = latitude,
      note = "Native saved values only: environmental sources are unavailable for profiles or new responses."))
  }
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
  detail_fields <- if (!is.null(focus)) paste0(c("effect_", "retained_", "effect_agreement_",
    "retained_agreement_", "context_direction_", "context_gap_", "support_effect_", "support_retained_"), focus)
  details <- if (length(detail_fields)) paired$data[, detail_fields, drop = FALSE] else NULL
  # Read the native analysis mask/results, not a smoothed display pixel.
  if (!is.null(result$raster_file)) {
    available <- intersect(detail_fields, names(terra::rast(result$raster_file)))
    native <- terra::extract(.cast_app_layer(result, c(names(values), available)), xy)
    values <- native[, names(values), drop = FALSE]
    if (length(available)) details[, available] <- native[, available, drop = FALSE]
    domain <- values$current
  } else domain <- terra::extract(terra::rast(file.path(project$case$preparation, "analysis_domain_mask.tif")), xy)[1L, 1L]
  if (!is.finite(domain) || (is.null(result$raster_file) && domain == 0))
    cli::cli_abort("\u8be5\u4f4d\u7f6e\u4e0d\u5728\u8fd9\u9879\u5206\u6790\u7684\u5171\u540c\u6709\u6548\u57df\u5185\u3002")
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
    driver_response = if (!is.null(details)) details[, seq_len(6L), drop = FALSE],
    driver_support = if (!is.null(details)) details[, 7:8, drop = FALSE],
    environment = data.frame(variable = vars, before = unlist(before[vars], use.names = FALSE),
      after = unlist(after[vars], use.names = FALSE), row.names = NULL),
    curve = data.frame(environment = grid, suitability = curve), driver = driver,
    longitude = longitude, latitude = latitude)
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
    job_busy <- shiny::reactiveVal(FALSE)
    status_kind <- shiny::reactiveVal("ready")
    job <- NULL; logfile <- NULL; pending <- NULL
    notify_error <- function(e) {
      log_text(paste(log_text(), conditionMessage(e), sep = "\n"))
      while (inherits(e$parent, "condition")) e <- e$parent
      message <- substr(strsplit(conditionMessage(e), "\n", fixed = TRUE)[[1L]][1L], 1L, 180L)
      status(message); status_kind("error")
      shiny::showNotification(message, type = "error", duration = 10)
    }
    start <- function(request) {
      if (!is.null(job)) {
        shiny::showNotification("\u5df2\u6709\u540e\u53f0\u4efb\u52a1\uff0c\u8bf7\u7b49\u5f85\u5b83\u5b8c\u6210\u3002", type = "message")
        return(invisible(FALSE))
      }
      logfile <<- tempfile("cast-app-", fileext = ".log")
      pending <<- request
      job <<- .cast_app_start(request, logfile)
      job_busy(TRUE); status_kind("busy")
      status("\u8ba1\u7b97\u4e2d\u2026")
      invisible(TRUE)
    }
    load_project <- function(view = TRUE, directory = NULL) {
      p <- .cast_explore_project(directory %||% input$case_dir %||% case_dir)
      project(p); active(NULL); point(NULL)
      scenarios <- names(p$comparison$results)
      shiny::updateSelectInput(session, "scenario", choices = stats::setNames(scenarios, .cast_app_scenario_label(scenarios)), selected = p$case$report_scenario)
      groups <- names(p$comparison$groups)
      choices <- c(stats::setNames(".all", "\u5168\u90e8\u60c5\u666f\u56e0\u5b50"),
        stats::setNames(paste0(".group:", groups), paste0("\u7ec4\uff1a", groups)),
        stats::setNames(p$fit$env_vars, .cast_app_driver_label(p, p$fit$env_vars)))
      shiny::updateSelectInput(session, "variables", choices = choices)
      shiny::updateSelectInput(session, "profile_driver", choices = stats::setNames(p$fit$env_vars, .cast_app_driver_label(p, p$fit$env_vars)))
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
      if (!is.null(job)) {
        shiny::showNotification("\u8bf7\u7b49\u5f85\u5f53\u524d\u4efb\u52a1\u5b8c\u6210\u540e\u518d\u5207\u6362\u6848\u4f8b\u3002", type = "message")
      } else tryCatch(load_project(), error = notify_error)
    })
    shiny::observeEvent(input$view, {
      shiny::req(project(), input$scenario)
      tryCatch(start(list(operation = "view", case_dir = project()$directory,
        scenario = input$scenario)), error = notify_error)
    })
    shiny::observeEvent(input$cancel_job, {
      if (!is.null(job) && job$is_alive()) job$kill_tree()
      job <<- NULL; pending <<- NULL; job_busy(FALSE); status_kind("ready")
      status("\u5df2\u505c\u6b62\uff0c\u5df2\u6709\u7ed3\u679c\u4fdd\u7559\u3002")
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
    shiny::observeEvent(input$driver_pick, {
      shiny::req(active()$evidence, input$driver_pick %in% active()$settings$variables)
      shiny::updateSelectInput(session, "focus_driver", selected = input$driver_pick)
    })
    setup_result <- shiny::reactiveVal(NULL)
    setup_request <- function(run) {
      list(operation = "setup", path = input$setup_path, data_dir = input$setup_data,
        occurrence_file = input$setup_occurrence, species = input$setup_species,
        settings = list(domain = input$setup_domain, input_type = input$setup_input,
          lulc_product = input$setup_lulc, models = input$setup_models, cv_k = input$setup_cv,
          resources = list(cores = input$setup_cores, workers = input$setup_workers, memory_mb = input$setup_memory)), run = run)
    }
    shiny::observeEvent(input$check_inputs, {
      setup_result(NULL)
      tryCatch(start(setup_request(FALSE)), error = notify_error)
    })
    shiny::observeEvent(input$start_analysis, {
      setup_result(NULL)
      tryCatch({
        if (isTRUE(start(setup_request(TRUE)))) {
          # A cancelled/failed rerun must not expose a stale in-memory case
          # while its result files are being replaced by the workflow.
          project(NULL); active(NULL); point(NULL); active_request(NULL)
        }
      }, error = notify_error)
    })
    output$setup_summary <- shiny::renderUI({
      x <- setup_result(); if (is.null(x)) return(NULL)
      shiny::tags$p(class = "cast-setup-receipt", if (isTRUE(x$run)) "Analysis completed." else
        "Input checks passed. No preparation, fitting or projection started.",
        shiny::tags$br(), x$case_dir)
    })
    shiny::observe({
      shiny::invalidateLater(750, session)
      if (is.null(job)) return()
      if (!is.null(logfile) && file.exists(logfile))
        log_text(paste(utils::tail(readLines(logfile, warn = FALSE, encoding = "UTF-8"), 60L), collapse = "\n"))
      if (job$is_alive()) return()
      finished <- job; job <<- NULL; job_busy(FALSE)
      tryCatch({
        value <- finished$get_result()
        if (identical(pending$operation, "setup")) {
          setup_result(value)
          if (isTRUE(value$run)) {
            shiny::updateTextInput(session, "case_dir", value = value$case_dir)
            load_project(directory = value$case_dir)
          } else {
            status_kind("ready"); status("Input checks passed; ready to run.")
          }
        } else {
          active(value); active_request(pending); point(NULL)
          if (!is.null(value$evidence)) shiny::updateSelectInput(session, "focus_driver",
            choices = stats::setNames(value$evidence$summary$driver,
              .cast_app_driver_label(project(), value$evidence$summary$driver)), selected = value$settings$focus)
          shiny::updateTabsetPanel(session, "result_tab", selected =
            if (!is.null(value$evidence) && !isTRUE(pending$maps)) "drivers" else "map")
          layers <- c("total_change", "current", "future", "support_total", "total_sd")
          labels <- c("\u603b\u9002\u5b9c\u5ea6\u53d8\u5316", "\u6539\u53d8\u524d", "\u6539\u53d8\u540e",
            "\u603b\u53d8\u5316\u8303\u56f4\u652f\u6301", "\u603b\u53d8\u5316\u6a21\u578b\u5206\u6b67")
          focus <- value$settings$focus
          if (!is.null(focus) && is.null(restore_state))
            shiny::updateSelectInput(session, "profile_driver", selected = focus)
          if (!is.null(focus) && !is.null(value$raster_file)) {
            layers <- c(layers, paste0(c("effect_", "retained_", "context_direction_",
              "context_gap_", "effect_agreement_", "retained_agreement_"), focus))
            labels <- c(labels, paste(.cast_app_driver_label(project(), focus), c("| \u6bd4\u8f83\u524d\u80cc\u666f\u54cd\u5e94", "| \u6bd4\u8f83\u540e\u80cc\u666f\u54cd\u5e94",
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
          status_kind("ready"); status("\u8ba1\u7b97\u5b8c\u6210")
        }
      }, error = notify_error)
    })
    session$onSessionEnded(function() {
      if (!is.null(job) && job$is_alive()) job$kill_tree()
    })
    shiny::observe({
      x <- active(); p <- project()
      figure <- input$figure
      figure_available <- !is.null(p) && length(figure) == 1L &&
        identical(basename(figure), figure) && file.exists(file.path(p$directory, "figures", figure))
      session$sendCustomMessage("cast_ui", list(busy = job_busy(), case_available = !is.null(p),
        map_available = !is.null(x$evidence), bundle_available = !is.null(x),
        points_available = !is.null(x$points_file), figure_available = figure_available,
        case_export_available = !is.null(p)))
    })
    output$status <- shiny::renderUI(shiny::tags$span(class = paste0("cast-status-", status_kind()), status()))
    output$case_label <- shiny::renderUI({
      p <- project()
      label <- if (is.null(p)) "\u672a\u6253\u5f00\u6848\u4f8b" else if (identical(p$case$species, "Schizothoracinae"))
        "\u88c2\u8179\u9c7c\u4e9a\u79d1" else gsub("_", " ", p$case$species)
      shiny::tags$span(class = "cast-case-label", label)
    })
    output$log <- shiny::renderUI(shiny::tags$pre(log_text()))
    output$map_caption <- shiny::renderUI({
      x <- active(); if (is.null(x)) return(shiny::tags$h2("\u9009\u62e9\u60c5\u666f\uff0c\u67e5\u770b\u73af\u5883\u54cd\u5e94"))
      vars <- x$settings$variables
      selection <- if (length(vars) == length(project()$fit$env_vars)) "\u5168\u90e8\u5165\u6a21\u56e0\u5b50" else
        paste(.cast_app_driver_label(project(), vars), collapse = "\u3001")
      shiny::tagList(shiny::tags$h2(.cast_app_scenario_label(x$settings$scenario)),
        shiny::tags$p(class = "cast-subtitle", paste(selection,
          if (is.null(x$settings$shifts)) "\u5b9e\u9645\u60c5\u666f" else "\u5047\u8bbe\u6539\u53d8", sep = " \u00b7 ")))
    })
    output$map_hint <- shiny::renderUI({
      x <- active()
      if (!is.null(x$preview)) return(NULL)
      shiny::div(class = "cast-map-hint", if (!is.null(x$evidence))
        "\u6837\u70b9\u6982\u89c8\u5df2\u5b8c\u6210\u3002\u5728\u73af\u5883\u54cd\u5e94\u9875\u751f\u6210\u56e0\u5b50\u5730\u56fe\u540e\uff0c\u5373\u53ef\u67e5\u770b\u7a7a\u95f4\u53d8\u5316\u3002" else
        "\u8bfb\u53d6\u5df2\u6709\u60c5\u666f\u5730\u56fe\u2026")
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
      session$sendCustomMessage("cast_point", !is.null(p))
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
    output$point_summary <- shiny::renderUI({
      p <- point()
      if (is.null(p)) return(shiny::tags$p(class = "cast-note", "\u70b9\u51fb\u5730\u56fe\u4e2d\u7684\u5206\u6790\u4f4d\u7f6e\u3002"))
      v <- .cast_app_point_contrast(p, input$layer, active()$settings$focus)
      labels <- if (v$kind == "gap") c("\u6bd4\u8f83\u524d\u54cd\u5e94", "\u6bd4\u8f83\u540e\u54cd\u5e94", "\u54cd\u5e94\u5dee") else
        c("\u6539\u53d8\u524d", "\u6539\u53d8\u540e", "\u53d8\u5316")
      shiny::tagList(shiny::tags$p(class = "cast-point-meta", paste(sprintf("%.5f, %.5f", p$longitude, p$latitude),
        if (length(v$selected)) paste("|", sprintf("%.4f", v$selected)) else "")),
        shiny::div(class = "cast-point-summary",
          .cast_app_metric(labels[1L], sprintf("%.3f", v$before)),
          .cast_app_metric(labels[2L], sprintf("%.3f", v$after)),
          .cast_app_metric(labels[3L], sprintf("%+.3f", v$change))),
        shiny::tags$p(class = "cast-note", if (isTRUE(v$supported)) "\u8be5\u4f4d\u7f6e\u901a\u8fc7\u8303\u56f4\u68c0\u67e5\u3002" else
          "\u8be5\u4f4d\u7f6e\u8d85\u51fa\u8bad\u7ec3\u8303\u56f4\uff0c\u4ec5\u4f9b\u67e5\u770b\u6a21\u578b\u8f93\u51fa\u3002"))
    })
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
      shiny::validate(shiny::need(!is.null(p$curve), p$note %||% "No environmental profile available."))
      old <- graphics::par(mar = c(3.8, 4, 1, 1), fg = "#74818A", col.axis = "#74818A", col.lab = "#74818A", cex = .85)
      on.exit(graphics::par(old))
      graphics::plot(p$curve$environment, p$curve$suitability, type = "l", col = "#206B61", lwd = 2,
        xlab = p$driver, ylab = "Relative suitability", bty = "l")
      value <- p$environment$before[p$environment$variable == p$driver]
      graphics::abline(v = value, lty = 3, col = "#A8B1B7")
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
      shiny::tags$p(class = "cast-note", paste("\u5171\u540c\u652f\u6301\u6837\u70b9", p$n_common, "/", p$n_input,
        "\u00b7 \u54cd\u5e94\u548c\u53cd\u53f7\u7387\u4f7f\u7528\u540c\u4e00\u6279\u4f4d\u7f6e\uff0c\u540c\u5411\u7387\u4ec5\u9488\u5bf9\u5b9e\u9645\u53d8\u5316\u6837\u70b9\u3002\u6bd4\u4f8b\u4e0d\u4ee3\u8868\u9762\u79ef\u3002"))
    })
    output$response_summary <- shiny::renderUI({
      x <- active()
      if (is.null(x$evidence)) return(shiny::div(class = "cast-empty", "\u70b9\u51fb\u5de6\u4fa7\u201c\u5206\u6790\u73af\u5883\u54cd\u5e94\u201d\uff0c\u4f7f\u7528\u5df2\u6709\u6a21\u578b\u751f\u6210\u6837\u70b9\u6982\u89c8\u3002"))
      p <- x$evidence$population; loc <- x$evidence$locations
      flip <- if (p$n_common) sum(loc$context_reversals > 0L, na.rm = TRUE) else NA_integer_
      changed <- if (p$n_common) sum(loc$same_largest == 0L, na.rm = TRUE) else NA_integer_
      # Exact ties and no-response cases stay unresolved, outside this fraction.
      resolved <- sum(!is.na(loc$same_largest))
      shiny::div(class = "cast-summary",
        .cast_app_metric("\u652f\u6301\u6837\u70b9", p$n_common, paste0("/ ", p$n_input)),
        .cast_app_metric("\u81f3\u5c11\u4e00\u4e2a\u56e0\u5b50\u53cd\u53f7", if (is.na(flip)) "\u2014" else flip, paste0("/ ", p$n_common)),
        .cast_app_metric("\u6700\u5927\u54cd\u5e94\u56e0\u5b50\u6539\u53d8", if (is.na(changed)) "\u2014" else changed, paste0("/ ", resolved),
          title = "\u4ec5\u6bd4\u8f83\u4e24\u79cd\u80cc\u666f\u4e0b\u5747\u80fd\u552f\u4e00\u5f52\u5c5e\u7684\u4f4d\u7f6e\uff1b\u5e76\u5217\u548c\u96f6\u54cd\u5e94\u4f4d\u7f6e\u4e0d\u8ba1\u5165\u6b64\u5206\u6bcd\u3002"))
    })
    output$drivers <- shiny::renderUI({
      shiny::req(active()$evidence)
      d <- driver_summary(); d <- d[!is.na(d$rank), , drop = FALSE]
      if (!nrow(d)) return(shiny::tags$p(class = "cast-note", "\u672a\u5f97\u5230\u53ef\u6392\u5e8f\u7684\u53d8\u5316\u56e0\u5b50\uff0c\u8bf7\u67e5\u770b\u672a\u6392\u5e8f\u56e0\u5b50\u548c\u8303\u56f4\u652f\u6301\u3002"))
      .cast_app_driver_table(project(), d)
    })
    output$drivers_inactive <- shiny::renderUI({
      shiny::req(active()$evidence)
      d <- driver_summary(); d <- d[is.na(d$rank), , drop = FALSE]
      if (!nrow(d)) return(shiny::tags$p(class = "cast-note", "\u5168\u90e8\u56e0\u5b50\u5df2\u53c2\u4e0e\u6392\u5e8f\u3002"))
      statuses <- c(unchanged_scenario = "\u672c\u6279\u6837\u70b9\u672a\u6539\u53d8", unchanged_on_common_domain = "\u5171\u540c\u4f4d\u7f6e\u672a\u6539\u53d8",
        no_common_support = "\u7f3a\u4e4f\u5171\u540c\u652f\u6301", zero_model_response = "\u6a21\u578b\u54cd\u5e94\u4e3a\u96f6",
        missing_predictor_pairs = "\u7f3a\u5c11\u73af\u5883\u5bf9", unchanged_known_pairs = "\u5df2\u77e5\u4f4d\u7f6e\u672a\u6539\u53d8")
      shiny::div(class = "cast-inactive", lapply(seq_len(nrow(d)), function(i)
        shiny::tags$span(title = paste(d$driver[i], statuses[d$status[i]]), .cast_app_driver_label(project(), d$driver[i]))))
    })
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
      old <- graphics::par(mar = c(4, 17, 1, 1), cex = .85, col.axis = "#74818A", col.lab = "#74818A"); on.exit(graphics::par(old))
      graphics::barplot(rbind(rev(d$mean_only), rev(d$mean_retained)), names.arg = rev(d$driver), beside = TRUE,
        horiz = TRUE, las = 1, cex.names = .7,
        border = NA, col = c("#B9C6CA", "#206B61"), legend.text = c("Before", "After"),
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
    output$download_case <- shiny::downloadHandler(filename = function() "castSDM_case.zip",
      content = function(file) {
        shiny::req(project()); directory <- project()$directory
        check_suggested("zip", "for the downloadable case")
        files <- list.files(directory, recursive = TRUE, full.names = FALSE)
        files <- files[!grepl("(^|/)(interactive|stages|figure_data)/|[.]lock$", files)]
        zip::zipr(file, files = file.path(directory, files), root = directory)
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
        state <- list(schema = "cast_app_state/1", package_version = as.character(utils::packageVersion("castSDM")), request = active_request(), layer = input$layer,
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
#' @param case_dir Optional completed single-species workflow directory.
#'   Omit it to start with the no-code setup form. It may be a
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
cast_app <- function(case_dir = NULL, workflow = NULL, run = TRUE,
                      launch.browser = interactive(), host = "127.0.0.1", port = NULL,
                      state = NULL) {
  for (pkg in c("shiny", "bslib", "leaflet", "callr", "terra", "rpart"))
    check_suggested(pkg, "for the environmental response application")
  if (is.null(case_dir)) case_dir <- file.path(getwd(), "castSDM_analysis", "output", "cast031_chelsa_ensmean", "Schizothoracinae")
  if (!is.character(case_dir) || length(case_dir) != 1L || is.na(case_dir) || !nzchar(case_dir))
    cli::cli_abort("Supply a case directory.")
  case_dir <- normalizePath(case_dir, winslash = "/", mustWork = FALSE)
  if (!is.null(workflow)) workflow <- normalizePath(workflow, winslash = "/", mustWork = TRUE)
  if (!is.null(state) && (!is.list(state) || !identical(state$schema, "cast_app_state/1") ||
      !is.list(state$request) || length(state$request$operation) != 1L ||
      !state$request$operation %in% c("view", "response") || length(state$layer) != 1L ||
      length(state$profile_driver) != 1L))
    cli::cli_abort("Supply a castSDM view/response state exported by the app.")
  app <- shiny::shinyApp(.cast_app_ui(case_dir, workflow), .cast_app_server(case_dir, workflow, state))
  if (!run) return(app)
  shiny::runApp(app, launch.browser = launch.browser, host = host, port = port)
}
