# Rscript tests/real-data/interface.R <completed_real_case> <separate_QA_dir>
# Real Shiny server + installed background jobs + exported numeric replay.
# No workflow fitting or national scenario prediction is submitted.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
library(castSDM)
case <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2L], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(tolower(qa), "/"), paste0(tolower(case), "/"))) stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
shadow <- tempfile("interface_case_", tmpdir = qa)
dir.create(shadow)
relative <- c("case.rds", "completion.rds", "training_pool.csv", "predictor_metadata.csv",
  "scenario_table.csv", "training/result.rds", "comparison/comparison.rds")
rasters <- list.files(file.path(case, "comparison/rasters"), pattern = "[.]tif$", full.names = FALSE)
relative <- c(relative, file.path("comparison/rasters", rasters))
for (file in relative) {
  target <- file.path(shadow, file)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  # Existing artifacts are immutable inputs. Links avoid duplicating large
  # result rasters; all new interactive products go into this QA shadow case.
  if (!file.link(file.path(case, file), target)) stopifnot(file.copy(file.path(case, file), target))
}
info <- readRDS(file.path(shadow, "case.rds"))
scenario <- info$report_scenario
server <- castSDM:::.cast_app_server(shadow, NULL)
shiny::testServer(server, {
  # Initialise observers before their first action-button event; the mock
  # session does not apply browser-side updateInput messages automatically.
  session$flushReact()
  session$setInputs(case_dir = shadow, scenario = scenario, layer = "total_change",
    supported_only = FALSE, profile_driver = "bio08", context = "current")
  session$flushReact()
  wait_job <- function() {
    deadline <- Sys.time()+240
    while (isTRUE(job_busy())) {
      if (Sys.time() > deadline) stop("Background interface job timed out: ", log_text())
      Sys.sleep(.2); session$elapse(1000); session$flushReact()
    }
    if (identical(status_kind(), "error")) stop(log_text())
    stopifnot(!is.null(active()))
  }
  if (is.null(active()) && !isTRUE(job_busy())) session$setInputs(load_case = 1L)
  cat("View dispatch:", !is.null(project()), job_busy(), status_kind(), "\n")
  wait_job()
  stopifnot(identical(active()$schema, "cast_app_view/1"), file.exists(active()$preview))
  session$setInputs(mode = "scenario", variables = ".all", maxdepth = 3L,
    maps = FALSE, grouping = FALSE, cores = 2L)
  session$setInputs(analyse = 1L)
  session$flushReact(); wait_job()
  value <- active()
  stopifnot(identical(value$schema, "cast_explore/1"), is.null(value$settings$sample_n),
    value$evidence$population$n_input == nrow(project()$run$cv$fold_info),
    is.null(value$raster_file), file.exists(value$points_file))
  xy <- project()$run$cv$fold_info[1L, c("lon", "lat")]
  source <- terra::rast(project()$pool$path[1L])
  point_sf <- sf::st_as_sf(xy, coords = c("lon", "lat"), crs = terra::crs(source))
  click <- sf::st_coordinates(sf::st_transform(point_sf, 4326))
  session$setInputs(map_click = list(lng = click[1L, 1L], lat = click[1L, 2L]))
  session$flushReact()
  stopifnot(!is.null(point()), is.finite(point()$values$total_change))
  state <- list(schema = "cast_app_state/1", package_version = as.character(packageVersion("castSDM")),
    request = active_request(), layer = "total_change", supported_only = FALSE,
    profile_driver = "bio08", click = list(lng = click[1L, 1L], lat = click[1L, 2L]))
  zip <- file.path(qa, "interface_export.zip")
  castSDM:::.cast_app_bundle(value$directory, state, zip)
  bundle <- tempfile("replay_", tmpdir = qa); dir.create(bundle)
  utils::unzip(zip, exdir = bundle)
  stopifnot(file.exists(file.path(bundle, "reproduce.R")),
    file.exists(file.path(bundle, "reproduce_view.R")))
  before <- tools::md5sum(value$points_file)
  replay <- new.env(parent = globalenv())
  sys.source(file.path(bundle, "reproduce.R"), envir = replay)
  stopifnot(identical(replay$result$key, value$key), identical(tools::md5sum(value$points_file), before))
  restored <- readRDS(file.path(bundle, "view_state.rds"))
  stopifnot(inherits(cast_app(shadow, state = restored, run = FALSE), "shiny.appobj"))
  saveRDS(list(case = shadow, state = restored, response_key = value$key,
    n_locations = value$evidence$population$n_input, point = point(),
    export = zip, replay_unchanged = TRUE), file.path(qa, "interface_receipt.rds"))
})
cat("Actual saved-case view, all-location response, point query, export and numeric replay passed.\n")
