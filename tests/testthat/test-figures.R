# AI-assisted restoration; scope and verification: inst/AI_ASSISTANCE.md.
test_that("underflowed zero p-values remain in screening plots", {
  skip_if_not_installed("ggplot2")
  x <- structure(list(method = "screen", selected = "zero",
    scores = data.frame(variable = c("zero", "positive", "unscored"),
      stage1_p = c(0, .2, NA_real_))), class = "cast_select")
  d <- plot(x, top = NULL)$data
  expect_setequal(d$variable, c("zero", "positive"))
  expect_equal(d$filter_score[d$variable == "zero"], -log10(.Machine$double.xmin))
  expect_identical(x$scores$stage1_p[1], 0)
})

test_that("sparse regular maps retain cell centres and declared spacing without shifting warnings", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("sf")
  d <- expand.grid(lon = c(100, 100.02, 100.06), lat = c(30, 30.02, 30.06))
  d$value <- seq_len(nrow(d))
  scale <- function(guide) ggplot2::scale_fill_gradient(na.value = NA)
  p <- .cast_map_surface(d, "value", scale, "Sparse grid", basemap = "none",
    raster = TRUE, cell_size = .02)
  built <- expect_warning(ggplot2::ggplot_build(p), NA)
  painted <- built$data[[1]][!is.na(built$data[[1]]$fill), ]
  expect_setequal(paste(painted$x, painted$y), paste(d$lon, d$lat))
  expect_equal(painted$xmax-painted$xmin, rep(.02, nrow(d)))
  expect_equal(painted$ymax-painted$ymin, rep(.02, nrow(d)))
  expect_warning(ggplot2::ggplotGrob(p), NA)
})

test_that("range maps classify on the native grid before display resampling", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 4, crs = "EPSG:4326")
  names(r) <- c("current", "joint", "support", "common_support")
  terra::values(r) <- cbind(c(.5, .2, .2, .6), c(.2, .2, .5, .7),
    c(1, 1, 0, 1), c(1, 0, 0, 1))
  source <- tempfile(fileext = ".tif"); target <- tempfile(fileext = ".tif")
  on.exit(unlink(c(source, target)), add = TRUE)
  terra::writeRaster(r, source, datatype = "FLT8S")
  .cast_figure_classes(list(raster_file = source), .5, target)
  values <- terra::values(terra::rast(target))
  expect_equal(unname(values[, 1]), c(-1, 0, 1, 2))
  expect_equal(unname(values[, 2]), c(-1, 0, NA, 2))
  expect_equal(unname(values[, 3]), c(-1, NA, NA, 2))
})

test_that("figure MESS evaluates raw inputs and retains the saved domain", {
  skip_if_not_installed("terra")
  reference <- data.frame(x = c(0, .3, .7, 1), z = c(0, 1, 0, 1))
  fit <- cast_as_fit(list(one = 1), reference, function(model, newdata) stop("No model prediction needed for MESS"))
  current <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2)
  names(current) <- names(reference)
  terra::values(current) <- cbind(c(0, .5, 1, 2), 0)
  future <- terra::rast(current); terra::values(future) <- cbind(c(1, 1.5, 2, 3), 0)
  domain <- current[[1]]; names(domain) <- "current"; terra::values(domain) <- c(NA, .5, .5, .5)
  source <- tempfile(fileext = ".tif"); target <- tempfile(fileext = ".tif")
  on.exit(unlink(c(source, target)), add = TRUE)
  terra::writeRaster(domain, source)
  .cast_figure_mess(fit, current, future, source, target, block_size = 2)
  values <- terra::values(terra::rast(target))
  expect_equal(unname(values[, "mess_current"]), c(NA, 0, 0, -100))
  expect_equal(unname(values[, "mess_future"]), c(NA, -50, -100, -200))
  expect_equal(unname(values[, "novel_future"]), c(NA, 1, 1, 1))
})

test_that("complete figure inventory includes the legacy set and render-only never predicts", {
  for (pkg in c("terra", "sf", "ggplot2", "patchwork", "ragg", "ranger", "pROC", "svglite"))
    skip_if_not_installed(pkg)
  set.seed(280)
  n <- 100
  d <- data.frame(lon = runif(n, 100, 102), lat = runif(n, 30, 32), x = rnorm(n), z = rnorm(n))
  d$presence <- rbinom(n, 1, plogis(d$x-.3*d$z))
  root <- tempfile("figure-export-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  p <- cast_protocol(models = "rf", select_method = "full", rf_ntree = 12,
    cv_k = 2, block_method = "grid", aoa = FALSE, ensemble_min_score = 0, seed = 3)
  run <- cast_run(d, protocol = p, output_dir = file.path(root, "training"), verbose = FALSE)
  r <- terra::rast(nrows = 2, ncols = 3, nlyrs = 2, xmin = 100, xmax = 101,
    ymin = 30, ymax = 31, crs = "EPSG:4326")
  names(r) <- c("x", "z"); terra::values(r) <- cbind(seq(-1, 1, length.out = 6), sin(1:6))
  f <- terra::rast(r); terra::values(f) <- terra::values(r)+.2
  cmp <- cast_compare(run, r, list(joint_ssp585_2071_2100 = f), list(first = "x", second = "z"),
    output_dir = file.path(root, "comparison"), verbose = FALSE)
  figures <- cast_figures(run, cmp, file.path(root, "figures"), current = r, future = f,
    response_grid = 3, display_resolution = .2, include = c("01", "02", "07"),
    basemap = "china", width_mm = 89, dpi = 72, verbose = FALSE)
  expect_true(all(sprintf("%02d", 1:10) %in% substr(figures$inventory$figure, 1, 2)))
  expect_equal(nrow(figures$status), 6L)
  expect_true(all(figures$status$status == "written_not_visually_checked"))
  expect_true(all(file.exists(file.path(figures$output_dir, figures$status$path))))
  tracker <- new.env(parent = emptyenv()); tracker$n <- 0L
  trace("predict.ranger", where = asNamespace("ranger"), tracer = function() {
    tracker$n <- tracker$n + 1L
    stop("Render-only must not call an engine")
  }, print = FALSE)
  on.exit(untrace("predict.ranger", where = asNamespace("ranger")), add = TRUE)
  withr::local_options(castSDM.figure_pdf_device = "pdf")
  again <- cast_figures(run, cmp, figures$output_dir, response_grid = 3,
    display_resolution = .2, include = "07", render_only = TRUE,
    basemap = "china", width_mm = 89, dpi = 72, verbose = FALSE,
    formats = c("png", "pdf", "svg"))
  expect_equal(tracker$n, 0L)
  expect_identical(again$pdf_device, "pdf")
  expect_identical(again$font_family, "Helvetica")
  expect_equal(nrow(again$status), 7L)
  expect_identical(again$requested, "07_model_performance")
  expect_true(all(figures$status$path %in% again$status$path))
  expect_true(any(endsWith(again$status$path, ".svg")))
  expect_true(all(again$status$status == "written_not_visually_checked"))
})
