test_that("sparse native reads preserve actual cell values, order and duplicates", {
  grid <- real_grids()$current
  i <- which(is.finite(terra::values(grid[[1L]], mat = FALSE)))
  indices <- c(rev(i), i)
  expected <- castSDM:::.cast_raster_block(grid, 1L, terra::nrow(grid))[indices, , drop = FALSE]
  actual <- castSDM:::.cast_raster_cells(grid, 1L, terra::nrow(grid), indices)
  expect_equal(actual, expected, ignore_attr = TRUE)
  expect_identical(names(actual), names(grid))
})

test_that("display requests prepare only the selected products", {
  core <- castSDM:::.cast_figure_display_request("core")
  expect_false(core$range)
  expect_length(core$other_layers, 0L)
  expect_true(all(c("current", "joint", "joint_change") %in% core$report_layers))
  expect_identical(castSDM:::.cast_figure_display_request("10")$other_layers, "joint")
  expect_identical(core$driver_mode, "focus")
  expect_null(castSDM:::.cast_figure_display_request("all"))
  expect_true(castSDM:::.cast_figure_display_request("08")$range)
})

test_that("saved real response maps require matching data and numerical policy", {
  skip_if_not_installed("filelock")
  fit <- real_fit(); g <- real_grids()
  root <- tempfile("real-response-reuse-"); dir.create(file.path(root,"figures",".data"),recursive=TRUE)
  directory <- file.path(root,"figures",".data")
  attr <- cast_attribution(fit,g$current,g$future,drivers=fit$env_vars,
    models="rf",weights=c(rf=1),filename=file.path(directory,"driver_attribution.tif"),verbose=FALSE)
  focus <- attr$drivers[1L]
  display <- castSDM:::.cast_figure_driver_display(attr, file.path(directory,"focus.tif"),
    resolution=.02, drivers=focus, include_diagnostics=FALSE)
  expect_identical(names(terra::rast(display)), paste0("effect_", focus))
  full <- castSDM:::.cast_figure_driver_display(attr, file.path(directory,"full.tif"),
    resolution=.02, drivers=focus)
  expect_true(all(paste0(c("effect_", "retained_"), focus) %in% names(terra::rast(full))))
  preview <- castSDM:::.cast_explore_preview(attr$raster_file, file.path(directory,"preview.tif"))
  expect_true(all(c("total_change", "supported_total_change", "support_total") %in% names(terra::rast(preview))))
  expect_true(terra::ncell(terra::rast(preview)) <= 250000L)
  project <- list(directory=root,fit=fit,comparison=list(weights=c(rf=1),clamp=FALSE,
    support_probs=c(.01,.99),results=list(report=list(raster_file=attr$raster_file))))
  memo <- new.env(parent=emptyenv())
  contract <- list(fit=castSDM:::.cast_fit_identity(fit),scenario="report",weights=c(rf=1),
    clamp=FALSE,support_probs=c(.01,.99),
    current=castSDM:::.cast_data_identity(g$current,memo),future=castSDM:::.cast_data_identity(g$future,memo),
    comparison=unname(tools::md5sum(attr$raster_file)))
  receipt <- file.path(directory,"driver_catalog.rds")
  castSDM:::.cast_cache_write(list(contract=contract,attribution=attr),receipt,
    castSDM:::.cast_digest(contract),attr$raster_file)
  expect_equal(castSDM:::.cast_explore_saved(project,"report",g$current,g$future,memo)$summary,attr$summary)
  project$comparison$clamp <- TRUE
  expect_null(castSDM:::.cast_explore_saved(project,"report",g$current,g$future,memo))
  expect_true(all(attr$model_summary_common$n == unique(attr$common_summary$n)))
})

test_that("comparison producers hand real blocks, GAM terms and predictions to explanations", {
  fit <- real_fit(); g <- real_grids()
  weights <- stats::setNames(rep(1/length(fit$models), length(fit$models)), names(fit$models))
  root <- tempfile("real-comparison-handoff-")
  comparison <- cast_compare(fit, g$current, list(actual = g$future),
    groups = list(environment = fit$env_vars), weights = weights, clamp = TRUE,
    output_dir = root, verbose = FALSE)
  # Equality alone would miss a fallback that recalculates everything.
  prepared <- castSDM:::.cast_compare_reuse(fit, g$current, g$future, comparison, weights, TRUE)
  expect_false(is.null(prepared$blocks))
  expect_false(is.null(prepared$gam_current))
  expect_true(isTRUE(attr(prepared$gam_current, "ready")))
  expect_true(file.exists(prepared$prediction))
  different_clamp <- castSDM:::.cast_compare_reuse(fit, g$current, g$future, comparison, weights, FALSE)
  expect_identical(different_clamp$blocks$key, prepared$blocks$key)
  expect_null(different_clamp$prediction)
  expect_null(different_clamp$gam_current)
  first <- prepared$blocks$starts[1L]
  count <- min(prepared$blocks$rows, terra::nrow(g$current)-first+1L)
  block <- castSDM:::.cast_compare_current_block(prepared$blocks, first, count)
  expect_equal(nrow(block$values), length(block$indices))
  seed <- castSDM:::.cast_reuse_predictions(prepared, first, count, block$indices, weights)
  expect_equal(nrow(seed$current), nrow(block$values))
  reused <- castSDM:::.cast_attribution(fit, g$current, g$future,
    drivers = fit$env_vars, weights = weights, clamp = TRUE, prepared = prepared,
    filename = file.path(root, "reused.tif"), verbose = FALSE)
  direct <- cast_attribution(fit, g$current, g$future,
    drivers = fit$env_vars, weights = weights, clamp = TRUE,
    filename = file.path(root, "direct.tif"), verbose = FALSE)
  expect_equal(terra::values(terra::rast(reused$raster_file)),
    terra::values(terra::rast(direct$raster_file)), tolerance = 1e-12)
  expect_null(castSDM:::.cast_compare_reuse(fit, g$future, g$current, comparison, weights, TRUE))
  receipt <- sub("[.]tif$", ".rds", prepared$prediction)
  envelope <- readRDS(receipt)
  envelope$key <- paste0(envelope$key, "_different_run")
  saveRDS(envelope, receipt)
  # A self-consistent overwritten receipt is still not this comparison's key.
  stale <- castSDM:::.cast_compare_reuse(fit, g$current, g$future, comparison, weights, TRUE)
  expect_null(stale$prediction)
  expect_false(is.null(stale$blocks))
})

test_that("ordered parallel baseline blocks preserve all real native predictions", {
  fit <- real_fit(); g <- real_grids()
  weights <- stats::setNames(rep(1/length(fit$models), length(fit$models)), names(fit$models))
  root <- tempfile("real-parallel-baseline-"); dir.create(root)
  plan <- list(current = g$current, vars = fit$env_vars, raster = TRUE, mask = NULL)
  block_size <- terra::ncol(g$current)
  blocks <- castSDM:::.cast_compare_current_blocks(plan, block_size,
    file.path(root, "current"), list(castSDM:::.cast_data_identity(g$current), NULL), FALSE)
  cluster <- castSDM:::.cast_start_compare_cluster(2L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  run <- function(name, pool) {
    gam <- castSDM:::.cast_gam_current_cache(fit, weights, fit$env_vars, blocks,
      file.path(root, paste0(name, "_gam")), TRUE)
    castSDM:::.cast_compare_baseline(fit, plan, weights, TRUE, block_size,
      file.path(root, paste0(name, ".tif")), blocks, gam, pool)
  }
  serial <- run("serial", NULL); concurrent <- run("parallel", cluster)
  expect_equal(terra::values(terra::rast(serial)), terra::values(terra::rast(concurrent)), tolerance = 1e-12)
  expect_true(all(unlist(parallel::clusterCall(cluster, function()
    !exists(".cast_baseline_state", envir = .GlobalEnv, inherits = FALSE)))))
})
