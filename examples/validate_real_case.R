# Real-data software verification; never sources 00_RUN_ALL.R or fits an SDM.
# Rscript examples/validate_real_case.R <completed-case> [check-output] [old-baseline.rds]
# Uses only saved models, actual scenario rasters, existing CV rows and a tiny
# crop of the actual grid. Model responses are not ecological causal truth.
args <- commandArgs(TRUE)
if (!length(args)) stop("Supply a completed real case directory; no fallback/generated dataset is used.")
library(castSDM)
directory <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
output <- if (length(args) >= 2L) args[2L] else tempfile("cast-real-check-")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
internal <- function(name) getFromNamespace(name, "castSDM")
project <- internal(".cast_explore_project")(directory)
original <- tools::md5sum(file.path(directory, c("case.rds", "training/result.rds", "run_claim.rds")))
scenario <- project$case$report_scenario
stacks <- internal(".cast_explore_stacks")(project, scenario)
info <- project$run$cv$fold_info
at <- unique(as.integer(round(seq(1L, nrow(info), length.out = min(256L, nrow(info))))))
xy <- as.matrix(info[at, c("lon", "lat")])
a <- terra::extract(stacks$current, xy)
b <- terra::extract(stacks$future, xy)
weights <- project$comparison$weights
clamp <- project$run$protocol$clamp
records <- character()
verified <- function(name, condition) {
  if (!isTRUE(condition)) stop("FAILED: ", name)
  records <<- c(records, name)
  cat("PASS: ", name, "\n", sep = "")
}
same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE))

full <- cast_attribution(project$fit, a, b, weights = weights, clamp = clamp,
  verbose = FALSE, resources = list(cores = 2L))
total <- cast_attribution(project$fit, a, b, weights = weights, clamp = clamp,
  drivers = character(), verbose = FALSE, resources = list(cores = 2L))
verified("total-only and full attribution agree on actual locations",
  same(total$data$total_change, full$data$total_change))
independent <- cast_compare(project$fit, a, setNames(list(b), scenario),
  groups = project$comparison$groups, weights = weights, clamp = clamp,
  domain = "scenario", verbose = FALSE, resources = list(cores = 2L))
verified("attribution agrees with the separate factorial comparison kernel",
  same(total$data$total_change, independent$results[[1L]]$data$joint_change))

diagnostic <- cast_attribution(project$fit, a, b, weights = weights, clamp = clamp,
  diagnostics = TRUE, verbose = FALSE, resources = list(cores = 2L))
verified("response diagnostics preserve every original prediction and support value",
  same(diagnostic$data[names(full$data)], full$data))
members <- lapply(names(weights), function(m) cast_attribution(project$fit, a, b,
  models = m, weights = setNames(1, m), clamp = clamp, verbose = FALSE)$data)
for (v in project$fit$env_vars) {
  for (prefix in c("effect_", "retained_")) {
    field <- paste0(prefix, v)
    predictions <- do.call(cbind, lapply(members, `[[`, field))
    expected <- rowMeans(sign(predictions) == sign(full$data[[field]]))
    verified(paste("model sign agreement from independent saved engines", field),
      same(diagnostic$data[[paste0(prefix, "agreement_", v)]], expected))
  }
}
evidence <- internal(".cast_response_evidence")(diagnostic$data, project$fit$env_vars, a, b)
bounds <- internal(".cast_support_bounds")(internal(".cast_reference")(project$fit, names(a)),
  project$comparison$support_probs)
eligible <- vapply(seq_len(nrow(a)), function(i) all(vapply(names(a), function(v)
  is.finite(a[[v]][i]) && is.finite(b[[v]][i]) &&
    a[[v]][i] >= bounds[[v]][1] && a[[v]][i] <= bounds[[v]][2] &&
    b[[v]][i] >= bounds[[v]][1] && b[[v]][i] <= bounds[[v]][2], logical(1))), logical(1))
verified("common ranking population agrees with actual training-reference ranges",
  identical(evidence$locations$common_support, as.integer(eligible)))
verified("all compared drivers have the same ranking denominator",
  all(evidence$summary$n_common == sum(eligible)))
verified("unchanged actual scenario predictors are explicitly excluded from ranking",
  all(is.na(evidence$summary$rank[evidence$summary$status == "unchanged_scenario"])))
for (v in names(a)) {
  row <- evidence$summary[evidence$summary$driver == v, ]
  verified(paste("common-domain signed response matches raw rows", v),
    same(row$mean_only, mean(full$data[[paste0("effect_", v)]][eligible])) &&
      same(row$mean_retained, mean(full$data[[paste0("retained_", v)]][eligible])))
}
write.csv(evidence$summary, file.path(output, "drivers_common.csv"), row.names = FALSE)
write.csv(evidence$population, file.path(output, "support_population.csv"), row.names = FALSE)

# Use genuinely missing cells from the real grid, never manufactured NAs.
edge_a <- as.data.frame(terra::values(stacks$current, row = 1L, nrows = 1L))
edge_b <- as.data.frame(terra::values(stacks$future, row = 1L, nrows = 1L))
missing <- which(rowSums(is.finite(as.matrix(edge_a))) == 0L &
  rowSums(is.finite(as.matrix(edge_b))) == 0L)
if (length(missing)) {
  edge_a <- edge_a[missing[1L], , drop = FALSE]; edge_b <- edge_b[missing[1L], , drop = FALSE]
  absent <- cast_attribution(project$fit, edge_a, edge_b, weights = weights, verbose = FALSE)
  absent <- internal(".cast_response_evidence")(absent$data, project$fit$env_vars, edge_a, edge_b)
  verified("actual NoData cells are labelled missing, never unchanged or ranked",
    all(absent$summary$status == "missing_predictor_pairs") && all(is.na(absent$summary$rank)))
}

# Reconstruct an actual projected future vector using its recorded differences.
# No fabricated covariate values or response labels enter this check.
i <- which(internal(".cast_finite_rows")(a) & internal(".cast_finite_rows")(b))[1L]
delta <- setNames(as.numeric(unlist(b[i, ] - a[i, ], use.names = FALSE)), names(a))
shifted <- cast_attribution(project$fit, a[i, , drop = FALSE], a[i, , drop = FALSE],
  shifts = delta, drivers = character(), weights = weights, clamp = clamp, verbose = FALSE)
verified("block shifts reproduce an actual current/future vector pair",
  same(shifted$data$total_change, total$data$total_change[i]))

groups <- cast_response_groups(a, total$data$total_change, info$fold[at],
  supported = total$data$support_total == 1)
keep <- !is.na(groups$group)
verified("group predictions reproduce fitted leaf means",
  same(predict(groups, a, type = "response")[keep], groups$fitted[keep]))
saveRDS(groups, file.path(output, "groups.rds"))
verified("response grouping survives serialisation",
  identical(predict(groups, a), predict(readRDS(file.path(output, "groups.rds")), a)))

# Exercise split traversal on the real responses even if spatial tuning chose
# a single group. This is a diagnostic tree, not a claimed ecological result.
x <- a[keep, , drop = FALSE]; names(x) <- paste0("v", seq_len(ncol(x)))
d <- x; d$.response <- total$data$total_change[keep]
tree <- rpart::rpart(.response ~ ., data = d, method = "anova",
  control = rpart::rpart.control(cp = 0, xval = 0L, maxdepth = 3L,
    maxsurrogate = 0L, usesurrogate = 0L))
nodes <- internal(".cast_response_nodes")(tree, x)
means <- tree$frame$yval[match(nodes, as.integer(row.names(tree$frame)))]
verified("numeric tree traversal matches rpart on actual rows",
  same(means, as.numeric(predict(tree, x))))
unchanged <- cast_attribution(project$fit, a, a, drivers = character(), weights = weights,
  clamp = clamp, verbose = FALSE)
root <- cast_response_groups(a, unchanged$data$total_change, info$fold[at],
  supported = unchanged$data$support_total == 1)
verified("unchanged real inputs yield one zero-response group",
  nrow(root$summary) == 1L && root$summary$mean_response == 0)

# Tiny actual raster crop, retaining its original geometry, values and NAs.
resolution <- terra::res(stacks$current)
e <- terra::ext(xy[i, 1L] - resolution[1L], xy[i, 1L] + resolution[1L],
  xy[i, 2L] - resolution[2L], xy[i, 2L] + resolution[2L])
ra <- terra::crop(stacks$current, e); rb <- terra::crop(stacks$future, e)
native <- cast_attribution(project$fit, ra, rb, drivers = character(), weights = weights,
  clamp = clamp, filename = file.path(output, "real_crop_response.tif"),
  verbose = FALSE, overwrite = TRUE, resources = list(cores = 2L))
table <- cast_attribution(project$fit,
  as.data.frame(terra::values(ra, mat = TRUE)), as.data.frame(terra::values(rb, mat = TRUE)),
  drivers = character(), weights = weights, clamp = clamp, verbose = FALSE)
verified("real native-grid crop and table responses agree",
  same(terra::values(terra::rast(native$raster_file)[["total_change"]], mat = FALSE), table$data$total_change))
focused <- cast_attribution(project$fit, ra, rb, drivers = names(a)[1L], weights = weights,
  clamp = clamp, filename = file.path(output, "real_crop_focused.tif"), diagnostics = TRUE,
  verbose = FALSE, overwrite = TRUE, resources = list(cores = 2L))
focused_table <- cast_attribution(project$fit,
  as.data.frame(terra::values(ra, mat = TRUE)), as.data.frame(terra::values(rb, mat = TRUE)),
  drivers = names(a)[1L], weights = weights, clamp = clamp, diagnostics = TRUE, verbose = FALSE)
verified("focused real raster and table agree for both backgrounds and all diagnostics",
  same(terra::values(terra::rast(focused$raster_file)), as.matrix(focused_table$data)))
focus_preview <- internal(".cast_explore_preview")(terra::rast(focused$raster_file),
  file.path(output, "real_crop_focused_preview.tif"), pixels = terra::ncell(ra))
discrete <- grep("context_direction_|agreement_", focused$layers, value = TRUE)
verified("real preview preserves direction classes and model-vote fractions",
  all(vapply(discrete, function(field) {
    a <- unique(na.omit(terra::values(terra::rast(focused$raster_file)[[field]], mat = FALSE)))
    b <- unique(na.omit(terra::values(terra::rast(focus_preview)[[field]], mat = FALSE)))
    length(a) > 0L && length(b) > 0L && all(vapply(b, function(x) any(abs(a - x) < 1e-6), logical(1)))
  }, logical(1))))
group_map <- internal(".cast_explore_group_map")(groups, ra,
  terra::rast(native$raster_file)[["support_total"]], file.path(output, "real_crop_groups.tif"))
supported <- terra::values(terra::rast(native$raster_file)[["support_total"]], mat = FALSE) == 1
membership <- predict(groups, as.data.frame(terra::values(ra, mat = TRUE)))
verified("group raster matches real-cell memberships on supported cells",
  same(terra::values(terra::rast(group_map), mat = FALSE)[supported %in% TRUE], membership[supported %in% TRUE]))
preview <- internal(".cast_explore_preview")(native$raster_file,
  file.path(output, "real_crop_preview.tif"), pixels = terra::ncell(ra))
verified("preview keeps support categorical",
  all(na.omit(terra::values(terra::rast(preview)[["support_total"]], mat = FALSE)) %in% c(0, 1)))
object_preview <- internal(".cast_explore_preview")(terra::rast(native$raster_file),
  file.path(output, "real_crop_object_preview.tif"), pixels = terra::ncell(ra))
verified("preview accepts raster objects without discarding their values",
  same(terra::values(terra::rast(preview)), terra::values(terra::rast(object_preview))))

explored <- cast_explore(directory, maps = FALSE, sample_n = length(at),
  output_dir = file.path(output, "exploration"), resources = list(cores = 2L))
verified("saved-case sample responses match independently extracted real inputs",
  same(read.csv(explored$points_file)$total_change, total$data$total_change))
stamp <- file.info(explored$points_file)$mtime
reused <- cast_explore(directory, maps = FALSE, sample_n = length(at),
  output_dir = file.path(output, "exploration"), resources = list(cores = 2L))
verified("content-verified cache reuses unchanged real outputs",
  identical(explored$key, reused$key) && identical(stamp, file.info(reused$points_file)$mtime))
replay <- new.env(parent = globalenv())
sys.source(file.path(explored$directory, "reproduce.R"), envir = replay)
verified("exported analysis script reuses the identical checked cache",
  identical(explored$key, replay$result$key) && identical(stamp, file.info(replay$result$points_file)$mtime))
verified("automatic map focus is the unique first-ranked changed driver",
  identical(explored$settings$focus,
    explored$evidence$summary$driver[explored$evidence$summary$rank %in% 1L]))
verified("real-case source-specific variable definitions are exported",
  file.exists(file.path(explored$directory, "predictor_metadata.csv")))

worker <- internal(".cast_app_start")(list(operation = "response", case_dir = directory,
  maps = FALSE, sample_n = length(at), output_dir = file.path(output, "exploration"),
  resources = list(cores = 2L)), file.path(output, "worker.log"))
worker$wait(timeout = 60000)
if (worker$is_alive()) { worker$kill(); stop("Background check timed out.") }
worker_result <- worker$get_result()
verified("isolated background worker uses the calling installation and same cache",
  identical(worker_result$key, explored$key))

app <- cast_app(directory, run = FALSE)
verified("Shiny object constructs without training or server startup", inherits(app, "shiny.appobj"))
rendered <- htmltools::renderTags(internal(".cast_app_ui")(directory))
verified("Shiny HTML contains all analysis controls", all(vapply(c("map", "analyse", "scenario", "profile_driver"),
  function(id) grepl(paste0('id="', id, '"'), rendered$html, fixed = TRUE), logical(1))))

# Query an actual location against the actual tiny native-grid result, using
# the same point/curve function as the app. No full-grid calculation is needed.
point_result <- explored
point_result$raster_file <- native$raster_file
ll <- terra::project(xy[i, , drop = FALSE], terra::crs(stacks$current), "EPSG:4326")
point <- internal(".cast_app_point")(project, point_result, ll[1L, 1L], ll[1L, 2L], names(a)[1L])
verified("map point query reads the correct native-grid response",
  same(point$values$total_change, total$data$total_change[i]))
verified("point environment values match actual raster inputs",
  same(point$environment$before, as.numeric(unlist(a[i, ], use.names = FALSE))) &&
    same(point$environment$after, as.numeric(unlist(b[i, ], use.names = FALSE))))

if (length(args) >= 3L) {
  baseline <- readRDS(args[3L])
  verified("0.28.3 legacy paired predictions are unchanged", same(full$data, baseline$baseline$data))
  verified("0.28.x frozen-model cache identity is retained",
    identical(internal(".cast_fit_identity")(project$fit), baseline$fit_identity))
  for (stage in names(baseline$stages)) verified(paste("unchanged training stage", stage),
    identical(internal(".cast_run_stage_identity")(stage, names(project$fit$models)), baseline$stages[[stage]]))
}
verified("original study model, case and run claim are unchanged",
  identical(original, tools::md5sum(names(original))))
write.csv(data.frame(check = records, status = "PASS"), file.path(output, "checks.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
cat(length(records), "real-data software checks passed; no SDM refitting or full-grid research run.\n")
