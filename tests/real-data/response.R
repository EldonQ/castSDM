# Rscript tests/real-data/response.R <real_case> <QA_dir> [scenario-io_bundle]
# All reference rows are retained. Optional plots use a verified actual crop.
args <- commandArgs(TRUE)
if (!length(args) %in% 2:3) stop("Supply a real completed case, QA directory and optional verified crop bundle.")
library(castSDM)
case <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[2L], winslash = "/", mustWork = FALSE)
if (startsWith(paste0(tolower(qa), "/"), paste0(tolower(case), "/"))) stop("Keep QA outside the source case.")
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
options(castSDM.predict_threads = 2L)
run <- readRDS(file.path(case, "training/result.rds"))
comparison <- readRDS(file.path(case, "comparison/comparison.rds"))
fit <- if (is.null(run$fit_full)) run$fit else run$fit_full
reference <- fit$scaling$reference
request <- castSDM:::.cast_figure_response_request("core", comparison, fit$env_vars, "ale")
stopifnot(identical(request$types, "ale"), setequal(request$variables, unlist(comparison$groups, use.names = FALSE)))
all_request <- castSDM:::.cast_figure_response_request("all", comparison, fit$env_vars, "ale")
stopifnot(identical(all_request$variables, fit$env_vars), setequal(all_request$types, c("ale", "fixed")),
  "mess" %in% castSDM:::.cast_figure_products(c("core", "novelty")))
variable <- if ("bio08" %in% request$variables) "bio08" else request$variables[1L]
stopifnot(length(variable) == 1L, !is.na(variable))
fast_time <- system.time(curve <- cast_response_curves(fit, variables = variable, type = "ale"))[["elapsed"]]
stopifnot(attr(curve, "reference_n") == nrow(reference))
errors <- list()
for (model in names(fit$models)) {
  info <- fit$models[[model]]
  predictor <- iml::Predictor$new(model = info, data = reference,
    predict.function = function(model, newdata) castSDM:::predict_single_model(model, newdata, clamp = TRUE))
  elapsed <- system.time(oracle <- iml::FeatureEffect$new(predictor, feature = variable,
    method = "ale", grid.size = 50)$results)[["elapsed"]]
  x <- curve[curve$model == model, ]
  difference <- max(abs(x$prediction-oracle$.value))
  stopifnot(identical(x$value, oracle[[variable]]), is.finite(difference), difference < 1e-12)
  errors[[model]] <- data.frame(model = model, n_reference = nrow(reference),
    max_ALE_difference = difference, ordinary_elapsed_s = elapsed)
}
utils::write.csv(do.call(rbind, errors), file.path(qa, "full_reference_ALE.csv"), row.names = FALSE)
cat("Published ALE parity on every reference row passed; total accelerated time:", fast_time, "s.\n")
fixed <- cast_response_curves(fit, variables = variable, type = "fixed")
for (model in names(fit$models)) {
  x <- fixed[fixed$model == model, ]
  profile <- reference[rep(1L, nrow(x)), , drop = FALSE]
  for (v in names(profile)) profile[[v]] <- stats::median(reference[[v]])
  profile[[variable]] <- x$value
  oracle <- castSDM:::predict_single_model(fit$models[[model]], profile, clamp = TRUE)
  stopifnot(max(abs(x$prediction-oracle)) < 1e-12)
}
cat("Median-reference fixed profiles equal ordinary prediction.\n")
if (length(args) == 3L) {
  saved <- readRDS(args[3L])
  out <- file.path(qa, "core_figures")
  result <- cast_figures(saved$run, saved$comparison, attribution = saved$attribution,
    output_dir = out, basemap = "china", dpi = 100, verbose = TRUE)
  stopifnot(is.data.frame(result$status), nrow(result$status) == 16L,
    length(result$requested) == 8L, all(result$status$status == "written_not_visually_checked"),
    all(result$status$warnings == ""))
  path <- result$data_file; receipt <- readRDS(path)
  data <- castSDM:::.cast_cache_read(path, receipt$key)
  stopifnot(identical(names(data$responses), "ale"),
    setequal(unique(data$responses$ale$variable), request$variables),
    attr(data$responses$ale, "reference_n") == nrow(reference),
    !file.exists(file.path(out, "source_data/fixed_profiles.csv")),
    !any(grepl("^R_", result$status$figure)))
  again <- cast_figures(saved$run, saved$comparison, output_dir = out,
    basemap = "china", dpi = 100, render_only = TRUE, verbose = FALSE)
  at <- match(result$status$path, again$status$path)
  stopifnot(!anyNA(at), identical(again$status$md5[at], result$status$md5), all(again$status$cache_hit))
  # Adding a fixed-profile request retains existing products without rewriting them.
  files <- data$artifact_files
  before <- file.info(files)$mtime
  extra <- cast_figures(saved$run, saved$comparison, output_dir = out,
    include = "14", basemap = "china", dpi = 100, verbose = FALSE)
  receipt <- readRDS(extra$data_file)
  expanded <- castSDM:::.cast_cache_read(extra$data_file, receipt$key)
  stopifnot(identical(before, file.info(files)$mtime),
    all(c("maps", "responses", "attribution") %in% expanded$products),
    setequal(names(expanded$responses), c("ale", "fixed")),
    identical(expanded$attribution_supported_display, data$attribution_supported_display))
  # The complete inventory still exists for the author's saved complete export.
  p <- file.path(case, "figures/.data/figure_data.rds"); receipt <- readRDS(p)
  complete <- castSDM:::.cast_cache_read(p, receipt$key)
  inventory <- names(castSDM:::.cast_figure_registry(complete, "china"))
  original <- utils::read.csv(file.path(case, "figures/figure_inventory.csv"))$figure
  stopifnot(setequal(inventory, original))
  cat("Eight core figures, full-reference scenario curves, render-only cache and complete inventory passed.\n")
}
