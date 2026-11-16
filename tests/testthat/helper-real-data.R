real_case <- local({
  fixture <- NULL
  function() {
    if (is.null(fixture)) fixture <<- readRDS(system.file("extdata", "schizothoracinae",
      "records.rds", package = "castSDM", mustWork = TRUE))
    fixture
  }
})

real_fit <- local({
  fit <- NULL
  function() {
    for (p in c("ranger", "maxnet", "mgcv")) skip_if_not_installed(p)
    if (is.null(fit)) {
      x <- real_case()
      fit <<- cast_fit(x$data[x$partition == "train", c("presence", x$variables)],
        models = x$fit_settings$models, rf_ntree = x$fit_settings$rf_ntree,
        maxent_classes = x$fit_settings$maxent_classes,
        gam_backend = x$fit_settings$gam_backend, seed = x$fit_settings$seed,
        num_threads = 1L, verbose = FALSE)
    }
    fit
  }
})

real_grids <- function() {
  skip_if_not_installed("terra")
  path <- system.file("extdata", "schizothoracinae", package = "castSDM", mustWork = TRUE)
  list(current = terra::rast(file.path(path, "current.tif")),
    future = terra::rast(file.path(path, "future.tif")))
}
