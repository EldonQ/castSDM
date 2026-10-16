# Direct interoperability check; equivalence, not a superiority benchmark.
.libPaths(c('E:/Package/.rlib_cast', .libPaths()))
pkgload::load_all('E:/Package/cast', quiet = TRUE)
root <- 'E:/Package/castSDM_workflows/core_method_validation'
output <- file.path(root, 'external_interoperability')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace('iml', quietly = TRUE), requireNamespace('blockCV', quietly = TRUE))
set.seed(1606)
reference <- data.frame(temperature = runif(200, -2, 2), moisture = runif(200, -1, 1))
y <- rbinom(200, 1, plogis(reference$temperature + reference$temperature * reference$moisture))
model <- stats::glm(y ~ temperature * moisture, reference, family = stats::binomial())
prediction <- function(model, newdata) as.numeric(stats::predict(model, newdata, type = 'response'))
explainer <- iml::Predictor$new(model, data = reference, y = y, predict.function = prediction)
adapter <- cast_as_fit(list(iml = explainer), reference,
  function(model, newdata) as.numeric(model$predict(newdata)[[1]]))
a <- reference[1:30, ]; b <- a
b$temperature <- b$temperature + seq(.05, .4, length.out = nrow(b))
b$moisture <- b$moisture + .1
z <- cast_attribution(adapter, a, b, verbose = FALSE)
only <- a; only$temperature <- b$temperature
held <- b; held$temperature <- a$temperature
oracle <- data.frame(total = prediction(model,b)-prediction(model,a),
  only = prediction(model,only)-prediction(model,a),
  retained = prediction(model,b)-prediction(model,held))
err <- c(total=max(abs(z$data$total_change-oracle$total)),
  only=max(abs(z$data$effect_temperature-oracle$only)),
  retained=max(abs(z$data$retained_temperature-oracle$retained)))
stopifnot(all(err < 1e-12))
utils::write.csv(data.frame(contrast=names(err),max_abs_error=err,tolerance=1e-12),
  file.path(output,'iml_prediction_equivalence.csv'),row.names=FALSE)
# Real blockCV design with a physical block width; no claim it is optimal for the fish case.
coords <- expand.grid(x = seq(1000, 19000, by = 2000), y = seq(1000, 19000, by = 2000))
coords$presence <- rep(c(0,1), length.out = nrow(coords))
pts <- sf::st_as_sf(coords, coords = c('x','y'), crs = 32648)
blocks <- blockCV::cv_spatial(x = pts, column = 'presence', size = 5000, k = 3,
  selection = 'random', iteration = 50, seed = 1606, plot = FALSE, progress = FALSE)
data <- data.frame(lon = coords$x, lat = coords$y, presence = coords$presence,
  temperature = seq(-1,1,length.out=nrow(coords)), moisture = rep(c(-1,0,1,0),25))
cv <- cast_cv(data, models='rf', select_method='full', k=3,
  fold_ids=blocks$folds_ids, rf_ntree=30, seed=1606, verbose=FALSE)
stopifnot(identical(cv$folds, as.integer(blocks$folds_ids)),cv$block_method=='external',
  all(cv$fold_status=='evaluated'))
utils::write.csv(data.frame(row=seq_len(nrow(data)),external_fold=blocks$folds_ids,
  cast_fold=cv$folds),file.path(output,'blockCV_fold_equivalence.csv'),row.names=FALSE)
utils::write.csv(cv$fold_metrics,file.path(output,'external_fold_metrics.csv'),row.names=FALSE)
utils::write.csv(data.frame(package=c('castSDM','iml','blockCV','sf'),
  version=vapply(c('castSDM','iml','blockCV','sf'),function(p) as.character(utils::packageVersion(p)),character(1))),
  file.path(output,'package_versions.csv'),row.names=FALSE)
utils::capture.output(sessionInfo(),file=file.path(output,'sessionInfo.txt'))
writeLines(c('Direct iml Predictor and blockCV interoperability passed.',
  'Paired contrasts reproduce the direct GLM predictions; exact external fold assignments preserved.',
  'This is implementation equivalence, not an advantage over iml or an optimal block-scale assessment.',
  'Callbacks inherit the external model prediction scale and preprocessing contract.'),
  file.path(output,'RESULTS.md'))
