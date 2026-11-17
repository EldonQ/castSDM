# Control-plane settings only: no occurrence records or environmental values
# are generated or altered by these selection/failure regressions.
test_that("an unconfigured entry never silently chooses a species", {
  root <- tempfile("unconfigured-workflow-"); dir.create(root)
  old <- Sys.getenv("CAST_SPECIES",unset=NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("CAST_SPECIES") else Sys.setenv(CAST_SPECIES=old))
  Sys.unsetenv("CAST_SPECIES")
  expect_null(castSDM:::.cast_workflow_selection(root,TRUE)$species)
  expect_error(castSDM:::.cast_workflow_selection(root),"No species is configured")
  expect_error(cast_workflow(root,root,root),"taxon identifier")
})

test_that("saved species selection survives fresh sessions and rejects conflicting overrides", {
  root <- tempfile("selected-workflow-"); dir.create(root)
  old <- Sys.getenv("CAST_SPECIES",unset=NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("CAST_SPECIES") else Sys.setenv(CAST_SPECIES=old))
  Sys.unsetenv("CAST_SPECIES")
  saveRDS(list(species="Callosciurus_erythraeus"),file.path(root,"workflow_settings.rds"))
  selection <- castSDM:::.cast_workflow_selection(root)
  expect_identical(selection$species,"Callosciurus_erythraeus")
  expect_identical(selection$source,"saved project settings")
  Sys.setenv(CAST_SPECIES="Schizothoracinae")
  expect_error(castSDM:::.cast_workflow_selection(root),"conflicts with the saved species")
  Sys.setenv(CAST_SPECIES="Callosciurus_erythraeus")
  expect_identical(castSDM:::.cast_workflow_selection(root)$species,"Callosciurus_erythraeus")
})

test_that("fresh no-code setup has no fish-specific target or domain", {
  skip_if_not_installed("shiny")
  root <- tempfile("neutral-setup-"); dir.create(root)
  old_dir <- getwd(); old <- Sys.getenv("CAST_SPECIES",unset=NA_character_)
  on.exit({setwd(old_dir); if(is.na(old)) Sys.unsetenv("CAST_SPECIES") else Sys.setenv(CAST_SPECIES=old)})
  Sys.unsetenv("CAST_SPECIES"); setwd(root)
  defaults <- castSDM:::.cast_app_setup_defaults(NULL,NULL)
  expect_identical(defaults$species,"")
  expect_identical(defaults$domain,"national")
  expect_false(grepl("Schizothoracinae",defaults$path,fixed=TRUE))
})

test_that("real training-cache identity depends on values and spatial row order", {
  x <- real_case(); original <- x$data
  a <- castSDM:::.cast_digest(original)
  reversed <- original[rev(seq_len(nrow(original))),,drop=FALSE]
  expect_false(identical(a,castSDM:::.cast_digest(reversed)))
  future <- original; future[x$variables] <- x$future
  expect_false(identical(a,castSDM:::.cast_digest(future)))
  expect_identical(x$data,original)
})


test_that("changing a declared real CRS representation cannot reuse the previous CV receipt", {
  for (pkg in c("terra","sf","blockCV","ranger","pROC")) skip_if_not_installed(pkg)
  x <- real_case(); data <- x$data
  alternative <- terra::crs(real_grids()$current,proj=TRUE)
  expect_true(isTRUE(sf::st_crs(x$crs) == sf::st_crs(alternative)))
  expect_false(identical(x$crs,alternative))
  prepared <- cast_prepare(data,coords_crs=x$crs,seed=x$fit_settings$seed,verbose=FALSE)
  protocol <- cast_protocol(models="rf",coords_crs=x$crs,cv_k=5L,ensemble_method="equal",
    rf_ntree=x$fit_settings$rf_ntree,num_threads=1L,seed=x$fit_settings$seed)
  cache <- tempfile("actual-training-cache-")
  first <- cast_run(data,prepared=prepared,protocol=protocol,cache_dir=cache,verbose=FALSE)
  protocol$coords_crs <- alternative
  second <- cast_run(data,prepared=prepared,protocol=protocol,cache_dir=cache,verbose=FALSE)
  hit <- function(run,stage) {
    i <- which(vapply(run$receipts,`[[`,character(1),"stage") == stage)
    run$receipts[[i]]$cache_hit
  }
  expect_false(hit(second,"04_cv"))
  expect_true(hit(second,"03_fit"))
  expect_equal(second$cv$oof,first$cv$oof,tolerance=1e-12)
  expect_identical(second$fit$screen,second$screen)
  expect_identical(second$fit_full$screen,second$screen)
  expect_identical(data,x$data)
})
