# Real occurrence/environment comparison, without national map projections.
# Rscript species-comparison.R <occurrence.csv> <environment_library> <new_QA_workspace>
# No observations, labels, future values or test fixtures are manufactured.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Supply occurrence CSV, environmental library and a separate new QA workspace.")
library(castSDM)
root <- normalizePath(args[3L], winslash="/", mustWork=FALSE)
if (dir.exists(root) && length(list.files(root, all.files=TRUE, no..=TRUE)))
  stop("Use an empty QA workspace; existing research results are never overwritten.")
cast_workflow(root, args[2L], args[1L], species="Callosciurus_erythraeus",
  settings=list(input_type="nsdm_occ", domain="national", lulc_product="glu"), run=FALSE)
e <- new.env(parent=globalenv())
for (f in c("runtime.R","data_adapter.R","inventory.R","prepare_data.R","01_config.R","scenario_case.R","02_pipeline.R"))
  sys.source(file.path(root,f), envir=e)
cfg <- e$castsdm_species_configs(root)$Callosciurus_erythraeus
preparation <- e$prepare_species(cfg)
d <- data.table::fread(file.path(preparation,paste0(cfg$species_name,"_model_data.csv")), data.table=FALSE)
p <- cfg$protocol
split <- cast_prepare(d,train_fraction=p$train_fraction,n_blocks=p$holdout_blocks,
  block_method=p$block_method,seed=p$seed,coords_crs=p$coords_crs,block_size=p$block_size)
meta <- e$case_pool_metadata(cfg,split$env_vars)
year <- e$case_training_year(preparation,split)
matched <- e$case_match_training_layers(d,split,meta,year$year,cfg$select_temporal_match)
meta <- matched$meta
response <- e$case_response_validation(cfg,matched$split,meta)
result <- cast_run(matched$data,prepared=matched$split,protocol=p,
  select_var_id=stats::setNames(meta$construct_id,meta$layer),
  select_category=stats::setNames(meta$group,meta$layer),select_roles=stats::setNames(meta$role,meta$layer),
  response_data=response,output_dir=file.path(root,"point_comparison"),cache_dir=cfg$cache_root)
cast_diagnose(result,output_dir=file.path(root,"diagnostics"))
cat("Real-data comparison complete; no national future maps were calculated.
")
