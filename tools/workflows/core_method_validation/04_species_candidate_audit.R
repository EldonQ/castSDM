# Candidate audit only; no ecological validation is inferred from record counts.
.libPaths(c("E:/Package/.rlib_cast", .libPaths()))
source("E:/Package/castSDM_workflows/R/nsdm_castsdm_adapter.R", encoding = "UTF-8")
input <- "F:/LLQ/24638538/SchiSOFT_occurrence_records.csv"
output <- "E:/Package/castSDM_workflows/core_method_validation/species_candidate_audit"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
raw <- data.table::fread(input, data.table = FALSE, encoding = "UTF-8")
targets <- c("Schizopygopsis pylzovi", "Schizopygopsis thermalis",
             "Schizothorax prenanti", "Schizothorax nukiangensis", "Ptychobarbus kaznakovi")
d <- raw[raw$scientificName %in% targets, , drop = FALSE]
valid <- is.finite(d$longitudeX) & is.finite(d$latitudeY) &
  abs(d$longitudeX) <= 180 & abs(d$latitudeY) < 90
if (!all(valid)) stop("Candidate input has invalid coordinates; inspect before continuing.")
mask <- terra::rast("E:/CHECO26/N-SDM/data/mask_china.tif")
xy <- nsdm_project_occurrences(d[, c("longitudeX", "latitudeY")], terra::crs(mask, proj = TRUE))
d$x_albers <- xy[, 1]; d$y_albers <- xy[, 2]
d$cell <- terra::cellFromXY(mask, xy)
inside <- terra::extract(mask, xy)[[1]]
d$inside_china_mask <- is.finite(inside) & inside > 0
d$record_year <- suppressWarnings(as.integer(substr(as.character(d$sampledDate), 1, 4)))
d$date_is_publication_proxy <- !is.na(d$remarksOfSampledDate) &
  grepl("Received date|Published date", d$remarksOfSampledDate, ignore.case = TRUE)
d$reference_missing <- is.na(d$referenceInEnglish) | d$referenceInEnglish %in% c("", "NA")
summary <- do.call(rbind, lapply(targets, function(s) {
  z <- d[d$scientificName == s, , drop = FALSE]
  kept <- z$inside_china_mask
  before <- unique(z$cell[kept & is.finite(z$record_year) &
    z$record_year <= 2017 & !z$date_is_publication_proxy])
  after <- unique(z$cell[kept & is.finite(z$record_year) &
    z$record_year > 2017 & !z$date_is_publication_proxy])
  data.frame(species = s, records = nrow(z),
    unique_coordinates = nrow(unique(z[, c("longitudeX", "latitudeY")])),
    records_inside_mask = sum(kept), unique_1km_cells = length(unique(z$cell[kept])),
    identifiable_references = length(unique(z$referenceInEnglish[!z$reference_missing])),
    records_missing_reference = sum(z$reference_missing),
    publication_date_proxies = sum(z$date_is_publication_proxy),
    true_dated_pre2018_cells = length(before), true_dated_post2017_cells = length(after),
    post2017_cells_not_in_pre2018 = length(setdiff(after, before)),
    unknown_year = sum(!is.finite(z$record_year)))
}))
utils::write.csv(summary, file.path(output, "candidate_summary.csv"), row.names = FALSE)
# Preserve original identifiers and references locally; do not redistribute records.
data.table::fwrite(d, file.path(output, "candidate_records.csv"))
utils::write.csv(data.frame(path = input, md5 = unname(tools::md5sum(input)),
  role = "original SchiSOFT records with date provenance"),
  file.path(output, "input_manifest.csv"), row.names = FALSE)
utils::capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
print(summary)
writeLines(c("Candidate audit completed; ecological validation has NOT been run.",
  "Provisional lead: Schizopygopsis pylzovi, based on species names and traceable references.",
  "Counts are records/cells, not independent sampling units.",
  "Publication/received dates cannot be treated as sampling years.",
  "Accessible aquatic domain, observation bias and independent evaluation remain to be defined."),
  file.path(output, "STATUS.txt"))
