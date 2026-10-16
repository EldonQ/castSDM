# 09_convert_schisoft.R ---------------------------------------------------------
# Convert the SchiSOFT literature-compilation occurrence dataset
# (F:/LLQ/24638538/SchiSOFT_occurrence_records.csv, 7333 records across all
# 12 Schizothoracinae genera; names validated via rfishbase/Eschmeyer and
# coordinates cleaned via CoordinateCleaner by the dataset authors) into the
# GBIF-compatible column schema that read_occurrences() input_type = "gbif"
# consumes.
#
# CASE DECISION (user, 2026-10-01): model the SUBFAMILY. All 7333 records are
# kept under the case taxon key "Schizothoracinae"; the output `genus` column
# carries that case key uniformly so the wrong_taxon rule
# (genus == cfg$species_name) passes every record. The TRUE identification
# stays in `species` (= scientificName, e.g. "Schizothorax davidi"), so the
# QC audit and ODMAP report preserve the underlying taxonomy. An earlier
# genus-level conversion (Schizopygopsis only, 1451 records ->
# Schizopygopsis_schisoft.csv, 586 cells) remains on disk for comparison runs.
#
# Field mapping and QC semantics under read_occurrences():
#   gbifID                    <- occurrenceID (integer, used for dedup order)
#   decimalLongitude/Latitude <- longitudeX / latitudeY
#   genus                     <- "Schizothoracinae" (case taxon key, uniform)
#   species                   <- scientificName (true name, kept in audit)
#   countryCode               <- "CN" when the WGS84 point falls inside the
#                                China bounding box (73-135E, 18-54N), else NA
#                                -> records outside are dropped by the
#                                outside_country rule; the mask_china extract
#                                (outside_mask) is the authoritative filter.
#   occurrenceStatus          <- "PRESENT"      (not_present check)
#   hasGeospatialIssues       <- FALSE          (dataset was CoordinateCleaner'd)
#   coordinateUncertaintyInMeters <- NA (literature coordinates carry no
#                                uncertainty) -> kept via
#                                occurrence_keep_unknown_uncertainty = TRUE
#   basisOfRecord             <- "LITERATURE" (sourceType Book/Journal/...)
#   year                      <- first 4 digits of sampledDate (YYYYMMDD)
#   locality                  <- referenceInEnglish (provenance trail)

src <- "F:/LLQ/24638538/SchiSOFT_occurrence_records.csv"
out_csv <- "E:/Package/castSDM_workflows/species_occ/Schizothoracinae_schisoft.csv"
case_taxon <- "Schizothoracinae"

d <- data.table::fread(src, data.table = FALSE, encoding = "UTF-8")
names(d) <- sub("^\ufeff", "", names(d))
message("SchiSOFT records read: ", nrow(d), " across ",
        length(unique(d$genus)), " genera")
print(table(d$genus))

lon <- suppressWarnings(as.numeric(d$longitudeX))
lat <- suppressWarnings(as.numeric(d$latitudeY))
valid <- is.finite(lon) & is.finite(lat) & abs(lon) <= 180 & abs(lat) < 90
if (any(!valid)) message("Dropping ", sum(!valid), " record(s) with non-finite coordinates.")
d <- d[valid, , drop = FALSE]
lon <- lon[valid]
lat <- lat[valid]

in_cn <- is.finite(lon) & is.finite(lat) & lon >= 73 & lon <= 135 & lat >= 18 & lat <= 54
message("Inside China bbox: ", sum(in_cn), " | outside: ", sum(!in_cn))

year <- suppressWarnings(as.integer(substr(as.character(d$sampledDate), 1L, 4L)))
year[!is.finite(year) | year < 1700 | year > 2100] <- NA_integer_

out <- data.frame(
  gbifID = as.character(d$occurrenceID),
  datasetKey = "SchiSOFT_24638538",
  species = d$scientificName,
  genus = case_taxon,
  taxonRank = "SUBFAMILY",
  basisOfRecord = "LITERATURE",
  decimalLongitude = lon,
  decimalLatitude = lat,
  year = year,
  countryCode = ifelse(in_cn, "CN", NA_character_),
  occurrenceStatus = "PRESENT",
  hasGeospatialIssues = "FALSE",
  coordinateUncertaintyInMeters = NA_real_,
  locality = d$referenceInEnglish,
  stringsAsFactors = FALSE)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out, out_csv)
cat("wrote ", nrow(out), " ", case_taxon, " records -> ", out_csv, "\n", sep = "")
cat("coordinate range: lon ", paste(round(range(lon), 2), collapse = ".."),
    " lat ", paste(round(range(lat), 2), collapse = ".."), "\n", sep = "")
