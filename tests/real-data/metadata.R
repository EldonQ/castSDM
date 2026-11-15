# Rscript tests/real-data/metadata.R <converted_SchiSOFT_csv> <original_source_csv>
# Both arguments must be real records; no simulated measurements or dates.
args <- commandArgs(TRUE)
if (length(args) != 2L || any(!file.exists(args))) stop("Supply both real source CSV files.")
library(castSDM)
env <- new.env(parent = asNamespace("castSDM"))
sys.source(system.file("workflows/single_species/prepare_data.R", package = "castSDM"), env)
raw <- data.table::fread(args[1L], data.table = FALSE, colClasses = "character")
joined <- env$case_source_metadata(raw, list(occurrence_metadata_file = args[2L]))
stopifnot(identical(raw$year, joined$year), identical(raw$decimalLongitude, joined$decimalLongitude),
  identical(raw$decimalLatitude, joined$decimalLatitude), nrow(raw) == nrow(joined))
trusted <- joined$sampling_date_status %in% c("reported_sampling_date", "reported_interval_within_one_year")
stopifnot(all(!is.finite(joined$sampling_year_verified[!trusted])),
  all(is.finite(joined$sampling_year_verified[trusted])))
print(table(joined$sampling_date_status))
cat("Source-verified metadata preserve original coordinates/years; ambiguous intervals and substitutes do not create precise sampling years.\n")
