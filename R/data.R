#' Ovis ammon (Argali Sheep) Occurrence Data
#'
#' Species occurrence and environmental data for Argali sheep (*Ovis ammon*)
#' on an ISEA3H Resolution 9 hexagonal grid covering China. Environmental
#' variables have been VIF-screened (threshold = 10) to 11 variables.
#'
#' @format A data frame with 3648 rows and 20 columns:
#' \describe{
#'   \item{HID}{Hexagon grid cell ID (ISEA3H Res9).}
#'   \item{lon, lat}{Geographic coordinates (WGS84).}
#'   \item{species}{Scientific name.}
#'   \item{sid}{IUCN species ID.}
#'   \item{family}{Taxonomic family (Bovidae).}
#'   \item{category}{IUCN Red List category (NT = Near Threatened).}
#'   \item{presence}{Binary presence/absence (1/0).}
#'   \item{fraction}{Occupancy fraction within hexagon (0-1).}
#'   \item{bio02}{Diurnal temperature range (BIO2, WorldClim).}
#'   \item{bio15}{Precipitation seasonality (BIO15, WorldClim).}
#'   \item{bio19}{Precipitation of coldest quarter (BIO19, WorldClim).}
#'   \item{elevation}{Elevation in meters (SRTM).}
#'   \item{aridityindexthornthwaite}{Thornthwaite aridity index (ENVIREM).}
#'   \item{maxtempcoldest}{Maximum temperature of coldest month (ENVIREM).}
#'   \item{tri}{Terrain Ruggedness Index (SRTM).}
#'   \item{topowet}{Topographic Wetness Index (ENVIREM).}
#'   \item{nontree}{Non-tree vegetation cover fraction.}
#'   \item{etccdi_cwd}{Consecutive Wet Days index (ETCCDI).}
#'   \item{landcover_igbp}{IGBP land cover class.}
#' }
#'
#' @source EcoISEA3H project. IUCN Red List spatial data, WorldClim v2,
#'   ENVIREM, SRTM, MODIS land cover.
#'
#' @examples
#' data(ovis_ammon)
#' head(ovis_ammon)
#' table(ovis_ammon$presence)
"ovis_ammon"


#' China Environmental Grid (ISEA3H Res9)
#'
#' Full environmental data grid for spatial prediction across China on the
#' ISEA3H Resolution 9 hexagonal grid. Contains the same VIF-screened
#' environmental variables as [ovis_ammon].
#'
#' @format A data frame with 3648 rows and 14 columns:
#' \describe{
#'   \item{HID}{Hexagon grid cell ID.}
#'   \item{lon, lat}{Geographic coordinates (WGS84).}
#'   \item{bio02, bio15, bio19, elevation, aridityindexthornthwaite,
#'     maxtempcoldest, tri, topowet, nontree, etccdi_cwd,
#'     landcover_igbp}{Environmental variables (see [ovis_ammon]).}
#' }
#'
#' @source EcoISEA3H project.
#'
#' @examples
#' data(china_env_grid)
#' head(china_env_grid)
"china_env_grid"
