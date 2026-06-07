#' Ovis ammon (Argali Sheep) Occurrence Data
#'
#' Species occurrence and environmental data for Argali sheep (*Ovis ammon*)
#' on an ISEA3H Resolution 9 hexagonal grid covering China. Contains 50
#' environmental variables (19 BioClim, ENVIREM, topographic, land cover,
#' ETCCDI climate extremes, and IGBP land cover fractions).
#'
#' @format A data frame with 3648 rows and 59 columns:
#' \describe{
#'   \item{HID}{Hexagon grid cell ID (ISEA3H Res9).}
#'   \item{lon}{Longitude (WGS84).}
#'   \item{lat}{Latitude (WGS84).}
#'   \item{species}{Scientific name.}
#'   \item{sid}{IUCN species ID.}
#'   \item{family}{Taxonomic family.}
#'   \item{category}{IUCN Red List category.}
#'   \item{presence}{Binary presence/absence (1/0).}
#'   \item{fraction}{Occupancy fraction within hexagon (0-1).}
#'   \item{bio01}{Annual Mean Temperature (WorldClim).}
#'   \item{bio02}{Mean Diurnal Range (WorldClim).}
#'   \item{bio03}{Isothermality (WorldClim).}
#'   \item{bio04}{Temperature Seasonality (WorldClim).}
#'   \item{bio05}{Max Temperature of Warmest Month (WorldClim).}
#'   \item{bio06}{Min Temperature of Coldest Month (WorldClim).}
#'   \item{bio07}{Temperature Annual Range (WorldClim).}
#'   \item{bio12}{Annual Precipitation (WorldClim).}
#'   \item{bio13}{Precipitation of Wettest Month (WorldClim).}
#'   \item{bio14}{Precipitation of Driest Month (WorldClim).}
#'   \item{bio15}{Precipitation Seasonality (WorldClim).}
#'   \item{bio16}{Precipitation of Wettest Quarter (WorldClim).}
#'   \item{bio17}{Precipitation of Driest Quarter (WorldClim).}
#'   \item{bio18}{Precipitation of Warmest Quarter (WorldClim).}
#'   \item{bio19}{Precipitation of Coldest Quarter (WorldClim).}
#'   \item{elevation}{Elevation in meters (SRTM).}
#'   \item{annualpet}{Annual Potential Evapotranspiration (ENVIREM).}
#'   \item{aridityindexthornthwaite}{Thornthwaite aridity index (ENVIREM).}
#'   \item{climaticmoistureindex}{Climatic Moisture Index (ENVIREM).}
#'   \item{embergerq}{Emberger's pluviothermic quotient (ENVIREM).}
#'   \item{maxtempcoldest}{Max temperature of coldest month (ENVIREM).}
#'   \item{petseasonality}{PET Seasonality (ENVIREM).}
#'   \item{tri}{Terrain Ruggedness Index (SRTM).}
#'   \item{topowet}{Topographic Wetness Index (SRTM).}
#'   \item{nontree}{Non-tree vegetation cover fraction.}
#'   \item{bare}{Bare ground cover fraction.}
#'   \item{tree}{Tree cover fraction.}
#'   \item{etccdi_cdd}{Consecutive Dry Days (ETCCDI).}
#'   \item{etccdi_cwd}{Consecutive Wet Days (ETCCDI).}
#'   \item{etccdi_fd}{Frost Days (ETCCDI).}
#'   \item{etccdi_gsl}{Growing Season Length (ETCCDI).}
#'   \item{etccdi_rx1day}{Max 1-day Precipitation (ETCCDI).}
#'   \item{etccdi_rx5day}{Max 5-day Precipitation (ETCCDI).}
#'   \item{etccdi_sdii}{Simple Daily Intensity Index (ETCCDI).}
#'   \item{igbp_01}{IGBP class 01 fraction.}
#'   \item{igbp_02}{IGBP class 02 fraction.}
#'   \item{igbp_03}{IGBP class 03 fraction.}
#'   \item{igbp_04}{IGBP class 04 fraction.}
#'   \item{igbp_05}{IGBP class 05 fraction.}
#'   \item{igbp_06}{IGBP class 06 fraction.}
#'   \item{igbp_07}{IGBP class 07 fraction.}
#'   \item{igbp_08}{IGBP class 08 fraction.}
#'   \item{igbp_09}{IGBP class 09 fraction.}
#'   \item{igbp_10}{IGBP class 10 fraction.}
#'   \item{igbp_11}{IGBP class 11 fraction.}
#'   \item{igbp_12}{IGBP class 12 fraction.}
#'   \item{igbp_13}{IGBP class 13 fraction.}
#'   \item{igbp_14}{IGBP class 14 fraction.}
#'   \item{igbp_15}{IGBP class 15 fraction.}
#'   \item{igbp_16}{IGBP class 16 fraction.}
#' }
#'
#' @source EcoISEA3H project. IUCN Red List spatial data, WorldClim v2,
#'   ENVIREM, SRTM, MODIS land cover.
#'
#' @examples
#' data(Ovis_ammon)
#' head(Ovis_ammon)
#' table(Ovis_ammon$presence)
"Ovis_ammon"


#' China Environmental Grid (ISEA3H Res9)
#'
#' Full environmental data grid for spatial prediction across China on the
#' ISEA3H Resolution 9 hexagonal grid. Contains the same 50 environmental
#' variables as the species datasets (see [Ovis_ammon]).
#'
#' @format A data frame with 3648 rows and 53 columns. Columns are: HID,
#'   lon, lat, bio01, bio02, bio03, bio04, bio05, bio06, bio07, bio12,
#'   bio13, bio14, bio15, bio16, bio17, bio18, bio19, elevation,
#'   annualpet, aridityindexthornthwaite, climaticmoistureindex,
#'   embergerq, maxtempcoldest, petseasonality, tri, topowet, nontree,
#'   bare, tree, etccdi_cdd, etccdi_cwd, etccdi_fd, etccdi_gsl,
#'   etccdi_rx1day, etccdi_rx5day, etccdi_sdii, igbp_01 through
#'   igbp_16. See [Ovis_ammon] for column descriptions.
#'
#' @source EcoISEA3H project.
#'
#' @examples
#' data(china_env_grid)
#' head(china_env_grid)
"china_env_grid"


# -- Bundled species datasets -------------------------------------------------
# All 32 species share the same ISEA3H Res9 grid and column structure as
# Ovis_ammon (3648 rows, 59 columns). Each is documented with a single
# roxygen block below.

#' @rdname species_datasets
#' @name species_datasets
#' @aliases Alces_alces Bos_mutus Budorcas_taxicolor Capra_sibirica
#'   Capreolus_pygargus Cervus_albirostris Cervus_canadensis
#'   Elaphodus_cephalophus Equus_kiang Gazella_subgutturosa
#'   Macaca_arctoides Macaca_assamensis Macaca_mulatta Macaca_thibetana
#'   Moschus_berezovskii Moschus_chrysogaster Moschus_fuscus
#'   Moschus_moschiferus Muntiacus_reevesi Muntiacus_vaginalis
#'   Naemorhedus_baileyi Nycticebus_bengalensis Pantholops_hodgsonii
#'   Procapra_picticaudata Pseudois_nayaur Rhinopithecus_roxellana
#'   Rusa_unicolor SID_14303 SID_3814 Sus_scrofa Trachypithecus_francoisi
#' @title Bundled Species Occurrence Datasets
#'
#' @description
#' Occurrence and environmental data for 31 additional mammal species on
#' the ISEA3H Resolution 9 hexagonal grid covering China. Each dataset
#' has the same structure as [Ovis_ammon]: 3648 rows and 59 columns
#' including `presence`, geographic coordinates, and 50 environmental
#' variables.
#'
#' Species: *Alces alces*, *Bos mutus*, *Budorcas taxicolor*,
#' *Capra sibirica*, *Capreolus pygargus*, *Cervus albirostris*,
#' *Cervus canadensis*, *Elaphodus cephalophus*, *Equus kiang*,
#' *Gazella subgutturosa*, *Macaca arctoides*, *Macaca assamensis*,
#' *Macaca mulatta*, *Macaca thibetana*, *Moschus berezovskii*,
#' *Moschus chrysogaster*, *Moschus fuscus*, *Moschus moschiferus*,
#' *Muntiacus reevesi*, *Muntiacus vaginalis*, *Naemorhedus baileyi*,
#' *Nycticebus bengalensis*, *Pantholops hodgsonii*,
#' *Procapra picticaudata*, *Pseudois nayaur*,
#' *Rhinopithecus roxellana*, *Rusa unicolor*, *Sus scrofa*,
#' *Trachypithecus francoisi*, plus two IUCN SID-referenced taxa
#' (SID 14303, SID 3814).
#'
#' @format See [Ovis_ammon] for the full column description.
#' @source EcoISEA3H project. IUCN Red List spatial data, WorldClim v2,
#'   ENVIREM, SRTM, MODIS land cover.
#' @seealso [Ovis_ammon], [china_env_grid]
#'
#' @examples
#' data(Alces_alces)
#' table(Alces_alces$presence)
NULL
