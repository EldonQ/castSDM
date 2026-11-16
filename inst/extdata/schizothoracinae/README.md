# Actual-data regression fixture

The 11,817 rows are the complete saved modelling sample from the author's
Schizothoracinae river/locality analysis: 1,817 occurrence records and 10,000
sampled background locations. `presence = 0` means background, not a surveyed
absence. The original training/hold-out partition and training spatial folds
are retained. There is no row subsampling, environmental imputation, coordinate
perturbation or manufactured response.

The six predictors are exactly the changing CHELSA climate variables selected
in that case: bio02, bio03, bio08, bio13, bio15 and bio19. This is an engineering
fixture for prediction contracts; fitting these six variables does not replace
the ecological analysis with its larger selected predictor set.

`future` contains actual extracted CHELSA ensemble-mean SSP585 2071-2100 values
at the same locations. The two TIFFs retain a 5-by-5 window of original cells
around an actual occurrence that encounters range extrapolation. No grid
resampling or new climate values are used. Neither the background sample nor
this window represents the entire geographic distribution of the subfamily.

## Sources and permissions

- Du, T., Ding, C., Yang, K., et al. SchiSOFT: A global dataset on species
  occurrences and functional traits of Schizothoracinae fish, version 1.
  https://doi.org/10.6084/m9.figshare.24638538.v1. Occurrence data: CC BY 4.0,
  https://creativecommons.org/licenses/by/4.0/. Records were taxonomically and
  spatially screened and thinned by the author's saved workflow.
- Karger, D. N., et al. (2021). Climatologies at high resolution for the earth's
  land surface areas. https://doi.org/10.16904/envidat.228.v2.1. CHELSA v2.1
  1981-2010 climatology; CC0 1.0.
- CHELSA v2.1 future bioclimatic data, ensemble-mean climate inputs supplied
  with the author's environmental library. Dataset description and CC0 1.0:
  https://www.chelsa-climate.org/datasets/chelsa_bioclim. Future values are
  scenario projections, not future observations or GCM uncertainty estimates.

`sources.csv` records original input basenames and content hashes.
`data-raw/build-real-fixture.R` reconstructs the fixture from the completed
case and physical source rasters. Original source files remain read-only.
