Analyze sub-district spatial heterogeneity in trachoma markers (Pgp3, TF, PCR)
in Amhara, Ethiopia using survey data and spatial modeling.

## Overview

This repository contains an R-based workflow that:

- Loads and cleans harmonized trachoma survey data for select studies in
  Amhara, Ethiopia.
- Builds cluster- and district-level summaries with geospatial coordinates.
- Computes spatial autocorrelation (Moran's I), variograms, and spatial models.
- Generates prediction grids and district-level spatial predictions.
- Produces tables, maps, and figures as Quarto reports.

## Repository structure

- `R/00-config.R`: Global configuration (libraries, paths, study list, markers,
  CRS, district metadata). Adjust this first to match your local data paths.
- `R/00-functions.R`: Shared helper functions (spatial statistics, variograms,
  model fitting, predictions).
- `R/01-download-hdx-data.R`: Optional script to download and unzip boundary
  shapefiles.
- `R/02-clean-data.R`: Data ingestion and cleaning; produces cluster-level
  summaries and processed shapefiles.
- `R/03-analysis.R`: Spatial analysis (Moran's I, variograms, spatial model
  fitting, prediction grid creation).
- `R/04-visualize-data.qmd`: Exploratory plots and maps.
- `R/05-results.qmd`: Result-focused tables and summaries.
- `R/06-figures.qmd`: Publication-ready figures.

## Inputs & outputs

**Inputs**
- Harmonized survey data RDS referenced in `R/00-config.R` (see
  `data_folder_path`).
- Ethiopia administrative boundary shapefiles (expected under
  `data/eth_admin_boundaries`).

**Outputs**
- Cleaned datasets and analysis artifacts under `data/clean/` (e.g.,
  `cluster_dat.rds`, `morans_i_results.rds`, `variograms.rds`,
  `spatial_models.rds`, `spatial_preds.rds`).
- HTML reports from the Quarto `.qmd` files.

## Workflow (by script)

Run scripts in order after updating paths in `R/00-config.R`:

1. **Configure paths and metadata**  
   Edit `R/00-config.R` to set `data_folder_path`, output locations, and
   review the markers and study lists.

2. **Optional: download boundary data**  
   `R/01-download-hdx-data.R` downloads an example shapefile bundle from HDX
   and unzips it into `data/clean/`.

3. **Clean and prepare datasets**  
   `R/02-clean-data.R`:
   - Loads harmonized survey data and Ethiopia admin boundaries.
   - Filters to the target studies and ages.
   - Computes cluster- and district-level summaries.
   - Writes cleaned data and shapefiles to `data/clean/`.

4. **Run spatial analyses**  
   `R/03-analysis.R`:
   - Calculates Moran's I over multiple neighbor definitions.
   - Builds and fits variograms.
   - Fits spatial models and creates prediction grids.
   - Saves model objects and prediction outputs to `data/clean/`.

5. **Generate reports and figures**  
   Render the Quarto files:
   - `R/04-visualize-data.qmd`: exploratory plots and maps.
   - `R/05-results.qmd`: summary results.
   - `R/06-figures.qmd`: final figures.

## Notes

- Prevalence data are publicly available at https://doi.org/10.5061/dryad.5qfttdzhx. The data used in
  this analysis was individual-level data (v5).
- GPS coordinates for cluster locations cannot be made public due to privacy considerations.

- Parallel processing is used in spatial modeling (see `n_cores` in
  `R/00-config.R` and `R/03-analysis.R`).
