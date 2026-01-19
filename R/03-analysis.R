# Setup ----

rm(list = ls())

library(here)

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))

# Load data ----

cluster_dat <- readRDS(here("data", "clean", "cluster_dat.rds"))
districts_tib <- readRDS(here("data", "clean", "districts_tib.rds"))
amhara_shp <- readRDS(here("data", "clean", "amhara_shp.rds"))

# Variograms ----

## Compute empirical variograms for all district + marker combos ----

combos <- cross_join(tibble(district = districts_tib$district), tibble(marker = unname(markers)))

variograms_meta <- parallel::mclapply(
  1:nrow(combos),
  compute_variogram,
  df = cluster_dat,
  combos = combos,
  proj_crs = proj_crs,
  mc.cores = n_cores
)
variograms_meta <- bind_rows(variograms_meta)

## Fit variogram models ----
fit_results <- variograms_meta %>%
  filter(!map_lgl(vgm, is.null)) %>%
  mutate(fit = map(vgm, ~ fit_variogram_models(.x))) %>%
  mutate(best_fit = map(fit, "best")) %>%
  mutate(district = factor(district, levels = districts_tib$district, ordered = TRUE))

fit_results <- variograms_meta %>%
  filter(!map_lgl(vgm, is.null)) %>%
  mutate(
    fit = map(vgm, fit_variogram_models),
    best_summary = map(fit, extract_best_name_and_range)
  ) %>%
  tidyr::unnest_wider(best_summary) %>%
  mutate(
    district = factor(district, levels = districts_tib$district, ordered = TRUE)
  )

# Save empirical variograms + metadata + warnings + fits
saveRDS(fit_results, here("data", "clean", "variograms.rds"))

# Spatial models ----

## Fit spatial models for all district + marker combos ----

spatial_models <- variograms_meta %>%
  mutate(models = pmap(list(district, marker), ~fit_spatial_model(..1, ..2, cluster_dat)),
         spatial_model = map(models, ~.x$spatial_model),
         simple_model = map(models, ~.x$simple_model),
         marg_aic_spatial = map(models, ~.x$aic_spatial[[1]]),
         marg_aic_simple = map(models, ~.x$aic_simple[[1]]),
         error_type = map_chr(models, ~if(!is.null(.x$error)) .x$error else NA_character_)) %>%
  mutate(lambda = spatial_model$lambda,
         nu = spatial_model$corrPars[[1]]$nu,
         rho = spatial_model$corrPars[[1]]$rho)

# Save spatial models

saveRDS(spatial_models, here("data", "clean", "spatial_models.rds"))

## Generate prediction grids and predictions ----

district_grids <- cluster_dat %>%
  left_join(districts_tib, by = "district") %>%
  filter(!is.na(district_shp)) %>%
  distinct(district_shp)

# Make a cluster for parallel processing
n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_cores)

parallel::clusterEvalQ(cl, {
  library(dplyr)
  library(purrr)
  library(sf)
})

parallel::clusterExport(
  cl,
  varlist = c("create_district_grid", "amhara_shp", "cluster_dat", "districts_tib"),
  envir = environment()
)

tictoc::tic("Parallel grid creation")

# Parallel compute grids
grids_list <- parallel::parLapply(
  cl,
  X = district_grids$district_shp,
  fun = function(shp_name) {
    create_district_grid(
      shp_name,
      shapefile = amhara_shp,
      data = cluster_dat,
      districts_lookup = districts_tib,
      n_grid = 1000
    )
  }
)

tictoc::toc()

parallel::stopCluster(cl)

# Attach back to tibble
district_grids <- district_grids %>%
  mutate(grid = grids_list) %>%
  left_join(districts_tib, by = "district_shp")

saveRDS(district_grids, here("data", "clean", "district_grids.rds"))

## Generate predictions for each district + marker combo ----

# Get predictions for all models

spatial_models <- spatial_models %>%
  mutate(
    predictions = pmap(
      list(district, marker, spatial_model),
      ~get_predictions_for_model(..1, ..2, ..3, district_grids)
    )
  )

saveRDS(spatial_models, here("output", "rds", "spatial_models_with_predictions.rds"))
