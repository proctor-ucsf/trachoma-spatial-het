# Setup ----

rm(list = ls())

library(here)

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))

# Load data ----

cluster_dat <- readRDS(here("data", "clean", "cluster_dat.rds"))
districts_tib <- readRDS(here("data", "clean", "districts_tib.rds"))
amhara_shp <- readRDS(here("data", "clean", "amhara_shp.rds"))

# Moran's I ----

# Define values of k (number of nearest neighbors) to test
k_grid <- c(4, 5, 6)

# Wrapper to pass k to compute_morans_i
compute_morans_i_k <- function(df, k) {
  compute_morans_i(df, k = k)
}

# Compute Moran's I for all district + marker + k combos
morans_results <- cluster_dat %>%
  tidyr::crossing(k = k_grid) %>%
  group_by(district, marker, k) %>%
  tidyr::nest() %>%
  mutate(res = purrr::map2(data, k, ~ compute_morans_i_k(.x, k = .y))) %>%
  tidyr::unnest(res) %>%
  dplyr::select(-data) %>%
  ungroup()

saveRDS(morans_results, here("data", "clean", "morans_i_results.rds"))

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

null_to_na_real <- function(x) {
  if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)
}

spatial_models_fit <- variograms_meta %>%
  mutate(models = pmap(list(district, marker), ~fit_spatial_model(..1, ..2, cluster_dat)),
         spatial_model = map(models, ~.x$spatial_model),
         binomial_model = map(models, ~.x$binomial_model),
         simple_model = map(models, ~.x$simple_model),
         marg_aic_spatial   = map(models, ~ null_to_na_real(.x$aic_spatial[[1]])),
         marg_aic_binomial  = map(models, ~ null_to_na_real(.x$aic_binomial[[1]])),
         marg_aic_simple    = map(models, ~ null_to_na_real(.x$aic_simple[[1]])),
         error_type = map_chr(models, ~if(!is.null(.x$error)) .x$error else NA_character_)) 

spatial_models <- spatial_models_fit %>%
  mutate(across(starts_with("marg_aic"), ~unlist(.x))) %>%
  mutate(spatial_fit_yn = ifelse(!map_lgl(spatial_model, is.null), "yes", "no"),
         binomial_fit_yn = ifelse(!map_lgl(binomial_model, is.null), "yes", "no"),
         simple_fit_yn = ifelse(!map_lgl(simple_model, is.null), "yes", "no")) %>%
  mutate(fit_yn = case_when(
    spatial_fit_yn == "yes" & simple_fit_yn == "yes" ~ "both fit",
    simple_fit_yn == "yes" ~ "only simple fit",
    spatial_fit_yn == "yes" ~ "only spatial fit",
    TRUE ~ "neither fit"))  %>%
  mutate(marg_aic_diff = marg_aic_spatial - marg_aic_simple,
         pref_model = case_when(is.na(marg_aic_diff) ~ "neither",
                                marg_aic_diff < -2 ~ "spatial",
                                marg_aic_diff > 2 ~ "simple",
                                 TRUE ~ "tie")) %>%
  mutate(lambda = spatial_model$lambda,
         nu = spatial_model$corrPars[[1]]$nu,
         rho = spatial_model$corrPars[[1]]$rho)

# Save spatial models

saveRDS(spatial_models, here("data", "clean", "spatial_models.rds"))

## Generate prediction grids and predictions ----

# Create prediction grids for each district
district_grids <- cluster_dat %>%
  left_join(districts_tib, by = "district") %>%
  filter(!is.na(district_shp)) %>%
  distinct(district_shp)

# Make and use a cluster for parallel processing of grid creation
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

## Generate predictions for each district + marker combo ----

# Get predictions for spatial model 

tictoc::tic("Generate predictions for spatial model")

pred_list <- mclapply(
  X = seq_len(nrow(spatial_models)),
  FUN = function(i) {
    get_predictions_for_model(
      spatial_models$district[[i]],
      spatial_models$marker[[i]],
      spatial_models$spatial_model[[i]],
      district_grids
    )
  },
  mc.cores = n_cores
)

spatial_models <- spatial_models %>%
  mutate(predictions = pred_list)

tictoc::toc()

# Get predictions for binomial model 

tictoc::tic("Generate predictions for binomial model")

binom_pred_list <- mclapply(
  X = seq_len(nrow(spatial_models)),
  FUN = function(i) {
    get_predictions_for_model(
      spatial_models$district[[i]],
      spatial_models$marker[[i]],
      spatial_models$binomial_model[[i]],
      district_grids
    )
  },
  mc.cores = n_cores
)

spatial_models <- spatial_models %>%
  mutate(binomial_predictions = binom_pred_list)

tictoc::toc()

saveRDS(spatial_models, here("data", "clean", "spatial_preds.rds"))
