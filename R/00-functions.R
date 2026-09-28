# collect_warnings
# Purpose: Evaluate an expression while capturing (and muffling) warnings.
# Inputs:
#   - expr: an expression to evaluate.
# Outputs:
#   - list with:
#     - value: evaluated result of expr.
#     - warnings: unique character vector of warning messages.

collect_warnings <- function(expr) {
  warns <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = unique(warns))
}

# Moran's I ----

# compute_morans_i
# Purpose: Compute Moran's I for a data frame of point prevalence values.
# Inputs:
#   - df: data frame with lon, lat, prev columns.
#   - k: number of nearest neighbors for spatial weights.
#   - nsim: number of Monte Carlo simulations.
# Outputs:
#   - tibble with morans_I, p_value, and n (number of points).

compute_morans_i <- function(df, k = 5, nsim = 1000) {
  
  # Need at least k + 1 points
  if (nrow(df) <= k) {
    return(tibble(
      morans_I = NA_real_,
      p_value  = NA_real_,
      n        = nrow(df)
    ))
  }
  
  coords <- as.matrix(df[, c("lon", "lat")])
  
  # k-nearest neighbors
  knn <- knearneigh(coords, k = k)
  nb  <- knn2nb(knn)
  
  # Row-standardized weights
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Monte Carlo Moran's I
  mi <- moran.mc(
    x = df$prev,
    listw = lw,
    nsim = nsim,
    zero.policy = TRUE
  )
  
  tibble(
    morans_I = unname(mi$statistic),
    p_value  = mi$p.value,
    n        = nrow(df)
  )
}


# Variograms ----

# empty_variogram_result
# Purpose: Standardize the empty/failed variogram result structure.
# Inputs:
#   - dist_i: district name.
#   - mrk: marker name.
#   - n_points: number of points in the subset.
#   - cutoff, width, n_bins, n_bins_populated, np_min, np_median: metadata.
#   - warnings: character vector of warning messages.
# Outputs:
#   - tibble with variogram metadata and NULL variogram object.
empty_variogram_result <- function(dist_i,
                                   mrk,
                                   n_points,
                                   cutoff = NA_real_,
                                   width = NA_real_,
                                   n_bins = 0L,
                                   n_bins_populated = 0L,
                                   np_min = NA_real_,
                                   np_median = NA_real_,
                                   warnings = character()) {
  tibble(
    district = dist_i,
    marker   = mrk,
    n_points = n_points,
    cutoff   = cutoff,
    width    = width,
    n_bins   = n_bins,
    n_bins_populated = n_bins_populated,
    np_min   = np_min,
    np_median = np_median,
    warnings = list(unique(warnings)),
    vgm      = list(NULL)
  )
}

# compute_variogram
# Purpose: Compute an empirical variogram for one district-marker combination.
# Inputs:
#   - i: row index into combos.
#   - df: cluster-level data frame with lon/lat/prev columns.
#   - combos: data frame with district and marker columns.
#   - proj_crs: target CRS for distance calculations.
#   - min_points: minimum points required to attempt a variogram.
#   - target_bins, min_bins_populated, min_pairs_per_bin: binning controls.
#   - max_dist_cap: maximum distance cap (meters).
# Outputs:
#   - tibble with variogram metadata and vgm object (or NULL if skipped).

compute_variogram <- function(i,
                              df,
                              combos,
                              proj_crs = proj_crs,
                              min_points = 10,
                              target_bins = 15,
                              min_bins_populated = 5,
                              min_pairs_per_bin = 1,
                              max_dist_cap = 50*1000) {
  dist_i <- combos$district[i]
  mrk    <- combos$marker[i]
  
  warnings_all <- character(0)
  
  subset_df <- df %>%
    dplyr::filter(district == dist_i, marker == mrk) %>%
    transmute(district, marker, lon, lat, prev) %>%
    dplyr::filter(!is.na(lon), !is.na(lat), !is.na(prev)) %>%
    mutate(lon  = as.numeric(lon),
           lat  = as.numeric(lat),
           prev = as.numeric(prev))
  
  n_points <- nrow(subset_df)
  if (n_points < min_points) {
    warnings_all <- c(warnings_all,
                      paste0("Too few points (", n_points, " < ", min_points, "). Skipping."))
    return(empty_variogram_result(dist_i, mrk, n_points, warnings = warnings_all))
  }
  
  # sf conversion + projection
  sf_obj <- st_as_sf(
    subset_df,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  sf_obj <- st_transform(sf_obj, proj_crs)
  
  # max pair distance
  dmat <- st_distance(sf_obj)
  max_pair_distance <- suppressWarnings(as.numeric(max(dmat, na.rm = TRUE))) # drop units to numeric
  
  if (!is.finite(max_pair_distance) || max_pair_distance <= 0) {
    warnings_all <- c(warnings_all, "Invalid max pairwise distance. Skipping.")
    return(empty_variogram_result(dist_i, mrk, n_points, warnings = warnings_all))
  }
  
  cutoff <- 0.5 * max_pair_distance
  if (is.finite(max_dist_cap))
    cutoff <- min(cutoff, max_dist_cap)
  
  # adaptive bin width: aim for ~target_bins
  width <- cutoff / target_bins
  if (!is.finite(width) || width <= 0) {
    warnings_all <- c(warnings_all, "Invalid adaptive bin width. Skipping.")
    return(
      empty_variogram_result(
        dist_i,
        mrk,
        n_points,
        cutoff = cutoff,
        warnings = warnings_all
      )
    )
  }
  
  # gstat::variogram expects sp-ish; convert sf -> sp safely
  sp_obj <- as(sf_obj, "Spatial")
  
  res <- collect_warnings(gstat::variogram(
    prev ~ 1,
    data   = sp_obj,
    cutoff = cutoff,
    width  = width
  ))
  warnings_all <- c(warnings_all, res$warnings)
  
  vgm_df <- res$value
  if (is.null(vgm_df) || nrow(vgm_df) == 0) {
    warnings_all <- c(warnings_all, "Variogram returned no bins. Skipping.")
    return(
      empty_variogram_result(
        dist_i,
        mrk,
        n_points,
        cutoff = cutoff,
        width = width,
        warnings = warnings_all
      )
    )
  }
  
  # store per-bin np and meta
  n_bins <- nrow(vgm_df)
  
  vgm_pop <- vgm_df[vgm_df$np >= min_pairs_per_bin, , drop = FALSE]
  n_bins_populated <- nrow(vgm_pop)

  if (n_bins_populated < min_bins_populated) {
    warnings_all <- c(
      warnings_all,
      paste0(
        "Too few populated bins (",
        n_bins_populated,
        " < ",
        min_bins_populated,
        ")."
      )
    )
    
    return(
      empty_variogram_result(
        dist_i,
        mrk,
        n_points,
        cutoff = cutoff,
        width  = width,
        n_bins = n_bins,
        n_bins_populated = n_bins_populated,
        np_min = if (n_bins_populated > 0)
          min(vgm_pop$np, na.rm = TRUE)
        else
          NA_real_,
        np_median = if (n_bins_populated > 0)
          median(vgm_pop$np, na.rm = TRUE)
        else
          NA_real_,
        warnings = warnings_all
      )
    )
  }
  
  tibble(
    district = dist_i,
    marker   = mrk,
    n_points = n_points,
    cutoff   = cutoff,
    width    = width,
    n_bins   = n_bins,
    n_bins_populated = n_bins_populated,
    np_min   = min(vgm_pop$np, na.rm = TRUE),
    np_median = median(vgm_pop$np, na.rm = TRUE),
    warnings = list(unique(warnings_all)),
    vgm      = list(vgm_df)
  )
}

# fit_variogram_models
# Purpose: Fit candidate variogram models and select the best by SSE.
# Inputs:
#   - vgm_df: empirical variogram data frame.
#   - fit_methods: gstat fit methods to try.
#   - min_bins_populated: minimum bins required to attempt fitting.
#   - min_pairs_per_bin: minimum pairs per bin.
#   - models: variogram model names to consider.
#   - kappa: Matern smoothness parameter.
# Outputs:
#   - list with fitted models, diagnostics, and best model.

fit_variogram_models <- function(vgm_df,
                                 fit_methods = c(1, 2, 6, 7),
                                 min_bins_populated = 5,
                                 min_pairs_per_bin = 1,
                                 models = c("Mat", "Sph"),
                                 kappa = 0.5) {
  if (is.null(vgm_df) || nrow(vgm_df) == 0) {
    return(list(
      fits = list(),
      best = NULL,
      diagnostics = tibble(),
      warnings = character(0)
    ))
  }
  
  warnings_all <- character(0)
  
  vgm_use <- as_tibble(vgm_df) %>%
    filter(np >= min_pairs_per_bin) %>%
    filter(is.finite(gamma), is.finite(dist))
  
  class(vgm_use) <- class(vgm_df)
  
  if (nrow(vgm_use) < min_bins_populated) {
    warnings_all <- c(warnings_all,
                      paste0(
                        "Insufficient populated bins for fitting (",
                        nrow(vgm_use),
                        ")."
                      ))
    return(list(
      fits = list(),
      best = NULL,
      diagnostics = tibble(),
      warnings = unique(warnings_all)
    ))
  }
  
  # More stable initial values
  init_nugget <- max(0, suppressWarnings(quantile(
    vgm_use$gamma, probs = 0.05, na.rm = TRUE
  )))
  init_sill_total <- suppressWarnings(quantile(vgm_use$gamma, probs = 0.95, na.rm = TRUE))
  init_psill <- max(0, init_sill_total - init_nugget)
  init_range <- suppressWarnings(quantile(vgm_use$dist, probs = 0.67, na.rm = TRUE))
  if (!is.finite(init_range) ||
      init_range <= 0)
    init_range <- max(vgm_use$dist, na.rm = TRUE) / 3
  
  # Build candidate "initial" models
  init_models <- list()
  if ("Mat" %in% models) {
    init_models$Mat <- gstat::vgm(
      psill = init_psill,
      model = "Mat",
      range = init_range,
      nugget = init_nugget,
      kappa = kappa
    )
  }
  if ("Sph" %in% models) {
    init_models$Sph <- gstat::vgm(
      psill = init_psill,
      model = "Sph",
      range = init_range,
      nugget = init_nugget
    )
  }
  
  # Fit across (model x fit.method) and compare SSErr
  fits <- list()
  diag_rows <- list()
  
  for (mname in names(init_models)) {
    for (fm in fit_methods) {
      res <- collect_warnings(tryCatch(
        gstat::fit.variogram(vgm_use, model = init_models[[mname]], fit.method = fm),
        error = function(e) {
          warning(paste0("Fit failed (", mname, ", method ", fm, "): ", e$message))
          NULL
        }
      ))
      warnings_all <- c(warnings_all, res$warnings)
      
      fit_obj <- res$value
      key <- ifelse(length(fit_methods) > 1, paste0(mname, "_method", fm), mname)
      
      if (is.null(fit_obj)) {
        fits[[key]] <- NULL
        diag_rows[[key]] <- tibble(
          model = mname,
          fit_method = fm,
          SSErr = NA_real_,
          ok = FALSE
        )
      } else {
        sse <- attr(fit_obj, "SSErr")
        fits[[key]] <- fit_obj
        diag_rows[[key]] <- tibble(
          model = mname,
          fit_method = fm,
          SSErr = as.numeric(sse),
          ok = TRUE
        )
      }
    }
  }
  
  diagnostics <- bind_rows(diag_rows) %>%
    arrange(!ok, SSErr, model, fit_method)
  
  best_key <- diagnostics %>%
    filter(ok, is.finite(SSErr)) %>%
    slice(1) %>%
    transmute(key = paste0(model, "_method", fit_method)) %>%
    pull(key)
  
  best_fit <- if (length(best_key) == 1)
    fits[[best_key]]
  else
    NULL
  
  list(
    fits = fits,
    best = best_fit,
    diagnostics = diagnostics,
    warnings = unique(warnings_all)
  )
}

# extract_best_name_and_range
# Purpose: Extract the best variogram model name and range from fit results.
# Inputs:
#   - fit: list returned by fit_variogram_models.
# Outputs:
#   - tibble with best_fit_name and best_range.

extract_best_name_and_range <- function(fit) {
  if (is.null(fit$best) || nrow(fit$diagnostics) == 0) {
    return(tibble(
      best_fit_name = NA_character_,
      best_range = NA_real_
    ))
  }
  
  # 1. Identify the winning row from diagnostics
  best_row <- fit$diagnostics %>%
    filter(ok, is.finite(SSErr)) %>%
    slice(1)
  
  best_fit_name <- paste0(
    best_row$model,
    "_method",
    best_row$fit_method
  )
  
  # 2. Extract range from the fitted variogram (non-nugget component)
  bf <- as_tibble(fit$best) %>%
    mutate(model = as.character(model))
  
  best_range <- bf %>%
    filter(model != "Nug") %>%
    pull(range) %>%
    { if (length(.) == 0) NA_real_ else as.numeric(.[1]) }
  
  tibble(
    best_fit_name = best_fit_name,
    best_range = best_range
  )
}


# Spatial models ----

# fit_spatial_model
# Purpose: Fit spatial (Matern) and non-spatial Gaussian models for a district-marker.
# Inputs:
#   - dist: district name.
#   - mrk: marker name.
#   - df: cluster-level data with prev, lon, lat.
# Outputs:
#   - list with spatial_model, simple_model, AICs, and optional error flag.

fit_spatial_model <- function(dist, mrk, df) {
  
  subset_data <- df %>%
    dplyr::filter(district == dist, marker == mrk) %>%
    dplyr::filter(!is.na(prev), !is.na(lat), !is.na(lon))
  
  if (nrow(subset_data) < 10) {
    warning(sprintf("Insufficient data for %s - %s: only %d observations", dist, mrk, nrow(subset_data)))
    return(list(spatial_model=NULL, simple_model=NULL, aic_spatial=NA_real_, aic_simple=NA_real_, error="insufficient_data"))
  }
  
  if (is.na(stats::var(subset_data$prev)) || stats::var(subset_data$prev) < 1e-6) {
    warning(sprintf("Zero/almost-zero variance in outcome for %s - %s", dist, mrk))
    return(list(spatial_model=NULL, simple_model=NULL, aic_spatial=NA_real_, aic_simple=NA_real_, error="zero_variance"))
  }
  
  if (stats::var(subset_data$lat) == 0 || stats::var(subset_data$lon) == 0) {
    warning(sprintf("No spatial variation for %s - %s (lat or lon constant)", dist, mrk))
    return(list(spatial_model=NULL, simple_model=NULL, aic_spatial=NA_real_, aic_simple=NA_real_, error="no_spatial_variation"))
  }
  
  # Keep warnings but don't treat as failure
  warn_spatial <- character(0)
  fit_spatial <- tryCatch(
      spaMM::fitme(
        prev ~ 1 + Matern(1 | lon + lat),
        data = subset_data,
        family = gaussian(link = "identity")
      ),
      warning = function(w) {
        # warn_spatial <<- c(warn_spatial, conditionMessage(w))
        # invokeRestart("muffleWarning")
        warning(sprintf("Spatial model warned for %s - %s: %s", dist, mrk, w$message))
        NULL
      },
    error = function(e) {
      warning(sprintf("Spatial model failed for %s - %s: %s", dist, mrk, e$message))
      NULL
    }
  )
  
  warn_binomial <- character(0)
  fit_binomial <- tryCatch(
      spaMM::fitme(
        cbind(pos, n - pos) ~ 1 + Matern(1 | lon + lat),
        data = subset_data,
        family = binomial(link = "logit")
      ),
      warning = function(w) {
        # warn_binomial <<- c(warn_binomial, conditionMessage(w))
        # invokeRestart("muffleWarning")
        warning(sprintf("Binomial spatial model warned for %s - %s: %s", dist, mrk, w$message))
        NULL
      },
    error = function(e) {
      warning(sprintf("Binomial spatial model failed for %s - %s: %s", dist, mrk, e$message))
      NULL
    }
  )
  
  warn_simple <- character(0)
  fit_simple <- tryCatch(
      spaMM::fitme(
        prev ~ 1,
        data = subset_data,
        family = gaussian(link = "identity")
      ),
      warning = function(w) {
        # warn_simple <<- c(warn_simple, conditionMessage(w))
        # invokeRestart("muffleWarning")
        warning(sprintf("Spatial model warned for %s - %s: %s", dist, mrk, w$message))
        NULL
      },
    error = function(e) {
      warning(sprintf("Nonspatial model failed for %s - %s: %s", dist, mrk, e$message))
      NULL
    }
  )
  
  aic_spatial <- if (!is.null(fit_spatial)) AIC(fit_spatial) else NA_real_
  aic_binomial <- if (!is.null(fit_binomial)) AIC(fit_binomial) else NA_real_
  aic_simple  <- if (!is.null(fit_simple))  AIC(fit_simple)  else NA_real_
  
  lambda_spatial <- if (!is.null(fit_spatial)) fit_spatial$lambda else NA_real_
  lambda_binomial <- if (!is.null(fit_binomial)) fit_binomial$lambda else NA_real_
  
  list(
    spatial_model = fit_spatial,
    binomial_model = fit_binomial,
    simple_model  = fit_simple,
    aic_spatial   = aic_spatial,
    aic_binomial  = aic_binomial,
    aic_simple    = aic_simple,
    lambda_spatial = lambda_spatial,
    lambda_binomial = lambda_binomial,
    warnings_spatial = warn_spatial,
    warnings_binomial = warn_binomial,
    warnings_simple  = warn_simple,
    error = NULL
  )
}

# create_district_grid
# Purpose: Create a grid of points within a district polygon.
# Inputs:
#   - dist_shp_name: district name as used in shapefile.
#   - shapefile: sf object with district polygons.
#   - data: cluster-level data.
#   - districts_lookup: lookup table with district name mappings.
#   - n_grid: number of grid cells per side.
# Outputs:
#   - sf geometry collection of grid points within the district polygon.

create_district_grid <- function(dist_shp_name, shapefile, data, districts_lookup, n_grid = 50) {
  
  # Get the district polygon from shapefile
  district_polygon <- shapefile %>%
    filter(adm3_name == dist_shp_name)
  
  # Filter data for this district
  district_data <- data %>%
    left_join(districts_lookup, by = "district") %>%
    filter(district_shp == dist_shp_name)
  
  # If no data for this district, return empty
  if(nrow(district_data) == 0) {
    return(st_sfc(crs = 4326))
  }
  
  # Create grid over the district polygon
  grid <- st_make_grid(district_polygon, n = c(n_grid, n_grid), what = "centers")
  
  # Keep only grid points that fall within the district boundary
  grid_within <- st_intersection(grid, district_polygon)
  
  return(grid_within)
}

# get_grid_preds
# Purpose: Predict prevalence over a grid of sf points using a spaMM model.
# Inputs:
#   - input_grid: sf points (grid).
#   - spamm_model_fit: fitted spaMM model.
# Outputs:
#   - sf object with predicted values in column `pred`.

get_grid_preds <- function(input_grid, spamm_model_fit) {
  # Convert grid to data frame
  grid_coords <- st_coordinates(input_grid)
  grid_df <- data.frame(
    lon = grid_coords[,1],
    lat = grid_coords[,2]
  )
  
  # Get predictions
  preds <- predict(spamm_model_fit, newdata = grid_df, type = "response")
  
  # Return as sf object
  grid_df$pred <- preds[,1]
  st_as_sf(grid_df, coords = c("lon", "lat"), crs = 4326)
}

# get_predictions_for_model
# Purpose: Run spatial predictions for a district/marker using a fitted model.
# Inputs:
#   - dist: district name.
#   - mrk: marker name.
#   - model: fitted spaMM model.
#   - grid_data: tibble with grid geometry per district.
# Outputs:
#   - tibble with x/y/value/district/marker columns, or NULL on failure.

get_predictions_for_model <- function(dist, mrk, model, grid_data) {
  if(is.null(model)) return(NULL)
  
  study_grid <- grid_data %>%
    filter(district == dist) %>%
    pull(grid) %>%
    .[[1]]
  
  preds_sf <- tryCatch(
    get_grid_preds(study_grid, model),
    error = function(e) {
      warning(paste0("Prediction failed for ", district, " - ", mrk, ": ", e$message))
      return(NULL)
    }
  )
  
  if(is.null(preds_sf)) return(NULL)
  
  # Convert to raster then tibble
  preds_coords <- st_coordinates(preds_sf)
  
  pred_tibble <- tibble(
    x = preds_coords[,1],
    y = preds_coords[,2],
    value = preds_sf$pred,  # Convert to percentage
    district = dist,
    marker = mrk
  )
  
  return(pred_tibble)
}

