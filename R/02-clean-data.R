# Setup ----

rm(list = ls())

source(here("R", "00-config.R"))

# Load raw data ----

# Load survey data
dat_raw <- readRDS(file.path(data_folder_path, "trachoma_serology_harmonized_data_indiv_v5.rds"))

# Load shapefiles
eth3_shp <- st_read(here("data", "eth_admin_boundaries", "eth_admin3.shp"))
# eth1_shp <- st_read(here("data", "eth_admin_boundaries", "eth_admin1.shp"))

# Pre-process survey data ----

dat <- dat_raw %>%
  filter(study_id %in% studies_of_interest) %>%
  # filter(district %in% districts_tib$districts) %>%
  rename(pgp3 = pgp3_pos,
         ct694 = ct694_pos)

## Cluster-level summaries ----

cluster_dat <- dat %>%
  group_by(study_id, cluster_id, district) %>%
  summarise(
    n_indiv = n(),
    mean_age = mean(age_years, na.rm = TRUE),
    across(c(tf, ti, pgp3, ct694, pcr),
           list(n = ~sum(!is.na(.)),
                pos = ~sum(., na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_pos"), 
                ~. / get(sub("_pos", "_n", cur_column())),
                .names = "{sub('_pos', '_prev', .col)}"))

### Cluster coordinates

cluster_gps <- dat %>%
  dplyr::select(study_id, country, district, eu, eu_name, cluster_id, lat, lon) %>%
  distinct() %>%
  group_by(study_id, country, district, eu, eu_name, cluster_id) %>%
  summarise(lat = mean(lat, na.rm = TRUE),
            lon = mean(lon, na.rm = TRUE),
            .groups = "drop") %>%
  distinct()

cluster_dat <- cluster_dat %>%
  left_join(cluster_gps, by = c("study_id", "district", "cluster_id")) %>%
  distinct()

## District-level summaries ----

cluster_dat_long <- cluster_dat %>%
  pivot_longer(starts_with("tf_") | starts_with("ti_") | starts_with("pgp3_") | starts_with("ct694_") | starts_with("pcr_"),
               names_to = c("marker", ".value"),
               names_sep = "_") %>%
  mutate(marker = recode(marker,
                         "tf" = "TF",
                         "ti" = "TI",
                         "pgp3" = "Pgp3",
                         "ct694" = "Ct694",
                         "pcr" = "PCR"))  %>%
  dplyr::filter(marker %in% markers) %>%
  mutate(district_n = gsub(" ", "\n", district)) %>%
  mutate(prev_perc = prev * 100) %>%
  mutate(study_id = factor(study_id, levels = studies_of_interest, ordered = TRUE),
         marker = factor(marker, levels = markers, ordered = TRUE))

### Order districts by median Pgp3 prevalence

district_order <- cluster_dat_long %>%
  dplyr::filter(marker == "Pgp3") %>%
  group_by(district) %>%
  summarise(median_prev = median(prev_perc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(median_prev)) %>%
  pull(district)

district_order_n <- gsub(" ", "\n", district_order)
district_order_n <- recode(district_order_n,
                           "Debre\nBirhan\nTown" = "Debre Birhan\nTown")

districts_tib <- districts_tib %>%
  mutate(district_n = ifelse(district == "Debre Birhan Town",
                             "Debre Birhan\nTown",
                             gsub(" ", "\n", district))) %>%
  mutate(district = factor(district, levels = district_order, ordered = TRUE),
         district_n = factor(district_n, levels = district_order_n, ordered = TRUE)) %>%
  arrange(district) %>%
  mutate(dist_color = viridis(n = nrow(.), option = "B", end = 0.9))

cluster_dat_long <- cluster_dat_long %>%
  mutate(district_n = ifelse(district == "Debre Birhan Town",
                             "Debre Birhan\nTown",
                             gsub(" ", "\n", district))) %>%
  mutate(district = factor(district, levels = district_order, ordered = TRUE),
         district_n = factor(district_n, levels = district_order_n, ordered = TRUE)) %>%
  arrange(district)

# Pre-process shapefiles ----

amhara_shp <- eth3_shp %>%
  filter(adm1_name == "Amhara") %>%
  mutate(survey_bin = ifelse(adm3_name %in% districts_tib$district_shp, TRUE, FALSE))

# Combine Goncha Siso Enebse and Sedae into Goncha 
# This is done because the survey data considers these two districts as one EU

combined_temp <- amhara_shp %>%
  dplyr::filter(adm3_name %in% c("Goncha Siso Enebse", "Sedae")) %>%
  summarise(
    area_sqkm = sum(area_sqkm),
    center_lat = mean(center_lat),
    center_lon = mean(center_lon),
    .groups = "drop"
  ) %>%
  mutate(adm3_name = "Goncha")

amhara_shp <- amhara_shp %>%
  dplyr::filter(!adm3_name %in% c("Goncha Siso Enebse", "Sedae")) %>%
  bind_rows(combined_temp)

# Save processed files ----

saveRDS(cluster_dat_long, here("data", "clean", "cluster_dat.rds"))

saveRDS(districts_tib, here("data", "clean", "districts_tib.rds"))
saveRDS(amhara_shp, here("data", "clean", "amhara_shp.rds"))

