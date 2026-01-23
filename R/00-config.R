# 8-spatial-analysis/R/00-config.R
# Spatial distribution of trachoma serology: config file

rm(list = ls())

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)
library(parallel)

## Spatial

library(sp)
library(sf)
library(spaMM)
library(gstat)
library(spdep)

## Plotting

library(ggrepel)
library(scales)
library(patchwork)
library(ggspatial)
library(viridis)
library(paletteer)

# Define global objects ---------------------------------------------------

data_folder_path <- "/Users/ariktha/Library/CloudStorage/Box-Box/trachoma-endgame/Data/final-v5"

# path_data_clean <- here("data", "clean")
n_cores <- max(1L, parallel::detectCores() - 1L)

proj_crs <- 32637  # WGS 84 / UTM zone 37N

markers <- c(
  "pgp3" = "Pgp3",
  # "ct694" = "Ct694",
  "tf" = "TF",
  # "ti" = "TI",
  "pcr" = "PCR"
)

studies_of_interest <- c("TCC-Ethiopia2019",
                         "TCC-Ethiopia2021",
                         "TCC-Amhara2022",
                         "TCC-Ethiopia2023")

# District labels and coordinates for maps -------------------------------

districts_tib <- tibble(
  district_shp = c(
    "Metema",
    "Ebenat",
    "Woreta town",
    "Fagta Lakoma",
    "Bibugn",
    "Goncha",
    "Michakel",
    "Debay Telatgen",
    "Albuko",
    "Debre Berhan town",
    "Fogera",
    "Tach Gayint"
  ),
  district = c(
    "Metema",
    "Ebinat",
    "Woreta Town",
    "Fagita Lekoma",
    "Bibugn",
    "Goncha",
    "Machakel",
    "Debay Tilatgin",
    "Albuko",
    "Debre Birhan Town",
    "Fogera",
    "Tach Gaynt"
  ),
  lat_lbl = c(
    12.95,  # Metema 
    12.55,  # Ebenat 
    11.90,  # Woreta town
    11.35,  # Fagta Lakoma
    10.85,  # Bibugn
    11.10,  # Goncha Siso Enbese
    10.25,  # Michaekel
    10.35,  # Debay Telatgen
    11.15,  # Albuko
    9.45,    # Debre Berhan town
    11.5,  # Fogera
    11.55   # Tach Gayint
  ),
  lon_lbl = c(
    35.75,  # Metema
    38.80,  # Ebenat
    37.05,  # Woreta town
    36.65,  # Fagta Lakoma
    37.30,  # Bibugn
    38.75,  # Goncha Siso Enbese
    37.35,  # Michaekel
    38.75,  # Debay Telatgen
    39.85,  # Albuko
    39.45,   # Debre Berhan town
    37.90,  # Fogera
    39.20   # Tach Gayint
  )
)
