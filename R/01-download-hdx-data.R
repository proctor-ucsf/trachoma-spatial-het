
# South Sudan admin boundaries
### Map downloaded from https://data.humdata.org/dataset/cod-ab-ssd

ssd_shp_url <- "https://data.humdata.org/dataset/cdd62bd9-e442-4eac-9b44-cfee8bf79153/resource/4d768ab6-c322-4af4-8fbf-6cbb442812d3/download/ssd_admbnda_imwg_nbs_20230829_shp.zip"

ssd_zip <- file.path(path_data_clean, "ssd_shapefiles.zip")
ssd_dir <- file.path(path_data_clean, "ssd_shapefiles")

if (!dir.exists(path_data_clean)) {
  dir.create(path_data_clean, recursive = TRUE)
}

download.file(ssd_shp_url, destfile = ssd_zip)

if (!dir.exists(ssd_dir)) {
  dir.create(ssd_dir, recursive = TRUE)
}

unzip(zipfile = ssd_zip, exdir = ssd_dir)

