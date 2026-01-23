library(here)

# Ethiopia admin boundaries
### Map downloaded from https://data.humdata.org/dataset/cod-ab-eth

eth_shp_url <- "https://data.humdata.org/dataset/cb58fa1f-687d-4cac-81a7-655ab1efb2d0/resource/274872ef-5add-48f4-95a4-5ce162af7f3f/download/eth_admin_boundaries.shp.zip"
eth_zip <- here("data", "eth_admin_boundaries.zip")
eth_dir <- here("data", "eth_admin_boundaries")

download.file(eth_shp_url, destfile = eth_zip)

unzip(zipfile = eth_zip, exdir = eth_dir)

