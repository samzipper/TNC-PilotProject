## Navarro_Habitat_CalculateWellStreamPairs.R
#' This script is intended to calculate analytical model inputs for each well-stream combination:
#'   -d = distance to stream [m]
#'   -S = effective storage coefficient [-]
#'   -Tr = effective transmissivity [L2/T]

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

## load data
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("modflow", "input", "iriv.shp"), stringsAsFactors=F)

# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# synthetic pumping wells
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), header=T)
sf.wel <- st_as_sf(df.wel, coords = c("lon", "lat"), crs = st_crs(sf.streams))

# rasters
r.lulc <- raster(paste0(dir.gis, "Navarro_Habitat_LULC_30m.tif"))
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Habitat_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Habitat_DTB_30m.tif"))
r.aquifers <- raster(paste0(dir.gis, "Navarro_Habitat_GroundwaterBasins_30m.tif"))

## extract some potentially relevant data
sf.wel$elev_m <- raster::extract(r.dem.30m, sf.wel)  # elevation of that grid cell [m]
sf.wel$dtb_m <- raster::extract(r.dtb, sf.wel)       # depth to bedrock [m]
sf.streams$elev_m <- raster::extract(r.dem.30m, sf.streams, fun='mean', na.rm=T)  # mean elevation of all grid cells stream touches [m]
sf.streams$dtb_m <- raster::extract(r.dtb, sf.streams, fun='mean', na.rm=T)       # mean DTB of all grid cells stream touches [m]

## predict stream width based on drainage area
sf.streams$width_m <- WidthFromDA(DA=sf.streams$TtDASKM, w.min=1, w.max=100)

## distance from each well to each stream segment
dist_all <- 
  st_distance(x=sf.streams, y=sf.wel)

df.all <- data.frame(
  SegNum = rep(sf.streams$SegNum, times = dim(dist_all)[2]),
  stream_elev_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  stream_dtb_m = rep(sf.streams$dtb_m, times = dim(dist_all)[2]),
  stream_width_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  WellNum = rep(sf.wel$WellNum, each = dim(dist_all)[1]),
  well_elev_m = rep(sf.wel$elev_m, times = dim(dist_all)[1]),
  well_dtb_m = rep(sf.wel$dtb_m, times = dim(dist_all)[1]),
  dist_wellToStream_m = as.numeric(dist_all)
)

## save output
df.all %>% 
  format(digits = 3) %>% 
  write.csv(file.path("results", "Navarro_Habitat_CalculateWellStreamPairs.csv"),
            row.names=F)
