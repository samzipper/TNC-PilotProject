## Navarro_Cannabis_CalculateWellStreamPairs.R
#' This script is intended to calculate analytical model inputs for each well-stream combination:
#'   -d = distance to stream [m]
#'   -S = effective storage coefficient [-]
#'   -Tr = effective transmissivity [L2/T]

source(file.path("src", "paths+packages.R"))

## load data
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F)

# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# synthetic pumping wells
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), header=T)
sf.wel <- st_as_sf(df.wel, coords = c("lon", "lat"), crs = st_crs(sf.streams))

# rasters
r.lulc <- raster(paste0(dir.gis, "Navarro_Cannabis_LULC_30m.tif"))
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))
r.aquifers <- raster(paste0(dir.gis, "Navarro_Cannabis_GroundwaterBasins_30m.tif"))
r.wte <- raster(paste0(dir.gis, "Navarro_Cannabis_WTE_30m.tif"))

## calculate saturated thickness of alluvial materials are saturated or not
r.wtd <- r.dem.30m - r.wte
r.alluvial.sat.thickness <- r.dtb - r.wtd
r.alluvial.sat.thickness[r.alluvial.sat.thickness < r.dtb] <- 0
r.alluvial.sat.thickness[r.alluvial.sat.thickness > r.dtb] <- r.dtb[r.alluvial.sat.thickness > r.dtb]

## extract some potentially relevant data
sf.wel$elev_m <- raster::extract(r.dem.30m, sf.wel)  # elevation of that grid cell [m]
sf.wel$dtb_m <- raster::extract(r.dtb, sf.wel)       # depth to bedrock 
sf.wel[,c("lon", "lat")] <- st_coordinates(sf.wel)
sf.streams$elev_m <- raster::extract(r.dem.30m, sf.streams, fun='mean', na.rm=T)  # mean elevation of all grid cells stream touches [m]
sf.streams$dtb_m <- raster::extract(r.dtb, sf.streams, fun='mean', na.rm=T)       # mean DTB of all grid cells stream touches [m]
#sf.streams$satThick_m <- raster::extract(r.alluvial.sat.thickness, sf.streams, fun='mean', na.rm=T)       # mean DTB of all grid cells stream touches [m]
#sf.streams$satThick_m[is.na(sf.streams$satThick_m)] <- 0
sf.streams[,c("lon", "lat")] <- 
  sf.streams %>% 
  st_centroid() %>% 
  st_coordinates()

sum(is.na(sf.streams$elev_m))
sum(is.na(sf.streams$dtb_m))

## predict stream width based on drainage area
sf.streams$width_m <- WidthFromDA(DA=sf.streams$TtDASKM, w.min=1, w.max=100)

## distance from each well to each stream segment
dist_all <- 
  st_distance(x=sf.streams, y=sf.wel)

df.all <- data.frame(
  SegNum = rep(sf.streams$SegNum, times = dim(dist_all)[2]),
#  stream_lon = rep(sf.streams$lon, times = dim(dist_all)[2]),
#  stream_lat = rep(sf.streams$lat, times = dim(dist_all)[2]),
  stream_elev_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  stream_dtb_m = rep(sf.streams$dtb_m, times = dim(dist_all)[2]),
  stream_width_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  WellNum = rep(sf.wel$WellNum, each = dim(dist_all)[1]),
#  well_lon = rep(sf.wel$lon, times = dim(dist_all)[1]),
#  well_lat = rep(sf.wel$lat, times = dim(dist_all)[1]),
  well_elev_m = rep(sf.wel$elev_m, times = dim(dist_all)[1]),
  well_dtb_m = rep(sf.wel$dtb_m, times = dim(dist_all)[1]),
  dist_wellToStream_m = as.numeric(dist_all)
)

## save output
df.all %>% 
  format(digits = 3) %>% 
  write.csv(file.path("results", "Navarro_Cannabis_CalculateWellStreamPairs.csv"),
            row.names=F)
