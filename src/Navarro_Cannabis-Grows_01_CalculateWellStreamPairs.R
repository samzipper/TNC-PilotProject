## Navarro_Cannabis-ExistingWells_01_CalculateWellStreamPairs+Depletion.R
#' This script is intended to calculate analytical model inputs for each well-stream combination:
#'   -d = distance to stream [m]
#'   -S = effective storage coefficient [-]
#'   -Tr = effective transmissivity [L2/T]

source(file.path("src", "paths+packages.R"))

## load data
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>% 
  sf::st_transform(crs.MODFLOW)

# grow locations shapefile
sf.grows <- 
  sf::st_read(file.path(dir.TNC, "nav_cannabis_ucbtnc", "nav_cannabis_ucbtnc.shp")) %>% 
  subset(year==16) %>%                  # 2016 data only
  st_zm(drop = TRUE, what = "ZM") %>%   # drop Z dimension from geometry
  sf::st_transform(crs.MODFLOW) %>% 
  dplyr::rename(GrowNum = FID_allgro)

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
sf.grows$elev_m <- raster::extract(r.dem.30m, sf.grows)  # elevation of that grid cell [m]
sf.grows$dtb_m <- raster::extract(r.dtb, sf.grows)       # depth to bedrock 
sf.grows[,c("lon", "lat")] <- sf::st_coordinates(sf.grows)
sf.streams$elev_m <- raster::extract(r.dem.30m, sf.streams, fun='mean', na.rm=T)  # mean elevation of all grid cells stream touches [m]
sf.streams$dtb_m <- raster::extract(r.dtb, sf.streams, fun='mean', na.rm=T)       # mean DTB of all grid cells stream touches [m]
sf.streams[,c("lon", "lat")] <- 
  sf.streams %>% 
  sf::st_centroid() %>% 
  sf::st_coordinates()

sum(is.na(sf.streams$elev_m))
sum(is.na(sf.streams$dtb_m))

## predict stream width based on drainage area
sf.streams$width_m <- WidthFromDA(DA=sf.streams$TtDASKM, w.min=1, w.max=100)

## distance from each well to each stream segment
dist_all <- 
  st_distance(x=sf.streams, y=sf.grows)

df.all <- data.frame(
  SegNum = rep(sf.streams$SegNum, times = dim(dist_all)[2]),
  stream_elev_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  stream_dtb_m = rep(sf.streams$dtb_m, times = dim(dist_all)[2]),
  stream_width_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  GrowNum = rep(sf.grows$GrowNum, each = dim(dist_all)[1]),
  well_elev_m = rep(sf.grows$elev_m, times = dim(dist_all)[1]),
  well_dtb_m = rep(sf.grows$dtb_m, times = dim(dist_all)[1]),
  dist_wellToStream_m = as.numeric(dist_all)
)

## save output
df.all %>% 
  format(digits = 3) %>% 
  write.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_CalculateWellStreamPairs.csv"),
            row.names=F)

sf.grows %>% 
  sf::st_write(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.shp"),
               delete_dsn=T, delete_layer=T)   # overwrite

## plot with everything
ggplot() +
  geom_sf(data=sf.basin) +
  geom_sf(data=sf.streams, color=col.cat.blu) +
  geom_sf(data=sf.grows, shape=21)
