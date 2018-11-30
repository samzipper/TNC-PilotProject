## Navarro_Cannabis_PrepGeospatialData.R
# This script is intended to load, plot, and save geospatial data needed for analysis of 
# pumping impacts on habitat. Everything should be at 30 m resolution and in the EPSG:26910 projection.

source("src/paths+packages.R")

## load data and reproject - goal is everything in 30 m resolution
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("modflow", "input", "iriv.shp"), stringsAsFactors=F) %>% 
  subset(lnLngt_ >= 100)
st_write(sf.streams, file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"))

# domain boundary shapefile
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# domain boundary including adjacent
shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent")

# relevant rasters (from script Navarro_Clip+PlotGeospatialData.R)
# # set up r.lulc to use as a boundary/mask
# r.lulc <- 
#   paste0(dir.nlcd, "NLCD2011_LC_California/NLCD2011_LC_California.tif") %>% 
#   raster() 
# r.lulc <- 
#   crop(r.lulc, extent(spTransform(shp.adj, crs(r.lulc)))) %>% 
#   projectRaster(from=., crs=crs(crs.MODFLOW), res=30, method="ngb") %>% 
#   mask(., spTransform(shp.adj, crs.MODFLOW),
#        filename=paste0(dir.gis, "Navarro_Habitat_LULC_30m.tif"), datatype="INT2S", overwrite=T)

# read in r.lulc raster
r.lulc <- raster(paste0(dir.gis, "Navarro_Cannabis_LULC_30m.tif"))

# # set up DEM
# r.dem <-
#   paste0(dir.gis, "Navarro_DEM_m_NED10m.tif") %>%
#   raster()
# r.dem.30m <- raster::aggregate(r.dem, fact=3)
# r.dem.30m[r.dem.30m < -20] <- NA
# r.dem.30m <- projectRaster(from = r.dem.30m,
#                            to = r.lulc,
#                            method="bilinear") %>%
#   mask(., spTransform(shp.adj, crs.MODFLOW),
#        filename=paste0(dir.gis, "Navarro_Habitat_DEM_30m.tif"), datatype="FLT4S", overwrite=T)

# read in DEM raster
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))

# # depth to bedrock
# r.dtb <-
#   paste0(dir.gis, "Navarro_DTB_m_SoilGrids1km.tif") %>%
#   raster() %>%
#   projectRaster(crs=crs.MODFLOW) %>%
#   resample(., r.lulc, method="bilinear") %>%
#   mask(., spTransform(shp.adj, crs.MODFLOW),
#        filename=paste0(dir.gis, "Navarro_Habitat_DTB_30m.tif"), overwrite=T)

# read in DTB raster
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))

# # alluvial aquifer extents - convert from shapefiles to rasters
# # this will be used to define spatial distribution of K and Sy
# shp.aquifers <-
#   readOGR(dsn="data/GroundwaterBasins_Bulletin118", layer="Navarro_GroundwaterBasins") %>%
#   spTransform(., crs.MODFLOW)
# r.aquifers <-
#   rasterize(shp.aquifers, r.lulc, field='OBJECTID', background=0) %>%
#   mask(., spTransform(shp.adj, crs.MODFLOW))
# writeRaster(r.aquifers, filename=paste0(dir.gis, "Navarro_Habitat_GroundwaterBasins_30m.tif"),
#             datatype="INT1U", overwrite=T)

# read in aquifers raster
r.aquifers <- raster(paste0(dir.gis, "Navarro_Cannabis_GroundwaterBasins_30m.tif"))
