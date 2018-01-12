## Navarro_GlobalGeospatialData.R
#' This script is intended to extract data from a variety of global
#' datasets using a shapefile of the Navarro River basin boundary.

source("src/paths+packages.R")

## open basin boundary and river network shapefiles
shp <- readOGR(dsn="data/HydroBASINS", layer="hybas_navarro_lev07")
shp.riv <- readOGR(dsn="data", layer="HydroSHEDS_Navarro_RiverNetwork_15s")

## HydroSHEDS data: dem and slope
# load raw data
r.dem <- raster(path.dem)
r.slope <- raster(path.slope)

# crop to basin boundary
r.dem.s <- crop.mask(r.dem, shp)
r.slope.s <- crop.mask(r.slope, shp)

# write rasters
writeRaster(r.dem.s, paste0(dir.gis, "DEM_m_HydroSHEDS_15s.tif"), datatype="INT2S", overwrite=T)
writeRaster(r.slope.s, paste0(dir.gis, "slope_deg_HydroSHEDS_15s.tif"), datatype="FLT4S", overwrite=T)

## soilgrids data: % sand (top 10 cm)
# load raw data
r.dtb <- raster(paste0(dir.SoilGrids, "1original/SoilGrids1km/BDTICM_M_1km_ll.tif"))
r.sand <- raster(paste0(dir.SoilGrids, "2derived/SoilGrids1km/SNDPPT_M_meanTop015cm_1km_ll.tif"))

# crop to basin boundary
r.dtb.s <- crop.mask(r.dtb, shp)
r.sand.s <- crop.mask(r.sand, shp)

# convert dtb from cm to m
r.dtb.s <- r.dtb.s/100

# write rasters
writeRaster(r.dtb.s, paste0(dir.gis, "DTB_m_SoilGrids1km.tif"), datatype="FLT4S", overwrite=T)
writeRaster(r.sand.s, paste0(dir.gis, "sand_prc_top15cm_SoilGrids1km.tif"), datatype="FLT4S", overwrite=T)

## GLHYMPS data: porosity and permeability
# load raw data
r.porosity <- raster(paste0(dir.GLHYMPS, "na_porosity_x100_EPSG4326.tif"))
r.logK <- raster(paste0(dir.GLHYMPS, "na_logK_Ice_x100_EPSG4326.tif"))

# crop to basin boundary
r.porosity.s <- crop.mask(r.porosity, shp)
r.logK.s <- crop.mask(r.logK, shp)

# anything that has a 0, set to NaN
r.porosity.s[r.porosity.s==0] <- NaN
r.logK.s[r.logK.s > -1000] <- NaN

# get rid of x100
r.porosity.s <- r.porosity.s/100
r.logK.s <- r.logK.s/100

# write rasters
writeRaster(r.porosity.s, paste0(dir.gis, "porosity_GLHYMPS.tif"), datatype="FLT4S", overwrite=T)
writeRaster(r.logK.s, paste0(dir.gis, "logK_GLHYMPS.tif"), datatype="FLT4S", overwrite=T)

## WTD from Fan et al. (2013)
# load raw data
r.wtd <- raster(path.wtd)

# crop to basin boundary
r.wtd.s <- crop.mask(r.wtd, shp)

# write rasters
writeRaster(r.wtd.s, paste0(dir.gis, "WTD_m_FanEtAl2013.tif"), datatype="FLT4S", overwrite=T)

## make plots
# tidy data
df.shp <- tidy(shp)
df.riv <- tidy(shp.riv)

df.dem <- as.data.frame(rasterToPoints(r.dem.s))
colnames(df.dem) <- c("lon", "lat", "dem")

df.slope <- as.data.frame(rasterToPoints(r.slope.s))
colnames(df.slope) <- c("lon", "lat", "slope")

df.dtb <- as.data.frame(rasterToPoints(r.dtb.s))
colnames(df.dtb) <- c("lon", "lat", "dtb")

df.sand <- as.data.frame(rasterToPoints(r.sand.s))
colnames(df.sand) <- c("lon", "lat", "sand")

df.porosity <- as.data.frame(rasterToPoints(r.porosity.s))
colnames(df.porosity) <- c("lon", "lat", "porosity")

df.logK <- as.data.frame(rasterToPoints(r.logK.s))
colnames(df.logK) <- c("lon", "lat", "logK")

df.wtd <- as.data.frame(rasterToPoints(r.wtd.s))
colnames(df.wtd) <- c("lon", "lat", "wtd")

# make plots
p.dem <- 
  ggplot() +
  geom_raster(data=df.dem, aes(x=lon, y=lat, fill=dem)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("USGS ", station, ": ", station.name)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Elevation [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "DEM_m_HydroSHEDS_15s.png"),
       p.dem, width=6.25, height=5.5, units="in")

p.slope <- 
  ggplot() +
  geom_raster(data=df.slope, aes(x=lon, y=lat, fill=slope)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("USGS ", station, ": ", station.name)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Slope [deg]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "slope_deg_HydroSHEDS_15s.png"),
       p.slope, width=6.25, height=5.5, units="in")

p.dtb <- 
  ggplot() +
  geom_raster(data=df.dtb, aes(x=lon, y=lat, fill=dtb)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("USGS ", station, ": ", station.name)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Depth to Bedrock [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "DTB_m_SoilGrids1km.png"),
       p.dtb, width=6.25, height=5.5, units="in")

p.sand <- 
  ggplot() +
  geom_raster(data=df.sand, aes(x=lon, y=lat, fill=sand)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("USGS ", station, ": ", station.name)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Sand, Top 15 cm [%]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "sand_prc_top15cm_SoilGrids1km.png"),
       p.sand, width=6.25, height=5.5, units="in")

p.porosity <- 
  ggplot() +
  geom_raster(data=df.porosity, aes(x=lon, y=lat, fill=porosity)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("USGS ", station, ": ", station.name)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Porosity [-]", direction=-1) +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "porosity_GLHYMPS.png"),
       p.porosity, width=6.25, height=5.5, units="in")

p.logK <- 
  ggplot() +
  geom_raster(data=df.logK, aes(x=lon, y=lat, fill=logK)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("USGS ", station, ": ", station.name)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="logK [m2]", direction=-1) +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "logK_GLHYMPS.png"),
       p.logK, width=6.25, height=5.5, units="in")

p.wtd <- 
  ggplot() +
  geom_raster(data=df.wtd, aes(x=lon, y=lat, fill=wtd)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("USGS ", station, ": ", station.name)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Water Table Depth [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "WTD_m_FanEtAl2013.png"),
       p.wtd, width=6.25, height=5.5, units="in")
