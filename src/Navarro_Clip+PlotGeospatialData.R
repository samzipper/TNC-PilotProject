## Navarro_GlobalGeospatialData.R
#' This script is intended to extract data from a variety of global
#' datasets using a shapefile of the Navarro River basin boundary.

source("src/paths+packages.R")

## open basin boundary and river network shapefiles
shp <- readOGR(dsn="data/NHD/WBD", layer="WBDHU10_Navarro")
shp.riv <- readOGR(dsn="data/NHD/HYD", layer="NHDFlowline_HU10_Navarro")

shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent")
shp.adj.riv <- readOGR(dsn="data/NHD/HYD", layer="NHDFlowline_HU12_Navarro+Adjacent")

# reproject to WGS to match global datasets
shp.WGS <- spTransform(shp, crs.WGS)
shp.riv.WGS <- spTransform(shp.riv, crs.WGS)
shp.adj.WGS <- spTransform(shp.adj, crs.WGS)
shp.adj.riv.WGS <- spTransform(shp.adj.riv, crs.WGS)

## NED data
# load raw data
r.dem.NED <- raster(paste0(dir.dem.NED, "Navarro_NED10m.vrt"))

# reproject boundary
r.dem.NED.crop <- crop(r.dem.NED, extent(spTransform(shp.adj, crs(r.dem.NED))))

# mask
r.dem.NED.mask <- mask(r.dem.NED.crop, spTransform(shp.adj, crs(r.dem.NED.crop)), 
                    filename=paste0(dir.gis, "Navarro_DEM_m_NED10m.tif"), datatype="FLT4S", overwrite=T)

## NLCD data: impervious, canopy, and lulc
# load raw data
r.can <- raster(paste0(dir.nlcd, "NLCD2011_CAN_California/NLCD2011_CAN_California.tif"))
r.imp <- raster(paste0(dir.nlcd, "NLCD2011_IMP_California/NLCD2011_IMP_California.tif"))
r.lulc <- raster(paste0(dir.nlcd, "NLCD2011_LC_California/NLCD2011_LC_California.tif"))

# crop to basin boundary
r.can.crop <- crop(r.can, extent(spTransform(shp.adj, crs(r.can))))
r.imp.crop <- crop(r.imp, extent(spTransform(shp.adj, crs(r.can))))
r.lulc.crop <- crop(r.lulc, extent(spTransform(shp.adj, crs(r.can))))

# mask
r.can.mask <- mask(r.can.crop, spTransform(shp.adj, crs(r.can)), 
                   filename=paste0(dir.gis, "Navarro_NLCD2011_CanopyCover_30m.tif"), datatype="INT2S", overwrite=T)
r.imp.mask <- mask(r.imp.crop, spTransform(shp.adj, crs(r.can)), 
                   filename=paste0(dir.gis, "Navarro_NLCD2011_ImpervCover_30m.tif"), datatype="INT2S", overwrite=T)
r.lulc.mask <- mask(r.lulc.crop, spTransform(shp.adj, crs(r.can)), 
                    filename=paste0(dir.gis, "Navarro_NLCD2011_LULC_30m.tif"), datatype="INT2S", overwrite=T)

## HydroSHEDS data: dem and slope
# load raw data
r.dem.hysheds <- raster(path.dem.hysheds)
r.slope.hysheds <- raster(path.slope.hysheds)

# crop to basin boundary
r.dem.hysheds.s <- crop.mask(r.dem.hysheds, shp.adj.WGS)
r.slope.hysheds.s <- crop.mask(r.slope.hysheds, shp.adj.WGS)

# write rasters
writeRaster(r.dem.hysheds.s, paste0(dir.gis, "Navarro_DEM_m_HydroSHEDS_15s.tif"), datatype="INT2S", overwrite=T)
writeRaster(r.slope.hysheds.s, paste0(dir.gis, "Navarro_slope_deg_HydroSHEDS_15s.tif"), datatype="FLT4S", overwrite=T)

## soilgrids data: % sand (top 10 cm)
# load raw data
r.dtb <- raster(paste0(dir.SoilGrids, "1original/SoilGrids1km/BDTICM_M_1km_ll.tif"))
r.sand <- raster(paste0(dir.SoilGrids, "2derived/SoilGrids1km/SNDPPT_M_meanTop015cm_1km_ll.tif"))

# crop to basin boundary
r.dtb.s <- crop.mask(r.dtb, shp.adj.WGS)
r.sand.s <- crop.mask(r.sand, shp.adj.WGS)

# convert dtb from cm to m
r.dtb.s <- r.dtb.s/100

# write rasters
writeRaster(r.dtb.s, paste0(dir.gis, "Navarro_DTB_m_SoilGrids1km.tif"), datatype="FLT4S", overwrite=T)
writeRaster(r.sand.s, paste0(dir.gis, "Navarro_sand_prc_top15cm_SoilGrids1km.tif"), datatype="FLT4S", overwrite=T)

## GLHYMPS data: porosity and permeability
# load raw data
r.porosity <- raster(paste0(dir.GLHYMPS, "na_porosity_x100_EPSG4326.tif"))
r.logK <- raster(paste0(dir.GLHYMPS, "na_logK_Ice_x100_EPSG4326.tif"))

# crop to basin boundary
r.porosity.s <- crop.mask(r.porosity, shp.adj.WGS)
r.logK.s <- crop.mask(r.logK, shp.adj.WGS)

# anything that has a 0, set to NaN
r.porosity.s[r.porosity.s==0] <- NaN
r.logK.s[r.logK.s > -1000] <- NaN

# get rid of x100
r.porosity.s <- r.porosity.s/100
r.logK.s <- r.logK.s/100

# write rasters
writeRaster(r.porosity.s, paste0(dir.gis, "Navarro_porosity_GLHYMPS.tif"), datatype="FLT4S", overwrite=T)
writeRaster(r.logK.s, paste0(dir.gis, "Navarro_logK_GLHYMPS.tif"), datatype="FLT4S", overwrite=T)

## WTD from Fan et al. (2013)
# load raw data
r.wtd <- raster(path.wtd)

# crop to basin boundary
r.wtd.s <- crop.mask(r.wtd, shp.adj.WGS)

# write rasters
writeRaster(r.wtd.s, paste0(dir.gis, "Navarro_WTD_m_FanEtAl2013.tif"), datatype="FLT4S", overwrite=T)

## make plots
# tidy data
df.shp <- tidy(shp.WGS)
df.riv <- tidy(shp.riv.WGS)

df.nlcd.shp <- tidy(spTransform(shp.adj, crs(r.can.mask)))
df.nlcd.riv <- tidy(spTransform(shp.adj.riv, crs(r.can.mask)))

df.ned.shp <- tidy(spTransform(shp.adj, crs(r.dem.NED.mask)))
df.ned.riv <- tidy(spTransform(shp.adj.riv, crs(r.dem.NED.mask)))

df.dem <- as.data.frame(rasterToPoints(r.dem.hysheds.s))
colnames(df.dem) <- c("lon", "lat", "dem")

df.slope <- as.data.frame(rasterToPoints(r.slope.hysheds.s))
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

df.nlcd.can <- as.data.frame(rasterToPoints(r.can.mask))
colnames(df.nlcd.can) <- c("lon", "lat", "canopy_prc")

df.nlcd.imp <- as.data.frame(rasterToPoints(r.imp.mask))
colnames(df.nlcd.imp) <- c("lon", "lat", "imp_prc")

df.nlcd.lulc <- as.data.frame(rasterToPoints(r.lulc.mask))
colnames(df.nlcd.lulc) <- c("lon", "lat", "lulc")

df.ned.dem <- as.data.frame(rasterToPoints(r.dem.NED.mask))
colnames(df.ned.dem) <- c("lon", "lat", "dem")

# make plots
p.dem <- 
  ggplot() +
  geom_raster(data=df.dem, aes(x=lon, y=lat, fill=dem)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Elevation [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_DEM_m_HydroSHEDS_15s.png"),
       p.dem, width=6.25, height=5.25, units="in")

p.slope <- 
  ggplot() +
  geom_raster(data=df.slope, aes(x=lon, y=lat, fill=slope)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Slope [deg]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_slope_deg_HydroSHEDS_15s.png"),
       p.slope, width=6.25, height=5.25, units="in")

p.dtb <- 
  ggplot() +
  geom_raster(data=df.dtb, aes(x=lon, y=lat, fill=dtb)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Depth to Bedrock [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_DTB_m_SoilGrids1km.png"),
       p.dtb, width=6.25, height=5.25, units="in")

p.sand <- 
  ggplot() +
  geom_raster(data=df.sand, aes(x=lon, y=lat, fill=sand)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Sand, Top 15 cm [%]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_sand_prc_top15cm_SoilGrids1km.png"),
       p.sand, width=6.25, height=5.25, units="in")

p.porosity <- 
  ggplot() +
  geom_raster(data=df.porosity, aes(x=lon, y=lat, fill=porosity)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Porosity [-]", direction=-1) +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_porosity_GLHYMPS.png"),
       p.porosity, width=6.25, height=5.25, units="in")

p.logK <- 
  ggplot() +
  geom_raster(data=df.logK, aes(x=lon, y=lat, fill=logK)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="logK [m2]", direction=-1) +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_logK_GLHYMPS.png"),
       p.logK, width=6.25, height=5.25, units="in")

p.wtd <- 
  ggplot() +
  geom_raster(data=df.wtd, aes(x=lon, y=lat, fill=wtd)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Water Table Depth [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_WTD_m_FanEtAl2013.png"),
       p.wtd, width=6.25, height=5.25, units="in")

## NLCD plots
p.nlcd.can <- 
  ggplot() +
  geom_raster(data=df.nlcd.can, aes(x=lon, y=lat, fill=canopy_prc)) +
  geom_path(data=df.nlcd.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.nlcd.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0)) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0)) +
  scale_fill_viridis(name="Canopy\nCover [%]", direction=-1, limits=c(0,100), breaks=seq(0,100,25)) +
  theme_scz() +
  theme(panel.border=element_blank())
ggsave(paste0(dir.gis, "Navarro_NLCD2011_CanopyCover_30m.png"),
       p.nlcd.can, width=5.5, height=6, units="in")

p.nlcd.imp <- 
  ggplot() +
  geom_raster(data=df.nlcd.imp, aes(x=lon, y=lat, fill=imp_prc)) +
  geom_path(data=df.nlcd.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.nlcd.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0)) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0)) +
  scale_fill_viridis(name="Impervious\nCover [%]", direction=1, limits=c(0,100), breaks=seq(0,100,25)) +
  theme_scz() +
  theme(panel.border=element_blank())
ggsave(paste0(dir.gis, "Navarro_NLCD2011_ImpervCover_30m.png"),
       p.nlcd.imp, width=5.5, height=6, units="in")

p.nlcd.lulc <- 
  ggplot() +
  geom_raster(data=df.nlcd.lulc, aes(x=lon, y=lat, fill=factor(lulc))) +
  geom_path(data=df.nlcd.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.nlcd.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0)) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0)) +
  scale_fill_manual(name="LULC, NLCD", values=pal.NLCD, labels=labels.NLCD) +
  theme_scz() +
  theme(panel.border=element_blank())
ggsave(paste0(dir.gis, "Navarro_NLCD2011_LULC_30m.png"),
       p.nlcd.lulc, width=6, height=6, units="in")

## NED DEM
p.ned.dem <- 
  ggplot() +
  geom_raster(data=df.ned.dem, aes(x=lon, y=lat, fill=dem)) +
  geom_path(data=df.ned.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.ned.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_equal() +
  labs(title=paste0("Navarro River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Elevation [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave(paste0(dir.gis, "Navarro_DEM_m_NED10m.png"),
       p.ned.dem, width=6.25, height=5.25, units="in")
