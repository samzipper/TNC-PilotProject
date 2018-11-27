## Figure_GeospatialData.R
#' This script is intended to make a plot of the geospatial data used in the analytical depletion
#' functions, which are created using the script Navarro_Cannabis_PrepGeospatialData.R

source(file.path("src", "paths+packages.R"))

## load data
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  subset(TermnlP == outlet.TerminalPa)

# domain boundary shapefile
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# rasters
r.lulc <- raster(paste0(dir.gis, "Navarro_Cannabis_LULC_30m.tif"))
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))
r.aquifers <- raster(paste0(dir.gis, "Navarro_Cannabis_GroundwaterBasins_30m.tif"))

## convert raster to points
df.dem <- 
  r.dem.30m %>% 
  rasterToPoints() %>% 
  as.data.frame() %>% 
  set_colnames(c("lon", "lat", "dem"))

df.dtb <- 
  r.dtb %>% 
  rasterToPoints() %>% 
  as.data.frame() %>% 
  set_colnames(c("lon", "lat", "dtb"))

df.aquifers <- 
  r.aquifers %>% 
  rasterToPoints() %>% 
  as.data.frame() %>% 
  set_colnames(c("lon", "lat", "aquifer"))

df.all <- 
  full_join(df.dem, df.dtb, by=c("lon", "lat")) %>% 
  full_join(df.aquifers, by=c("lon", "lat"))

# get rid of locations that have NA in any dataset
# this is primarily the coastal plain in aquifer and dem
df.all <- df.all[complete.cases(df.all), ]

p.elev <-
  ggplot() +
  geom_raster(data=df.all, aes(x=lon, y=lat, fill=dem)) +
  geom_sf(data=sf.basin, color=col.cat.red, fill=NA) +
  geom_sf(data=sf.streams, color="white") +
  scale_fill_viridis(name="Land Surface\nElevation [m]") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank())

p.DTB <- 
  ggplot() +
  geom_raster(data=df.all, aes(x=lon, y=lat, fill=dtb)) +
  geom_sf(data=sf.basin, color=col.cat.red, fill=NA) +
  geom_sf(data=sf.streams, color="white") +
  scale_fill_viridis(name="Depth To\nBedrock [m]", direction=-1) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank())

p1 <- ggplotGrob(p.elev)
p2 <- ggplotGrob(p.DTB)
p <- cbind(p1, p2, size="first")
p$heights <- unit.pmax(p1$heights, p1$heights)

ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_GeospatialData.pdf"),
       p, width=190, height=110, units="mm", device=cairo_pdf)

ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_GeospatialData.png"),
       p, width=190, height=110, units="mm")
