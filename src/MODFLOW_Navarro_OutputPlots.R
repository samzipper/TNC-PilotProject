## MODFLOW_Navarro_OutputPlots.R
#' This script is intended to plot output data from MODFLOW.
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source("src/paths+packages.R")

# X and Y resolution
DELR <- 250  # x spacing - column width
DELC <- 250  # y spacing - row height

## load common data
# domain boundary shapefile
shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

# river boundary shapefile
shp.riv.adj <- readOGR(dsn="data/NHDPlusV2/HYD", layer="NHDPlusV21_National_Seamless NHDFlowline_Network_WBDHU12_Navarro+Adjacent",
                       stringsAsFactors=F)
shp.riv.adj@data$StreamOrde <- as.numeric(shp.riv.adj@data$StreamOrde)
shp.riv.adj.streams.UTM <- subset(shp.riv.adj, StreamOrde >= 2)

## load ibound raster
r.ibound <- raster("modflow/input/ibound.tif")

# make data frame
df <- data.frame(rasterToPoints(r.ibound))
colnames(df) <- c("lon", "lat", "ibound")

# Head and WTD ---------------------------------------

# load matrices
m.head <- as.matrix(read.table("modflow/output/head.txt"))
m.wtd <- as.matrix(read.table("modflow/output/wtd.txt"))

# add to data frame
df$head <- as.vector(t(m.head))
df$wtd <- as.vector(t(m.wtd))

# Make plots --------------------------------------------------------------

## prep polygon boundaries
df.basin <- tidy(spTransform(readOGR(dsn="data/NHD/WBD", layer="WBDHU10_Navarro"), crs.MODFLOW))
df.riv <- tidy(shp.riv.adj.streams.UTM)

## head and WTD
p.head <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=head)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_viridis(name="Head [m]", na.value="white") +
  coord_equal() +
  theme_scz() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

p.wtd <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=wtd)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_viridis(name="Water Table\nDepth [m]", na.value="white") +
  coord_equal() +
  theme_scz() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

## save plots
# align
p1 <- ggplotGrob(p.head+theme(legend.position="bottom"))
p2 <- ggplotGrob(p.wtd+theme(legend.position="bottom"))
p <- cbind(p1, p2, size="first")
p$heights <- unit.pmax(p1$heights, p2$heights)

# save
ggsave("modflow/output/head+wtd.png",
       p, width=190, height=100, units="mm")
