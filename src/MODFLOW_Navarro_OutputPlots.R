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
shp <- readOGR(dsn="data/NHD/WBD", layer="WBDHU10_Navarro")
shp.UTM <- spTransform(shp, crs.MODFLOW)

shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

# river boundary shapefile
shp.riv.adj <- readOGR(dsn="data/NHDPlusV2/HYD", layer="NHDPlusV21_National_Seamless NHDFlowline_Network_WBDHU12_Navarro+Adjacent",
                       stringsAsFactors=F)
shp.riv.adj@data$StreamOrde <- as.numeric(shp.riv.adj@data$StreamOrde)
shp.riv.adj.streams <- subset(shp.riv.adj, StreamOrde >= 2)
shp.riv.adj.streams.UTM <- spTransform(shp.riv.adj.streams, crs.MODFLOW)

## prep polygon boundaries for plots
df.basin <- tidy(shp.UTM)
df.basin.adj <- tidy(shp.adj.UTM)
df.riv <- tidy(shp.riv.adj.streams.UTM)

## load ibound raster
r.ibound <- raster("modflow/input/ibound.tif")

# make data frame
df <- data.frame(rasterToPoints(r.ibound))
colnames(df) <- c("lon", "lat", "ibound")

# Head and WTD ---------------------------------------

# load matrices
m.head <- as.matrix(read.table("modflow/Navarro-SteadyState/output/head.txt"))
m.wtd <- as.matrix(read.table("modflow/Navarro-SteadyState/output/wtd.txt"))

# add to data frame
df$head <- as.vector(t(m.head))
df$wtd <- as.vector(t(m.wtd))

## make plots
p.head <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=head)) +
  geom_contour(data=df, aes(x=lon, y=lat, z=head), color="white", binwidth=50) +
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
ggsave("modflow/Navarro-SteadyState/output/head+wtd.png",
       p, width=190, height=100, units="mm")

# Pumping results ---------------------------------------------------------

# load output dataframe
df.wel <- read.table("modflow/Navarro-SteadyState-WithPumping/output/iwel_out.txt")

# column names
colnames(df.wel) <- c("row", "col", "ztop_m", "lay", "leakage", "depletion")

# Python indices are 0-based; add 1 for R
df.wel[,c("row", "col", "lay")] <- df.wel[,c("row", "col", "lay")]+1

# get x/y coordinates
df.wel$lon <- xFromCol(r.ibound, df.wel["col"])[,1]
df.wel$lat <- yFromRow(r.ibound, df.wel["row"])[,1]

# convert to spatialpoints
sp.wel <- SpatialPoints(df.wel[,c("lon", "lat")], proj4string=crs(r.ibound))
sp.wel <- SpatialPointsDataFrame(sp.wel, df.wel)

# build IDW model
gs <- gstat(formula=as.formula("depletion~1"), locations=sp.wel, nmax=8)

# interpolate to raster
r.wel.depletion <- interpolate(r.ibound, gs)

# mask with shapefile
r.wel.depletion.HUC <- mask(r.wel.depletion, shp.UTM)

# convert to data frame
df.wel.depletion <- as.data.frame(rasterToPoints(r.wel.depletion.HUC))
colnames(df.wel.depletion) <- c("lon", "lat", "depletion")

p.wel.depletion <- 
  ggplot() +
  geom_raster(data=df.wel.depletion, aes(x=lon, y=lat, fill=100*depletion)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_viridis(name="Depletion [%]", na.value="white", limits=c(0,100), breaks=seq(0,100,25)) +
  coord_equal() +
  theme_scz() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))
ggsave("modflow/Navarro-SteadyState-WithPumping/output/p.wel.depletion.png",
       p.wel.depletion, width=130, height=100, units="mm")
