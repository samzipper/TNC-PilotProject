## MODFLOW_Navarro_OutputPlots.R
#' This script is intended to plot output data from MODFLOW.
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source("src/paths+packages.R")

## load common data
# domain boundary shapefile
shp <- readOGR(dsn="data/NHD/WBD", layer="WBDHU10_Navarro")
shp.UTM <- spTransform(shp, crs.MODFLOW)

shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

# streamlines
shp.streams <- readOGR(dsn="modflow/input", layer="iriv", stringsAsFactors=F)

## prep polygon boundaries for plots
df.basin <- tidy(shp.UTM)
df.basin.adj <- tidy(shp.adj.UTM)
df.riv <- tidy(shp.streams)

## load ibound raster
r.ibound <- raster("modflow/input/ibound.tif")
DELR <- res(r.ibound)[1]
DELC <- res(r.ibound)[2]

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
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

p.wtd <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=wtd)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_viridis(name="Water Table\nDepth [m]", na.value="white") +
  coord_equal() +
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

# SFR output --------------------------------------------------------------

## load SFR data
df.sfr <- read.csv("modflow/Navarro-SteadyState/output/sfr.csv", stringsAsFactors=F)

# get coordinates
df.sfr$cellNum <- cellFromRowCol(r.ibound, rownr=df.sfr$row, colnr=df.sfr$column)
df.sfr[,c("lon", "lat")] <- xyFromCell(r.ibound, df.sfr$cellNum)

# for plotting, set a lower threshold for color scaling
Q.min <- 1
df.sfr$Qout[df.sfr$Qout<Q.min] <- Q.min

# aquifer exchange
df.sfr$StreamType <- "No Exchange"
df.sfr$StreamType[df.sfr$Qaquifer<0] <- "Gaining"
df.sfr$StreamType[df.sfr$Qaquifer>0] <- "Losing"
df.sfr$StreamType <- factor(df.sfr$StreamType, levels=c("Losing", "No Exchange", "Gaining"))

df.sfr <- left_join(df.sfr, df)

## going down a segment
seg.plot <- 325
df.seg.melt <-
  df.sfr[,c("segment", "reach", "Qin", "Qaquifer", "Qout", "stage", "head")] %>% 
  subset(segment==seg.plot) %>% 
  melt(., id=c("segment", "reach", "Qin"))

p.reach.stage <-
  ggplot(subset(df.seg.melt, variable %in% c("stage", "head")), 
         aes(x=reach, y=value, color=variable)) +
  geom_line() +
  scale_x_continuous(name="Reach", expand=c(0,0)) +
  scale_y_continuous(name="Elevation [m]", expand=c(0,0)) +
  labs(title=paste0("Segment ", seg.plot)) +
  scale_color_manual(name=NULL, labels=c("stage"="Stream\nStage", "head"="Head"),
                     values=c("stage"=col.cat.red, "head"=col.cat.blu))

p.reach.Qaquifer <-
  ggplot(subset(df.seg.melt, variable %in% c("Qaquifer")), 
         aes(x=reach, y=value*1000/(DELR*DELC))) +
  geom_line(color=col.cat.blu) +
  scale_x_continuous(name="Reach", expand=c(0,0)) +
  scale_y_continuous(name="Flux into Aquifer [mm/d]", expand=c(0,0))

p.reach.Qout <-
  ggplot(subset(df.seg.melt, variable %in% c("Qout")), 
         aes(x=reach, y=value)) +
  geom_line(color=col.cat.blu) +
  scale_x_continuous(name="Reach", expand=c(0,0)) +
  scale_y_continuous(name="Stream Discharge [m3/d]", expand=c(0,0))

# save plots
# align
p1 <- ggplotGrob(p.reach.stage+theme(legend.position=c(0.01,0.01), legend.justification=c(0,0)))
p2 <- ggplotGrob(p.reach.Qaquifer)
p3 <- ggplotGrob(p.reach.Qout)
p <- rbind(p1, p2, p3, size="first")
p$widths <- unit.pmax(p1$widths, p2$widths, p3$widths)

# save
ggsave("modflow/Navarro-SteadyState/output/SFR_transects.png",
       p, width=100, height=190, units="mm")

## maps
p.sfr.Qout <-
  ggplot() + 
  geom_raster(data=df.sfr, aes(x=lon, y=lat, fill=Qout)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_viridis(name="Discharge [m3/d]", na.value="white", trans="log10", direction=-1) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

p.sfr.StreamType <-
  ggplot(df.sfr, aes(x=lon, y=lat, fill=StreamType)) +
  geom_raster() +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_manual(values=c("Losing"=col.cat.red, "No Exchange"=col.cat.yel, "Gaining"=col.cat.blu)) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

## save plots
# align
p1 <- ggplotGrob(p.sfr.Qout+theme(legend.position="bottom"))
p2 <- ggplotGrob(p.sfr.StreamType+theme(legend.position="bottom"))
p <- cbind(p1, p2, size="first")
p$heights <- unit.pmax(p1$heights, p2$heights)

# save
ggsave("modflow/Navarro-SteadyState/output/SFR_maps.png",
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
  theme(axis.text.y=element_text(angle=90, hjust=0.5))
ggsave("modflow/Navarro-SteadyState-WithPumping/output/p.wel.depletion.png",
       p.wel.depletion, width=130, height=100, units="mm")
