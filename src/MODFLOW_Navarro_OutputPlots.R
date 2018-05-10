## MODFLOW_Navarro_OutputPlots.R
#' This script is intended to plot output data from MODFLOW:
#'  WTE = water table elevation [m]
#'  WTD = water table depth [m] = ELEV - WTE
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source(file.path("src", "paths+packages.R"))

## choose stream boundary condition and modflow version
stream_BC <- "SFR"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"

## load common data
# domain boundary shapefile
shp <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro")
shp.UTM <- spTransform(shp, crs.MODFLOW)

shp.adj <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

# streamlines
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv", stringsAsFactors=F)

## prep polygon boundaries for plots
df.basin <- tidy(shp.UTM)
df.basin.adj <- tidy(shp.adj.UTM)
df.riv <- tidy(shp.streams)

## load ibound raster
r.ibound <- raster(file.path("modflow", "input", "ibound.tif"))
DELR <- res(r.ibound)[1]
DELC <- res(r.ibound)[2]

# make data frame
df <- data.frame(rasterToPoints(r.ibound))
colnames(df) <- c("lon", "lat", "ibound")

# WTE and WTD ---------------------------------------

# load matrices
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wte.csv"), header=F))
m.wtd <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wtd.csv"), header=F))

# add to data frame
df$wte <- as.vector(t(m.wte))
df$wtd <- as.vector(t(m.wtd))

## make plots
p.wte <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=wte)) +
  geom_contour(data=df, aes(x=lon, y=lat, z=wte), color="white", binwidth=50) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  labs(subtitle=paste0("Range: ", round(min(m.wte, na.rm=T), 2), " to ", round(max(m.wte, na.rm=T), 2), " m")) +
  scale_fill_viridis(name="WTE [m]", na.value="white") +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

p.wtd <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=wtd)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  labs(subtitle=paste0("Range: ", round(min(m.wtd, na.rm=T), 2), " to ", round(max(m.wtd, na.rm=T), 2), " m")) +
  scale_fill_viridis(name="Water Table\nDepth [m]", na.value="white") +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

## save plots
# align
p1 <- ggplotGrob(p.wte+theme(legend.position="bottom"))
p2 <- ggplotGrob(p.wtd+theme(legend.position="bottom"))
p <- cbind(p1, p2, size="first")
p$heights <- unit.pmax(p1$heights, p2$heights)

# save
ggsave(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wte+wtd.png"),
       p, width=190, height=110, units="mm")

# SFR output --------------------------------------------------------------

if (modflow_v=="SFR"){
  ## load SFR data
  df.sfr.ReachData <- read.table(file.path("modflow", "input", "isfr_ReachData.txt"), stringsAsFactors=F, header=T)
  df.sfr <- read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "sfr.csv"), stringsAsFactors=F)
  
  # combine output with selected ReachData
  df.sfr <- left_join(df.sfr, df.sfr.ReachData[,c("OBJECTID", "SegNum", "StreamOrde", "SFR_NSEG", "SFR_IREACH")],
                      by=c("segment"="SFR_NSEG", "reach"="SFR_IREACH"))
  
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
  seg.plot <- 172
  df.seg.melt <-
    df.sfr[,c("SegNum", "reach", "Qin", "Qaquifer", "Qout", "stage", "head")] %>% 
    subset(SegNum==seg.plot) %>% 
    melt(., id=c("SegNum", "reach", "Qin"))
  
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
  ggsave(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "SFR_transects.png"),
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
  ggsave(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "SFR_maps.png"),
         p, width=190, height=100, units="mm")
  
}
