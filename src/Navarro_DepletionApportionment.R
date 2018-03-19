## Navarro_DepletionApportionment.R
#' This script is intended to calculate depletion apportionment
#' among different well reaches for a bunch of stream reaches.
#' 
#' The well locations are created using the script MODFLOW_Navarro_InputPrepData.R
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source("src/paths+packages.R")

# Prep input data ---------------------------------------------------------

## load well locations
df.wel <- read.table("modflow/input/iwel.txt", sep=" ", header=T)

# make a spatial points data frame
xy <- df.wel[,c("lon", "lat")]
spdf.wel <- SpatialPointsDataFrame(coords = xy, data = df.wel,
                                   proj4string = CRS(crs.MODFLOW))


## load stream locations
# input files
shp.riv.adj <- readOGR(dsn="data/NHDPlusV2/HYD", layer="NHDPlusV21_National_Seamless NHDFlowline_Network_WBDHU12_Navarro+Adjacent",
                       stringsAsFactors=F)
shp.riv.adj@data$StreamOrde <- as.numeric(shp.riv.adj@data$StreamOrde)

# get rid of coastline (StreamOrde = -9)
shp.riv.adj.streams <- subset(shp.riv.adj, StreamOrde>0)

# eliminate tiny streams - some useful options might be:
shp.riv.adj.streams <- subset(shp.riv.adj.streams, StreamOrde >= 2)  # MAKE SURE THIS MATCHES MODFLOW_Navarro_InputPrepData.R !!!

# reproject to UTM
shp.riv.adj.streams.UTM <- spTransform(shp.riv.adj.streams, crs.MODFLOW)

# make a reach number column
shp.riv.adj.streams.UTM@data$ReachNum <- seq(1, dim(shp.riv.adj.streams.UTM@data)[1], 1)

# only retain some useful columns
shp.streams <- subset(shp.riv.adj.streams.UTM, select=c("ReachNum", "OBJECTID", "REACHCODE", "FromNode", "ToNode", "TerminalPa"))

# data frame for ggplots
df.streams <- tidy(shp.streams, id=ReachNum)

# combine into 1 line per stream feature (which is defined by TerminalPa)
shp.streams.dissolve <- gLineMerge(shp.streams)

# Calculate local area ----------------------------------------------------

# load raster
r.ibound <- raster("modflow/input/ibound.tif")

# get x/y/z data frame
df.ibound <- as.data.frame(rasterToPoints(r.ibound))
df.ibound <- subset(df.ibound, ibound != 0)
names(df.ibound) <- c("lon", "lat", "ibound")

# convert to spatial points
spdf.ibound <- SpatialPointsDataFrame(coords = df.ibound[,c("lon", "lat")], data = df.ibound,
                                      proj4string = CRS(crs.MODFLOW))

# calculate distance to nearest stream feature
spdf.ibound$distToStream <- gDistance(spdf.ibound, shp.streams.dissolve, byid=T)[1,]

# calculate local area (95th percentile) [m]
local.area.m <- quantile(spdf.ibound$distToStream, 0.95)

# double the local area because it seems too small...
local.area.m <- local.area.m*2

# # diagnostic plots to make sure things worked...
# ggplot() + 
#   geom_point(data=as.data.frame(spdf.ibound), aes(x=lon, y=lat, color=distToStream)) +
#   geom_path(data=df.streams, aes(x=long, y=lat, group=group))
# 
# ggplot() +
#   geom_histogram(data=as.data.frame(spdf.ibound), aes(x=distToStream), binwidth=100) +
#   geom_vline(xintercept=local.area.m, color="red")

# Convert stream lines to points ------------------------------------------

# define point spacing and figure out how many points to make
pt.spacing <- 10  # [m]
n.pts <- round(gLength(shp.streams)/pt.spacing)

# sample points
shp.streams.pts <- spsample(shp.streams, n=n.pts, type="regular")
df.streams.pts <- as.data.frame(shp.streams.pts)
colnames(df.streams.pts) <- c("lon", "lat")
  
# figure out what ReachNum each point corresponds to
shp.streams.buffer <- buffer(shp.streams, 0.1, dissolve=F)
int <- intersect(shp.streams.pts, shp.streams.buffer)
df.streams.pts <- cbind(df.streams.pts, int@data)

# Calculate depletion apportionment ---------------------------------------

## loop through wells
start.flag <- T
for (wel in df.wel$WellNum){
  ## for a given well, calculate the distance to each point
  i.wel <- which(df.wel$WellNum==wel)
  
  # get distance to all stream points
  df.wel.dist <- data.frame(ReachNum = df.streams.pts$ReachNum,
                            distToWell.m = round(sqrt((df.streams.pts$lon-df.wel$lon[i.wel])^2 + (df.streams.pts$lat-df.wel$lat[i.wel])^2), 2))
  
  # grab the lat/lon for these points
  df.wel.dist$lon <- df.streams.pts$lon
  df.wel.dist$lat <- df.streams.pts$lat
  
  # subset to only points within local area
  df.wel.dist.local <- subset(df.wel.dist, 
                              lon >= (df.wel$lon[i.wel]-local.area.m) &
                                lon <= (df.wel$lon[i.wel]+local.area.m) &
                                lat >= (df.wel$lat[i.wel]-local.area.m) &
                                lat <= (df.wel$lat[i.wel]+local.area.m))
  
  # calculate depletion apportionment fractions for different methods
  df.id <- apportion.inv.dist(reach=df.wel.dist.local$ReachNum, 
                              dist=df.wel.dist.local$distToWell.m, 
                              w=1, col.names=c("ReachNum", "f.InvDist"))
  
  df.idsq <- apportion.inv.dist(reach=df.wel.dist.local$ReachNum, 
                                dist=df.wel.dist.local$distToWell.m, 
                                w=2, col.names=c("ReachNum", "f.InvDistSq"))
  
  df.web <- apportion.web.dist(reach=df.wel.dist.local$ReachNum, 
                               dist=df.wel.dist.local$distToWell.m, 
                               w=1, col.names=c("ReachNum", "f.Web"))
  
  df.websq <- apportion.web.dist(reach=df.wel.dist.local$ReachNum, 
                                 dist=df.wel.dist.local$distToWell.m, 
                                 w=2, col.names=c("ReachNum", "f.WebSq"))
  
  df.tpoly <- apportion.tpoly(reach=df.wel.dist.local$ReachNum, 
                              dist=df.wel.dist.local$distToWell.m, 
                              lon=df.wel.dist.local$lon, 
                              lat=df.wel.dist.local$lat, 
                              wel.lon=df.wel$lon[i.wel],
                              wel.lat=df.wel$lat[i.wel],
                              wel.num=wel,
                              local.area.m=local.area.m,
                              coord.ref=CRS(crs.MODFLOW),
                              col.names=c("ReachNum", "f.TPoly"))
  
  # combine into single data frame
  df.apportion <- 
    full_join(df.id, df.idsq, by="ReachNum") %>% 
    full_join(x=., y=df.web, by="ReachNum") %>% 
    full_join(x=., y=df.websq, by="ReachNum") %>% 
    full_join(x=., y=df.tpoly, by="ReachNum")
  df.apportion$WellNum <- wel
  
  # add column for minimum distance to well from anywhere on this reach
  df.apportion <- 
    group_by(df.wel.dist, ReachNum) %>% 
    summarize(distToWell.min.m = min(distToWell.m)) %>% 
    left_join(x=df.apportion, y=., by=c("ReachNum"))
  
  if (start.flag){
    df.apportion.all <- df.apportion
    start.flag <- F
  } else {
    df.apportion.all <- rbind(df.apportion.all, df.apportion)
  }
  
  # status update
  print(paste0("Well ", wel, " of ", max(df.wel$WellNum), " complete @ ", Sys.time()))
  
}

# any NA in the 'f' columns should be converted to a 0
df.apportion.all[is.na(df.apportion.all)] <- 0

# Save output -------------------------------------------------------------

write.csv(df.apportion.all, 
          "results/Navarro_DepletionApportionment_AllMethods+Wells+Reaches.csv", 
          row.names=F)