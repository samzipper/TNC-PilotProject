## Navarro_DepletionApportionment_AdjacentOnly.R
#' This script is intended to calculate depletion apportionment
#' among different well reaches for a bunch of stream reaches
#' using only adjacent catchments to the well.
#' 
#' The well locations are created using the script MODFLOW_Navarro_InputPrepData.R
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

# Prep input data ---------------------------------------------------------

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

# make a spatial points data frame
xy <- df.wel[,c("lon", "lat")]
spdf.wel <- SpatialPointsDataFrame(coords = xy, data = df.wel,
                                   proj4string = CRS(crs.MODFLOW))


## load stream data - created in MODFLOW_Navarro_InputPrepData.R
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# data frame for ggplots
df.streams <- tidy(shp.streams, id=SegNum)

# combine into 1 line per stream feature (which is defined by TerminalPa)
shp.streams.dissolve <- gLineMerge(shp.streams)

# Calculate local area ----------------------------------------------------

### don't need to redo this commented out section unless ibound changes
# # load raster
# r.ibound <- raster(file.path("modflow", "input", "ibound.tif"))
# 
# # get x/y/z data frame
# df.ibound <- as.data.frame(rasterToPoints(r.ibound))
# df.ibound <- subset(df.ibound, ibound != 0)
# names(df.ibound) <- c("lon", "lat", "ibound")
# 
# # convert to spatial points
# spdf.ibound <- SpatialPointsDataFrame(coords = df.ibound[,c("lon", "lat")], data = df.ibound,
#                                       proj4string = CRS(crs.MODFLOW))
# 
# # calculate distance to nearest stream feature
# spdf.ibound$distToStream <- gDistance(spdf.ibound, shp.streams.dissolve, byid=T)[1,]
# 
# # calculate local area (95th percentile) [m]
# local.area.m <- quantile(spdf.ibound$distToStream, 0.95)
local.area.m <- 4616.931  # for 100 m resolution

# increase local area by 5x to include streams several catchments distant
local.area.m <- local.area.m*5

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
  
# figure out what SegNum each point corresponds to
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
  df.wel.dist <- data.frame(SegNum = df.streams.pts$SegNum,
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
  
  ## calculate depletion apportionment fractions for different methods
  # for now: we will use TPoly to define the adjacent catchments
  df.tpoly <- apportion.tpoly(reach=df.wel.dist.local$SegNum, 
                              dist=df.wel.dist.local$distToWell.m, 
                              lon=df.wel.dist.local$lon, 
                              lat=df.wel.dist.local$lat, 
                              wel.lon=df.wel$lon[i.wel],
                              wel.lat=df.wel$lat[i.wel],
                              wel.num=wel,
                              local.area.m=local.area.m,
                              coord.ref=CRS(crs.MODFLOW),
                              col.names=c("SegNum", "f.TPoly"))
  
  # create reach_dist data frame with adjacent catchments only
  SegNum.adj <- df.tpoly$SegNum
  reach_dist <- data.frame(reach = df.wel.dist.local$SegNum[df.wel.dist.local$SegNum %in% SegNum.adj], 
                           dist = df.wel.dist.local$distToWell.m[df.wel.dist.local$SegNum %in% SegNum.adj])
  
  df.id <- 
    apportion_inverse(reach_dist = reach_dist, w=1) %>% 
    set_colnames(c("SegNum", "f.InvDist"))
  
  df.idsq <- 
    apportion_inverse(reach_dist = reach_dist, w=2) %>% 
    set_colnames(c("SegNum", "f.InvDistSq"))
  
  df.web <-
    apportion_web(reach_dist = reach_dist, w=1) %>% 
    set_colnames(c("SegNum", "f.Web"))
  
  df.websq <-
    apportion_web(reach_dist = reach_dist, w=2) %>% 
    set_colnames(c("SegNum", "f.WebSq"))
  
  # combine into single data frame
  df.apportion <- 
    full_join(df.id, df.idsq, by="SegNum") %>% 
    full_join(x=., y=df.web, by="SegNum") %>% 
    full_join(x=., y=df.websq, by="SegNum") %>% 
    full_join(x=., y=df.tpoly, by="SegNum")
  df.apportion$WellNum <- wel
  
  # add column for minimum distance to well from anywhere on this reach
  df.apportion <- 
    group_by(df.wel.dist, SegNum) %>% 
    summarize(distToWell.min.m = min(distToWell.m)) %>% 
    left_join(x=df.apportion, y=., by=c("SegNum"))
  
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

# check: does it sum to 1 for each WellNum and method?
df.sum <-
  df.apportion.all %>% 
  group_by(WellNum) %>% 
  summarize(Qf.InvDist = sum(f.InvDist),
            Qf.InvDistSq = sum(f.InvDistSq),
            Qf.Web = sum(f.Web),
            Qf.WebSq = sum(f.WebSq),
            Qf.TPoly = sum(f.TPoly))

table(df.sum$Qf.InvDist)
table(df.sum$Qf.InvDistSq)
table(df.sum$Qf.Web)
table(df.sum$Qf.WebSq)
table(df.sum$Qf.TPoly)

# Save output -------------------------------------------------------------

# round to 5 significant digits to reduce file size
df.apportion.all$f.InvDist <- round(df.apportion.all$f.InvDist, 5)
df.apportion.all$f.InvDistSq <- round(df.apportion.all$f.InvDistSq, 5)
df.apportion.all$f.Web <- round(df.apportion.all$f.Web, 5)
df.apportion.all$f.WebSq <- round(df.apportion.all$f.WebSq, 5)
df.apportion.all$f.TPoly <- round(df.apportion.all$f.TPoly, 5)

write.csv(df.apportion.all, 
          file.path("results","Navarro_DepletionApportionment_AdjacentOnly_AllMethods+Wells+Reaches.csv"), 
          row.names=F)
