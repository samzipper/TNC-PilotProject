## Navarro_Cannabis_04_DepletionApportionment.R
#' This script is intended to calculate depletion apportionment for all stream segments.
#' 
#' It requires output from Navarro_Cannabis_03_CalculateWellStreamPairs.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

### what dates do you want output for? (units: number of days since start of pumping)
# convert years and dates to DOY since pumping started
yrs.model <- c(1, 10, 50) 
DOYs.model <- yday(c("2017-09-15"))
DOYs.all <- rep(DOYs.model, times=length(yrs.model)) + rep(365*(yrs.model-1), each=length(DOYs.model))

### load and pre-process data
# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>%
  sf::st_transform(crs.MODFLOW)

# load stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

# synthetic pumping wells
sf.wel <- 
  file.path("modflow", "input", "iwel.txt") %>% 
  read.table(header=T) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(sf.streams))

## well-stream geometry output from Navarro_Cannabis_03_CalculateWellStreamPairs.R
df.all <- 
  file.path("results", "Navarro_Cannabis_03_CalculateWellStreamPairs.csv") %>% 
  read.csv(stringsAsFactors=F)

# # compare different Transmissivity approaches and pick your poison
# ggplot(df.all, aes(x=Tr_bulk_m2d/Tr_well_m2d)) +
#   geom_histogram() + 
#   geom_vline(xintercept=1, color="red") +
#   scale_x_log10()
# 
# ggplot(df.all, aes(x=Tr_bulk_m2d/Tr_stream_m2d)) +
#   geom_histogram() + 
#   geom_vline(xintercept=1, color="red") +
#   scale_x_log10()

# subset to method of choice
df.all <-
  df.all %>% 
  dplyr::select(SegNum, WellNum, Tr_bulk_m2d, S_bulk, lmda_m2d) %>% 
  transform(SegNum_WellNum = paste0(SegNum, "_", WellNum))

## subdivide sf.streams into points for web squared calculation
# define point spacing and figure out how many points to make
pt.spacing <- 20  # [m]
shp.streams <- as(sf.streams, "Spatial")
n.pts <- round(gLength(shp.streams)/pt.spacing)
set.seed(1)
shp.streams.pts <- spsample(shp.streams, n=n.pts, type="regular")
df.streams.pts <- as.data.frame(shp.streams.pts)
colnames(df.streams.pts) <- c("lon", "lat")
shp.streams.buffer <- buffer(shp.streams, 0.1, dissolve=F) # figure out what SegNum each point corresponds to
int <- intersect(shp.streams.pts, shp.streams.buffer)
df.streams.pts <- cbind(df.streams.pts, int@data)
sf.streams.pts <- sf::st_as_sf(df.streams.pts, coords=c("lon", "lat"), crs=st_crs(sf.streams))

## find distance from each well point to each stream segment
pt_coords <- st_coordinates(sf.streams.pts)
dist_all_pts <- 
  st_distance(x=sf.streams.pts, y=sf.wel) %>% 
  as.numeric() %>% 
  data.frame(
    SegNum = rep(sf.streams.pts$SegNum, times = dim(sf.wel)[1]),
    WellNum = rep(sf.wel$WellNum, each = dim(sf.streams.pts)[1]),
    dist_wellToStream_m = .,
    lon = rep(pt_coords[,"X"], times = dim(sf.wel)[1]),
    lat = rep(pt_coords[,"Y"], times = dim(sf.wel)[1])
  ) %>% 
  transform(SegNum_WellNum = paste0(SegNum, "_", WellNum)) %>% 
  subset(SegNum_WellNum %in% df.all$SegNum_WellNum) %>% 
  left_join(df.all, by=c("SegNum", "WellNum", "SegNum_WellNum"))

### depletion calculations
## loop through times [d] for depletion calculations
min_frac <- 0.01  # minimum depletion considered

w.start.flag <- T
counter <- 0
for (w in unique(df.all$WellNum)){
  
  # get lat/lon of well
  wel_coord <- 
    sf.wel %>% 
    subset(WellNum == w) %>% 
    st_coordinates()
  
  ## retain only relevant points
  # all points
  dist_w_pts <- subset(dist_all_pts, WellNum==w)
  
  # closest point only
  dist_w_pts_closest <- 
    dist_w_pts %>% 
    group_by(SegNum) %>% 
    filter(dist_wellToStream_m == min(dist_wellToStream_m))
  
  # first: Theissen polygons to figure out adjacent catchments
  df.apportion.t <- 
    dist_w_pts_closest %>% 
    apportion_polygon(., 
                      wel_lon = wel_coord[1,"X"], 
                      wel_lat = wel_coord[1,"Y"],
                      crs = CRS(st_crs(sf.streams)[["proj4string"]]),
                      reach_name = "SegNum",
                      dist_name = "dist_wellToStream_m") %>% 
    set_colnames(c("SegNum", "frac_depletion"))
  
  max_dist_prev <- max(c(500, min(dist_w_pts$dist_wellToStream_m)))
  for (time_days in DOYs.all){
    
    # find maximum distance, based on maximum observed S, Tr, lmda (inclusive estimate)
    max_dist <- depletion_max_distance(Qf_thres = min_frac,
                                       d_interval = 250,
                                       d_min = max_dist_prev,
                                       d_max = max(dist_w_pts$dist_wellToStream_m),
                                       method="hunt",
                                       t = time_days,
                                       S = min(dist_w_pts$S_bulk),
                                       Tr = max(dist_w_pts$Tr_bulk_m2d),
                                       lmda = max(dist_w_pts$lmda_m2d))
    
    # second: use web^2 to apportion to any stream segment that is within max_dist
    #         OR that is adjacent to well based on Theissen Polygon
    df.apportion.w.t <- 
      dist_w_pts %>% 
      subset((dist_wellToStream_m <= max_dist) | (SegNum %in% df.apportion.t$SegNum)) %>% 
      apportion_web(reach_dist = ., 
                    w = 2, 
                    min_frac = min_frac,
                    reach_name = "SegNum",
                    dist_name = "dist_wellToStream_m") %>% 
      set_colnames(c("SegNum", "frac_depletion")) %>% 
      transform(WellNum = w,
                time_days = time_days) %>% 
      left_join(dist_w_pts_closest, by=c("SegNum", "WellNum"))
    
    if (w.start.flag){
      df.out <- df.apportion.w.t[ , !(names(df.apportion.w.t) %in% c("lon", "lat", "SegNum_WellNum"))]
      w.start.flag <- F
    } else {
      df.out <- rbind(df.out, df.apportion.w.t[ , !(names(df.apportion.w.t) %in% c("lon", "lat", "SegNum_WellNum"))])
    }
    
    # update max_dist starting value
    max_dist_prev <- max_dist
    
  }
  
  # status update
  counter <- counter+1
  print(paste0("Well ", counter, " of ", length(unique(df.all$WellNum)), " depletion apportionment complete, ", Sys.time()))
}

# save output
df.out %>% 
  dfDigits(x=., digits=3) %>% 
  write.csv(file.path("results", "Navarro_Cannabis_04_DepletionApportionment.csv"),
            row.names=F)
