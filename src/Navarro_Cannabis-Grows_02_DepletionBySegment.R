## Navarro_Cannabis-Grows_02_DepletionBySegment.R
#' This script is intended to calculate depletion for all stream segments.
#' 
#' It requires output from Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

### what dates do you want output for? (units: number of days since start of pumping)
# convert years and dates to DOY since pumping started
yrs.model <- c(1, 2, 3, 4, 5, 10, 15, 20, 30, 50) 
DOYs.model <- yday(c("2017-01-15", "2017-02-14", "2017-03-15", "2017-04-15", "2017-05-15", "2017-06-15",
                     "2017-07-15", "2017-08-15", "2017-09-15", "2017-10-15", "2017-11-15", "2017-12-15"))
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

# grow locations from TNC (filtered and transformed in Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R)
sf.grows <-
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.shp"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

## load pumping rates - this is proprietary from TNC, cannot be shared (until Wilson et al paper published)
df.pump <- read.csv(file.path(dir.TNC, "CannabisMonthlyWaterUse_WilsonEtAl.csv"), stringsAsFactors=F)
df.pump$Month <- factor(df.pump$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump$Setting <- factor(df.pump$Setting, levels=c("Outdoor", "Greenhouse"))
df.pump$MonthNum <- match(df.pump$Month, month.abb)
df.pump$MonthLengthDays <- lubridate::days_in_month(df.pump$MonthNum)
df.pump$m3PlantDay <- df.pump$MeanWaterUse_GalPlantDay*gal.to.m3

# set up long data frames for intermittent_pumping script; separate for outdoor and greenhouse
t.max.yrs <- max(yrs.model)
df.pump.outdoor <- 
  subset(df.pump, Setting=="Outdoor") %>% 
  dplyr::select(MonthNum, MonthLengthDays, m3PlantDay) %>% 
  replicate(t.max.yrs, ., simplify = FALSE) %>% 
  dplyr::bind_rows() %>% 
  transform(EndOfMonthDays = cumsum(MonthLengthDays))
df.pump.outdoor$StartOfMonthDays <- c(1, df.pump.outdoor$EndOfMonthDays[1:(t.max.yrs*12)-1]+1)

df.pump.greenhouse <- 
  subset(df.pump, Setting=="Greenhouse") %>% 
  dplyr::select(MonthNum, MonthLengthDays, m3PlantDay) %>% 
  replicate(t.max.yrs, ., simplify = FALSE) %>% 
  dplyr::bind_rows() %>% 
  transform(EndOfMonthDays = cumsum(MonthLengthDays))
df.pump.greenhouse$StartOfMonthDays <- c(1, df.pump.greenhouse$EndOfMonthDays[1:(t.max.yrs*12)-1]+1)

## subdivide sf.streams into points for web squared calculation
# define point spacing and figure out how many points to make
pt.spacing <- 25  # [m]
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
  st_distance(x=sf.streams.pts, y=sf.grows) %>% 
  as.numeric() %>% 
  data.frame(
    SegNum = rep(sf.streams.pts$SegNum, times = dim(sf.grows)[1]),
    GrowNum = rep(sf.grows$GrowNum, each = dim(sf.streams.pts)[1]),
    dist_wellToStream_m = .,
    lon = rep(pt_coords[,"X"], times = dim(sf.grows)[1]),
    lat = rep(pt_coords[,"Y"], times = dim(sf.grows)[1])
  )

## well-stream geometry output from Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R
df.all <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_CalculateWellStreamPairs.csv") %>% 
  read.csv(stringsAsFactors=F)

# predict top of screen based on linear relationship with dtb_m
# see script Navarro_WellCompletionReports.R
coef_m <- -3.02
coef_b <- 101.11
df.all$well_screenTop_depth_m <- coef_m * df.all$well_dtb_m + coef_b   # depth to top of well screen [m]
df.all$well_screenTop_elev_m <- df.all$well_elev_m - df.all$well_screenTop_depth_m  # elevation of top of well screen [m]

## calculate aquifer and stream hydrostratigraphic values
# define K and S for unconsolidated sediment (K_u) and fractured bedrock (K_f)
K_u <- 86400*10^(-3.52)  # [m/d]
K_f <- 86400*10^(-8.2)   # [m/d]

S_u <- 0.27  # [-]
S_f <- 0.19  # [-]

# calculate transmissivity
# Reeves et al (2009) associate transmissivity with each stream segment,
# such that a stream segment has the same transmissivity regardless of where
# the pumping well is located.
#
# preliminary version: for each stream segment, define the transmissivity as:
# The weighted transmissivity of the alluvial aquifer and any bedrock between the
# top of the well screen and the stream. If the top of the well screen is higher
# than the stream, the transmissivity will be equal to just the transmissivity of
# the alluvial aquifer.
#
# this is following Kollet & Zlotnik

df.all$m_u <- df.all$stream_dtb_m  # unconsolutated sediment thickness [m]
df.all$m_f <- abs((df.all$stream_elev_m - df.all$stream_dtb_m) - df.all$well_screenTop_elev_m)  # fractured bedrock thickness [m]

df.all$K_eff <- (K_u*df.all$m_u + K_f*df.all$m_f)/(df.all$m_u + df.all$m_f)  # effective K [m/d]
df.all$Tr_eff <- df.all$K_eff*(df.all$m_u + df.all$m_f)  # effective Tr [m2/d]

df.all$S_eff <- (S_u*df.all$m_u + S_f*df.all$m_f)/(df.all$m_u + df.all$m_f)  # effective S [-]

# calculate streambed conductance following Reeves
df.all$lmda <- 0.1*K_u*df.all$stream_width_m/df.all$stream_dtb_m  # streambed conductance [m2/d]

# add information about grow conditions
df.all <- left_join(df.all, as.data.frame(sf.grows)[,c("GrowNum", "greenhouse", "growsize", "plants")], by="GrowNum")

### depletion calculations
## loop through times [d] for depletion calculations
min_frac <- 0.01  # minimum depletion considered

w.start.flag <- T
counter <- 0
for (w in unique(df.all$GrowNum)){
  
  wel_coord <- 
    sf.grows %>% 
    subset(GrowNum == w) %>% 
    st_coordinates()
  
  # retain only relevant points
  dist_w_pts <- subset(dist_all_pts, GrowNum==w)
  
  # first: Theissen polygons to figure out adjacent catchments
  df.apportion.t <- 
    dist_w_pts %>% 
    group_by(SegNum) %>% 
    filter(dist_wellToStream_m == max(dist_wellToStream_m)) %>% 
    apportion_polygon(., 
                      wel_lon = wel_coord[1,"X"], 
                      wel_lat = wel_coord[1,"Y"],
                      crs = CRS(st_crs(sf.streams)[["proj4string"]]),
                      reach_name = "SegNum",
                      dist_name = "dist_wellToStream_m") %>% 
    set_colnames(c("SegNum", "frac_depletion"))
  
  
  t.start.flag <- T
  max_dist_prev <- 500
  for (time_days in DOYs.all){
    
    # find maximum distance, based on maximum observed S, Tr, lmda (inclusive estimate)
    max_dist <- depletion_max_distance(Qf_thres = min_frac,
                                       d_interval = 250,
                                       d_min = max_dist_prev,
                                       d_max = max(df.all$dist_wellToStream_m[df.all$GrowNum==w]),
                                       method="hunt",
                                       t = time_days,
                                       S = max(df.all$S_eff[df.all$GrowNum==w]),
                                       Tr = max(df.all$Tr_eff[df.all$GrowNum==w]),
                                       lmda = max(df.all$lmda[df.all$GrowNum==w]))
    
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
      transform(GrowNum = w,
                time_days = time_days)
    
    if (t.start.flag){
      df.apportion <- df.apportion.w.t
      t.start.flag <- F
    } else {
      df.apportion <- rbind(df.apportion, df.apportion.w.t)
    }
    
    # update max_dist starting value
    max_dist_prev <- max_dist
    
  }
  
  # join with df_all
  df.all.t <- 
    left_join(df.apportion, 
              df.all[,c("SegNum", "GrowNum", "dist_wellToStream_m", "S_eff", "Tr_eff", "lmda", "greenhouse", "growsize", "plants")], 
              by=c("SegNum", "GrowNum"))
  
  if (w.start.flag){
    df.out <- df.all.t
    w.start.flag <- F
  } else {
    df.out <- rbind(df.out, df.all.t)
  }
  
  # status update
  counter <- counter+1
  print(paste0("Well ", counter, " of ", length(unique(df.all$GrowNum)), " depletion apportionment complete, ", Sys.time()))
}

## unique well-seg combos
df.combos <- 
  df.out %>% 
  dplyr::select(SegNum, GrowNum, dist_wellToStream_m, S_eff, Tr_eff, lmda, greenhouse, growsize, plants) %>% 
  unique()

## calculate depletion for all combos
start.flag.Qs <- T
for (i in 1:dim(df.combos)[1]){
  # identify well-seg combo
  seg <- df.combos$SegNum[i]
  grow <- df.combos$GrowNum[i]
  
  # get times
  output_t_days <- df.out$time_days[df.out$SegNum==seg & df.out$GrowNum==grow]
  output_frac <- df.out$frac_depletion[df.out$SegNum==seg & df.out$GrowNum==grow]
  
  # different calculation depending on outdoor vs indoor plants
  if (df.combos$greenhouse[i] == 1){
    Qs <- intermittent_pumping(t = output_t_days,
                               starts = df.pump.greenhouse$StartOfMonthDays,
                               stops  = df.pump.greenhouse$EndOfMonthDays,
                               rates  = df.pump.greenhouse$m3PlantDay*df.combos$plants[i],
                               method = "hunt",
                               d = df.combos$dist_wellToStream_m[i],
                               S = df.combos$S_eff[i],
                               Tr = df.combos$Tr_eff[i],
                               lmda = df.combos$lmda[i])
  } else {
    Qs <- intermittent_pumping(t = output_t_days,
                               starts = df.pump.outdoor$StartOfMonthDays,
                               stops  = df.pump.outdoor$EndOfMonthDays,
                               rates  = df.pump.outdoor$m3PlantDay*df.combos$plants[i],
                               method = "hunt",
                               d = df.combos$dist_wellToStream_m[i],
                               S = df.combos$S_eff[i],
                               Tr = df.combos$Tr_eff[i],
                               lmda = df.combos$lmda[i])
  }
  
  # compile output
  df.depletion <- data.frame(SegNum = seg,
                             GrowNum = grow,
                             time_days = output_t_days,
                             Qs = Qs)
  
  if (start.flag.Qs){
    df.Qs <- df.depletion
    start.flag.Qs <- F
  } else {
    df.Qs <- rbind(df.Qs, df.depletion)
  }
  
  print(paste0("Depletion ", i, " of ", dim(df.combos)[1], " complete, ", Sys.time()))
  
}

# combine and save
df.Qs %>% 
  left_join(df.out, by=c("SegNum", "GrowNum", "time_days")) %>% 
  transform(depletion_m3d = Qs*frac_depletion) %>% 
  format(digits = 3) %>% 
  write.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_DepletionBySegment.csv"),
            row.names=F)
