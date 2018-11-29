## Navarro_Cannabis_DepletionBySegment.R
#' This script is intended to calculate depletion for stream segments 
#' with high intrinsic habitat potential.
#' 
#' It requires output from Navarro_Cannabis_CalculateWellStreamPairs.R 
#' and Navarro_Cannabis_HabitatIntrinsicPotential.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

### load and pre-process data
# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

## output from Navarro_Cannabis_HabitatIntrinsicPotential.R
df.habitat <- 
  file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  melt(id=c("SegNum"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1), 
                           labels=c("Low", "High"),
                           include.lowest=T)) %>% 
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3])

# summarize each segment to max value for any species
df.habitat.summary <-
  df.habitat %>% 
  subset(metric=="mean") %>% 
  group_by(SegNum) %>% 
  summarize(IP = max(IP)) %>% 
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1, 5),   # nothing will get Outside Navarro category
                           labels=c("Low", "High", "Outside Navarro"),
                           include.lowest=T))

# join habitat data with stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  left_join(df.habitat.summary, by=c("SegNum")) %>% 
  replace_na(list("IP_class" = "Outside Navarro"))

# subdivide sf.streams into points for web squared calculation
# define point spacing and figure out how many points to make
pt.spacing <- 10  # [m]
shp.streams <- as(sf.streams, "Spatial")
n.pts <- round(gLength(shp.streams)/pt.spacing)
set.seed(1)
shp.streams.pts <- spsample(shp.streams, n=n.pts, type="regular")
df.streams.pts <- as.data.frame(shp.streams.pts)
colnames(df.streams.pts) <- c("lon", "lat")
shp.streams.buffer <- buffer(shp.streams, 0.1, dissolve=F) # figure out what SegNum each point corresponds to
int <- intersect(shp.streams.pts, shp.streams.buffer)
df.streams.pts <- cbind(df.streams.pts, int@data)
sf.streams.pts <- st_as_sf(df.streams.pts, coords=c("lon", "lat"), crs=st_crs(sf.streams))

## well locations
# synthetic pumping wells
sf.wel <- 
  file.path("modflow", "input", "iwel.txt") %>% 
  read.table(header=T) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(sf.streams))

## find distance from each stream point to each stream segment
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
  )

## well-stream geometry output from Navarro_Cannabis_CalculateWellStreamPairs.R
df.all <- 
  file.path("results", "Navarro_Cannabis_CalculateWellStreamPairs.csv") %>% 
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

### depletion calculations
## loop through times [d] for depletion calculations
min_frac <- 0.01  # minimum depletion considered

start.flag <- T
for (time_days in c(1*365, 5*365, 10*365, 25*365)){

  # find maximum distance, based on maximm observed S, Tr, lmda (inclusive estimate)
  max_dist <- depletion_max_distance(Qf_thres = min_frac,
                                     d_interval = 100,
                                     d_min = NULL,
                                     d_max = Inf,
                                     method="hunt",
                                     t = time_days,
                                     S = max(df.all$S_eff),
                                     Tr = max(df.all$Tr_eff),
                                     lmda = max(df.all$lmda))
  
  ## apportion for each WellNum
  w.start.flag <- T
  for (w in unique(df.all$WellNum)){
    wel_coord <- 
      sf.wel %>% 
      subset(WellNum == w) %>% 
      st_coordinates()
    
    # first: Theissen polygons to figure out adjacent catchments
    df.apportion.t <- 
      subset(dist_all_pts, WellNum==w) %>% 
      group_by(SegNum) %>% 
      filter(dist_wellToStream_m == max(dist_wellToStream_m)) %>% 
      apportion_polygon(., 
                        wel_lon = wel_coord[1,"X"], 
                        wel_lat = wel_coord[1,"Y"],
                        crs = CRS(st_crs(sf.streams)[["proj4string"]]),
                        reach_name = "SegNum",
                        dist_name = "dist_wellToStream_m") %>% 
      set_colnames(c("SegNum", "frac_depletion"))
    
    # second: use web^2 to apportion to any stream segment that is within max_dist
    #         OR that is adjacent to well based on Theissen Polygon
    df.apportion.w <- 
      dist_all_pts %>% 
      subset(WellNum==w & ((dist_wellToStream_m <= max_dist) | (SegNum %in% df.apportion.t$SegNum))) %>% 
      apportion_web(reach_dist = ., 
                    w = 2, 
                    min_frac = min_frac,
                    reach_name = "SegNum",
                    dist_name = "dist_wellToStream_m") %>% 
      set_colnames(c("SegNum", "frac_depletion")) %>% 
      transform(WellNum = w,
                time_d = time_days)
    
    if (w.start.flag){
      df.apportion <- df.apportion.w
      w.start.flag <- F
    } else {
      df.apportion <- rbind(df.apportion, df.apportion.w)
    }
    
    print(paste0("Time ", time_days, ", well ", w, " depletion apportionment complete, ", Sys.time()))
    
  }
  
  # join with df_all
  df.all.t <- 
    left_join(df.apportion, df.all, by=c("SegNum", "WellNum"))
  
  if (start.flag){
    df.out <- df.all.t
    start.flag <- F
  } else {
    df.out <- rbind(df.out, df.all.t)
  }
}

## calculate depletion
df.out$Qf <- hunt(t = df.out$time_d,
                  d = df.out$dist_wellToStream_m,
                  S = df.out$S_eff,
                  Tr = df.out$Tr_eff,
                  lmda = df.out$lmda) * 
  df.out$frac_depletion

## save output
df.out %>% 
  format(digits = 3) %>% 
  write.csv(file.path("results", "Navarro_Cannabis_DepletionBySegment.csv"),
            row.names=F)
