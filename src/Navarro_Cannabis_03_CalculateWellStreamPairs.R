## Navarro_Cannabis_CalculateWellStreamPairs.R
#' This script is intended to calculate analytical model inputs for each well-stream combination:
#'   -d = distance to stream [m]
#'   -S = effective storage coefficient [-]
#'   -Tr = effective transmissivity [L2/T]

source(file.path("src", "paths+packages.R"))

## load data
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F)

# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# synthetic pumping wells
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), header=T)
sf.wel <- st_as_sf(df.wel, coords = c("lon", "lat"), crs = st_crs(sf.streams))

# well completion reports
sf.wcr <- 
  sf::st_read(file.path("data", "WellCompletionReports", "Navarro_WellCompletionReports.shp"), stringsAsFactors=F) %>% 
  sf::st_transform(crs.MODFLOW) %>% 
  subset(is.finite(Top.Of.Per) &
           is.finite(Bottom.of) &
           !(Planned.Us %in% c("Injection", "Monitoring"))) %>% 
  transform(screenTopDepth_m = Top.Of.Per*0.3048,
            screenBotDepth_m = Bottom.of*0.3048) %>% 
  dplyr::select(WCR.Number, Planned.Us, screenTopDepth_m, screenBotDepth_m, geometry) %>% 
  sf::st_as_sf()

avg.screen.length <- round(mean(sf.wcr$screenBotDepth_m - sf.wcr$screenTopDepth_m), 1)

# rasters
r.lulc <- raster(paste0(dir.gis, "Navarro_Cannabis_LULC_30m.tif"))
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))
r.aquifers <- raster(paste0(dir.gis, "Navarro_Cannabis_GroundwaterBasins_30m.tif"))
r.wte <- raster(paste0(dir.gis, "Navarro_Cannabis_WTE_30m.tif"))

## calculate saturated thickness of alluvial materials are saturated or not
r.wtd <- r.dem.30m - r.wte
r.alluvial.sat.thickness <- r.dtb - r.wtd
r.alluvial.sat.thickness[r.alluvial.sat.thickness < r.dtb] <- 0
r.alluvial.sat.thickness[r.alluvial.sat.thickness > r.dtb] <- r.dtb[r.alluvial.sat.thickness > r.dtb]

## extract some potentially relevant data
sf.wel$elev_m <- raster::extract(r.dem.30m, sf.wel)  # elevation of that grid cell [m]
sf.wel$dtb_m <- raster::extract(r.dtb, sf.wel)       # depth to bedrock 
sf.wel$wtd_m <- raster::extract(r.wtd, sf.wel)       # water table depth
sf.wel$wte_m <- sf.wel$elev_m - sf.wel$wtd_m         # water table elevation
sf.wel$screenTopDepth_m <- sf.wel$wtd_m        # depth of top of well screen [m]
sf.wel$screenTopDepth_m[sf.wel$screenTopDepth_m < 0] <- 0
sf.wel$screenTopElev_m <- sf.wel$elev_m - sf.wel$screenTopDepth_m  # elevation of top of well screen [m]
sf.wel$screenBotDepth_m <- sf.wel$screenTopDepth_m + avg.screen.length  # depth of bottom of well screen [m]
sf.wel$screenBotElev_m <- sf.wel$elev_m - sf.wel$screenBotDepth_m
sf.wel[,c("lon", "lat")] <- st_coordinates(sf.wel)

sf.streams$elev_m <- raster::extract(r.dem.30m, sf.streams, fun='mean', na.rm=T)  # mean elevation of all grid cells stream touches [m]
sf.streams$dtb_m <- raster::extract(r.dtb, sf.streams, fun='mean', na.rm=T)       # mean DTB of all grid cells stream touches [m]
#sf.streams$satThick_m <- raster::extract(r.alluvial.sat.thickness, sf.streams, fun='mean', na.rm=T)       # mean DTB of all grid cells stream touches [m]
#sf.streams$satThick_m[is.na(sf.streams$satThick_m)] <- 0
sf.streams[,c("lon", "lat")] <- 
  sf.streams %>% 
  st_centroid() %>% 
  st_coordinates()

sum(is.na(sf.streams$elev_m))
sum(is.na(sf.streams$dtb_m))

## predict stream width based on drainage area
sf.streams$width_m <- WidthFromDA(DA=sf.streams$TtDASKM, w.min=1, w.max=100)

## distance from each well to each stream segment
dist_all <- 
  st_distance(x=sf.streams, y=sf.wel)

df.all <- data.frame(
  SegNum = rep(sf.streams$SegNum, times = dim(dist_all)[2]),
#  stream_lon = rep(sf.streams$lon, times = dim(dist_all)[2]),
#  stream_lat = rep(sf.streams$lat, times = dim(dist_all)[2]),
  stream_elev_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  stream_dtb_m = rep(sf.streams$dtb_m, times = dim(dist_all)[2]),
  stream_width_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
  WellNum = rep(sf.wel$WellNum, each = dim(dist_all)[1]),
#  well_lon = rep(sf.wel$lon, times = dim(dist_all)[1]),
#  well_lat = rep(sf.wel$lat, times = dim(dist_all)[1]),
  well_elev_m = rep(sf.wel$elev_m, times = dim(dist_all)[1]),
  well_dtb_m = rep(sf.wel$dtb_m, times = dim(dist_all)[1]),
  dist_wellToStream_m = as.numeric(dist_all)
)

# max number of segments to consider as potential options for each well
#   selected this number iteratively - did a preliminary analysis considering all well-stream combos
#   then looked at the max number of segments affected by a given well (which was 29)
max.well.segs <- 30

df.all <- 
  data.frame(
    SegNum = rep(sf.streams$SegNum, times = dim(dist_all)[2]),
    stream_elev_m = rep(sf.streams$elev_m, times = dim(dist_all)[2]),
    stream_dtb_m = rep(sf.streams$dtb_m, times = dim(dist_all)[2]),
    stream_width_m = rep(sf.streams$width_m, times = dim(dist_all)[2]),
    WellNum = rep(sf.wel$WellNum, each = dim(dist_all)[1]),
    well_elev_m = rep(sf.wel$elev_m, times = dim(dist_all)[1]),
    well_screenTopDepth_m = rep(sf.wel$screenTopDepth_m, times = dim(dist_all)[1]),
    well_screenBotDepth_m = rep(sf.wel$screenBotDepth_m, times = dim(dist_all)[1]),
    well_wte_m = rep(sf.wel$wte_m, times = dim(dist_all)[1]),
    well_dtb_m = rep(sf.wel$dtb_m, times = dim(dist_all)[1]),
    dist_wellToStream_m = as.numeric(dist_all),
    mu_well_m = NaN,
    mf_well_m = NaN,
    mu_stream_m = NaN,
    mf_stream_m = NaN,
    Tr_well_m2d = NaN,
    Tr_stream_m2d = NaN,
    Tr_bulk_m2d = NaN,
    S_well = NaN,
    S_stream = NaN,
    S_bulk = NaN,
    lmda_m2d = NaN
  ) %>% 
  group_by(WellNum) %>% 
  top_n(-max.well.segs, dist_wellToStream_m)  # - value to select smallest distances for each WellNum

## calculate aquifer and stream hydrostratigraphic values
# coarsen data for extraction
r.wte.coarse <- raster::aggregate(r.wte, fact=5)
r.elev.coarse <- raster::aggregate(r.dem.30m, fact=5)
r.dtb.coarse <- raster::aggregate(r.dtb, fact=5)

## for each well-stream pair...
start.time <- Sys.time()
for (i in 1:dim(df.all)[1]){
  ## T and S calculations (single points)
  # extract useful values
  stream_elev_m <- df.all$stream_elev_m[i]
  stream_bedrock_elev_m <- stream_elev_m - df.all$stream_dtb_m[i]
  stream_width_m <- df.all$stream_width_m[i]
  well_elev_m <- df.all$well_elev_m[i]
  well_wte_m <- df.all$well_wte_m[i]
  well_bedrock_elev_m <- well_elev_m - df.all$well_dtb_m[i]
  screenTop_elev_m <- well_elev_m - df.all$well_screenTopDepth_m[i]
  screenBot_elev_m <- well_elev_m - df.all$well_screenBotDepth_m[i]
  
  # minimum allowed aquifer thickness [m] - needed because some wells only have ~3 m screen itnervals which is not reasonable
  min_thickness_m <- avg.screen.length  # don't allow < this value
  min_aq_thick_m <- 
    max(c(
      (max(c(screenTop_elev_m, screenBot_elev_m, stream_elev_m)) - 
         min(c(screenTop_elev_m, screenBot_elev_m, stream_elev_m))),
      min_thickness_m
    ))
  
  ## calculate transmissivity at well
  # set top and bottom of flow domain
  # top = max of well screen top and WTE, constrained by land surface 
  well_T_top_m <- min(c(well_elev_m, max(c(well_wte_m, screenTop_elev_m))))
  # bottom = min of well screen bottom and stream elevation, constrained by min allowed thickness
  well_T_bot_m <- min(c((well_T_top_m - min_aq_thick_m), screenBot_elev_m, stream_elev_m))
  
  if (well_T_top_m > well_bedrock_elev_m){
    # top of transmissivity calculation is in alluvial
    mu_well <- well_T_top_m - well_bedrock_elev_m
  } else {
    # top of transmissivity calculation is in bedrock
    mu_well <- 0
  }
  
  if (well_T_bot_m > well_bedrock_elev_m){
    # bottom of well screen is in alluvial
    mf_well <- 0
  } else {
    # bottom of well screen is in bedrock
    mf_well <- well_T_top_m - well_T_bot_m - mu_well
  }
  
  # weight based on flow parallel to layering (horizontal)
  Tr_well <- (mu_well + mf_well)*(K_u*mu_well + K_f*mf_well)/(mu_well + mf_well)
  S_well <- (S_u*mu_well + S_f*mf_well)/(mu_well + mf_well)
  
  ## calculate transmissivity at stream - bottom of stream to bottom of well screen
  mu_stream <- max(c(0, (stream_elev_m - stream_bedrock_elev_m)))
  mf_stream <- max(c(0, (abs(stream_elev_m - screenBot_elev_m) - mu_stream)))
  
  # weight based on flow parallel to layering (horizontal)
  Tr_stream <- (mu_stream + mf_stream)*(K_u*mu_stream + K_f*mf_stream)/(mu_stream + mf_stream)
  S_stream <- (S_u*mu_stream + S_f*mf_stream)/(mu_stream + mf_stream)
  
  # calculate lmda based on 1/10 of Tr of alluvial only (similar to Reeves)
  lmda <- 0.1*K_u*stream_width_m/mu_stream  # streambed conductance [m2/d]
  
  ## now the hard part... calculate based on all points intersecting a line
  # make a line connecting well to closest point on stream segment
  sf.connector <- st_nearest_points(subset(sf.wel, WellNum==df.all$WellNum[i]), 
                                    subset(sf.streams, SegNum==df.all$SegNum[i])) %>% 
    st_sf() %>% 
    as(., "Spatial") %>% 
    as(., "SpatialLines")
  
  # at each point along segment extract data
  df.connector <- 
    data.frame(
      wte_m = raster::extract(r.wte.coarse, sf.connector)[[1]],
      dtb_m = raster::extract(r.dtb.coarse, sf.connector)[[1]],
      elev_m = raster::extract(r.elev.coarse, sf.connector)[[1]]
    )
  df.connector$bedrock_elev_m <- df.connector$elev_m - df.connector$dtb_m
  df.connector <- df.connector[complete.cases(df.connector), ]  # remove NAs - happens if connector cuts across an area outside domain
  
  # elevation of top of aquifer at each location
  df.connector$aq_top_elev_m <- pmax(df.connector$wte_m, stream_elev_m)
  df.connector$aq_bot_elev_m <- min(c(stream_elev_m, screenBot_elev_m))
  df.connector$aq_thickness_m <- df.connector$aq_top_elev_m - df.connector$aq_bot_elev_m
  
  # check for too-thin aquifers
  df.connector$aq_bot_elev_m[df.connector$aq_thickness_m < min_aq_thick_m] <-
    df.connector$aq_top_elev_m[df.connector$aq_thickness_m < min_aq_thick_m] - min_aq_thick_m
  df.connector$aq_thickness_m <- df.connector$aq_top_elev_m - df.connector$aq_bot_elev_m
  
  if (sum(df.connector$aq_thickness_m < 0) > 0) stop("negative aquifer thickness")
  
  df.connector$mu <- df.connector$aq_top_elev_m - df.connector$bedrock_elev_m
  df.connector$mu[df.connector$mu < 0] <- 0
  df.connector$mu[df.connector$mu > df.connector$aq_thickness_m] <- 
    df.connector$aq_thickness_m[df.connector$mu > df.connector$aq_thickness_m]
  
  df.connector$mf <- df.connector$aq_thickness_m - df.connector$mu
  df.connector$mf[df.connector$mf < 0] <- 0
  
  if (sum((df.connector$mf + df.connector$mu) - (df.connector$aq_top_elev_m - df.connector$aq_bot_elev_m) > 1e-2) > 0) stop("Miscalculated thicknesses")
  
  # calculate Tr in each cell based on flow perpendicular to layers since each layer = 1 raster cell 
  # which already has bulk transmissivity averaged for parallel flow
  df.connector$Tr_m2d <- 
    (df.connector$mu + df.connector$mf)*
    (K_u*df.connector$mu + K_f*df.connector$mf)/(df.connector$mu + df.connector$mf)
  df.connector$S <- 
    (S_u*df.connector$mu + S_f*df.connector$mf)/(df.connector$mu + df.connector$mf)
  cell.size <- mean(res(r.wte.coarse))  # cell size in m
  Tr_bulk_m2d <- (cell.size*dim(df.connector)[1])/sum(cell.size/df.connector$Tr_m2d)
  S_bulk <- mean(df.connector$S)
  
  # fill in data frame
  df.all$mu_well_m[i] <- mu_well
  df.all$mf_well_m[i] <- mf_well
  df.all$mu_stream_m[i] <- mu_stream
  df.all$mf_stream_m[i] <- mf_stream
  df.all$Tr_well_m2d[i] <- Tr_well
  df.all$S_well[i] <- S_well
  df.all$Tr_stream_m2d[i] <- Tr_stream
  df.all$S_stream[i] <- S_stream
  df.all$lmda_m2d[i] <- lmda
  df.all$Tr_bulk_m2d[i] <- Tr_bulk_m2d
  df.all$S_bulk[i] <- S_bulk
  
  # status update
  print(paste0(i, " of ", dim(df.all)[1], " complete, ", round(difftime(Sys.time(), start.time, units="min"), 2), " min"))
  
}

## save output
df.all %>% 
  dfDigits(x=., digits=3) %>% 
  write.csv(file.path("results", "Navarro_Cannabis_03_CalculateWellStreamPairs.csv"),
            row.names=F)