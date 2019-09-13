## Figure_ScaleSegment_Depletion.R
#' This script is intended to calculate depletion for stream segments 
#' with high intrinsic habitat potential.
#' 
#' It requires output from Navarro_Cannabis_05_DepletionBySegment.R 
#' and Navarro_Cannabis_HabitatIntrinsicPotential.R

source(file.path("src", "paths+packages.R"))

### load and pre-process data
## read in data from Navarro_Cannabis_DepletionBySegment.R
df.depletion <- read.csv(file.path("results", "Navarro_Cannabis_05_DepletionBySegment.csv"), stringsAsFactors=F)

## output from Navarro_Cannabis_02_HabitatIntrinsicPotential.R
df.habitat <- 
  file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  melt(id=c("SegNum"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1.01,100),     # nothing will get Outside Navarro category
                           labels=c("Low", "High", "Outside Navarro"),
                           include.lowest=T)) %>% 
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3]) %>% 
  subset(species == "Coho" & metric=="max")  # focus on Coho based on conversation with Jen - most sensitive to habitat conditions

# join habitat data with stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  left_join(df.habitat, by=c("SegNum")) %>% 
  replace_na(list("IP_class" = "Outside Navarro"))

## well locations
# synthetic pumping wells
sf.wel <- 
  file.path("modflow", "input", "iwel.txt") %>% 
  read.table(header=T) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(sf.streams))

## domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>% 
  sf::st_transform(crs=st_crs(sf.streams))

### combine data
# add IP_class to depletion
df.depletion <- 
  left_join(df.depletion, df.habitat, by="SegNum") %>% 
  replace_na(list("IP_class"="Outside Navarro"))

# summarize depletion in high ecological value segments by time_days and WellNum
df.value <-
  df.depletion %>% 
  group_by(time_days, WellNum, IP_class, pump_factor) %>% 
  summarize(depletion_m3d.HighValue = sum(depletion_m3d))

sf.value <- 
  full_join(sf.wel, df.value, by="WellNum")

### plots
# times to plot
pumps_plot <- unique(df.value$pump_factor)
times_plot <- unique(df.value$time_days)

## for facets, need to duplicate basin and stream for all timesteps
## also create rasters of interpolated depletion
r.empty <- raster(extent(sf.basin), crs=crs(sf.basin))
res(r.empty) <- c(150,150)

#### Analysis of streamflow depletion at different buffer distances
# also include depth to bedrock and water table depth
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))
r.wte <- raster(paste0(dir.gis, "Navarro_Cannabis_WTE_30m.tif"))
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
r.wtd <- r.dem.30m - r.wte

# scroll through pump factors
for (p in pumps_plot){
  ## inverse distance interpolation to make rasters
  # not doing this as a for loop because I want to save individual rasters
  df.interp.1 <- data.frame(Z = subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[1])$depletion_m3d.HighValue,
                            X = st_coordinates(subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[1]))[,"X"],
                            Y = st_coordinates(subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[1]))[,"Y"])
  sp.interp.1 <- SpatialPoints(df.interp.1[,c("X", "Y")], proj4string=CRS(st_crs(sf.basin)[["proj4string"]]))
  sp.interp.1 <- SpatialPointsDataFrame(sp.interp.1, df.interp.1)
  gs.1 <- gstat(formula = Z~1, locations = sp.interp.1, nmax = 4)  # nmax controls weighting
  r.t.1 <- interpolate(r.empty, gs.1)
  r.t.1.mask <- mask(r.t.1, sf.basin)
  df.Qf.t.1 <- 
    r.t.1.mask %>% 
    raster::rasterToPoints() %>% 
    raster::as.data.frame() %>% 
    set_colnames(c("lon", "lat", "depletion_m3d.HighValue")) %>% 
    transform(time_days = times_plot[1])
  
  df.interp.2 <- data.frame(Z = subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[2])$depletion_m3d.HighValue,
                            X = st_coordinates(subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[2]))[,"X"],
                            Y = st_coordinates(subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[2]))[,"Y"])
  sp.interp.2 <- SpatialPoints(df.interp.2[,c("X", "Y")], proj4string=CRS(st_crs(sf.basin)[["proj4string"]]))
  sp.interp.2 <- SpatialPointsDataFrame(sp.interp.2, df.interp.2)
  gs.2 <- gstat(formula = Z~1, locations = sp.interp.2, nmax = 4)  # nmax controls weighting
  r.t.2 <- interpolate(r.empty, gs.2)
  r.t.2.mask <- mask(r.t.2, sf.basin)
  df.Qf.t.2 <- 
    r.t.2.mask %>% 
    raster::rasterToPoints() %>% 
    raster::as.data.frame() %>% 
    set_colnames(c("lon", "lat", "depletion_m3d.HighValue")) %>% 
    transform(time_days = times_plot[2])
  
  df.interp.3 <- data.frame(Z = subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[3])$depletion_m3d.HighValue,
                            X = st_coordinates(subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[3]))[,"X"],
                            Y = st_coordinates(subset(sf.value, IP_class=="High" & pump_factor == p & time_days == times_plot[3]))[,"Y"])
  sp.interp.3 <- SpatialPoints(df.interp.3[,c("X", "Y")], proj4string=CRS(st_crs(sf.basin)[["proj4string"]]))
  sp.interp.3 <- SpatialPointsDataFrame(sp.interp.3, df.interp.3)
  gs.3 <- gstat(formula = Z~1, locations = sp.interp.3, nmax = 4)  # nmax controls weighting
  r.t.3 <- interpolate(r.empty, gs.3)
  r.t.3.mask <- mask(r.t.3, sf.basin)
  df.Qf.t.3 <- 
    r.t.3.mask %>% 
    raster::rasterToPoints() %>% 
    raster::as.data.frame() %>% 
    set_colnames(c("lon", "lat", "depletion_m3d.HighValue")) %>% 
    transform(time_days = times_plot[3])
  
  ## combine all rasters for ggplot
  df.Qf <- rbind(df.Qf.t.1, df.Qf.t.2, df.Qf.t.3)
  
  # make stream and basin boundaries facet-compatible
  sf.streams.facet <- 
    sf.streams %>% 
    dplyr::mutate(time_days = times_plot[1])
  sf.streams.facet <- 
    sf.streams %>% 
    dplyr::mutate(time_days = times_plot[2]) %>% 
    rbind(sf.streams.facet, .)
  sf.streams.facet <- 
    sf.streams %>% 
    dplyr::mutate(time_days = times_plot[3]) %>% 
    rbind(sf.streams.facet, .)
  
  sf.basin.facet <- 
    sf.basin %>% 
    dplyr::mutate(time_days = times_plot[1])
  sf.basin.facet <- 
    sf.basin %>% 
    dplyr::mutate(time_days = times_plot[2]) %>% 
    rbind(sf.basin.facet, .)
  sf.basin.facet <- 
    sf.basin %>% 
    dplyr::mutate(time_days = times_plot[3]) %>% 
    rbind(sf.basin.facet, .)
  
  ## cycle through buffer intervals
  buff_interval <- 100
  buffers <- seq(buff_interval, 3000, buff_interval)
  r.wtd[r.wtd < 0] <- 0
  
  r.dtb.coarse <- raster::aggregate(r.dtb, fact=res(r.t.3.mask)[1]/res(r.dtb)[1])
  r.wtd.coarse <- raster::aggregate(r.wtd, fact=res(r.t.3.mask)[1]/res(r.wtd)[1])
  for (b in 1:length(buffers)){
    
    # bookkeeping
    buff_dist <- buffers[b]
    if (b > 1) buff_previous <- buffers[b-1] else buff_previous <- 0
    if (b > 1) sf_previous <- sf.streams.buffer else sf_previous <- NA
    
    # create buffer around stream segments
    sf.streams.buffer <-
      sf.streams %>%
      subset(IP_class=="High") %>%
      sf::st_buffer(dist=buff_dist, endCapStyle="FLAT", joinStyle="MITRE")
    
    # remove part of polygon already previously extracted
    if (b > 1){
      sf.streams.donut <-
        sf::st_difference(sf.streams.buffer, sf_previous) %>%
        subset(SegNum==SegNum.1)
    } else {
      sf.streams.donut <- sf.streams.buffer
    }
    
    # extract mean from raster within that buffer for all stream segments
    df.buff.dist <-
      rbind(data.frame(buff_dist = buff_dist,
                       buff_range = paste0(buff_previous, "-", buff_dist, " m"),
                       buff_range_center = (buff_previous + buff_dist)/2,
                       time_days = times_plot[1],
                       pump_factor = p,
                       SegNum = sf.streams.buffer$SegNum,
                       IP_class = sf.streams.buffer$IP_class,
                       depletion_m3d_HighValue = raster::extract(r.t.1.mask, sf.streams.donut, fun=mean, na.rm=T, weight=T)),
            data.frame(buff_dist = buff_dist,
                       buff_range = paste0(buff_previous, "-", buff_dist, " m"),
                       buff_range_center = (buff_previous + buff_dist)/2,
                       time_days = times_plot[2],
                       pump_factor = p,
                       SegNum = sf.streams.buffer$SegNum,
                       IP_class = sf.streams.buffer$IP_class,
                       depletion_m3d_HighValue = raster::extract(r.t.2.mask, sf.streams.donut, fun=mean, na.rm=T, weight=T)),
            data.frame(buff_dist = buff_dist,
                       buff_range = paste0(buff_previous, "-", buff_dist, " m"),
                       buff_range_center = (buff_previous + buff_dist)/2,
                       time_days = times_plot[3],
                       pump_factor = p,
                       SegNum = sf.streams.buffer$SegNum,
                       IP_class = sf.streams.buffer$IP_class,
                       depletion_m3d_HighValue = raster::extract(r.t.3.mask, sf.streams.donut, fun=mean, na.rm=T, weight=T))
      )
    
    # extract bedrock depth
    if (p == pumps_plot[1]) {
      df.dtb.dist <- data.frame(buff_dist = buff_dist,
                                buff_range = paste0(buff_previous, "-", buff_dist, " m"),
                                buff_range_center = (buff_previous + buff_dist)/2,
                                SegNum = sf.streams.buffer$SegNum,
                                IP_class = sf.streams.buffer$IP_class,
                                DTB_m_HighValue = raster::extract(r.dtb.coarse, sf.streams.donut, fun=mean, na.rm=T, weight=T),
                                WTD_m_HighValue = raster::extract(r.wtd.coarse, sf.streams.donut, fun=mean, na.rm=T, weight=T))
      
    }
    
    if (p == pumps_plot[1] & b == 1){
      df.buff <- df.buff.dist
      df.dtb <- df.dtb.dist
      buff_list <- paste0(buff_previous, "-", buff_dist, " m")
    } else {
      df.buff <- rbind(df.buff, df.buff.dist)
      if (p == pumps_plot[1]) df.dtb <- rbind(df.dtb, df.dtb.dist)
      buff_list <- c(buff_list, paste0(buff_previous, "-", buff_dist, " m"))
    }
    
    # status update
    print(paste0("p ", p, " ", buff_dist, " complete"))
  }
  
}

write.csv(df.buff, file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleSegment_Depletion-SensitivityAnalysis_df.buff.csv"), quote=F, row.names=F)
write.csv(df.dtb, file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleSegment_Depletion-SensitivityAnalysis_df.dtb.csv"), quote=F, row.names=F)

df.buff <- read.csv(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleSegment_Depletion-SensitivityAnalysis_df.buff.csv"), stringsAsFactors=F)
df.dtb <- read.csv(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleSegment_Depletion-SensitivityAnalysis_df.dtb.csv"), stringsAsFactors=F)

# join df.buff and df.dtb
df.buff <-
  dplyr::left_join(df.buff, df.dtb) %>% 
  transform(well_in_alluvium = DTB_m_HighValue > WTD_m_HighValue)

# mean for each buffer distance
df.buff.mean <-
  df.buff %>% 
  group_by(time_days, buff_range_center, pump_factor) %>% 
  summarize(depletion_m3d_HighValue_mean = mean(depletion_m3d_HighValue)) %>% 
  transform(deriv1 = c(NA, diff(depletion_m3d_HighValue_mean, lag=1))) %>% 
  transform(deriv2 = c(NA, diff(deriv1, lag=1)))

## Figure: depletion as a function of buffer distance
p.dist <- 
  ggplot() + 
  geom_line(data=df.buff, aes(x=buff_range_center, y=depletion_m3d_HighValue, group=SegNum), alpha=0.3, color=col.cat.red) +
  geom_line(data=subset(df.buff, well_in_alluvium), aes(x=buff_range_center, y=depletion_m3d_HighValue, group=SegNum), color=col.cat.red) +
  geom_line(data=df.buff.mean, aes(x=buff_range_center, y=depletion_m3d_HighValue_mean), size=2, color=col.cat.blu) +
  facet_wrap(~time_days,
             nrow = 1,
             labeller = as_labeller(time_labels_df)) +
  scale_x_continuous(name=paste0("Distance from Stream [center of ", buff_interval, " m bin]")) +
  scale_y_continuous(name="Depletion from High\nPotential Streams [m\u00b3 d\u207b\u00b9]") +
  NULL