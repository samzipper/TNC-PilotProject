## Figure_ScaleWellt_Depletion.R
#' This script is intended to relate segment-level and well-level depletion
#' Very similar to Figure_Habitat+GrowLocations.R

source(file.path("src", "paths+packages.R"))

## load data from Navarro_Cannabis_HabitatIntrinsicPotential.R
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  subset(TermnlP == outlet.TerminalPa) %>% 
  left_join(read.csv(file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv")),
            by=c("SegNum"))

# grow locations shapefile
sf.grows <- 
  sf::st_read(file.path(dir.TNC, "nav_cannabis_ucbtnc", "nav_cannabis_ucbtnc.shp")) %>% 
  subset(year==16) %>%                  # 2016 data only
  st_zm(drop = TRUE, what = "ZM") %>%   # drop Z dimension from geometry
  sf::st_transform(crs.MODFLOW) %>% 
  dplyr::rename(GrowNum = FID_allgro)

# rasters
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))
r.wte <- raster(paste0(dir.gis, "Navarro_Cannabis_WTE_30m.tif"))

## calculate saturated thickness of alluvial materials are saturated or not
r.wtd <- r.dem.30m - r.wte
r.alluvial.sat.thickness <- r.dtb - r.wtd
r.alluvial.sat.thickness[r.alluvial.sat.thickness < r.dtb] <- 0
r.alluvial.sat.thickness[r.alluvial.sat.thickness > r.dtb] <- r.dtb[r.alluvial.sat.thickness > r.dtb]

## extract some potentially relevant data
sf.grows$elev_m <- raster::extract(r.dem.30m, sf.grows)  # elevation of that grid cell [m]
sf.grows$dtb_m <- raster::extract(r.dtb, sf.grows)       # depth to bedrock 

# predict top of screen based on linear relationship with dtb_m
# see script Navarro_WellCompletionReports.R
coef_m <- -3.02
coef_b <- 101.11
sf.grows$well_screenTop_depth_m <- coef_m * sf.grows$dtb_m + coef_b   # depth to top of well screen [m]
sf.grows$well_screenTop_elev_m <- sf.grows$elev_m - sf.grows$well_screenTop_depth_m  # elevation of top of well screen [m]


## read in cannabis water use data - this is proprietary from TNC, cannot be shared (until Wilson et al paper published)
df.pump <- read.csv(file.path(dir.TNC, "CannabisMonthlyWaterUse_WilsonEtAl.csv"), stringsAsFactors=F)

## set monthly factor
df.pump$Month <- factor(df.pump$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump$Setting <- factor(df.pump$Setting, levels=c("Outdoor", "Greenhouse"))

# calculate total water use
df.pump.yr <- 
  df.pump %>% 
  transform(WaterUse_m3PlantDay = MeanWaterUse_GalPlantDay*gal.to.m3,
            MonthDays = days_in_month(as.numeric(Month))) %>% 
  group_by(Setting) %>% 
  summarize(WaterUse_m3PlantYr = sum(WaterUse_m3PlantDay*MonthDays))

# calculate water use for each grow
sf.grows$WaterUse_m3y <- NaN
sf.grows$WaterUse_m3y[sf.grows$greenhouse==0] <- df.pump.yr$WaterUse_m3PlantYr[df.pump.yr$Setting=="Outdoor"]*sf.grows$plants[sf.grows$greenhouse==0]
sf.grows$WaterUse_m3y[sf.grows$greenhouse==1] <- df.pump.yr$WaterUse_m3PlantYr[df.pump.yr$Setting=="Greenhouse"]*sf.grows$plants[sf.grows$greenhouse==1]

# domain boundary shapefile
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>% 
  sf::st_transform(crs.MODFLOW)

sf.subbasin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU12_Navarro.shp")) %>% 
  sf::st_transform(crs.MODFLOW)
sf.subbasin$area_m2 <- as.numeric(sf::st_area(sf.subbasin))

#### load and process cannabis data
# depletion by segment associated with each well - output from Navarro_Cannabis-Grows_02_DepletionBySegment.R
df <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(time_yrs = time_days/365,
            depletion_m3s = depletion_m3d*86400)

# calculate distance to closest stream for each grow
df.grow.closest <- 
  df %>% 
  group_by(GrowNum) %>% 
  summarize(dist_m_closestStream = min(dist_wellToStream_m))

quantile(df.grow.closest$dist_m_closestStream, (1-0.73))
sum(df.grow.closest$dist_m_closestStream > 500)/length(df.grow.closest$dist_m_closestStream)
sum(df.grow.closest$dist_m_closestStream < 100)/length(df.grow.closest$dist_m_closestStream)

# category for water source
dist.thres.sw <- 100
dist.thres.gw <- 500
df.grow.closest$WaterSource <- "Surface Water Only"
df.grow.closest$WaterSource[df.grow.closest$dist_m_closestStream > dist.thres.sw & df.grow.closest$dist_m_closestStream < dist.thres.gw] <- "Surface Water or Groundwater"
df.grow.closest$WaterSource[df.grow.closest$dist_m_closestStream > dist.thres.gw] <- "Groundwater Only"
df.grow.closest$WaterSource <- factor(df.grow.closest$WaterSource, levels=c("Surface Water Only", 
                                                                            "Surface Water or Groundwater",
                                                                            "Groundwater Only"))

## first: do some data trimming to get df down to a more manageable size
# habitat - output from Navarro_Cannabis_02_HabitatIntrinsicPotential.R - 
#   used to figure out what stream reaches are in Navarro
df.habitat <- 
  file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  melt(id=c("SegNum"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3]) %>% 
  subset(species == "Coho" & metric == "mean") %>%  # focus on Coho based on conversation with Jen - most sensitive to habitat conditions
  transform(IP_class = addNA(cut(IP, 
                                 breaks=c(0,0.7,1), 
                                 labels=c("Low", "High"),
                                 include.lowest=T)))

# calculate year and month
df$year.dec <- df$time_days/365 + 1
df$year <- floor(df$year.dec)
df$month <- round((df$year.dec - df$year)*12)+1

# what years to plot? options: 1  2  3  4  5 10 15 20 30 50
yrs.plot <- c(1, 10, 50)

# trim data frame
df.analysis <-
  df %>% 
  subset(SegNum %in% df.habitat$SegNum) %>% 
  subset(year %in% yrs.plot) %>% 
  left_join(df.grow.closest, by="GrowNum") %>% 
  subset(WaterSource != "Surface Water Only")

## iterate many times to calculate depletion by year and month 
## for different combinations of pumping wells
n.iter <- 100

# figure out how many to sample
n.grows <- length(unique(df.grow.closest$GrowNum))
prc.gw <- 0.73  # 73% of grow sites in Mendocino County use GW (from Dillis)
n.gw <- round(n.grows*prc.gw)
grows.gw.only <- df.grow.closest$GrowNum[df.grow.closest$WaterSource=="Groundwater Only"]
grows.sw.or.gw <- df.grow.closest$GrowNum[df.grow.closest$WaterSource=="Surface Water or Groundwater"]
n.sample <- n.gw - length(grows.gw.only)  # need to randomly select n.sample from the "Surface Water or Groundwater" category

# what month(s) to plot?
mo.plot <- 9

set.seed(1)
for (iter in 1:n.iter){
  grows.gw <- c(grows.gw.only, base::sample(grows.sw.or.gw, size=n.sample))
  
  # summarize total depletion by year, month, and grow
  df.depletion.total.i <- 
    df.analysis %>% 
    subset(GrowNum %in% grows.gw & month %in% mo.plot) %>% 
    group_by(year, month) %>% 
    summarize(depletion_m3d_Navarro = sum(depletion_m3d))
  
  df.depletion.i <-
    df.analysis %>% 
    subset(GrowNum %in% grows.gw & month %in% mo.plot) %>% 
    group_by(year, month, GrowNum) %>% 
    summarize(depletion_m3d_grow = sum(depletion_m3d)) %>% 
    left_join(df.depletion.total.i, by=c("year", "month")) %>% 
    transform(iter = iter,
              depletion_prc_grow = depletion_m3d_grow/depletion_m3d_Navarro)
  
  df.depletion.i <- 
    df.depletion.i[order(df.depletion.i$year, df.depletion.i$depletion_prc_grow), ] %>% 
    # rank from smallest contribution to depletion (1) to largest contribution to depletion (n.gw) within each year
    transform(depletion_rank = rep(seq(n.gw, 1, -1), length(yrs.plot))) %>% 
    # calculate rank and depletion normalized 0-1
    transform(depletion_rank_prc = depletion_rank/n.gw,
              depletion_prc_cum = ave(depletion_prc_grow, year, FUN=cumsum))

  if (iter==1){
    df.depletion <- df.depletion.i
  } else {
    df.depletion <- rbind(df.depletion, df.depletion.i)
  }
  
  print(paste0("iter ", iter, " complete"))
  
}

## summarize to mean and std based on rank
df.depletion.rank.summary <-
  df.depletion %>% 
  group_by(year, month, depletion_rank) %>% 
  summarize(depletion_rank_prc = mean(depletion_rank_prc),
            depletion_prc_cum_mean = mean(depletion_prc_cum),
            depletion_prc_cum_std = sd(depletion_prc_cum))

# for each year, figure out what rank corresponds to 50% depletion
df.lines <- 
  df.depletion.rank.summary %>% 
  filter(abs(depletion_prc_cum_mean-0.5) == min(abs(depletion_prc_cum_mean-0.5)))

## summarize to mean and std based on GrowNum
df.depletion.rank.summary <-
  df.depletion %>% 
  group_by(year, month, GrowNum) %>% 
  summarize(depletion_m3d_mean = mean(depletion_m3d_grow),
            depletion_m3d_std = sd(depletion_m3d_grow),
            depletion_prc_mean = mean(depletion_prc_grow),
            depletion_prc_std = sd(depletion_prc_grow),
            depletion_rank_mean = mean(depletion_rank),
            depletion_rank_std = mean(depletion_rank))

sf.grows.depletion <- 
  sf.grows %>% 
  dplyr::select(geometry, GrowNum, greenhouse, growsize, plants, WaterUse_m3y, 
                elev_m, dtb_m, well_screenTop_depth_m, well_screenTop_elev_m) %>% 
  left_join(df.grow.closest, by="GrowNum")

for (yr in yrs.plot){
  sf.grows.depletion.yr <- 
    left_join(sf.grows.depletion, subset(df.depletion.rank.summary, year==yr), by="GrowNum") %>% 
    replace_na(list(year = yr, 
                    depletion_m3d_mean = 0,
                    depletion_m3d_std = 0,
                    depletion_prc_mean = 0,
                    depletion_prc_std = 0,
                    depletion_rank_mean = 0,
                    depletion_rank_std = 0))
  
  if (yr==yrs.plot[1]){
    sf.grows.depletion.all <- sf.grows.depletion.yr
  } else {
    sf.grows.depletion.all <- rbind(sf.grows.depletion.all, 
                                    sf.grows.depletion.yr)
  }
}


# maps with specific grow sites - these are for exploration and 
#    cannot be included in repository or publication
ggplot(sf.grows.depletion.all) +
  geom_sf(aes(color=depletion_prc_mean, size=depletion_prc_mean), alpha=0.5) +
  facet_wrap(~year) + 
  scale_color_viridis() +
  theme(legend.position="bottom")


ggplot() + 
  geom_sf(data=sf.basin, fill=NA, color=col.gray) +
  geom_sf(data=sf.grows, aes(color=WaterSource), alpha=0.5) +
  geom_sf(data=sf.streams, color="blue")

## calculate importance of various well attributes
df.grows.depletion <- 
  sf.grows.depletion.all %>% 
  subset(WaterSource != "Surface Water Only") %>% 
  as.data.frame()

ggplot(df.grows.depletion, aes(y=depletion_prc_mean, x=dist_m_closestStream, color=factor(year))) +
  geom_point(alpha=0.5) +
  stat_smooth(method="loess")

set.seed(1)
for (iter in 1:n.iter){
  grows.gw <- c(grows.gw.only, base::sample(grows.sw.or.gw, size=n.sample))
  
  # aggregate data
  df.depletion.grows.i <-
    df.analysis %>% 
    subset(GrowNum %in% grows.gw & month %in% mo.plot) %>% 
    group_by(year, month, GrowNum) %>% 
    summarize(depletion_m3d_grow = sum(depletion_m3d)) %>% 
    left_join(df.grows.depletion, by=c("year", "month", "GrowNum"))
  
  # build statistical models and extract R2
  for (yr in yrs.plot){
    lm.full <- lm(depletion_m3d_mean ~ dist_m_closestStream + WaterUse_m3y + well_screenTop_depth_m, 
                  data=subset(df.depletion.grows.i, year==yr & depletion_prc_mean > 0.001))
    
    lm.drop.dist <- lm(depletion_m3d_mean ~ WaterUse_m3y + well_screenTop_depth_m, 
                       data=subset(df.depletion.grows.i, year==yr & depletion_prc_mean > 0.001))
    
    lm.drop.water <- lm(depletion_m3d_mean ~ dist_m_closestStream + well_screenTop_depth_m, 
                        data=subset(df.depletion.grows.i, year==yr & depletion_prc_mean > 0.001))
    
    lm.drop.screendepth <- lm(depletion_m3d_mean ~ dist_m_closestStream + WaterUse_m3y, 
                              data=subset(df.depletion.grows.i, year==yr & depletion_prc_mean > 0.001))
    
    if (iter==1 & yr==yrs.plot[1]){
      df.fit <- data.frame(iter = iter, 
                           year = yr,
                           adj.R2.full = summary(lm.full)$adj.r.squared,
                           adj.R2.drop.dist = summary(lm.drop.dist)$adj.r.squared,
                           adj.R2.drop.water = summary(lm.drop.water)$adj.r.squared,
                           adj.R2.drop.screendepth = summary(lm.drop.screendepth)$adj.r.squared)
    } else {
      df.fit <- rbind(df.fit, 
                      data.frame(iter = iter, 
                                 year = yr,
                                 adj.R2.full = summary(lm.full)$adj.r.squared,
                                 adj.R2.drop.dist = summary(lm.drop.dist)$adj.r.squared,
                                 adj.R2.drop.water = summary(lm.drop.water)$adj.r.squared,
                                 adj.R2.drop.screendepth = summary(lm.drop.screendepth)$adj.r.squared))
    }
  }
  print(paste0("iter ", iter, " complete"))
}

# calculate change in R2 for each variable
df.fit$del.adj.R2.drop.dist <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.dist)/df.fit$adj.R2.full
df.fit$del.adj.R2.drop.water <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.water)/df.fit$adj.R2.full
df.fit$del.adj.R2.drop.screendepth <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.screendepth)/df.fit$adj.R2.full

df.fit.summary <- 
  df.fit %>% 
  group_by(year) %>% 
  summarize(adj.R2.mean = mean(adj.R2.full, na.rm=T),
            adj.R2.std = sd(adj.R2.full, na.rm=T),
            change.dist.mean = mean(del.adj.R2.drop.dist, na.rm=T),
            change.dist.std = sd(del.adj.R2.drop.dist, na.rm=T),
            change.water.mean = mean(del.adj.R2.drop.water, na.rm=T),
            change.water.std = sd(del.adj.R2.drop.water, na.rm=T),
            change.screen.mean = mean(del.adj.R2.drop.screendepth, na.rm=T),
            change.screen.std = sd(del.adj.R2.drop.screendepth, na.rm=T)) %>% 
  melt(id=c("year", "adj.R2.mean", "adj.R2.std")) %>% 
  transform(var = str_split_fixed(variable, pattern="[.]", n=3)[,2],
            metric = str_split_fixed(variable, pattern="[.]", n=3)[,3]) %>% 
  dplyr::select(year, value, var, metric) %>% 
  dcast(year + var ~ metric) %>% 
  transform(lab = factor(var, levels=c("water", "dist", "screen"), 
                         labels=c("Water Use", "Distance to\nClosest Stream", "Screen Depth")))

## Figure: cumulative distribution function
ggplot(df.depletion.rank.summary) +
  geom_line(aes(x=depletion_rank_prc, y=depletion_prc_cum_mean, color=factor(year))) + 
  geom_ribbon(aes(x=depletion_rank_prc, 
                  ymin=(depletion_prc_cum_mean-depletion_prc_cum_std),
                  ymax=(depletion_prc_cum_mean+depletion_prc_cum_std),
                  fill = factor(year)),
              alpha = 0.25) +
  annotate("segment", x=-Inf, xend=df.lines$depletion_rank_prc, 
           y=0.5, yend=0.5,
           color=c(col.cat.grn, col.cat.org, col.cat.red),
           linetype = "dashed") + 
  annotate("segment", x=df.lines$depletion_rank_prc, xend=df.lines$depletion_rank_prc, 
           y=-Inf, yend=0.5,
           color=c(col.cat.grn, col.cat.org, col.cat.red),
           linetype = "dashed") + 
  scale_x_continuous(name = "Percent of Wells", labels = scales::percent, 
                     limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(name = "Percent of Total Streamflow Depletion", labels = scales::percent, 
                     limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name = "Years of\nPumping", 
                     values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  scale_fill_manual(name = "Years of\nPumping", 
                    values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background=element_blank(),
        plot.margin = unit(c(1.5,4.5,0,0), "mm")) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_PercentCDF.png"),
         width=95, height=85, units="mm")

## Figure: importance of different variables
ggplot(df.fit.summary, aes(x=lab, y=mean, fill=factor(year))) + 
  geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0, color=col.gray) +
  geom_errorbar(aes(ymin=(mean-std), 
                    ymax=(mean+std)),
                width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(name="% Reduction in Adjusted R2", 
                     labels=scales::percent) +
  scale_x_discrete(name="Variable Dropped from Model") +
  scale_fill_manual(name = "Years of\nPumping", 
                    values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background=element_blank()) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_VariableImportance.png"),
         width=95, height=85, units="mm")

df.fit %>% 
  group_by(year) %>% 
  summarize(adj.R2.mean = mean(adj.R2.full, na.rm=T),
            adj.R2.std = sd(adj.R2.full, na.rm=T))

# ### analysis by HUC12 subbasin
# # aggregate number of plants by subbasin
# sf.subbasins.plants <-
#   sf::st_join(sf.grows, sf.subbasin, join=st_within) %>% 
#   subset(HUC12 != 180101080407) %>%   # grows within this HUC12 were not mapped, but 2 dots are inside it anyways; due to projection rounding?
#   group_by(HUC12) %>% 
#   summarize(TotalPlants = sum(plants),
#             TotalGrowsize = sum(growsize),
#             NumberGrows = sum(is.finite(plants)),
#             TotalWaterUse_m3y = sum(WaterUse_m3y)) %>% 
#   as.data.frame() %>% 
#   left_join(sf.subbasin, ., by="HUC12") %>% 
#   transform(TotalWaterUse_m3y.m2 = TotalWaterUse_m3y/area_m2)
# 
# # aggregate streams by subbasin
# sf.subbasins.streams <-
#   sf.streams %>% 
#   sf::st_transform(st_crs(sf.subbasin)) %>% 
#   sf::st_join(sf.subbasin, join=st_intersects) %>% 
#   group_by(HUC12) %>% 
#   summarize(TotalStreamLength_m = sum(segmentLength_m)) %>% 
#   as.data.frame() %>% 
#   dplyr::select(HUC12, TotalStreamLength_m) %>% 
#   left_join(sf.subbasin, ., by="HUC12")
# 
# sf.subbasins.streams <-
#   sf.streams %>% 
#   sf::st_transform(st_crs(sf.subbasin)) %>% 
#   sf::st_join(sf.subbasin, join=st_intersects) %>% 
#   subset(IP_class=="High") %>%  
#   group_by(HUC12) %>% 
#   summarize(HighValueStreamLength_m = sum(segmentLength_m)) %>% 
#   as.data.frame() %>% 
#   dplyr::select(HUC12, HighValueStreamLength_m) %>% 
#   left_join(sf.subbasins.streams, ., by="HUC12") %>% 
#   replace_na(list(HighValueStreamLength_m = 0)) %>% 
#   transform(HighValueStreamLength_prc = HighValueStreamLength_m/TotalStreamLength_m) 
# 
# ## aggregate streamflow depletion by subbasin
# # summarize by segment
# df.depletion.seg <- 
#   df.depletion %>% 
#   subset(SegNum %in% sf.streams$SegNum) %>% 
#   group_by(time_days, SegNum) %>% 
#   summarize(TotalDepletion_m3d = sum(depletion_m3d)) %>% 
#   right_join(sf.streams, ., by="SegNum") %>% 
#   sf::st_transform(st_crs(sf.subbasin)) %>% 
#   sf::st_join(sf.subbasin, join=st_intersects) %>% 
#   subset(HUC12 != 180101080407) %>%   # grows within this HUC12 were not mapped, but 2 dots are inside it anyways; due to projection rounding?
#   group_by(HUC12, time_days) %>% 
#   summarize(TotalDepletion_m3d_HUC12 = sum(TotalDepletion_m3d)) %>% 
#   as.data.frame()
# 
# ### collect data to plot
# df.subbasin.plot <-
#   sf.subbasins.plants[,c("HUC12", "TotalWaterUse_m3y", "TotalWaterUse_m3y.m2")] %>% 
#   as.data.frame() %>% 
#   left_join(sf.subbasins.streams[,c("HUC12", "HighValueStreamLength_prc")], by="HUC12") %>% 
#   right_join(df.depletion.seg[,c("HUC12", "time_days", "TotalDepletion_m3d_HUC12")], by="HUC12")
# 
# 
# p.WaterUse <- 
#   ggplot(subset(df.subbasin.plot, time_days==365), aes(y=HighValueStreamLength_prc, x=TotalWaterUse_m3y.m2)) +
#   geom_point() +
#   stat_smooth(method="lm") +
#   scale_x_continuous(name="Total Annual Cannabis Water Use [m3/m2]") +
#   scale_y_continuous(name="Percent of Stream Segment Length\nwith High Habitat Potential")
# 
# p.Depletion <- 
#   ggplot(df.subbasin.plot, aes(y=HighValueStreamLength_prc, x=TotalDepletion_m3d_HUC12, color=factor(time_days))) +
#   geom_point() +
#   stat_smooth(method="lm") +
#   scale_x_continuous(name="Total Streamflow Depletion [m3/d]") +
#   scale_y_continuous(name="Percent of Stream Segment Length\nwith High Habitat Potential")
# 
# plot_grid(p.WaterUse,
#           p.Depletion,
#           ncol=2,
#           align="tb") %>% 
#   save_plot(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_HUC12.png"),
#             plot = .,
#             ncol = 2,
#             base_width=95/25.4,
#             base_height=95/25.4)
