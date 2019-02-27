## Navarro_Cannabis-Grows_03_SummarizeDepletion.R
#' This script is intended to take estimated depletion associated with existing cannabis grow locations
#' and summarize it along different levels: by stream segment, by ecological sensitivity, by watershed, etc.
#'
#' It requires output from Navarro_Cannabis-Grows_02_DepletionBySegment.R

source(file.path("src", "paths+packages.R"))

## load data

# habitat - output from Navarro_Cannabis_02_HabitatIntrinsicPotential.R
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

# depletion by segment associated with each well
df <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  left_join(df.habitat[,c("SegNum", "IP", "IP_class")], by="SegNum") %>% 
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

ggplot(df.grow.closest, aes(x=dist_m_closestStream, fill=WaterSource)) +
  geom_histogram(breaks=seq(0, 3000, 100)) +
  geom_vline(xintercept=c(dist.thres.sw,dist.thres.gw), color=col.gray) +
  scale_x_continuous(name="Distance to Closest Stream [m]") +
  scale_y_continuous(name="Number of Cultivation Sites") +
  scale_fill_manual(name="Assumed Water Source", 
                    values=c("Surface Water Only"=col.cat.blu, 
                             "Surface Water or Groundwater"=col.cat.org,
                             "Groundwater Only"=col.cat.grn)) +
  theme(legend.position=c(0.99,0.99),
        legend.justification=c(1,1))

# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>%
  sf::st_transform(crs.MODFLOW)

# grow locations from TNC (filtered and transformed in Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R)
sf.grows <-
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.shp"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

# calculate annual water use by grow location
df.pump <- read.csv(file.path(dir.TNC, "CannabisMonthlyWaterUse_WilsonEtAl.csv"), stringsAsFactors=F)
df.pump$Month <- factor(df.pump$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump$Setting <- factor(df.pump$Setting, levels=c("Outdoor", "Greenhouse"))
df.pump$MonthNum <- match(df.pump$Month, month.abb)
df.pump$MonthLengthDays <- lubridate::days_in_month(df.pump$MonthNum)
df.pump$m3PlantDay <- df.pump$MeanWaterUse_GalPlantDay*gal.to.m3
df.pump$m3PlantMonth <- df.pump$m3PlantDay*df.pump$MonthLengthDays
df.pump.yr <- 
  df.pump %>% 
  group_by(Setting) %>% 
  summarize(m3PlantYear = sum(m3PlantMonth))

sf.grows$waterUse_m3yr[sf.grows$greenhouse==1] <- 
  sf.grows$plants[sf.grows$greenhouse==1]*df.pump.yr$m3PlantYear[df.pump.yr$Setting=="Greenhouse"]
sf.grows$waterUse_m3yr[sf.grows$greenhouse==0] <- 
  sf.grows$plants[sf.grows$greenhouse==0]*df.pump.yr$m3PlantYear[df.pump.yr$Setting=="Outdoor"]

# load stream shapefile and join habitat data (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  left_join(df.habitat[,c("SegNum", "IP", "IP_class")], by=c("SegNum")) %>%
  sf::st_transform(crs.MODFLOW)

## summarize data
# depletion by segment
df.segment <-
  df %>% 
  group_by(time_yrs, SegNum) %>% 
  summarize(depletion_m3d_sum = sum(depletion_m3d))

# depletion by ecological sensitivity
df.habitat <- 
  df %>% 
  group_by(time_yrs, IP_class) %>% 
  summarize(depletion_m3d_sum = sum(depletion_m3d))

# depletion by Navarro or not
df.navarro <-
  df %>% 
  group_by(time_yrs, !is.na(IP_class)) %>% 
  summarize(depletion_m3d_sum = sum(depletion_m3d)) %>% 
  set_colnames(c("time_yrs", "navarro", "depletion_m3d_sum"))

# depletion by GrowNum
df.grow <-
  df %>% 
  group_by(time_yrs, GrowNum) %>% 
  summarize(depletion_m3d_sum = sum(depletion_m3d))

## plots
# depletion by segment
qplot(df.segment$depletion_m3d)

# depletion by ecological sensitivity
ggplot(df.habitat, aes(x=time_yrs, y=depletion_m3d_sum, color=IP_class)) +
  geom_point() +
  geom_line()

df.habitat %>% 
  group_by(time_yrs) %>% 
  summarize(depletion_m3d_allStreams = sum(depletion_m3d_sum)) %>% 
  left_join(df.habitat, ., by="time_yrs") %>% 
  subset(IP_class=="High") %>% 
  ggplot(aes(x=time_yrs, y=depletion_m3d_sum/depletion_m3d_allStreams)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(name="% of Streamflow Depletion from High Ecological Value Streams", 
                     labels=scales::percent)

# depletion within Navarro
ggplot(df.navarro, aes(x=time_yrs, y=depletion_m3d_sum, color=navarro)) +
  geom_point() +
  geom_line()

# depletion by GrowNum
df.grow %>% 
  left_join(sf.grows, by="GrowNum") %>% 
  ggplot() +
  geom_sf(aes(color=depletion_m3d_sum)) +
  facet_wrap(~time_yrs)

ggplot(sf.grows, aes(size=waterUse_m3yr)) +
  geom_sf()

df.grow %>% 
  left_join(sf.grows, by="GrowNum") %>% 
  ggplot() +
  geom_histogram(aes(x=depletion_m3d_sum)) +
  facet_wrap(~time_yrs)
