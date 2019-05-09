## Figure_ScaleWatershed_Depletion.R

source(file.path("src", "paths+packages.R"))

## what years to plot? options: 1  2  3  4  5 10 15 20 30 50
yrs.plot <- c(1, 10, 50)

#### load and process streamflow data
## get data from USGS - output from script Navarro_StreamflowData.R
df.Q <- read.csv(file.path("results", "Navarro_StreamflowData.csv"), stringsAsFactors=F)
df.Q.info <- read.csv(file.path("results", "Navarro_StreamflowData_siteInfo.csv"), stringsAsFactors=F)
# df.Q colnames:
#   staid = station (char)
#   val = discharge [cfs] (numeric)
#   dates = date (Date)
#   qualcode = quality code

## data inspection and cleanup
if (sum(is.na(df.Q$val))>0) stop(paste0('no data: ', paste(df.Q$dates[is.na(df.Q$val)], collapse=", ")))

# gap-fill with linear interpolation
df.Q$val <- na.approx(df.Q$val, maxgap=14)

## unit conversions
# cfs to m3/d
df.Q$discharge_m3d <- df.Q$val*cfs.to.m3d

# calculate baseflow using Nathan & McMahon digital filter
bf <- BaseflowSeparation(df.Q$discharge_m3d)
df.Q$baseflow_m3d <- bf[,1]
df.Q$quickflow_m3d <- bf[,2]

# year and water year
df.Q$date <- ymd(df.Q$date)
df.Q$year <- year(df.Q$date)
df.Q$water.year <- year(df.Q$date+days(sum(days_in_month(c(10,11,12)))))
df.Q$month <- month(df.Q$date)
df.Q$DOY <- yday(df.Q$date)

# summarize by year and month
df.Q.yr.mo <-
  summarize(group_by(df.Q, water.year, month),
            discharge_m3d.mean = mean(discharge_m3d),
            discharge_m3d.cum = sum(discharge_m3d),
            baseflow_m3d.mean = mean(baseflow_m3d),
            baseflow_m3d.cum = sum(baseflow_m3d))

# summarize by month
# what years to average for comparison?
yrs.average <- seq(1999, 2018)  # last 20 years
df.Q.mo <-
  df.Q.yr.mo %>% 
  subset(water.year %in% yrs.average) %>% 
  group_by(month) %>% 
  summarize(baseflow_m3d = mean(baseflow_m3d.mean),
            baseflow_m3d.std = sd(baseflow_m3d.mean),
            baseflow_m3d.min = min(baseflow_m3d.mean)) %>% 
  transform(baseflow_m3d.ribbon.min = baseflow_m3d - baseflow_m3d.std,
            baseflow_m3d.ribbon.max = baseflow_m3d + baseflow_m3d.std)
df.Q.mo$baseflow_m3d.ribbon.min[df.Q.mo$baseflow_m3d.ribbon.min < df.Q.mo$baseflow_m3d.min] <- 
  df.Q.mo$baseflow_m3d.min[df.Q.mo$baseflow_m3d.ribbon.min < df.Q.mo$baseflow_m3d.min]

## habitat - output from Navarro_Cannabis_02_HabitatIntrinsicPotential.R - 
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

#### load and process cannabis data
# depletion by segment associated with each well - output from Navarro_Cannabis-Grows_02_DepletionBySegment.R
df.grow <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_03_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(time_yrs = time_days/365,
            depletion_m3s = depletion_m3d*86400)

# load well-stream pairs which has distance
df.pairs <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  dplyr::select(SegNum, GrowNum, dist_wellToStream_m)

df.grow <-
  left_join(df.grow, df.pairs, by=c("SegNum", "GrowNum"))

# calculate year and month
df.grow$year.dec <- df.grow$time_days/365 + 1
df.grow$year <- floor(df.grow$year.dec)
df.grow$month <- round((df.grow$year.dec - df.grow$year)*12)+1

# determine what grows use groundwater, from Chris' model
# grow locations shapefile
sf.grows <- 
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg")) %>% 
  sf::st_transform(crs.MODFLOW)

# summarize to year/month
df.grow.depletion.summary <-
  # trim data frame
  df.grow %>% 
  subset(SegNum %in% df.habitat$SegNum) %>% 
  subset(year %in% yrs.plot) %>% 
  subset(GrowNum %in% sf.grows$GrowNum[sf.grows$Well.rf.pred=="Yes"]) %>% 
  # summarize
  group_by(year, month) %>% 
  summarize(depletion_m3d_Navarro = sum(depletion_m3d)) %>% 
  # join with streamflow
  left_join(df.Q.mo, by="month") %>% 
  transform(WaterUser="Cannabis")

#### load and process residential data
# depletion by segment associated with each well - output from Navarro_Residential_03_DepletionBySegment.R
df.res <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Residential_03_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F)

# calculate year and month
df.res$year.dec <- df.res$time_days/365 + 1
df.res$year <- floor(df.res$year.dec)
df.res$month <- round((df.res$year.dec - df.res$year)*12)+1

# trim and summarize data frame
df.res.depletion.summary <-
  df.res %>% 
  subset(SegNum %in% df.habitat$SegNum) %>% 
  subset(year %in% yrs.plot) %>% 
  # sum for all segments in Navarro
  group_by(year, month) %>% 
  summarize(depletion_m3d_Navarro = sum(depletion_m3d)) %>% 
  left_join(df.Q.mo, by="month") %>% 
  transform(WaterUser="Residential")

## combine residential and cannabis
df.depletion.summary <- 
  rbind(df.grow.depletion.summary, df.res.depletion.summary)

### Figure: monthly baseflow and streamflow depletion
p.vol <- 
  ggplot(df.depletion.summary, aes(x=month, y=depletion_m3d_Navarro, 
                                   color=factor(year), shape=WaterUser, linetype=WaterUser)) +
  geom_line() + 
  geom_point() +
  scale_x_continuous(name = "Month", breaks=seq(1,12)) +
  scale_y_continuous(name = "Streamflow Depletion [m\u00b3 d\u207b\u00b9]") +
  scale_color_manual(name = "Years of\nPumping", 
                     values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  scale_linetype_manual(name = "Water User", 
                        values=c("Cannabis"="solid", "Residential"="dashed")) +
  scale_shape_manual(name = "Water User", 
                     values=c("Cannabis"=16, "Residential"=1)) +
  theme(legend.position=c(1,0),
        legend.justification=c(1,0),
        legend.background=element_blank()) +
  guides(color=F, fill=F, shape=F, linetype=F) +
  NULL

p.prc <- 
  ggplot(df.depletion.summary, aes(x=month, y=100*depletion_m3d_Navarro/baseflow_m3d, 
                                   color=factor(year), shape=WaterUser, linetype=WaterUser)) +
  geom_line() + 
  geom_point() +
  scale_x_continuous(name = "Month", breaks=seq(1,12)) +
  scale_y_continuous(name = "Streamflow Depletion [% of Baseflow]") +
  scale_color_manual(name = "Years of Pumping", 
                     values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  scale_linetype_manual(name = "Water User", 
                        values=c("Cannabis"="solid", "Residential"="dashed")) +
  scale_shape_manual(name = "Water User", 
                     values=c("Cannabis"=16, "Residential"=1)) +
  theme(legend.position=c(0,0.93),
        legend.justification=c(0,1),
        legend.box="vertical",
        legend.background=element_blank()) +
  guides(color = guide_legend(title.position="top", title.hjust=0.5, direction="horizontal"), 
         shape = guide_legend(title.position="top", title.hjust=0.5, direction="vertical"),
         linetype = guide_legend(title.position="top", title.hjust=0.5, direction="vertical")) +
  NULL

plot_grid(p.vol, 
          p.prc, 
          labels = c("(a)", "(b)"), 
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = 0.16,
          label_y = 0.99,
          rel_heights=c(1,1),
          align="v", nrow=2) %>% 
  save_plot(filename = file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWatershed_Depletion_Depletion+Baseflow.png"), 
            plot = .,
            ncol = 1,
            nrow = 2,
            base_width = 95/25.4,
            base_height = 70/25.4)