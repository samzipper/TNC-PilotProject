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
df.Q.mo <-
  summarize(group_by(df.Q.yr.mo, month),
            baseflow_m3d = mean(baseflow_m3d.mean),
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
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(time_yrs = time_days/365,
            depletion_m3s = depletion_m3d*86400)

# calculate distance to closest stream for each grow
df.grow.closest <- 
  df.grow %>% 
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
# calculate year and month
df.grow$year.dec <- df.grow$time_days/365 + 1
df.grow$year <- floor(df.grow$year.dec)
df.grow$month <- round((df.grow$year.dec - df.grow$year)*12)+1

# trim data frame
df.grow.analysis <-
  df.grow %>% 
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

set.seed(1)
for (iter in 1:n.iter){
  grows.gw <- c(grows.gw.only, base::sample(grows.sw.or.gw, size=n.sample))
  
  # summarize total depletion by year and month
  df.grow.depletion.i <-
    df.grow.analysis %>% 
    subset(GrowNum %in% grows.gw) %>% 
    group_by(year, month) %>% 
    summarize(depletion_m3d_Navarro = sum(depletion_m3d)) %>% 
    transform(iter = iter)
  
  if (iter==1){
    df.grow.depletion <- df.grow.depletion.i
  } else {
    df.grow.depletion <- rbind(df.grow.depletion, df.grow.depletion.i)
  }
  
}

# summarize to mean and std
df.grow.depletion.summary <-
  df.grow.depletion %>% 
  group_by(year, month) %>% 
  summarize(depletion_m3d_Navarro_mean = mean(depletion_m3d_Navarro),
            depletion_m3d_Navarro_std = sd(depletion_m3d_Navarro)) %>% 
  left_join(df.Q.mo, by="month") %>% 
  transform(WaterUser="Cannabis")

#### load and process residential data
# depletion by segment associated with each well - output from Navarro_Residential_03_DepletionBySegment.R
df.res <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Residential_03_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(time_yrs = time_days/365,
            depletion_m3s = depletion_m3d*86400)

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

### Figure: monthly baseflow and streamflow depletion
p.vol <- 
  ggplot() +
  geom_ribbon(data=df.grow.depletion.summary, aes(x=month, 
                                                  ymin=(depletion_m3d_Navarro_mean-depletion_m3d_Navarro_std), 
                                                  ymax=(depletion_m3d_Navarro_mean+depletion_m3d_Navarro_std),
                                                  fill=factor(year)),
              alpha = 0.25) +
  geom_line(data=df.grow.depletion.summary, aes(x=month, y=depletion_m3d_Navarro_mean, color=factor(year))) +
  geom_point(data=df.grow.depletion.summary, aes(x=month, y=depletion_m3d_Navarro_mean, color=factor(year))) +
  scale_x_continuous(name = "Month", breaks=seq(1,12)) +
  scale_color_manual(name = "Years of\nPumping", 
                     values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  scale_fill_manual(name = "Years of\nPumping", 
                    values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  scale_y_continuous(name = "Streamflow Depletion [m3/d]") +
  theme(legend.position=c(1,0),
        legend.justification=c(1,0),
        legend.background=element_blank()) +
  guides(color=F, fill=F) +
  NULL

p.prc <-
  ggplot(df.grow.depletion.summary, aes(x=month)) +
  geom_ribbon(aes(ymin=100*(depletion_m3d_Navarro_mean-depletion_m3d_Navarro_std)/baseflow_m3d, 
                  ymax=100*(depletion_m3d_Navarro_mean+depletion_m3d_Navarro_std)/baseflow_m3d,
                  fill=factor(year)), alpha=0.5) +
  geom_line(aes(y=100*depletion_m3d_Navarro_mean/baseflow_m3d, color=factor(year))) +
  geom_point(aes(y=100*depletion_m3d_Navarro_mean/baseflow_m3d, color=factor(year))) +
  scale_x_continuous(name = "Month", breaks=seq(1,12)) +
  scale_y_continuous(name="Streamflow Depletion\n[% of Monthly Baseflow]") +
  scale_fill_manual(name="Years of\nPumping", values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  scale_color_manual(name="Years of\nPumping", values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  theme(legend.position=c(0,0.93),
        legend.justification=c(0,1),
        legend.background=element_blank()) +
  #  guides(color=F, fill=F) +
  NULL

plot_grid(p.vol, 
          p.prc, 
          labels = c("(a)", "(b)"), 
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = 0.16,
          label_y = 0.99,
          align="v", nrow=2) %>% 
  save_plot(filename = file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWatershed_Depletion_Depletion+Baseflow.png"), 
            plot = .,
            ncol = 1,
            nrow = 2,
            base_width = 95/25.4,
            base_height = 60/25.4)

## SI Figure: how water source was selected
ggplot(df.grow.closest, aes(x=dist_m_closestStream, fill=WaterSource)) +
  geom_histogram(breaks=seq(0, 3000, 100)) +
  #  geom_vline(xintercept=c(dist.thres.sw,dist.thres.gw), color=col.gray) +
  scale_x_continuous(name="Distance to Closest Stream [m]") +
  scale_y_continuous(name="Number of Cultivation Sites") +
  scale_fill_manual(name="Assumed Water Source", 
                    values=c("Surface Water Only"=col.cat.blu, 
                             "Surface Water or Groundwater"=col.cat.org,
                             "Groundwater Only"=col.cat.grn)) +
  theme(legend.position=c(0.99,0.99),
        legend.justification=c(1,1)) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWatershed_Depletion_AllocatingWaterUse.png"),
         width = 95, height = 95, units="mm") +
  NULL

## SI Figure: Comparison of residential and cannabis depletion
df.comparison <-
  rbind(df.res.depletion.summary[,c("year", "month", "depletion_m3d_Navarro", "WaterUser")],
        setNames(df.grow.depletion.summary[,c("year", "month", "depletion_m3d_Navarro_mean", "WaterUser")], c("year", "month", "depletion_m3d_Navarro", "WaterUser")))

time_labels <- 
  setNames(paste0("Year ", yrs.plot),
           yrs.plot)

ggplot() +
  geom_line(data=df.comparison, aes(x=month, y=depletion_m3d_Navarro, color=factor(year), linetype=WaterUser)) +
  scale_x_continuous(name = "Month", breaks=seq(1,12)) +
  facet_wrap(~year, ncol=3, labeller=as_labeller(time_labels), scales="free_y") +
  scale_color_manual(name = "Years of Pumping", 
                     values=c(col.cat.grn, col.cat.org, col.cat.red), guide=F) +
  scale_linetype_manual(name = "Water User", 
                     values=c("Residential"="dashed", "Cannabis"="solid")) +
  scale_y_continuous(name = "Streamflow Depletion [m3/d]", 
                     limits=c(min(df.comparison$depletion_m3d_Navarro), max(df.comparison$depletion_m3d_Navarro))) +
  theme(legend.position="bottom",
        legend.background=element_blank()) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWatershed_Depletion_CompareResidential.png"),
         width = 190, height = 95, units="mm") +
  NULL
