## Figure_CannabisWaterUse.R
#' This script is intended to make a figure showing monthly cannabis water use for indoor and outdoor growers.

source(file.path("src", "paths+packages.R"))

## read in cannabis water use data - this is proprietary from TNC/Dillis, cannot be shared
# created using script Navarro_Cannabis-Grows_03_DepletionBySegment.R
df.pump <- read.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_03_PumpingRate.csv"), stringsAsFactors=F)
df.pump$WaterUseMean_m3d[df.pump$WaterUseMean_m3d < 0] <- 0

## set monthly factor
df.pump$Month <- factor(df.pump$MonthNum, labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

## create max and min ribbons for each site
df.pump$WaterUseMax_m3d <- df.pump$WaterUseMean_m3d + df.pump$WaterUseStd_m3d
df.pump$WaterUseMin_m3d <- df.pump$WaterUseMean_m3d - df.pump$WaterUseStd_m3d
df.pump$WaterUseMin_m3d[df.pump$WaterUseMin_m3d < 0] <- 0
df.pump$WaterUseSum_m3mo <- df.pump$WaterUseMean_m3d*lubridate::days_in_month(df.pump$MonthNum)

#### SI figure: monthly water use at each site
ggplot(df.pump, aes(x=Month)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_ribbon(aes(ymin=WaterUseMin_m3d, ymax=WaterUseMax_m3d, group=GrowNum), alpha=0.1) +
#  geom_line(aes(y=WaterUseMean_m3d, group=GrowNum)) +
  theme(legend.position=c(0.01, 0.99),
        legend.justification=c(0,1)) +
  scale_y_continuous(name="Estimated Water Use [m3/day]") +
  scale_x_discrete() +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_CannabisWaterUse.png"),
         width = 95, height = 95, units="mm") +
  NULL

#### SI figure: annual water use histogram
# summarize annual sum
df.pump %>% 
  dplyr::group_by(GrowNum) %>% 
  dplyr::summarize(WaterUseSum_m3yr = sum(WaterUseSum_m3mo)) %>% 
  transform(WaterUseCut = base::cut(WaterUseSum_m3yr, 
                                    breaks=c(seq(0,1000,100), Inf), 
                                    labels = c("< 100", "100-200", "200-300", "300-400", "400-500", "500-600",
                                               "600-700", "700-800", "800-900", "900-1000", "> 1000"),
                                    include.lowest=T)) %>% 
  # plot
  ggplot(aes(x=WaterUseCut)) +
  geom_bar() +
  scale_x_discrete(name="Estimated Annual Water Use [m3]") +
  scale_y_continuous(name = "Number of Parcels") +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_CannabisWaterUse-Annual.png"),
         width = 190, height = 95, units="mm") +
  NULL

#### figure: monthly streamflow and total cannabis water use within Navarro
df.pump.month <-
  df.pump %>% 
  group_by(Month, MonthNum) %>% 
  summarize(WaterUseMean_m3d_sum = sum(WaterUseMean_m3d),
            WaterUseMin_m3d_sum = sum(WaterUseMin_m3d),
            WaterUseMax_m3d_sum = sum(WaterUseMax_m3d))

# load streamflow - same as script 'Figure_StreamflowData.R'
df <- read.csv(file.path("results", "Navarro_StreamflowData.csv"), stringsAsFactors=F)
if (sum(is.na(df$val))>0) stop(paste0('no data: ', paste(df$dates[is.na(df$val)], collapse=", ")))

# 4 missing days, gap-fill with linear interpolation
df$val <- na.approx(df$val, maxgap=14)

# convert from cfs to m3d
cfs.to.m3d <- (0.3048^3)*86400
df$discharge_m3d <- df$val*cfs.to.m3d

# calculate baseflow using Nathan & McMahon digital filter
bf <- BaseflowSeparation(df$discharge_m3d)
df$baseflow_m3d <- bf[,1]
df$quickflow_m3d <- bf[,2]

# year and water year
df$date <- ymd(df$date)
df$year <- year(df$date)
df$water.year <- year(df$date+days(sum(days_in_month(c(10,11,12)))))
df$month <- month(df$date)
df$DOY <- yday(df$date)

# summarize to monthly means - last 20 years only
df.streamflow.month <-
  df %>% 
  subset(water.year >= 1999) %>% 
  group_by(month) %>% 
  summarize(discharge_m3d_mean = mean(discharge_m3d),
            baseflow_m3d_mean = mean(baseflow_m3d),
            baseflow_m3d_std = sd(baseflow_m3d)) %>% 
  transform(baseflow_m3d_max = baseflow_m3d_mean+baseflow_m3d_std,
            baseflow_m3d_min = baseflow_m3d_mean-baseflow_m3d_std)

# join and plot
df.plot <- left_join(df.pump.month, df.streamflow.month, by=c("MonthNum"="month"))

p.waterUse.abs <-
  ggplot(df.plot, aes(x=Month)) +
  geom_ribbon(aes(ymin=WaterUseMin_m3d_sum, ymax=WaterUseMax_m3d_sum, group=1), fill=col.cat.grn, alpha=0.25) +
  geom_line(aes(y=WaterUseMean_m3d_sum, group=1), color=col.cat.grn) +
  geom_point(aes(y=WaterUseMean_m3d_sum, group=1), color=col.cat.grn) +
  # plot streamflow - ribbon shows 10% GW presumptive standard from Gleeson & Richter
  geom_ribbon(aes(ymin=baseflow_m3d_mean*0.9, ymax=baseflow_m3d_mean*1.1, group=1), fill=col.cat.blu, alpha=0.25) +
  geom_line(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  geom_point(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  scale_x_discrete(labels=seq(1,12)) +
  scale_y_log10(name="Volumetric Flux [m3/d]")

p.waterUse.prc <-
  ggplot(df.plot, aes(x=Month)) +
  geom_ribbon(aes(ymin=100*WaterUseMin_m3d_sum/baseflow_m3d_mean, ymax=100*WaterUseMax_m3d_sum/baseflow_m3d_mean, group=1), 
              fill=col.cat.grn, alpha=0.25) +
  geom_line(aes(y=100*WaterUseMean_m3d_sum/baseflow_m3d_mean, group=1), color=col.cat.grn) +
  geom_point(aes(y=100*WaterUseMean_m3d_sum/baseflow_m3d_mean, group=1), color=col.cat.grn) +
  scale_x_discrete(labels=seq(1,12)) +
  scale_y_continuous(name="Cannabis Water Use\n[% of Monthly Baseflow]")

plot_grid(p.waterUse.abs + 
            annotate("text", x=8, y=5e2, label="Cannabis\nWater Use", color=col.cat.grn, hjust=0.5, vjust=1) +
            annotate("text", x=1, y=3e5, label="Baseflow", color=col.cat.blu, hjust=0), 
          p.waterUse.prc, 
          labels = c("(a)", "(b)"),
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = c(0.02, 0.1),
          align="v", nrow=2) %>% 
  save_plot(filename = file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_CannabisWaterUse+Streamflow.png"), 
            plot = .,
            ncol = 1,
            nrow = 2,
            base_width = 95/25.4,
            base_height = 60/25.4)
