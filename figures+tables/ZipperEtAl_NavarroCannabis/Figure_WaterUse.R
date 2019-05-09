## Figure_WaterUse.R
#' This script is intended to make a figure showing monthly cannabis water use.

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

## figure out if predicted groundwater user
# grow locations shapefile
sf.grows <- 
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg")) %>% 
  sf::st_transform(crs.MODFLOW)

df.pump <- left_join(df.pump, sf.grows[,c("GrowNum", "Well.rf.pred")], by="GrowNum")

### residential water use
# define pumping rates: 250 gpm in winter, 500 gpm in summer
df.pump.res <- data.frame(
  Month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  MeanWaterUse_GalHouseDay = c(250, 250, 250, 250, 500, 500, 500, 500, 500, 500, 250, 250)
)
df.pump.res$Month <- factor(df.pump.res$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump.res$MonthNum <- match(df.pump.res$Month, month.abb)
df.pump.res$MonthLengthDays <- lubridate::days_in_month(df.pump.res$MonthNum)
df.pump.res$m3HouseDay <- df.pump.res$MeanWaterUse_GalHouseDay*gal.to.m3

# load house locations to figure out total
sf.houses <- 
  sf::st_read(file.path(dir.TNC, "Structures_Navarro_NAIP_2016", "Structures_Navarro_NAIP_2016.shp")) %>% 
  subset(Structure=="Res H")
df.pump.res$ResWaterUseSum_m3d <- df.pump.res$m3HouseDay*dim(sf.houses)[1]

#### SI figure: monthly water use at each site with a well
ggplot(subset(df.pump, Well.rf.pred=="Yes"), aes(x=Month)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_ribbon(aes(ymin=WaterUseMin_m3d, ymax=WaterUseMax_m3d, group=GrowNum), alpha=0.1) +
#  geom_line(aes(y=WaterUseMean_m3d, group=GrowNum)) +
  theme(legend.position=c(0.01, 0.99),
        legend.justification=c(0,1)) +
  scale_y_continuous(name="Estimated Water Use [m\u00b3 d\u207b\u00b9]") +
  scale_x_discrete() +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_WaterUse-IndividualParcels.png"),
         width = 95, height = 95, units="mm") +
  NULL

#### SI figure: annual groundwater use histogram
# summarize annual sum
df.pump %>% 
  subset(Well.rf.pred=="Yes") %>% 
  dplyr::group_by(GrowNum) %>% 
  dplyr::summarize(WaterUseSum_m3yr = sum(WaterUseSum_m3mo)) %>% 
  transform(WaterUseCut = base::cut(WaterUseSum_m3yr, 
                                    breaks=c(seq(0,1000,100), Inf), 
                                    labels = c("< 100", "100-200", "200-300", "300-400", "400-500", "500-600",
                                               "600-700", "700-800", "800-900", "900-1000", "> 1000"),
                                    include.lowest=T)) %>% 
  # plot
  ggplot(aes(x=WaterUseCut)) +
  geom_bar(fill=col.cat.grn) +
  scale_x_discrete(name="Estimated Groundwater Use [m\u00b3 yr\u207b\u00b9]") +
  scale_y_continuous(name = "Number of Parcels", limits=c(0, 103), expand=c(0,0)) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_WaterUse-AnnualHistogram.png"),
         width = 190, height = 95, units="mm") +
  NULL

#### figure: monthly streamflow and total cannabis groundwater use within Navarro
df.pump.month <-
  df.pump %>% 
  subset(Well.rf.pred=="Yes") %>% 
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
df.plot <- left_join(df.pump.month, df.streamflow.month, by=c("MonthNum"="month")) %>% 
  left_join(df.pump.res[,c("MonthNum", "ResWaterUseSum_m3d")], by="MonthNum")

p.waterUse.abs <-
  ggplot(df.plot, aes(x=Month)) +
  # plot cannabis
  geom_ribbon(aes(ymin=WaterUseMin_m3d_sum, ymax=WaterUseMax_m3d_sum, group=1), fill="black", alpha=0.15) +
  geom_line(aes(y=WaterUseMean_m3d_sum, group=1)) +
  geom_point(aes(y=WaterUseMean_m3d_sum, group=1)) +
  # plot residential
  geom_line(aes(y=ResWaterUseSum_m3d, group=1), linetype="dashed") +
  geom_point(aes(y=ResWaterUseSum_m3d, group=1), shape=1) +
  # plot streamflow - ribbon shows 10% GW presumptive standard from Gleeson & Richter
  geom_ribbon(aes(ymin=baseflow_m3d_mean*0.9, ymax=baseflow_m3d_mean*1.1, group=1), fill=col.cat.blu, alpha=0.25) +
  geom_line(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  geom_point(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  scale_x_discrete(labels=seq(1,12)) +
  scale_y_log10(name="Volumetric Flux [m\u00b3 d\u207b\u00b9]",
                expand=c(0.06,0.06))

p.waterUse.prc <-
  ggplot(df.plot, aes(x=Month)) +
  # plot cannabis
  geom_ribbon(aes(ymin=100*WaterUseMin_m3d_sum/baseflow_m3d_mean, ymax=100*WaterUseMax_m3d_sum/baseflow_m3d_mean, group=1), 
              fill="black", alpha=0.15) +
  geom_line(aes(y=100*WaterUseMean_m3d_sum/baseflow_m3d_mean, group=1)) +
  geom_point(aes(y=100*WaterUseMean_m3d_sum/baseflow_m3d_mean, group=1)) +
  # plot residential
  geom_line(aes(y=100*ResWaterUseSum_m3d/baseflow_m3d_mean, group=1), linetype="dashed") +
  geom_point(aes(y=100*ResWaterUseSum_m3d/baseflow_m3d_mean, group=1), shape=1) +
  scale_x_discrete(labels=seq(1,12)) +
  scale_y_continuous(name="Groundwater Use [% of Baseflow]")

plot_grid(p.waterUse.abs + 
            annotate("text", x=4.1, y=1.2e2, label="Cannabis", color="black", hjust=0, vjust=1) +
            annotate("text", x=1, y=4e3, label="Residential", color="black", hjust=0, vjust=1) +
            annotate("text", x=1, y=3e5, label="Baseflow", color=col.cat.blu, hjust=0), 
          p.waterUse.prc, 
          labels = c("(a)", "(b)"),
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = c(0.17, 0.17),
          label_y = c(1, 0.99),
          align="v", nrow=2) %>% 
  save_plot(filename = file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_WaterUse-BaseflowComparison.png"), 
            plot = .,
            ncol = 1,
            nrow = 2,
            base_width = 95/25.4,
            base_height = 60/25.4)
