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

### residential water use - output from Navarro_Residential_03_DepletionBySegment.R
# define pumping rates
df.pump.res <- 
  file.path("results", "Residential_WaterUse.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(domain == "Healdsburg")

# load house locations to figure out total
sf.houses <- 
  file.path(dir.TNC, "Structures_Navarro_NAIP_2016", "Structures_Navarro_NAIP_2016.shp") %>% 
  sf::st_read() %>% 
  subset(Structure=="Res H")
sf.houses$HouseNum <- seq(1, dim(sf.houses)[1])

## load point of diversion to screen out surface water users
# points of diversion
sf.diversions <- 
  file.path(dir.TNC, "nav_pointsofdiversion_domestic", "nav_pointsofdiversion_domestic.shp") %>% 
  sf::st_read() %>% 
  subset(Beneficial == "Domestic") %>%  # domestic users only
  subset(POD_Status %in% c("Active", "Certified", "Claimed", "Licensed", "Permitted", "Registered"))  # remove cancelled, closed, inactive, rejected, or revoked

# find nearest house to each point
nearest_house <- 
  sf::st_nearest_feature(sf.diversions, sf.houses)
sf.houses$groundwater <- T
sf.houses$groundwater[nearest_house] <- F

sum(sf.houses$groundwater==F)

# there are some houses that are nearest to multiple POD; need to remove and re-do until we have 1 house per POD
sf.houses.2 <- subset(sf.houses, groundwater)
sf.diversion.2 <- sf.diversions[which(duplicated(nearest_house)), ]
nearest_house.2 <- 
  sf::st_nearest_feature(sf.diversion.2, sf.houses.2)
sf.houses$groundwater[sf.houses$HouseNum %in% sf.houses.2$HouseNum[nearest_house.2]] <- F

sum(sf.houses$groundwater==F)

sf.houses.3 <- subset(sf.houses, groundwater)
sf.diversion.3 <- sf.diversion.2[which(duplicated(nearest_house.2)), ]
nearest_house.3 <- 
  sf::st_nearest_feature(sf.diversion.3, sf.houses.3)
sf.houses$groundwater[sf.houses$HouseNum %in% sf.houses.3$HouseNum[nearest_house.3]] <- F

sum(sf.houses$groundwater==F) == dim(sf.diversions)[1]

## convert to cubic meters
df.pump.res$m3HouseDay <- df.pump.res$WaterUseMean_GalHouseDay*gal.to.m3
df.pump.res$m3HouseDay_std <- df.pump.res$WaterUseStd_GalHouseDay*gal.to.m3
df.pump.res$ResWaterUseSum_m3d <- df.pump.res$m3HouseDay*sum(sf.houses$groundwater)
df.pump.res$ResWaterUseMin_m3d <- (df.pump.res$m3HouseDay-df.pump.res$m3HouseDay_std)*sum(sf.houses$groundwater)
df.pump.res$ResWaterUseMax_m3d <- (df.pump.res$m3HouseDay+df.pump.res$m3HouseDay_std)*sum(sf.houses$groundwater)

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
  left_join(df.pump.res[,c("MonthNum", "ResWaterUseSum_m3d", "ResWaterUseMin_m3d", "ResWaterUseMax_m3d")], by="MonthNum")

p.waterUse.abs <-
  ggplot(df.plot, aes(x=Month)) +
  # plot cannabis
  geom_ribbon(aes(ymin=WaterUseMin_m3d_sum, ymax=WaterUseMax_m3d_sum, group=1), fill="black", alpha=0.15) +
  geom_line(aes(y=WaterUseMean_m3d_sum, group=1)) +
  geom_point(aes(y=WaterUseMean_m3d_sum, group=1)) +
  # plot residential
  geom_ribbon(aes(ymin=ResWaterUseMin_m3d, ymax=ResWaterUseMax_m3d, group=1), fill="black", alpha=0.15) +
  geom_line(aes(y=ResWaterUseSum_m3d, group=1), linetype="dashed") +
  geom_point(aes(y=ResWaterUseSum_m3d, group=1), shape=1) +
  # plot baseflow - ribbon shows 10% GW presumptive standard from Gleeson & Richter
  geom_ribbon(aes(ymin=baseflow_m3d_mean*0.9, ymax=baseflow_m3d_mean*1.1, group=1), fill=col.cat.blu, alpha=0.25) +
  geom_line(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  geom_point(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  # aesthetics
  annotation_logticks(sides = "l", color=col.gray) +
  scale_x_discrete(labels=seq(1,12), expand=c(0.025,0.025)) +
  scale_y_log10(name="Volumetric Flux [m\u00b3 d\u207b\u00b9]",
                expand=c(0.06,0.06),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))

p.waterUse.prc <-
  ggplot(df.plot, aes(x=Month)) +
  # plot cannabis
  geom_ribbon(aes(ymin=100*ResWaterUseMin_m3d/baseflow_m3d_mean, ymax=100*ResWaterUseMax_m3d/baseflow_m3d_mean, group=1), 
              fill="black", alpha=0.15) +
  geom_line(aes(y=100*WaterUseMean_m3d_sum/baseflow_m3d_mean, group=1)) +
  geom_point(aes(y=100*WaterUseMean_m3d_sum/baseflow_m3d_mean, group=1)) +
  # plot residential
  geom_ribbon(aes(ymin=100*WaterUseMin_m3d_sum/baseflow_m3d_mean, ymax=100*WaterUseMax_m3d_sum/baseflow_m3d_mean, group=1), 
              fill="black", alpha=0.15) +
  geom_line(aes(y=100*ResWaterUseSum_m3d/baseflow_m3d_mean, group=1), linetype="dashed") +
  geom_point(aes(y=100*ResWaterUseSum_m3d/baseflow_m3d_mean, group=1), shape=1) +
  scale_x_discrete(labels=seq(1,12), expand=c(0.025,0.025)) +
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
          label_x = 0.91,
          label_y = 0.99,
          align="v", nrow=2) %>% 
  save_plot(filename = file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_WaterUse-BaseflowComparison.png"), 
            plot = .,
            ncol = 1,
            nrow = 2,
            base_width = 95/25.4,
            base_height = 60/25.4)


# useful stats
df.plot[which.max(df.plot$WaterUseMean_m3d_sum), c("MonthNum", "WaterUseMean_m3d_sum")]
df.plot[which.max(df.plot$ResWaterUseSum_m3d), c("MonthNum", "ResWaterUseSum_m3d")]

df.plot$WaterUseMean_m3d_sum[9]/df.plot$baseflow_m3d_mean[9]
df.plot$ResWaterUseSum_m3d[9]/df.plot$baseflow_m3d_mean[9]

# total water use [m3/yr]
sum(df.plot$ResWaterUseSum_m3d*lubridate::days_in_month(df.plot$MonthNum))
sum(df.plot$WaterUseMean_m3d_sum*lubridate::days_in_month(df.plot$MonthNum))

sum(df.plot$ResWaterUseSum_m3d*lubridate::days_in_month(df.plot$MonthNum))+sum(df.plot$WaterUseMean_m3d_sum*lubridate::days_in_month(df.plot$MonthNum))

# traditional ag estimated to use 1825 acre-feet (=2251101 m3), of which 97% from surface water
2251101*0.97

sum(df.plot$ResWaterUseSum_m3d*lubridate::days_in_month(df.plot$MonthNum))/(2251101*0.97)
sum(df.plot$WaterUseMean_m3d_sum*lubridate::days_in_month(df.plot$MonthNum))/(2251101*0.97)
