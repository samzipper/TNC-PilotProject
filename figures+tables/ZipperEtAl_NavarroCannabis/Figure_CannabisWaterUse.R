## Figure_CannabisWaterUse.R
#' This script is intended to make a figure showing monthly cannabis water use for indoor and outdoor growers.

source(file.path("src", "paths+packages.R"))

## read in cannabis water use data - this is proprietary from TNC, cannot be shared (until Wilson et al paper published)
df.pump <- read.csv(file.path(dir.TNC, "CannabisMonthlyWaterUse_WilsonEtAl.csv"), stringsAsFactors=F)

## set monthly factor
df.pump$Month <- factor(df.pump$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump$Setting <- factor(df.pump$Setting, levels=c("Outdoor", "Greenhouse"))

## gallons to liters conversion factor
gal.to.L <- 3.78541
gal.to.m3 <- gal.to.L/1000

#### SI figure: monthly water use per plant
ggplot(df.pump, aes(x=Month)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_ribbon(aes(ymin=MinWaterUse_GalPlantDay*gal.to.L, ymax=MaxWaterUse_GalPlantDay*gal.to.L, 
                  fill=Setting, group=Setting), alpha=0.25) +
  geom_point(aes(y=MeanWaterUse_GalPlantDay*gal.to.L, color=Setting, group=Setting)) +
  geom_line(aes(y=MeanWaterUse_GalPlantDay*gal.to.L, color=Setting, group=Setting)) +
  theme(legend.position=c(0.01, 0.99),
        legend.justification=c(0,1)) +
  scale_y_continuous(name="Reported Water Use [L/plant/day]") +
  scale_x_discrete(labels=seq(1,12)) +
  scale_color_manual(name="Cultivation Setting", values=c("Greenhouse"=col.cat.red, "Outdoor"=col.cat.blu)) +
  scale_fill_manual(name="Cultivation Setting", values=c("Greenhouse"=col.cat.red, "Outdoor"=col.cat.blu)) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_CannabisWaterUse.png"),
         width = 95, height = 95, units="mm") +
  NULL

#### figure: monthly streamflow and total cannabis water use within Navarro
# grow locations from TNC (filtered and transformed in Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R)
sf.grows <-
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.shp"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

# calculate total number of plants: greenhouse and outdoor
df.plants <-
  data.frame(Setting = factor(c("Outdoor", "Greenhouse")),
             TotalPlants = c(sum(sf.grows$plants[sf.grows$greenhouse==0]),
                             sum(sf.grows$plants[sf.grows$greenhouse==1])))

# calculate total water use
df.pump <- 
  df.pump %>% 
  left_join(df.plants, by="Setting") %>% 
  transform(WaterUse_m3d = MeanWaterUse_GalPlantDay*TotalPlants*gal.to.m3)

df.pump.month <-
  df.pump %>% 
  group_by(Month) %>% 
  summarize(WaterUse_m3d_sum = sum(WaterUse_m3d))

# load streamflow - same as script 'Figure_StreamflowData.R'
df <- read.csv(file.path("results", "Navarro_StreamflowData.csv"), stringsAsFactors=F)
if (sum(is.na(df$val))>0) stop(paste0('no data: ', paste(df$dates[is.na(df$val)], collapse=", ")))

# 3 missing days, gap-fill with linear interpolation
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

# summarize to monthly means
df.streamflow.month <-
  df %>% 
  group_by(month) %>% 
  summarize(discharge_m3d_mean = mean(discharge_m3d),
            baseflow_m3d_mean = mean(baseflow_m3d))

# join and plot
df.pump.month$month <- match(df.pump.month$Month, month.abb)
df.pump.month <- left_join(df.pump.month, df.streamflow.month, by=c("month"))

p.waterUse.abs <-
  ggplot(df.pump.month, aes(x=Month)) +
  geom_line(aes(y=WaterUse_m3d_sum, group=1), color=col.cat.grn) +
  geom_point(aes(y=WaterUse_m3d_sum, group=1), color=col.cat.grn) +
  geom_line(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  geom_point(aes(y=baseflow_m3d_mean, group=2), color=col.cat.blu) +
  scale_x_discrete(labels=seq(1,12)) +
  scale_y_log10(name="Volumetric Flux [m3/d]")

p.waterUse.prc <-
  ggplot(df.pump.month, aes(x=Month, y=100*WaterUse_m3d_sum/baseflow_m3d_mean, group=1)) +
  geom_line(color=col.cat.grn) +
  geom_point(color=col.cat.grn) +
  scale_x_discrete(labels=seq(1,12)) +
  scale_y_continuous(name="Cannabis Water Use\n[% of Monthly Baseflow]")

plot_grid(p.waterUse.abs + 
            annotate("text", x=4, y=7e1, label="Cannabis\nWater Use", color=col.cat.grn, hjust=0.5, vjust=1) +
            annotate("text", x=1, y=3e5, label="Baseflow", color=col.cat.blu, hjust=0), 
          p.waterUse.prc, 
          labels = c("(a)", "(b)"), 
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = 0.1,
          align="v", nrow=2) %>% 
  save_plot(filename = file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_CannabisWaterUse+Streamflow.png"), 
            plot = .,
            ncol = 1,
            nrow = 2,
            base_width = 95/25.4,
            base_height = 60/25.4)
