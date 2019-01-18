## Figure_CannabisWaterUse.R
#' This script is intended to make a figure showing monthly cannabis water use for indoor and outdoor growers.

source("src/paths+packages.R")

## read in cannabis water use data - this is proprietary from TNC, cannot be shared (until Wilson et al paper published)
df <- read.csv(file.path(dir.TNC, "CannabisMonthlyWaterUse_WilsonEtAl.csv"), stringsAsFactors=F)

## set monthly factor
df$Month <- factor(df$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df$Setting <- factor(df$Setting, levels=c("Outdoor", "Greenhouse"))

## gallons to liters conversion factor
gal.to.L <- 3.78541

## plot
ggplot(df, aes(x=Month)) +
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
