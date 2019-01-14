## Navarro_WellCompletionReports.R
#' This script is intended to investigate the CA well completion reports dataset
#' for the Navarro River Watershed. It is at township & range resolution.

source(file.path("src", "paths+packages.R"))

## load data
# well completion reports
sf.wells <- 
  sf::st_read(file.path("data", "WellCompletionReports", "Navarro_WellCompletionReports.shp"), stringsAsFactors=F) %>% 
  subset(is.finite(Top.Of.Per))

# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("modflow", "input", "iriv.shp"), stringsAsFactors=F) %>% 
  subset(TermnlP == outlet.TerminalPa)

# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# rasters
r.lulc <- raster(paste0(dir.gis, "Navarro_Cannabis_LULC_30m.tif"))
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))
r.aquifers <- raster(paste0(dir.gis, "Navarro_Cannabis_GroundwaterBasins_30m.tif"))

# extract elevation
sf.wells$elev_m <- raster::extract(r.dem.30m, sf.wells)
sf.wells$dtb_m <- raster::extract(r.dtb, sf.wells)

# load in hand-coded bedrock info
df.bedrock <- 
  file.path("data", "WellCompletionReports", "Navarro_WellCompletionReports_BedrockInfo.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  left_join(sf.wells[,c("WCR.Number", "elev_m", "dtb_m")], by="WCR.Number")

## make some simple plots
ggplot(sf.wells) + geom_sf()

ggplot(sf.wells, aes(x=Top.Of.Per)) + geom_histogram(binwidth=10)

ggplot(sf.wells, aes(x=Well.Yield)) + geom_histogram(binwidth=1)

ggplot(sf.wells, aes(y=Top.Of.Per, x=elev_m)) +
  geom_point() +
  stat_smooth(method="lm")

ggplot(sf.wells, aes(y=Top.Of.Per, x=dtb_m)) +
  geom_point() +
  stat_smooth(method="lm")

## compare global DTB dataset to hand-coded DTB from well logs
# wells that are in bedrock - should plot along 1:1 line
#  above 1:1 line means global dataset underpredicts bedrock depth; below 1:1 line means global dataset overpredicts
df.bedrock %>% 
  subset(BedrockDepth_ft != 9999) %>% 
  ggplot(aes(x=dtb_m, y=BedrockDepth_ft*0.3048)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red") +
  stat_smooth(method="lm")
lm(BedrockDepth_ft*0.3048 ~ dtb_m, data=subset(df.bedrock, BedrockDepth_ft != 9999)) %>% 
  summary()

# same but including weathered bedrock
df.bedrock %>% 
  subset(WeatheredBedrockDepth_ft != 9999) %>% 
  ggplot(aes(x=dtb_m, y=WeatheredBedrockDepth_ft*0.3048)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red") +
  stat_smooth(method="lm")
lm(WeatheredBedrockDepth_ft*0.3048 ~ dtb_m, data=subset(df.bedrock, WeatheredBedrockDepth_ft != 9999)) %>% 
  summary()

# positive value means global dataset underpredicts; negative value means global dataset overpredicts
df.bedrock %>% 
  subset(WeatheredBedrockDepth_ft != 9999) %>% 
  ggplot(aes(color=WeatheredBedrockDepth_ft*0.3048-dtb_m)) +
  geom_sf() +
  scale_color_gradient2()

# wells that do not reach bedrock - should plot below 1:1 line
#  gets 17/22 correct (~77%)
df.bedrock %>% 
  subset(BedrockDepth_ft == 9999) %>% 
  ggplot(aes(y=ScreenBottomDepth_ft*0.3048, x=dtb_m)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

## look at statistical relationships with dtb and elev
lm(Top.Of.Per ~ dtb_m, data=sf.wells) %>% 
  summary()

lm(Top.Of.Per ~ elev_m, data=sf.wells) %>% 
  summary()

lm(Top.Of.Per ~ elev_m + dtb_m, data=sf.wells) %>% 
  summary()

# group by township, range, section
sf.wells.summary <-
  sf.wells %>%  
  group_by(Township, Range, Section) %>% 
  summarize(n.wells = n(),
            lat = mean(Decimal.La),
            lon = mean(Decimal.Lo),
            ScreenTopDepth_m_mean = mean(Top.Of.Per)*0.3048,
            ScreenTopDepth_m_median = median(Top.Of.Per)*0.3048,
            ScreenTopDepth_m_sd = sd(Top.Of.Per)*0.3048)

ggplot() +
  geom_sf(data=sf.wells.summary, aes(color=ScreenTopDepth_m_mean))

ggplot() +
  geom_sf(data=sf.wells.summary, aes(color=n.wells))
