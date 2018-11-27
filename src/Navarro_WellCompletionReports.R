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
r.lulc <- raster(paste0(dir.gis, "Navarro_Habitat_LULC_30m.tif"))
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Habitat_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Habitat_DTB_30m.tif"))
r.aquifers <- raster(paste0(dir.gis, "Navarro_Habitat_GroundwaterBasins_30m.tif"))

# extract elevation
sf.wells$elev_m <- raster::extract(r.dem.30m, sf.wells)
sf.wells$dtb_m <- raster::extract(r.dtb, sf.wells)

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
