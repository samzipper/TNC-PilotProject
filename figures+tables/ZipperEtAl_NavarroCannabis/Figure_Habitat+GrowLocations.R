## Figure_Habitat+GrowLocations.R
#' Map of Navarro River Watershed with intrinsic potential habitat suitability.
#' Requires output from Navarro_Cannabis_HabitatIntrinsicPotential.R

source(file.path("src", "paths+packages.R"))

## load data from Navarro_Cannabis_HabitatIntrinsicPotential.R
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  subset(TermnlP == outlet.TerminalPa) %>% 
  left_join(read.csv(file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv")),
            by=c("SegNum"))

# grow locations shapefile
sf.grows <- 
  sf::st_read(file.path(dir.TNC, "nav_cannabis_ucbtnc", "nav_cannabis_ucbtnc.shp")) %>% 
  subset(year==16) %>%                  # 2016 data only
  st_zm(drop = TRUE, what = "ZM") %>%   # drop Z dimension from geometry
  sf::st_transform(crs.MODFLOW)

## read in cannabis water use data - this is proprietary from TNC, cannot be shared (until Wilson et al paper published)
df.pump <- read.csv(file.path(dir.TNC, "CannabisMonthlyWaterUse_WilsonEtAl.csv"), stringsAsFactors=F)

## set monthly factor
df.pump$Month <- factor(df.pump$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump$Setting <- factor(df.pump$Setting, levels=c("Outdoor", "Greenhouse"))

# calculate total water use
df.pump.yr <- 
  df.pump %>% 
  transform(WaterUse_m3PlantDay = MeanWaterUse_GalPlantDay*gal.to.m3,
            MonthDays = days_in_month(as.numeric(Month))) %>% 
  group_by(Setting) %>% 
  summarize(WaterUse_m3PlantYr = sum(WaterUse_m3PlantDay*MonthDays))

# calculate water use for each grow
sf.grows$WaterUse_m3y <- NaN
sf.grows$WaterUse_m3y[sf.grows$greenhouse==0] <- df.pump.yr$WaterUse_m3PlantYr[df.pump.yr$Setting=="Outdoor"]*sf.grows$plants[sf.grows$greenhouse==0]
sf.grows$WaterUse_m3y[sf.grows$greenhouse==1] <- df.pump.yr$WaterUse_m3PlantYr[df.pump.yr$Setting=="Greenhouse"]*sf.grows$plants[sf.grows$greenhouse==1]

# domain boundary shapefile
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>% 
  sf::st_transform(crs.MODFLOW)

sf.subbasin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU12_Navarro.shp")) %>% 
  sf::st_transform(crs.MODFLOW)

# aggregate number of plants by subbasin
sf.subbasins.plants <-
  sf::st_join(sf.grows, sf.subbasin, join=st_within) %>% 
  subset(HUC12 != 180101080407) %>%   # grows within this HUC12 were not mapped, but 2 dots are inside it anyways; due to projection rounding?
  group_by(HUC12) %>% 
  summarize(TotalPlants = sum(plants),
            TotalGrowsize = sum(growsize),
            NumberGrows = sum(is.finite(plants)),
            TotalWaterUse_m3y = sum(WaterUse_m3y)) %>% 
  as.data.frame() %>% 
  left_join(sf.subbasin, ., by="HUC12")

# combine habitat sf into one
sf.all <-
  sf.streams %>% 
  dplyr::select(SegNum, Coho_IP_mean, Coho_IP_max, Chinook_IP_mean, Chinook_IP_max, Steel_IP_mean, Steel_IP_max) %>% 
  melt(id=c("SegNum", "geometry"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1), 
                           labels=c("Low", "High"),
                           include.lowest=T)) %>% 
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3]) %>% 
  subset(species == "Coho" & metric == "max")


# plot
ggplot() +
  geom_sf(data=sf.subbasins.plants, aes(fill=TotalWaterUse_m3y/10000), color="white") +
  geom_sf(data=sf.streams, color="white") +
  geom_sf(data=sf.all, aes(color=IP_class)) +
  scale_fill_viridis(name="Estimated Cannabis Water\nUse [10,000 cu. m/yr]", 
                     trans="log10", limits=c(1, max(sf.subbasins.plants$TotalWaterUse_m3y/10000))) +
  scale_color_manual(name="Intrinsic Habitat Potential", values=c("Low"="black", "High"=col.cat.red)) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.box="vertical",
        legend.margin = margin(0,0,0,0, unit="mm")) +
  guides(color=guide_legend(override.aes = list(fill=c("black", col.cat.red)), keyheight=0.05, keywidth=3,
                            title.position="top", title.hjust=0.5, title.vjust=0.5),
         fill=guide_colorbar(title.hjust=0.5, title.vjust=1)) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat+GrowLocations.png"),
         width=95, height=120, units="mm") +
  NULL
  