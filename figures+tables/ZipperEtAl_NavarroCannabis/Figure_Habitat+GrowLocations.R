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
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg")) %>% 
  subset(Well.rf.pred=="Yes") %>%
  sf::st_transform(crs.MODFLOW)

## read in cannabis water use data - this is proprietary from TNC/Dillis, cannot be shared
# created using script Navarro_Cannabis-Grows_03_DepletionBySegment.R
df.pump <- read.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_03_PumpingRate.csv"), stringsAsFactors=F) %>% 
  subset(GrowNum %in% sf.grows$GrowNum)
df.pump$WaterUseMean_m3d[df.pump$WaterUseMean_m3d < 0] <- 0
df.pump$WaterUseSum_m3mo <- df.pump$WaterUseMean_m3d*lubridate::days_in_month(df.pump$MonthNum)

## set monthly factor
df.pump$Month <- factor(df.pump$MonthNum, labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

# calculate total water use
df.pump.yr <- 
  df.pump %>% 
  dplyr::group_by(GrowNum) %>% 
  dplyr::summarize(WaterUseSum_m3yr = sum(WaterUseSum_m3mo))

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

# combine habitat sf into one
sf.all <-
  sf.streams %>% 
  dplyr::select(SegNum, Coho_IP_mean, Coho_IP_max, Chinook_IP_mean, Chinook_IP_max, Steel_IP_mean, Steel_IP_max, geometry) %>% 
  reshape2::melt(id=c("SegNum", "geometry"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1), 
                           labels=c("Low", "High"),
                           include.lowest=T)) %>% 
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3]) %>% 
  subset(species == "Coho" & metric == "max") %>% 
  sf::st_sf()

## plot as a raster
sf.grows.wateruse <-
  dplyr::left_join(sf.grows, df.pump.yr, by="GrowNum")

r.wateruse <-
  sf.grows.wateruse %>% 
  rasterize(., raster(crs=crs.MODFLOW, ext=extent(sf.subbasin), resolution=2000), 
            field = "WaterUseSum_m3yr", fun = "sum")
values(r.wateruse)[is.na(values(r.wateruse))] <- 0

df.wateruse <- 
  r.wateruse %>% 
  mask(sf.subbasin) %>% 
  as.data.frame(xy=T) %>% 
  set_colnames(c("x", "y", "can_m3yr"))

## map of residential water use
# define pumping rates: 250 gpm in winter, 500 gpm in summer
df.pump.res <- data.frame(
  Month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  MeanWaterUse_GalHouseDay = c(250, 250, 250, 250, 500, 500, 500, 500, 500, 500, 250, 250)
)
df.pump.res$Month <- factor(df.pump.res$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump.res$MonthNum <- match(df.pump.res$Month, month.abb)
df.pump.res$MonthLengthDays <- lubridate::days_in_month(df.pump.res$MonthNum)
df.pump.res$m3HouseDay <- df.pump.res$MeanWaterUse_GalHouseDay*gal.to.m3
house.pump <- sum(df.pump.res$m3HouseDay*df.pump.res$MonthLengthDays)

# load point locations
sf.houses <- 
  sf::st_read(file.path(dir.TNC, "Structures_Navarro_NAIP_2016", "Structures_Navarro_NAIP_2016.shp")) %>% 
  subset(Structure=="Res H") %>%                  # houses only
  st_zm(drop = TRUE, what = "ZM") %>%   # drop Z dimension from geometry
  sf::st_transform(crs.MODFLOW)
sf.houses$WaterUseSum_m3yr <- house.pump

r.res.wateruse <-
  sf.houses %>% 
  rasterize(., raster(crs=crs.MODFLOW, ext=extent(sf.subbasin), resolution=2000), 
            field = "WaterUseSum_m3yr", fun = "sum")
values(r.res.wateruse)[is.na(values(r.res.wateruse))] <- 0

df.res.wateruse <- 
  r.res.wateruse %>% 
  mask(sf.subbasin) %>% 
  as.data.frame(xy=T) %>% 
  set_colnames(c("x", "y", "res_m3yr"))

## raster, cannabis and residential
dplyr::left_join(df.wateruse, df.res.wateruse, by=c("x", "y")) %>% 
  melt(id=c("x", "y")) %>% 
  ggplot() +
  geom_raster(aes(x=x, y=y, fill=value/1000)) +
  geom_sf(data=sf.basin, color=col.gray, fill=NA, aes(geometry = geometry)) +
  geom_sf(data=sf.streams, color="white", aes(geometry = geometry)) +
  geom_sf(data=sf.all, aes(color=IP_class, geometry = geometry)) +
  facet_wrap(~variable, labeller=as_labeller(c("can_m3yr"="(a) Cannabis", "res_m3yr"="(b) Residential"))) +
  scale_fill_viridis(name="Annual Groundwater Use [x1000 m\u00b3]", na.value=NA, trans="log10",
                     limits=c(0.1, max(df.res.wateruse$res_m3yr/1000, na.rm=T)),
                     breaks=c(0.1,1,10,70), labels=c("0.1", "1", "10", "70")) +
  scale_color_manual(name="Stream Intrinsic Habitat Potential", values=c("Low"="black", "High"=col.cat.red)) +
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
  guides(color=guide_legend(override.aes = list(fill=c("black", col.cat.red)), 
                            keyheight=0.01, keywidth=4, title.hjust=0.5, title.vjust=0.6),
         fill=guide_colorbar(title.hjust=0.5, title.vjust=0.8, order=1)) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat+GrowLocations-raster.png"),
         width=190, height=128, units="mm") +
  NULL

# ## raster, residential only
# ggplot() +
#   geom_raster(data=df.res.wateruse, aes(x=x, y=y, fill=res_m3yr)) +
#   geom_sf(data=sf.basin, color=col.gray, fill=NA) +
#   geom_sf(data=sf.streams, color="white") +
#   geom_sf(data=sf.all, aes(color=IP_class)) +
#   scale_fill_viridis(name="Annual Residential \nGroundwater Use [m\u00b3]", na.value=NA, 
#                      trans="log10") +
#   scale_color_manual(name="Stream Intrinsic Habitat Potential", values=c("Low"="black", "High"=col.cat.red)) +
#   scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
#   scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
#   coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
#   theme(panel.grid=element_line(color="transparent")) +
#   theme(axis.text.y=element_text(angle=90, hjust=0.5),
#         legend.position="bottom",
#         legend.background=element_blank(),
#         legend.box.background=element_blank(),
#         legend.box="vertical",
#         legend.margin = margin(0,0,0,0, unit="mm")) +
#   guides(color=guide_legend(override.aes = list(fill=c("black", col.cat.red)), keyheight=0.05, keywidth=3,
#                             title.position="top", title.hjust=0.5, title.vjust=0.5, show.legend=F),
#          fill=guide_colorbar(title.hjust=0.5, title.vjust=1)) +
#   # ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat+GrowLocations-raster.png"),
#   #         width=95, height=120, units="mm") +
#   NULL
# 
# ## raster, cannabis only
# p.cannabis.raster <- 
#   ggplot() +
#   geom_raster(data=df.wateruse, aes(x=x, y=y, fill=can_m3yr)) +
#   geom_sf(data=sf.basin, color=col.gray, fill=NA) +
#   geom_sf(data=sf.streams, color="white") +
#   geom_sf(data=sf.all, aes(color=IP_class)) +
#   scale_fill_viridis(name="Annual Cannabis \nGroundwater Use [m\u00b3]", na.value=NA, 
#                      trans="log10", limits=c(100, max(df.wateruse$layer, na.rm=T)),
#                      breaks = c(100, 1000, 5000)) +
#   scale_color_manual(name="Stream Intrinsic Habitat Potential", values=c("Low"="black", "High"=col.cat.red)) +
#   scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
#   scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
#   coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
#   theme(panel.grid=element_line(color="transparent")) +
#   theme(axis.text.y=element_text(angle=90, hjust=0.5),
#         legend.position="bottom",
#         legend.background=element_blank(),
#         legend.box.background=element_blank(),
#         legend.box="vertical",
#         legend.margin = margin(0,0,0,0, unit="mm")) +
#   guides(color=guide_legend(override.aes = list(fill=c("black", col.cat.red)), keyheight=0.05, keywidth=3,
#                             title.position="top", title.hjust=0.5, title.vjust=0.5),
#          fill=guide_colorbar(title.hjust=0.5, title.vjust=1)) +
# #  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat+GrowLocations-raster.png"),
# #         width=95, height=120, units="mm") +
#   NULL
# 
# ## plot cannabis use by subbasin
# # aggregate number of plants by subbasin
# sf.subbasins.plants <-
#   dplyr::left_join(sf.grows, df.pump.yr, by="GrowNum") %>% 
#   sf::st_join(sf.subbasin, join=st_within) %>% 
#   subset(HUC12 != 180101080407) %>%   # grows within this HUC12 were not mapped, but 2 dots are inside it anyways; due to projection rounding?
#   group_by(HUC12) %>% 
#   summarize(TotalPlants = sum(plants),
#             TotalGrowsize = sum(growsize),
#             NumberGrows = sum(is.finite(plants)),
#             TotalWaterUse_m3y = sum(WaterUseSum_m3yr)) %>% 
#   as.data.frame() %>% 
#   left_join(sf.subbasin, ., by="HUC12")
# 
# ggplot() +
#   geom_sf(data=sf.subbasins.plants, aes(fill=TotalWaterUse_m3y/1000), color="white") +
#   geom_sf(data=sf.streams, color="white") +
#   geom_sf(data=sf.all, aes(color=IP_class)) +
#   scale_fill_viridis(name="Annual Cannabis\nGroundwater Use [1000 m\u00b3]", trans="log10") +
#   scale_color_manual(name="Stream Intrinsic Habitat Potential", values=c("Low"="black", "High"=col.cat.red)) +
#   scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
#   scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
#   coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
#   theme(panel.grid=element_line(color="transparent")) +
#   theme(axis.text.y=element_text(angle=90, hjust=0.5),
#         legend.position="bottom",
#         legend.background=element_blank(),
#         legend.box.background=element_blank(),
#         legend.box="vertical",
#         legend.margin = margin(0,0,0,0, unit="mm")) +
#   guides(color=guide_legend(override.aes = list(fill=c("black", col.cat.red)), keyheight=0.05, keywidth=3,
#                             title.position="top", title.hjust=0.5, title.vjust=0.5),
#          fill=guide_colorbar(title.hjust=0.5, title.vjust=1)) +
#   ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat+GrowLocations-subbasins.png"),
#          width=95, height=120, units="mm") +
#   NULL
