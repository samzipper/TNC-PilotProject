## Figure_SummaryResults.R
#' This script is intended to summarize depletion caused by wells located in 
#' different buffer distances from streams.

source(file.path("src", "paths+packages.R"))

## prep stream data
## output from Navarro_Cannabis_02_HabitatIntrinsicPotential.R
df.habitat <- 
  file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  melt(id=c("SegNum"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1.01,100),     # nothing will get Outside Navarro category
                           labels=c("Low", "High", "Outside Navarro"),
                           include.lowest=T)) %>% 
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3]) %>% 
  subset(species == "Coho" & metric=="max")  # focus on Coho based on conversation with Jen - most sensitive to habitat conditions

# join habitat data with stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  left_join(df.habitat, by=c("SegNum")) %>% 
  replace_na(list("IP_class" = "Outside Navarro")) %>% 
  subset(IP_class != "Outside Navarro")

# make 1.2 km buffer around stream segments (From Figure_ScaleSegment_Depletino-SensitivityAnalysis.R)
buff_dist <- 1200  # [m]
sf.streams.high.buffer <-
  sf.streams %>%
  subset(IP_class=="High") %>%
  sf::st_buffer(dist=buff_dist, endCapStyle="ROUND", joinStyle="MITRE") %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs.MODFLOW) %>% 
  smoothr::fill_holes(threshold=1e8)  # fill in little holes due to stream geometry and segmentation

sf.streams.all.buffer <-
  sf.streams %>%
  sf::st_buffer(dist=buff_dist, endCapStyle="ROUND", joinStyle="MITRE") %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs.MODFLOW) %>% 
  smoothr::fill_holes(threshold=1e8)

## prep cannabis data
# grow locations
sf.grows <- 
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg")) %>% 
  subset(Well.rf.pred=="Yes") %>%  # subset to GW only
  sf::st_transform(crs.MODFLOW)

# pumping rates
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

## intersect with stream buffer
sf.grows.high <- sf::st_intersection(sf.streams.high.buffer, sf.grows)  # within 1 km of high-value stream segment
sf.grows.low <- 
  sf::st_intersection(sf.streams.all.buffer, sf.grows) %>%     # within 1 km of any stream segment
  subset(!(GrowNum %in% sf.grows.high$GrowNum))

GrowNum.high <- unique(sf.grows.high$GrowNum)
GrowNum.low <- unique(sf.grows.low$GrowNum)
GrowNum.far <- unique(sf.grows$GrowNum)
GrowNum.far <- GrowNum.far[!(GrowNum.far %in% GrowNum.high) & !(GrowNum.far %in% GrowNum.low)]

# summarize pumping by proximity group
df.pump.summary <-
  data.frame(grow_class = c("< 1.2 km to\nHigh Potential", "< 1.2 km to\nLow Potential", "> 1.2 km\nto Stream"),
             WaterUse_m3yr = c(sum(df.pump.yr$WaterUseSum_m3yr[df.pump.yr$GrowNum %in% GrowNum.high]),
                               sum(df.pump.yr$WaterUseSum_m3yr[df.pump.yr$GrowNum %in% GrowNum.low]),
                               sum(df.pump.yr$WaterUseSum_m3yr[df.pump.yr$GrowNum %in% GrowNum.far])),
             n.grows = c(length(GrowNum.high), length(GrowNum.low), length(GrowNum.far)))

## prep residential data
# load house locations
sf.houses <- 
  file.path(dir.TNC, "Structures_Navarro_NAIP_2016", "Structures_Navarro_NAIP_2016.shp") %>% 
  sf::st_read() %>% 
  subset(Structure=="Res H")
sf.houses$HouseNum <- seq(1, dim(sf.houses)[1])

# load point of diversion to screen out surface water users
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

## intersect with stream buffer
sf.houses.GW <- 
  sf::st_transform(sf.houses, sf::st_crs(sf.streams.high.buffer)) %>% 
  subset(groundwater)
sf.houses.high <- sf::st_intersection(sf.streams.high.buffer, sf.houses.GW)  # within 1 km of high-value stream segment
sf.houses.low <- 
  sf::st_intersection(sf.streams.all.buffer, sf.houses.GW) %>%     # within 1 km of any stream segment
  subset(!(HouseNum %in% sf.houses.high$HouseNum))

n.res.gw <- dim(sf.houses.GW)[1]
n.res.high <- dim(sf.houses.high)[1]
n.res.low <- dim(sf.houses.low)[1]

# define pumping rates per house
df.pump.res <- 
  file.path("results", "Residential_WaterUse.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(domain == "Healdsburg")
df.pump.res$m3HouseDay <- df.pump.res$WaterUseMean_GalHouseDay*gal.to.m3
df.pump.res$m3HouseDay_std <- df.pump.res$WaterUseStd_GalHouseDay*gal.to.m3

m3.house.yr <- sum(df.pump.res$m3HouseDay*lubridate::days_in_month(df.pump.res$MonthNum))

# summarize residential pumping by proximity group
df.pump.res.summary <-
  data.frame(grow_class = c("< 1.2 km to\nHigh Potential", "< 1.2 km to\nLow Potential", "> 1.2 km\nto Stream"),
             WaterUse_m3yr = c(n.res.high*m3.house.yr,
                               n.res.low*m3.house.yr,
                               (n.res.gw - n.res.high - n.res.low)*m3.house.yr),
             n.grows = c(n.res.high, n.res.low, (n.res.gw - n.res.high - n.res.low)))

## prep map background
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>%
  sf::st_transform(crs.MODFLOW)

r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
df.dem <- 
  r.dem.30m %>% 
  mask(sf.basin) %>% 
  rasterToPoints() %>% 
  as.data.frame() %>% 
  set_colnames(c("lon", "lat", "dem"))

## make plots
p.bar <-
  ggplot(df.pump.summary, aes(x=grow_class)) +
  geom_bar(aes(y=WaterUse_m3yr/1e4, fill=grow_class), stat="identity") +
  geom_text(aes(y=WaterUse_m3yr/1e4+median(df.pump.summary$WaterUse_m3yr/1e4)*.075, label=paste0(n.grows, " parcels"))) +
  scale_fill_manual(guide=F, values=c(col.cat.red, "black", col.cat.grn)) +
  scale_y_continuous(name="Cannabis Groundwater Use\n[x 10,000 m\u00b3/yr]",
                     limits=c(0, max(df.pump.summary$WaterUse_m3yr/1e4)*1.1), expand=c(0,0)) +
  scale_x_discrete(name="Proximity to Stream")

p.bar.res <-
  ggplot(df.pump.res.summary, aes(x=grow_class)) +
  geom_bar(aes(y=WaterUse_m3yr/1e4, fill=grow_class), stat="identity") +
  geom_text(aes(y=WaterUse_m3yr/1e4+median(df.pump.res.summary$WaterUse_m3yr/1e4)*.14, label=paste0(n.grows, " houses"))) +
  scale_fill_manual(guide=F, values=c(col.cat.red, "black", col.cat.grn)) +
  scale_y_continuous(name="Residential Groundwater Use\n[x 10,000 m\u00b3/yr]",
                     limits=c(0, max(df.pump.res.summary$WaterUse_m3yr/1e4)*1.1), expand=c(0,0)) +
  scale_x_discrete(name="Proximity to Stream")

p.map <-
  ggplot() +
  geom_raster(data=df.dem, aes(x=lon, y=lat, fill=dem)) +
  geom_sf(data=sf.streams.all.buffer, fill=col.gray, color=NA, alpha=0.75) +
  geom_sf(data=sf.streams.high.buffer, fill=col.cat.red, alpha=0.25, color=NA) +
  geom_sf(data=sf.streams, aes(color=IP_class)) +
  scale_color_manual(name="Intrinsic Habitat Potential", values=c("Low"="black", "High"=col.cat.red)) +
  scale_fill_viridis(name="Land Surface\nElevation [m]", limits=c(0, max(df.dem$dem, na.rm=T)), breaks=c(0,500,1000)) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position=c(0,0),
        legend.justification=c(0,0),
        legend.background=element_blank(),
        legend.box.background=element_blank()) +
    guides(color=guide_legend(override.aes = list(fill=c("black", col.cat.red)), 
                              keyheight=0.05, keywidth=3, title.hjust=0, title.vjust=0.6, order=1),
           fill=guide_colorbar(title.hjust=0.5, title.vjust=0.8, order=2, barwidth=2, barheight=5))

## save plots
# top row: just map
# bottom row: both bar plots
plot_grid(p.bar, p.bar.res,
          nrow = 1,
          align = "b",
          labels = c("(b)", "(c)"),
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = c(0.05, 0.07),
          label_y = c(0.99, 0.92)) %>% 
  plot_grid(p.map, 
            .,
            nrow=2,
            rel_heights = c(1,0.5),
            labels = c("(a)", " "),
            label_size = 10,
            label_fontfamily = "Arial",
            label_fontface = "plain",
            label_x = c(0.2, 0.07),
            label_y = c(0.99, 0.99)) %>% 
  save_plot(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_SummaryResults.png"),
            plot = .,
            nrow=1,
            base_width=190/25.4,
            base_height=190/25.4)

plot_grid(p.bar, p.bar.res,
          nrow = 1,
          align = "b",
          labels = c("(b)", "(c)"),
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = c(0.05, 0.07),
          label_y = c(0.99, 0.92)) %>% 
  plot_grid(p.map, 
            .,
            nrow=2,
            rel_heights = c(1,0.5),
            labels = c("(a)", " "),
            label_size = 10,
            label_fontfamily = "Arial",
            label_fontface = "plain",
            label_x = c(0.2, 0.07),
            label_y = c(0.99, 0.99)) %>% 
  save_plot(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_SummaryResults.pdf"),
            plot = .,
            nrow=1,
            base_width=190/25.4,
            base_height=190/25.4,
            device=cairo_pdf)

# ## save plots (no residential)
# plot_grid(p.map, 
#           p.bar, 
#           nrow=1,
#           rel_widths=c(1,0.95),
#           align="t",
#           labels = c("(a)", "(b)"),
#           label_size = 10,
#           label_fontfamily = "Arial",
#           label_fontface = "plain",
#           label_x = c(0.1, 0.08),
#           label_y = c(0.99, 0.99)) %>% 
#   save_plot(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_SummaryResults.png"),
#             plot = .,
#             nrow=1,
#             base_width=190/25.4,
#             base_height=95/25.4)
