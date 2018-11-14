## Figure_Habitat_IntrinsicPotential.R
#' Map of Navarro River Watershed with intrinsic potential habitat suitability.

source(file.path("src", "paths+packages.R"))

## load data
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("modflow", "input", "iriv.shp"), stringsAsFactors=F) %>% 
  subset(TermnlP == outlet.TerminalPa)

# domain boundary shapefile
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# habitat suitability
sf.coho <- 
  sf::st_read(file.path("data", "HabitatSuitability", "Navarro_NOAA_IP_Coho.shp")) %>% 
  dplyr::select(LLID, ORDER_, CO_IP215_C)
sf.coho$CO_IP215_C[sf.coho$CO_IP215_C < 0] <- 0

sf.chinook <- 
  sf::st_read(file.path("data", "HabitatSuitability", "Navarro_NOAA_IP_Chinook.shp")) %>% 
  dplyr::select(LLID, ORDER_, CHK_IP_CUR)

sf.steelhead <- 
  sf::st_read(file.path("data", "HabitatSuitability", "Navarro_NOAA_IP_Steelhead.shp")) %>% 
  dplyr::select(LLID, ORDER_, ST_IP_CURV)

## create length column
sf.streams$length_m <- as.numeric(sf::st_length(sf.streams))
sf.chinook$length_m <- as.numeric(sf::st_length(sf.chinook))
sf.coho$length_m <- as.numeric(sf::st_length(sf.coho))
sf.steelhead$length_m <- as.numeric(sf::st_length(sf.steelhead))

qplot(sf.steelhead$length_m)

mean(sf.streams$length_m)
mean(sf.chinook$length_m)
mean(sf.coho$length_m)
mean(sf.steelhead$length_m)

## the NOAA data has super short stream segments - intersect with sf.streams to aggregate
sf.streams.buffer <- sf::st_buffer(sf.streams, dist=2)  # lines don't exactly match, so create buffer and intersect

# summarize by SegNum
df.chinook.summarize <- 
  sf::st_intersection(sf.chinook, sf.streams.buffer) %>% 
  as.data.frame() %>% 
  dplyr::group_by(SegNum) %>% 
  dplyr::summarize(Chinook_IP_mean = weighted.mean(CHK_IP_CUR, length_m),
                   Chinook_IP_max = max(CHK_IP_CUR))

df.coho.summarize <- 
  sf::st_intersection(sf.coho, sf.streams.buffer) %>% 
  as.data.frame() %>% 
  dplyr::group_by(SegNum) %>% 
  dplyr::summarize(Coho_IP_mean = weighted.mean(CO_IP215_C, length_m),
                   Coho_IP_max = max(CO_IP215_C))

df.steelhead.summarize <- 
  sf::st_intersection(sf.steelhead, sf.streams.buffer) %>% 
  as.data.frame() %>% 
  dplyr::group_by(SegNum) %>% 
  dplyr::summarize(Steelhead_IP_mean = weighted.mean(ST_IP_CURV, length_m),
                   Steelhead_IP_max = max(ST_IP_CURV))

# make overall output sf object
sf.streams.habitat <- 
  sf.streams %>% 
  left_join(df.chinook.summarize, by=c("SegNum")) %>% 
  left_join(df.coho.summarize, by=c("SegNum")) %>% 
  left_join(df.steelhead.summarize, by=c("SegNum"))

# combine sf into one
sf.all <-
  sf.streams.habitat %>% 
  dplyr::select(SegNum, Coho_IP_mean, Coho_IP_max, Chinook_IP_mean, Chinook_IP_max, Steelhead_IP_mean, Steelhead_IP_max) %>% 
  melt(id=c("SegNum", "geometry"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1), 
                           labels=c("Low", "High"),
                           include.lowest=T)) %>% 
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3])

sf.n.high <-
  sf.all %>% 
  as.data.frame() %>% 
  group_by(SegNum, metric) %>% 
  summarize(n.species.high = sum(IP > 0.7)) %>% 
  subset(metric=="mean") %>% 
  left_join(sf.streams, ., by=c("SegNum"))

## make plot
sf.all %>% 
  subset(metric=="mean") %>% 
  ggplot() +
  facet_wrap(~species) +
  geom_sf(aes(color=IP)) +
  geom_sf(data=sf.basin, color=col.cat.red, fill=NA) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_color_viridis(name="Intrinsic\nHabitat\nPotential", limits=c(0,1), breaks=seq(0,1,0.25),
                      guide = guide_colourbar(barheight=4)) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position=c(0,0),
        legend.justification=c(0,0),
        legend.background=element_blank(),
        legend.box.background=element_blank()) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat_IntrinsicPotential_Continuous.png"),
         width=190, height=74, units="mm") +
  NULL

sf.all %>% 
  subset(metric=="mean") %>% 
  ggplot() +
  facet_wrap(~species) +
  geom_sf(aes(color=IP_class)) +
  geom_sf(data=sf.basin, color=col.cat.red, fill=NA) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_color_manual(name="Intrinsic\nHabitat\nPotential", values=c("Low"=col.cat.blu, "High"=col.cat.org)) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position=c(0,0),
        legend.justification=c(0,0),
        legend.background=element_blank(),
        legend.box.background=element_blank())  +
  guides(color=guide_legend(override.aes = list(fill=c(col.cat.blu, col.cat.org)),
                            keyheight=0.1, keywidth=2)) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat_IntrinsicPotential_Discrete.png"),
         width=190, height=74, units="mm") +
  NULL

sf.n.high %>% 
  ggplot() +
  geom_sf(aes(color=n.species.high>0)) +
  geom_sf(data=sf.basin, color=col.cat.red, fill=NA) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_color_manual(name="Ecological Value", 
                     values=c("FALSE"=col.cat.blu, "TRUE"=col.cat.org),
                     labels=c("FALSE"="Low", "TRUE"="High")) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position=c(0,0),
        legend.justification=c(0,0),
        legend.background=element_blank(),
        legend.box.background=element_blank())  +
  guides(color=guide_legend(override.aes = list(fill=c(col.cat.blu, col.cat.org)),
                            keyheight=0.1, keywidth=2)) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_Habitat_IntrinsicPotential_Discrete-AnySpecies.png"),
         width=95, height=95, units="mm") +
  NULL

sum(sf.n.high$n.species.high >= 1)
sum(is.finite(sf.n.high$n.species.high))

sum(subset(sf.n.high, n.species.high >= 1)$length_m)/1000
sum(sf.n.high$length_m)/1000

