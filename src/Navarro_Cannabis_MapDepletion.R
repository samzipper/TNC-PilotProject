## Navarro_Cannabis_DepletionBySegment.R
#' This script is intended to calculate depletion for stream segments 
#' with high intrinsic habitat potential.
#' 
#' It requires output from Navarro_Cannabis_CalculateWellStreamPairs.R 
#' and Navarro_Cannabis_HabitatIntrinsicPotential.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

### load and pre-process data
## read in data from Navarro_Cannabis_DepletionBySegment.R
df.depletion <- read.csv(file.path("results", "Navarro_Cannabis_DepletionBySegment.csv"), stringsAsFactors=F)

## output from Navarro_Cannabis_HabitatIntrinsicPotential.R
df.habitat <- 
  file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  melt(id=c("SegNum"), value.name="IP", variable.name="Species_IP_metric") %>% 
  replace_na(list(IP=0)) %>%   # stream segments with no IP data indicates not suitable habitat
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1), 
                           labels=c("Low", "High"),
                           include.lowest=T)) %>% 
  transform(species = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,1],
            metric = str_split_fixed(Species_IP_metric, pattern="_", n=3)[,3])

# summarize each segment to max value for any species
df.habitat.summary <-
  df.habitat %>% 
  subset(metric=="mean") %>% 
  group_by(SegNum) %>% 
  summarize(IP = max(IP)) %>% 
  transform(IP_class = cut(IP, 
                           breaks=c(0,0.7,1, 5),   # nothing will get Outside Navarro category
                           labels=c("Low", "High", "Outside Navarro"),
                           include.lowest=T))

# join habitat data with stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  left_join(df.habitat.summary, by=c("SegNum")) %>% 
  replace_na(list("IP_class" = "Outside Navarro"))

## well locations
# synthetic pumping wells
sf.wel <- 
  file.path("modflow", "input", "iwel.txt") %>% 
  read.table(header=T) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(sf.streams))

## domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>% 
  sf::st_transform(crs=st_crs(sf.streams))

### combine data
# add IP_class to depletion
df.depletion <- 
  left_join(df.depletion, df.habitat.summary, by="SegNum") %>% 
  replace_na(list("IP_class"="Outside Navarro"))

# summarize depletion in high ecological value segments by time_d and WellNum
df.value <-
  df.depletion %>% 
  group_by(time_d, WellNum, IP_class) %>% 
  summarize(Qf.sum = sum(Qf))

sf.value <- 
  full_join(sf.wel, df.value, by="WellNum")

### plots
# times to plot
times_plot <- unique(df.value$time_d)

# for facets, need to duplicate basin and stream for all timesteps
# also create rasters of interpolated depletion
r.empty <- raster(extent(sf.basin), crs=crs(sf.basin))
res(r.empty) <- c(200,200)

for (t in times_plot){
  ## inverse distance interpolation to make rasters
  df.interp <- data.frame(Z = subset(sf.value, IP_class=="High" & time_d == t)$Qf.sum,
                          X = st_coordinates(subset(sf.value, IP_class=="High" & time_d == t))[,"X"],
                          Y = st_coordinates(subset(sf.value, IP_class=="High" & time_d == t))[,"Y"])
  sp.interp <- SpatialPoints(df.interp[,c("X", "Y")], proj4string=CRS(st_crs(sf.basin)[["proj4string"]]))
  sp.interp <- SpatialPointsDataFrame(sp.interp, df.interp)
  gs <- gstat(formula = Z~1, locations = sp.interp, nmax = 6)  # nmax controls weighting
  r.t <- interpolate(r.empty, gs)
  r.t.mask <- mask(r.t, sf.basin)
  df.Qf.t <- 
    r.t.mask %>% 
    rasterToPoints() %>% 
    as.data.frame() %>% 
    set_colnames(c("lon", "lat", "Qf.sum")) %>% 
    transform(time_d = t)
  
  if (t == times_plot[1]){
    sf.streams.facet <- 
      sf.streams %>% 
      transform(time_d = t)
    
    sf.basin.facet <- 
      sf.basin %>% 
      transform(time_d = t)
    
    df.Qf <- df.Qf.t
    
  } else {
    sf.streams.facet <- 
      sf.streams %>% 
      transform(time_d = t) %>% 
      rbind(sf.streams.facet, .)
    
    sf.basin.facet <- 
      sf.basin %>% 
      transform(time_d = t) %>% 
      rbind(sf.basin.facet, .)
    
    df.Qf <- rbind(df.Qf, df.Qf.t)
    
  }
}

# build labeller
time_labels <- 
  setNames(paste0("Year ", unique(df.Qf$time_d)/365),
           unique(df.Qf$time_d))

ggplot(data=df.Qf) +
  geom_raster(aes(x = lon, y=lat, fill=Qf.sum)) +
  geom_sf(data=subset(sf.streams.facet, IP_class != "Outside Navarro"), aes(color=IP_class)) +
#  geom_sf(data=sf.basin.facet, color=col.cat.red, fill=NA) +
  facet_wrap(~(time_d),
             nrow = 1,
             labeller = as_labeller(time_labels)) +
  scale_fill_viridis(limits=c(0,1), breaks=seq(0,1,0.5),
                     labels = scales::percent) +
  scale_color_manual(values=c("Low"=col.cat.blu, "High"=col.cat.org, "Outside Navarro"=col.gray)) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank())  +
  guides(color=guide_legend(title="Stream Ecological Value", 
                            override.aes = list(fill=c(col.cat.blu, col.cat.org)),
                            keyheight=0.1, 
                            keywidth=2, 
                            title.position="top", 
                            title.hjust=0.5,
                            order=2),
         fill = guide_colorbar(title = "Depletion from High-Value\nStreams [% of pumping rate]",
                               #title.position="top", 
                               title.hjust=0.5,
                               order=1)) +
  ggsave(file.path("results", "Navarro_Cannabis_MapDepletion_MapQf.png"),
         width=190, height=80, units="mm")
 