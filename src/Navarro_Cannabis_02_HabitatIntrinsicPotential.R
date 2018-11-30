## Navarro_Cannabis_HabitatIntrinsicPotential.R
#' Process intrinsic potential habitat suitability.

source(file.path("src", "paths+packages.R"))

## load data
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
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
  dplyr::summarize(Steel_IP_mean = weighted.mean(ST_IP_CURV, length_m),
                   Steel_IP_max = max(ST_IP_CURV))

# make overall output sf object
sf.streams.habitat <- 
  sf.streams %>% 
  left_join(df.chinook.summarize, by=c("SegNum")) %>% 
  left_join(df.coho.summarize, by=c("SegNum")) %>% 
  left_join(df.steelhead.summarize, by=c("SegNum"))

head(sf.streams.habitat)

## write to results
sf.streams.habitat %>% 
  as.data.frame %>% 
  dplyr::select(SegNum, Chinook_IP_mean, Chinook_IP_max, Coho_IP_mean, Coho_IP_max, Steel_IP_mean, Steel_IP_max) %>% 
  write.csv(file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv"),
            row.names=F)
