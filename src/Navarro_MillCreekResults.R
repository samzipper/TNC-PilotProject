## Navarro_MillCreekResults.R
#' Looking at the results from Mill Creek subwatershed per Jen's request.

source(file.path("src", "paths+packages.R"))

## load mill creek subwatershed boundary
sf_streams <- 
  file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp") %>% 
  sf::st_read()

sf_mill <- 
  file.path(dir.TNC, "nav_watersheds_mill", "nav_watersheds_mill.shp") %>% 
  sf::st_read() %>% 
  sf::st_transform(crs = sf::st_crs(sf_streams))

# figure out which SegNum is Mill Creek
sf_streams$mill_creek <- as.numeric(sf::st_intersects(sf_streams, sf_mill))

ggplot() + 
  geom_sf(data = sf_streams, aes(geometry = geometry, color = factor(mill_creek))) + 
  geom_sf(data = sf_mill, aes(geometry = geometry), color = "red", fill = NA)

segs_mill <- sf_streams$SegNum[is.finite(sf_streams$mill_creek)]

#### load and process cannabis data
# depletion by segment associated with each well - output from Navarro_Cannabis-Grows_02_DepletionBySegment.R
df.grow <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_03_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F)

# load well-stream pairs which has distance
df.pairs <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  dplyr::select(SegNum, GrowNum, dist_wellToStream_m)

df.grow <-
  left_join(df.grow, df.pairs, by=c("SegNum", "GrowNum"))

# calculate year and month
df.grow$year.dec <- df.grow$time_days/365 + 1
df.grow$year <- floor(df.grow$year.dec)
df.grow$month <- round((df.grow$year.dec - df.grow$year)*12)+1

# determine what grows use groundwater, from Chris' model
# grow locations shapefile
sf.grows <- 
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg")) %>% 
  sf::st_transform(crs.MODFLOW)

# summarize to year/month
df.grow.depletion.summary <-
  # trim data frame
  df.grow %>% 
  subset(SegNum %in% segs_mill) %>% 
  subset(GrowNum %in% sf.grows$GrowNum[sf.grows$Well.rf.pred=="Yes"]) %>% 
  # summarize
  group_by(year, month) %>% 
  summarize(depletion_m3d_Mill = sum(depletion_m3d)) %>% 
  # join with streamflow
  transform(WaterUser="Cannabis")

#### load and process residential data
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

# depletion by segment associated with each well - output from Navarro_Residential_03_DepletionBySegment.R
df.res <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Residential_03_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) 

# calculate year and month
df.res$year.dec <- df.res$time_days/365 + 1
df.res$year <- floor(df.res$year.dec)
df.res$month <- round((df.res$year.dec - df.res$year)*12)+1

# trim and summarize data frame
df.res.depletion.summary <-
  df.res %>% 
  subset(HouseNum %in% sf.houses$HouseNum[sf.houses$groundwater]) %>%  # groundwater users only
  subset(SegNum %in% segs_mill) %>%   # navarro only
  # sum for all segments in Navarro
  group_by(year, month) %>% 
  summarize(depletion_m3d_Mill = sum(depletion_m3d)) %>% 
  transform(WaterUser="Residential")

## combine residential and cannabis
df.depletion.summary <- 
  rbind(df.grow.depletion.summary, df.res.depletion.summary)

## save output
df.depletion.summary %>% 
  readr::write_csv(file.path("results", "Navarro_MillCreekResults.csv"))
