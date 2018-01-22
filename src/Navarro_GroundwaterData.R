## Navarro_GroundwaterData.R
#' Compile groundwater level data for basin.
#' 
#' Not included: USGS NWIS data. There are 3 wells which have 1-2 readings each in the 1980s.

source("src/paths+packages.R")

# # CA-DWR groundwater data -------------------------------------------------
# 
# ## note: this only has the recent data! Skipping it in favor of using 
# ##       recent + historic data which was manually downloaded
# 
# ## open basin boundary and river network shapefiles
# shp <- readOGR(dsn="data/NHD/WBD", layer="WBDHU10_Navarro")
# shp.riv <- readOGR(dsn="data/NHD/HYD", layer="NHDFlowline_HU10_Navarro")
# 
# shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU10_Navarro+Adjacent")
# shp.adj.riv <- readOGR(dsn="data/NHD/HYD", layer="NHDFlowline_HU10_Navarro+Adjacent")
# 
# # reproject to WGS to match global datasets
# shp.WGS <- spTransform(shp, crs.WGS)
# shp.riv.WGS <- spTransform(shp.riv, crs.WGS)
# shp.adj.WGS <- spTransform(shp.adj, crs.WGS)
# shp.adj.riv.WGS <- spTransform(shp.adj.riv, crs.WGS)
# 
# ## read in list of stations
# df.dwr.stations <- read.csv(paste0(dir.gw.dwr, "gst_file.csv"))
# 
# ## convert to spatialpoints
# # get longitude/latitude (in that order)
# xy <- df.dwr.stations[,c("LONGITUDE","LATITUDE")]
# 
# # convert to spatial points
# spdf.dwr.stations <- SpatialPointsDataFrame(coords = xy, data = df.dwr.stations,
#                                             proj4string = CRS(crs.WGS))
# 
# ## find wells inside basin boundary
# inside <- !is.na(over(spdf.dwr.stations, as(shp.WGS, "SpatialPolygons")))
# inside.adj <- !is.na(over(spdf.dwr.stations, as(shp.adj.WGS, "SpatialPolygons")))
# df.dwr.stations <- df.dwr.stations[inside,]
# df.dwr.stations.adj <- df.dwr.stations[inside.adj,]
# spdf.dwr.stations.adj <- spdf.dwr.stations[inside.adj,]
# 
# ## load groundwater level data
# df.dwr.gwl <- read.csv(paste0(dir.gw.dwr, "gwl_file.csv"))
# df.dwr.gwl.adj <- subset(df.dwr.gwl, SITE_CODE %in% df.dwr.stations.adj$SITE_CODE)
# 
# ## add in site info
# df.dwr.gwl.adj <- left_join(df.dwr.gwl.adj, df.dwr.stations.adj, by="SITE_CODE")
# 
# ## make output data frame
# df.dwr.out <- data.frame(agency = "CA-DWR",
#                          site = df.dwr.gwl.adj$SITE_CODE,
#                          datetime = mdy_hms(df.dwr.gwl.adj$MEASUREMENT_DATE),
#                          ground_elev_ft = df.dwr.gwl.adj$GS_ELEVATION,
#                          water_elev_ft = df.dwr.gwl.adj$RP_ELEVATION-df.dwr.gwl.adj$RP_READING,
#                          wtd_ft = df.dwr.gwl.adj$GS_ELEVATION - (df.dwr.gwl.adj$RP_ELEVATION-df.dwr.gwl.adj$RP_READING),
#                          lon = df.dwr.gwl.adj$LONGITUDE,
#                          lat = df.dwr.gwl.adj$LATITUDE,
#                          well_depth_ft = df.dwr.gwl.adj$TOTAL_DEPTH_FT,
#                          well_use = df.dwr.gwl.adj$CASGEM_STATION_USE_DESC,
#                          catchment = "Adjacent",
#                          stringsAsFactors=F)
# df.dwr.out$catchment[df.dwr.out$site %in% df.dwr.stations$SITE_CODE] <- "Navarro"
# 
# ## save data
# write.csv(df.dwr.out, "results/GroundwaterLevels_CA-DWR_Navarro.csv", row.names=F)

# CA-DWR data including historic ------------------------------------------

## get list of all files
files <- list.files("data/GroundwaterLevels/Navarro/", pattern="*.csv")

## loop and open files
start.flag <- T
for (f in files){
  # open file
  df.f.in <- read.csv(paste0("data/GroundwaterLevels/Navarro/", f), stringsAsFactors=F)
  
  # extract well name
  well <- str_split_fixed(f, "_", 3)[2]
  era <- str_split(str_split_fixed(f, "_", 3)[3], "[.]")[[1]][1]
  
  if (era=="Historic"){
    # get rid of last 7, which is metadata about well (lat/lon, etc)
    df.f.in <- df.f.in[1:(dim(df.f.in)[1]-7), ]
    
    # extract common data for historic and recent
    df.f <- data.frame(agency = "CA-DWR",
                       site = well,
                       datetime = mdy_hm(df.f.in$Measurement_Date),
                       ground_elev_ft = as.numeric(df.f.in$GS_Elevation),
                       water_elev_ft = as.numeric(df.f.in$WSE),
                       wtd_ft = as.numeric(df.f.in$GSWS),
                       era = era,
                       catchment = "Navarro",
                       stringsAsFactors=F)
  } else {
    df.f <- data.frame(agency = "CA-DWR",
                       site = well,
                       datetime = mdy_hm(df.f.in$Measurement_Date),
                       ground_elev_ft = as.numeric(df.f.in$GS_Elevation),
                       water_elev_ft = as.numeric(df.f.in$WSE),
                       wtd_ft = as.numeric(df.f.in$GSWS),
                       era = era,
                       catchment = "Navarro",
                       stringsAsFactors=F)
  }

  # add to overall data frame
  if (start.flag){
    df <- df.f
    start.flag <- F
  } else {
    df <- rbind(df.f, df)
  }
}

# correct for datum shift between historic and recent
df.era.comp <- dplyr::summarize(group_by(df, site, era),
                                ground_elev_max = max(ground_elev_ft, na.rm=T),
                                ground_elev_mean = mean(ground_elev_ft, na.rm=T),
                                ground_elev_min = min(ground_elev_ft, na.rm=T))
# for each site/era, max-mean-min are all equal to each other, meaning simple arithmetic shift
df.era.shift <- data.frame(site = subset(df.era.comp, era=="Historic")$site,
                           shift = subset(df.era.comp, era=="Historic")$ground_elev_mean - subset(df.era.comp, era=="Recent")$ground_elev_mean)

# shift elevation accordingly
df <- left_join(df, df.era.shift)
df$shift[df$era=="Recent"] <- 0
df$ground_elev_ft <- df$ground_elev_ft - df$shift
df$water_elev_ft <- df$water_elev_ft - df$shift

# Make plots --------------------------------------------------------------

# timeseries plots for sites within the basin
df.dwr.navarro <- subset(df, catchment=="Navarro")

# annual maximum water table depth, only for years with frequent data points
df.dwr.navarro$year <- year(df.dwr.navarro$datetime)
df.dwr.navarro.max <- dplyr::summarize(group_by(subset(df.dwr.navarro, year<=2013), site, year),
                                       max.wtd_ft = max(wtd_ft, na.rm=T))

# color palette
pal.wells <- 
  c("390047N1233685W001"="#2ca25f", 
    "390047N1233685W002"="#99d8c9", 
    "390106N1233775W001"="#88419d", 
    "390106N1233775W002"="#9ebcda", 
    "390254N1233823W001"="#d7301f", 
    "390254N1233823W002"="#fc8d59", 
    "390254N1233823W003"="#fdcc8a")

# GW elevation plot
p.dwr.wse <-
  ggplot(df.dwr.navarro, aes(x=datetime, y=water_elev_ft*0.3048, color=site)) +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_datetime(name="Date") +
  scale_y_continuous(name="Groundwater Elevation [m]") +
  theme_scz()

# WTD plot
p.dwr.wtd <-
  ggplot(df.dwr.navarro, aes(x=datetime, y=wtd_ft*0.3048, color=site)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_datetime(name="Date") +
  scale_y_reverse(name="Water Table Depth [m]") +
  scale_color_manual(values=pal.wells) +
  theme_scz()

## statistics: significant trend?
df.trend <- data.frame(site = unique(df.dwr.navarro$site),
                       trend.mm_yr = NaN,
                       p = NaN,
                       R2 = NaN,
                       max.trend.mm_yr = NaN,
                       max.p = NaN,
                       max.R2 = NaN)
for (i in 1:length(df.trend$site)){
  # get site
  s <- df.trend$site[i]
  
  # linear fit
  fit.s <- lm(wtd_ft*0.3048*1000 ~ datetime, data=subset(df.dwr.navarro, site==s))
  fit.max.s <- lm(max.wtd_ft*0.3048*1000 ~ year, data=subset(df.dwr.navarro.max, site==s))
  
  # extract data
  df.trend$trend.mm_yr[i] <- coef(fit.s)[2]*86400*365.25
  df.trend$p[i] <- lmp(fit.s)
  df.trend$R2[i] <- summary(fit.s)$r.squared
  df.trend$max.trend.mm_yr[i] <- coef(fit.max.s)[2]
  df.trend$max.p[i] <- lmp(fit.max.s)
  df.trend$max.R2[i] <- summary(fit.max.s)$r.squared
  
}

p.dwr.wtd.max <- 
  ggplot(df.dwr.navarro.max, aes(x=year, y=max.wtd_ft, color=site)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Year") +
  scale_y_reverse(name="Lowest Water Table Depth [m]") +
  scale_color_manual(values=pal.wells) +
  theme_scz()

## save plot
site.legend <- g_legend(p.dwr.wtd.max)

# align scatter and density plots
p1 <- ggplotGrob(p.dwr.wtd + guides(color="none"))
p2 <- ggplotGrob(p.dwr.wtd.max + guides(color="none"))
p <- cbind(p1, p2, size="first")
p$heights <- unit.pmax(p1$heights, p2$heights)

ggsave("results/GroundwaterLevels_CA-DWR_Navarro.png", 
       arrangeGrob(p, site.legend, ncol=2, widths=c(0.75,0.25)), width=12, height=6, units="in")
