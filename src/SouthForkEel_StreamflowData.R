## SouthForkEel_GetStreamflowData.R
#' This script is intended to download streamflow data for the SouthForkEel River Watershed
#' and do some simple analysis and plotting.

source("src/paths+packages.R")

## first: make a plot showing gauging station locations
# load geospatial data
shp <- readOGR(dsn="data/NHD/WBD", layer="WBDHU8_SouthForkEel")
shp.riv <- readOGR(dsn="data/NHD/HYD", layer="NHDFlowline_SouthForkEel")
r.dem <- raster("results/GIS/SouthForkEel_DEM_m_HydroSHEDS_15s.tif")

# transform for plot
df.shp <- tidy(spTransform(shp, crs.WGS))
df.riv <- tidy(spTransform(shp.riv, crs.WGS))
df.dem <- as.data.frame(rasterToPoints(r.dem))
colnames(df.dem) <- c("lon", "lat", "dem")

# plot
p.map.stations <- 
  ggplot() +
  geom_raster(data=df.dem, aes(x=lon, y=lat, fill=dem)) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="white") +
  geom_polygon(data=df.shp, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_point(data=df.info, aes(x=lng, y=lat), color="red") +
  geom_label(data=df.info, aes(x=lng, y=lat, label=staid), color="red", hjust=0.5, vjust=1.25) +
  coord_equal() +
  labs(title=paste0("South Fork Eel River, HUC " , HUC)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_fill_viridis(name="Elevation [m]") +
  theme_scz() +
  theme(panel.border=element_blank(),
        legend.position=c(0.05,0.05),
        legend.justification=c(0,0))
ggsave("results/streamflow/SouthForkEel_Map_GagingStations.png",
       p.map.stations, width=5.25, height=6, units="in")

## second: scroll through each station, plot
for (s in 1:length(df.info$staid)){
  
  # station
  station <- df.info$staid[s]
  station.name <- df.info$staname[s]
  baseflow.mm_d <- df.info$baseflow.mm_d[s]
  
  ## get data from USGS
  df <- importDVs(station, code="00060", stat="00003", sdate="1900-01-01", edate="2017-12-31")
  # df colnames:
  #   staid = station (char)
  #   val = discharge [cfs] (numeric)
  #   dates = date (Date)
  #   qualcode = quality code
  
  ## data inspection and cleanup
  plotParam(df)
  cleanUp(df, task = "view")  # appear to be no bad values
  if (df$date[dim(df)[1]] - df$date[1] != dim(df)[1]-1) stop('missing dates')
  #if (sum(is.na(df$val))>0) stop(paste0('no data: ', paste(df$dates[is.na(df$val)], collapse=", ")))
  
  ## unit conversions
  # cfs to mm/d
  cfs.to.mm <- (0.3048^3)*(0.001^3)*(1/area.km2)*86400*1000*1000
  df$discharge.mm_d <- df$val*cfs.to.mm
  
  # year and water year
  df$year <- year(df$date)
  df$water.year <- year(df$date+days(sum(days_in_month(c(10,11,12)))))
  df$month <- month(df$date)
  df$DOY <- yday(df$date)
  
  # summarize by year and month
  df.yr.mo <-
    summarize(group_by(df, year, month),
              discharge.mm_d.mean = mean(discharge.mm_d),
              discharge.mm_d.cum = sum(discharge.mm_d))
  df.yr.mo$date.mid <- ymd(paste0(df.yr.mo$year, "-", df.yr.mo$month, "-", round(days_in_month(df.yr.mo$month)/2)))
  
  # summarize by water year
  df.yr.water <-
    summarize(group_by(df, water.year),
              discharge.mm = sum(discharge.mm_d))
  
  # summarize by DOY
  df.DOY <- summarize(group_by(df, DOY),
                      discharge.mm_d.mean = mean(discharge.mm_d, na.rm=T))
  
  # flow-duration curve
  df.fdc <- df[order(-df$discharge.mm_d), ]
  df.fdc <- subset(df.fdc, is.finite(discharge.mm_d))
  df.fdc$rank <- seq(1,dim(df.fdc)[1])
  df.fdc$exceed.prob <- 100*df.fdc$rank/(dim(df.fdc)[1]+1)
  
  ## plots
  # daily timeseries: log(discharge) vs date
  p.logQ.date <- 
    ggplot(df, aes(x=dates, y=discharge.mm_d)) + 
    geom_hline(yintercept=baseflow.mm_d, color="red") +
    geom_line(color="blue") +
    labs(title=paste0("USGS ", station, ": ", station.name)) +
    scale_y_log10(name="Discharge [mm/d]") +
    scale_x_date(name="Date", date_breaks="10 years", date_labels="%Y", expand=c(0,0)) +
    theme_scz() +
    theme(panel.grid=element_line(color="grey65"))
  ggsave(paste0("results/streamflow/SouthForkEel-", station, "_StreamflowData_p.logQ.date.png"),
         p.logQ.date, width=6.25, height=6, units="in")
  
  # monthly timeseries: log(discharge) vs date
  p.logQ.date.mo <-
    ggplot(df.yr.mo, aes(x=date.mid, y=discharge.mm_d.mean)) + 
    geom_hline(yintercept=baseflow.mm_d, color="red") +
    geom_line(color="blue") +
    labs(title=paste0("USGS ", station, ": ", station.name)) +
    scale_y_log10(name="Mean Monthly Discharge [mm/d]") +
    scale_x_date(name="Date", date_breaks="10 years", date_labels="%Y", expand=c(0,0)) +
    theme_scz() +
    theme(panel.grid=element_line(color="grey65"))
  ggsave(paste0("results/streamflow/SouthForkEel-", station, "_StreamflowData_p.logQ.date.mo.png"),
         p.logQ.date.mo, width=6.25, height=6, units="in")
  
  # water year timeseries
  p.Q.yr.water <-
    ggplot(df.yr.water, aes(x=water.year, y=discharge.mm)) + 
    geom_line(color="blue") +
    geom_point(color="blue") +
    labs(title=paste0("USGS ", station, ": ", station.name)) +
    scale_y_continuous(name="Cumulative Discharge [mm]") +
    scale_x_continuous(name="Water Year", expand=c(0,0)) +
    theme_scz() +
    theme(panel.grid=element_line(color="grey65"))
  ggsave(paste0("results/streamflow/SouthForkEel-", station, "_StreamflowData_p.Q.yr.water.png"),
         p.Q.yr.water, width=6.25, height=6, units="in")
  
  # mean annual hydrograph
  p.Q.DOY <-
    ggplot() +
    geom_hline(yintercept=baseflow.mm_d, color="red") +
    geom_line(data=df, aes(x=DOY, y=discharge.mm_d, group=year), color="grey65", alpha=0.25) +
    geom_line(data=df.DOY, aes(x=DOY, y=discharge.mm_d.mean), color="blue") +
    labs(title=paste0("USGS ", station, ": ", station.name)) +
    scale_y_log10(name="Discharge [mm/d]") +
    scale_x_continuous(name="Day of Year", expand=c(0,0)) +
    theme_scz()
  ggsave(paste0("results/streamflow/SouthForkEel-", station, "_StreamflowData_p.Q.DOY.png"),
         p.Q.DOY, width=6.25, height=6, units="in")
  
  # flow duration curve
  p.fdc <-
    ggplot(df.fdc, aes(y=discharge.mm_d, x=exceed.prob)) +
    geom_hline(yintercept=baseflow.mm_d, color="red") +
    geom_point(color="blue", shape=21) +
    scale_y_log10(name="Discharge [mm/d]") + 
    scale_x_continuous("Percent Exceedance", limits=c(0,100), breaks=seq(0,100,25)) +
    theme_scz() +
    theme(panel.grid=element_line(color="grey65"))
  ggsave(paste0("results/streamflow/SouthForkEel-", station, "_StreamflowData_p.fdc.png"),
         p.fdc, width=6.25, height=6, units="in")
  
  ## mean monthly trends
  # mean month trends table
  df.mo.trend <- data.frame(month=seq(1,12),
                            R2 = NaN,
                            p = NaN)
  for (mo in seq(1,12)){
    df.mo.trend$R2[mo] <- summary(lm(discharge.mm_d.mean ~ year, subset(df.yr.mo, month==mo)))$r.squared
    df.mo.trend$p[mo] <- lmp(lm(discharge.mm_d.mean ~ year, subset(df.yr.mo, month==mo)))
  }
  
  # significance label
  df.mo.trend$sig.label <- ""
  df.mo.trend$sig.label[df.mo.trend$p<0.05] <- "*"
  df.mo.trend$sig.label[df.mo.trend$p<0.01] <- "**"
  df.mo.trend$sig.label[df.mo.trend$p<0.001] <- "***"
  labs.mo.sig <- labs.mo
  for (mo in seq(1,12)){
    labs.mo.sig[mo] <- paste0(labs.mo[mo], df.mo.trend$sig.label[mo])
  }
  
  p.Q.mo.trend <-
    ggplot(df.yr.mo, aes(x=year, y=discharge.mm_d.mean)) +
    geom_hline(yintercept=0, color="gray65") +
    geom_hline(yintercept=baseflow.mm_d, color="red") +
    geom_point() +
    stat_smooth(method="lm") +
    facet_wrap(~month, scales="free_y", labeller=as_labeller(labs.mo.sig)) +
    labs(title=paste0("USGS ", station, ": ", station.name)) +
    scale_x_continuous(name="Year", expand=c(0,0)) +
    scale_y_continuous(name="Mean Monthly Discharge [mm/d]") +
    theme_scz()
  ggsave(paste0("results/streamflow/SouthForkEel-", station, "_StreamflowData_p.Q.mo.trend.png"),
         p.Q.mo.trend, width=8, height=6, units="in")
  
}