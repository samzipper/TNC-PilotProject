## Navarro_GetStreamflowData.R
#' This script is intended to download streamflow data for the Navarro River Watershed
#' and do some simple analysis and plotting.

source(file.path("src", "paths+packages.R"))

## get data from USGS - only have to run once
#df <- importDVs(station.outlet, code="00060", stat="00003", sdate="1950-10-01", edate="2018-09-30")
#df.info <- siteInfo(station.outlet)
#write.csv(df, file.path("results", "Navarro_StreamflowData.csv"), row.names=F)
#write.csv(df.info, file.path("results", "Navarro_StreamflowData_siteInfo.csv"), row.names=F)
df <- read.csv(file.path("results", "Navarro_StreamflowData.csv"))
df.info <- read.csv(file.path("results", "Navarro_StreamflowData_siteInfo.csv"))
df$date <- ymd(df$date)
# df colnames:
#   staid = station (char)
#   val = discharge [cfs] (numeric)
#   dates = date (Date)
#   qualcode = quality code

## data inspection and cleanup
if (df$date[dim(df)[1]] - df$date[1] != dim(df)[1]-1) stop('missing dates')
if (sum(is.na(df$val))>0) stop(paste0('no data: ', paste(df$dates[is.na(df$val)], collapse=", ")))

# gap-fill with linear interpolation
df$val <- na.approx(df$val, maxgap=14)

## unit conversions
# cfs to mm/d
cfs.to.mm <- (0.3048^3)*(0.001^3)*(1/area.km2)*86400*1000*1000
df$discharge.mm_d <- df$val*cfs.to.mm

# calculate baseflow using Nathan & McMahon digital filter
bf <- BaseflowSeparation(df$discharge.mm_d)
df$baseflow.mm_d <- bf[,1]
df$quickflow.mm_d <- bf[,2]

# year and water year
df$year <- year(df$date)
df$water.year <- year(df$date+days(sum(days_in_month(c(10,11,12)))))
df$month <- month(df$date)
df$DOY <- yday(df$date)

# summarize by year and month
df.yr.mo <-
  summarize(group_by(df, year, month),
            discharge.mm_d.mean = mean(discharge.mm_d),
            discharge.mm_d.cum = sum(discharge.mm_d),
            quickflow.mm_d.mean = mean(quickflow.mm_d),
            quickflow.mm_d.cum = sum(quickflow.mm_d),
            baseflow.mm_d.mean = mean(baseflow.mm_d),
            baseflow.mm_d.cum = sum(baseflow.mm_d))
df.yr.mo$date.mid <- ymd(paste0(df.yr.mo$year, "-", df.yr.mo$month, "-", round(days_in_month(df.yr.mo$month)/2)))

# summarize by water year
df.yr.water <-
  summarize(group_by(df, water.year),
            discharge.mm = sum(discharge.mm_d),
            baseflow.mm = sum(baseflow.mm_d))

# summarize by DOY
df.DOY <- summarize(group_by(df, DOY),
                    discharge.mm_d.mean = mean(discharge.mm_d, na.rm=T),
                    quickflow.mm_d.mean = mean(quickflow.mm_d, na.rm=T),
                    baseflow.mm_d.mean = mean(baseflow.mm_d, na.rm=T))

# summarize to monthly means
df.mo <- 
  group_by(df.yr.mo, month) %>% 
  summarize(discharge.mm_mo = mean(discharge.mm_d.cum),
            quickflow.mm_mo = mean(quickflow.mm_d.cum),
            baseflow.mm_mo = mean(baseflow.mm_d.cum))

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
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  scale_y_log10(name="Discharge [mm/d]") +
  scale_x_date(name="Date", date_breaks="10 years", date_labels="%Y", expand=c(0,0)) +
  theme_scz() +
  theme(panel.grid=element_line(color="grey65"))
ggsave("results/streamflow/Navarro_StreamflowData_p.logQ.date.png",
       p.logQ.date, width=6.25, height=6, units="in")

# monthly timeseries: log(discharge) vs date
p.logQ.date.mo <-
  ggplot(df.yr.mo, aes(x=date.mid, y=discharge.mm_d.mean)) + 
  geom_hline(yintercept=baseflow.mm_d, color="red") +
  geom_line(color="blue") +
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  scale_y_log10(name="Mean Monthly Discharge [mm/d]") +
  scale_x_date(name="Date", date_breaks="10 years", date_labels="%Y", expand=c(0,0)) +
  theme_scz() +
  theme(panel.grid=element_line(color="grey65"))
ggsave("results/streamflow/Navarro_StreamflowData_p.logQ.date.mo.png",
       p.logQ.date.mo, width=6.25, height=6, units="in")
  
# water year timeseries
p.Q.yr.water <-
  ggplot(df.yr.water, aes(x=water.year, y=discharge.mm)) + 
  geom_line(color="blue") +
  geom_point(color="blue") +
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  scale_y_continuous(name="Cumulative Discharge [mm]") +
  scale_x_continuous(name="Water Year", expand=c(0,0)) +
  theme_scz() +
  theme(panel.grid=element_line(color="grey65"))
ggsave("results/streamflow/Navarro_StreamflowData_p.Q.yr.water.png",
       p.Q.yr.water, width=6.25, height=6, units="in")

# mean annual hydrograph
p.Q.DOY <-
  ggplot() +
  geom_hline(yintercept=baseflow.mm_d, color="red") +
  geom_line(data=df, aes(x=DOY, y=discharge.mm_d, group=year), color="grey65", alpha=0.25) +
  geom_line(data=df.DOY, aes(x=DOY, y=discharge.mm_d.mean), color="black", size=2) +
  geom_line(data=df.DOY, aes(x=DOY, y=baseflow.mm_d.mean), color="blue", size=2) +
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  scale_y_log10(name="Discharge [mm/d]") +
  scale_x_continuous(name="Day of Year", expand=c(0,0)) +
  theme_scz()
ggsave("results/streamflow/Navarro_StreamflowData_p.Q.DOY.png",
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
ggsave("results/streamflow/Navarro_StreamflowData_p.fdc.png",
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
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Mean Monthly Discharge [mm/d]") +
  theme_scz()
ggsave(paste0("results/streamflow/Navarro_StreamflowData_p.Q.mo.trend.png"),
       p.Q.mo.trend, width=8, height=6, units="in")
