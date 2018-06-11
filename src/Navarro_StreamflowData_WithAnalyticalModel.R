## Navarro_GetStreamflowData.R
#' This script is intended to download streamflow data for the Navarro River Watershed
#' and do some simple analysis and plotting.

source(file.path("src", "paths+packages.R"))

## get data from USGS - only have to run once
#df <- importDVs(station.outlet, code="00060", stat="00003", sdate="1900-01-01", edate="2017-12-31")
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

# summarize by DOY
df.DOY <- summarize(group_by(df, DOY),
                    discharge.mm_d.mean = mean(discharge.mm_d, na.rm=T),
                    quickflow.mm_d.mean = mean(quickflow.mm_d, na.rm=T),
                    baseflow.mm_d.mean = mean(baseflow.mm_d, na.rm=T))

# mean annual hydrograph
ggplot() +
  geom_hline(yintercept=baseflow.mm_d, color="red") +
  geom_line(data=df, aes(x=DOY, y=discharge.mm_d, group=year), color="grey65", alpha=0.25) +
  geom_line(data=df.DOY, aes(x=DOY, y=discharge.mm_d.mean), color="black", size=2) +
  geom_line(data=df.DOY, aes(x=DOY, y=baseflow.mm_d.mean), color="blue", size=2) +
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  scale_y_log10(name="Discharge [mm/d]") +
  scale_x_continuous(name="Day of Year", expand=c(0,0)) +
  theme_scz() + 
  ggsave("results/streamflow/Navarro_StreamflowData_p.Q.DOY.png",
         width=6.25, height=6, units="in")

# mean annual baseflow hydrograph with presumptive standard
ggplot() +
#  geom_hline(yintercept=baseflow.mm_d, color="red") +
  geom_line(data=df.DOY, aes(x=DOY, y=discharge.mm_d.mean), color="black", size=2) +
  geom_line(data=df.DOY, aes(x=DOY, y=baseflow.mm_d.mean), color="blue", size=2) +
  geom_ribbon(data=df.DOY, aes(x=DOY, ymin=baseflow.mm_d.mean*0.8, ymax=baseflow.mm_d.mean*1.2), fill="blue", alpha=0.5) +
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  #scale_y_continuous(name="Baseflow [mm/d]") +
  scale_y_log10(name="Baseflow [mm/d]") +
  scale_x_continuous(name="Day of Year", expand=c(0,0)) +
  theme_scz()

# calculate depletion through time
n.years <- 3
mo <- c(rep(1,31), rep(2,28), rep(3,31), rep(4,30), rep(5,31), rep(6,30), 
        rep(7,31), rep(8,31), rep(9,30), rep(10,31), rep(11,30), rep(12,31))
Qf <- glover(t=seq(1,365*n.years), d=400, S=0.1, Tr=100*0.1)  # depletion fraction
Qw <- 0.02  # pump rate [mm/d]
Qs <- Qf*Qw  # streamflow depletion [mm/d]
df.depl <- data.frame(year = rep(seq(10, by=10, length.out=n.years), each=365),
                      month = rep(mo, times=n.years),
                      DOY = rep(seq(1,365), times=n.years),
                      discharge = rep(df.DOY$discharge.mm_d.mean[1:365], times=n.years),
                      baseflow = rep(df.DOY$baseflow.mm_d.mean[1:365], times=n.years),
                      Qf = Qf,
                      Qs = Qs)

df.depl.mo <- 
  df.depl %>% 
  group_by(year, month) %>% 
  summarize(discharge.mo = mean(discharge),
            baseflow.mo = mean(baseflow),
            Qs.mo = mean(Qs))

# monthly
ggplot() +
  #  geom_hline(yintercept=baseflow.mm_d, color="red") +
  geom_line(data=df.depl.mo[1:12, ], aes(x=month, y=discharge.mo), color="black", size=2) +
  geom_ribbon(data=df.depl.mo[1:12, ], aes(x=month, ymin=baseflow.mo*0.8, ymax=baseflow.mo), fill="blue", alpha=0.5) +
  geom_line(data=df.depl.mo, aes(x=month, y=baseflow.mo-Qs.mo, color=factor(year)), size=1) +
  geom_line(data=df.depl.mo[1:12, ], aes(x=month, y=baseflow.mo), color="blue", size=2) +
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  #scale_y_continuous(name="Baseflow [mm/d]") +
  scale_y_log10(name="Baseflow [mm/d]") +
  scale_x_continuous(name="Month", expand=c(0,0), breaks=seq(1,12)) +
  scale_color_manual(name="Year Since Start of Pumping",
                     values=c("yellow", "orange", "red")) +
  theme_scz() +
  theme(legend.position="bottom") +
  ggsave("results/streamflow/Navarro_StreamflowData_WithAnalyticalModel_Monthly.png",
         width=6.25, height=6, units="in")

# daily
ggplot() +
  #  geom_hline(yintercept=baseflow.mm_d, color="red") +
  geom_line(data=df.DOY, aes(x=DOY, y=discharge.mm_d.mean), color="black", size=2) +
  geom_line(data=df.DOY, aes(x=DOY, y=baseflow.mm_d.mean), color="blue", size=2) +
  geom_ribbon(data=df.DOY, aes(x=DOY, ymin=baseflow.mm_d.mean*0.8, ymax=baseflow.mm_d.mean*1.2), fill="blue", alpha=0.5) +
  geom_line(data=df.depl, aes(x=DOY, y=baseflow-Qs, color=factor(year))) +
  labs(title=paste0("USGS ", station.outlet, ": ", station.outlet.name)) +
  #scale_y_continuous(name="Baseflow [mm/d]") +
  scale_y_log10(name="Baseflow [mm/d]") +
  scale_x_continuous(name="Day of Year", expand=c(0,0)) +
  scale_color_discrete(name="Year Since Start of Pumping") +
  theme_scz() +
  theme(legend.position="bottom")
