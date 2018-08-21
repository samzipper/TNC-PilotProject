## Figure_StreamflowData.R
#' This script is intended to download streamflow data for the Navarro River Watershed
#' and make a figure.
#' Very similar to script Navarro_StreamflowtData.R

source("src/paths+packages.R")

## get data from USGS - output from script Navarro_StreamflowData.R
df <- read.csv(file.path("results", "Navarro_StreamflowData.csv"))
df.info <- read.csv(file.path("results", "Navarro_StreamflowData_siteInfo.csv"))
# df colnames:
#   staid = station (char)
#   val = discharge [cfs] (numeric)
#   dates = date (Date)
#   qualcode = quality code

## data inspection and cleanup
plotParam(df)
cleanUp(df, task = "view")  # appear to be no bad values
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
  summarize(group_by(df, water.year, month),
            discharge.mm_d.mean = mean(discharge.mm_d),
            discharge.mm_d.cum = sum(discharge.mm_d),
            baseflow.mm_d.mean = mean(baseflow.mm_d),
            baseflow.mm_d.cum = sum(baseflow.mm_d))

# summarize by month
df.mo <-
  summarize(group_by(df.yr.mo, month),
            discharge.mm = mean(discharge.mm_d.cum),
            baseflow.mm = mean(baseflow.mm_d.cum),
            quickflow.mm = discharge.mm - baseflow.mm)

# summarize by DOY
df.DOY <- summarize(group_by(df, DOY),
                    discharge.mm_d.mean = mean(discharge.mm_d, na.rm=T),
                    baseflow.mm_d.mean = mean(baseflow.mm_d, na.rm=T))
df.DOY$baseflow.cfs.mean <- df.DOY$baseflow.mm_d.mean/cfs.to.mm

## plots
# mean annual hydrograph
p.Q.DOY <-
  ggplot() +
  geom_line(data=df, aes(x=DOY, y=discharge.mm_d, group=year), color=col.gray, alpha=0.25) +
  geom_line(data=df.DOY, aes(x=DOY, y=discharge.mm_d.mean), color="black", size=2) +
  geom_line(data=df.DOY, aes(x=DOY, y=baseflow.mm_d.mean), color=col.cat.blu, size=2) +
  annotate("text", x=60, y=90, label="Individual Years", color=col.gray) +
  annotate("text", x=60, y=20, label="Mean Annual Hydrograph", color="black") +
  annotate("text", x=60, y=0.3, label="Mean Annual\nBaseflow Hydrograph", color=col.cat.blu) +
  scale_y_continuous(name="Discharge [mm/d]", trans="log10", breaks=(10^seq(-3,2)), labels=c("0.001", "0.01", "0.1", "1", "10", "100")) +
  scale_x_continuous(name="Day of Year", limits=c(0,366), breaks=seq(0,360,90), expand=c(0,0))
ggsave(file.path("figures+tables", "Figure_StreamflowData.png"),
       p.Q.DOY, width=190, height=65, units="mm")

### mean monthly trends - don't need for first paper
## mean monthly trends
# mean month trends table
df.mo.trend <- data.frame(month=seq(1,12),
                          R2 = NaN,
                          p = NaN)
for (mo in seq(1,12)){
  df.mo.trend$R2[mo] <- summary(lm(baseflow.mm_d.mean ~ water.year, subset(df.yr.mo, month==mo)))$r.squared
  df.mo.trend$p[mo] <- lmp(lm(baseflow.mm_d.mean ~ water.year, subset(df.yr.mo, month==mo)))
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

p.baseflow.mo.trend <-
  ggplot(df.yr.mo, aes(x=water.year, y=baseflow.mm_d.mean, group=month)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_point(color=col.cat.blu) +
  stat_smooth(method="lm", color="black") +
  facet_wrap(~month, scales="free_y", labeller=as_labeller(labs.mo.sig)) +
  scale_x_continuous(name="Water Year", breaks=seq(1950, 2010, 20)) +
  scale_y_continuous(name="Mean Monthly Baseflow [mm/d]")

## save plots
ggsave(file.path("figures+tables", "Figure_StreamflowData.png"),
       grid.arrange(p.Q.DOY+labs(title="(a)") + theme(plot.title=element_text(hjust=0.01, vjust=-7)), 
                    p.baseflow.mo.trend+labs(title="(b)") + theme(plot.title=element_text(hjust=0, vjust=-7)), 
                    ncol=1, heights=c(0.75,1)),
       width=190, height=190, units="mm")
