## Figure_WatershedScale_Depletion.R

source(file.path("src", "paths+packages.R"))

## get data from USGS - output from script Navarro_StreamflowData.R
df.Q <- read.csv(file.path("results", "Navarro_StreamflowData.csv"), stringsAsFactors=F)
df.Q.info <- read.csv(file.path("results", "Navarro_StreamflowData_siteInfo.csv"), stringsAsFactors=F)
# df.Q colnames:
#   staid = station (char)
#   val = discharge [cfs] (numeric)
#   dates = date (Date)
#   qualcode = quality code

## data inspection and cleanup
if (sum(is.na(df.Q$val))>0) stop(paste0('no data: ', paste(df.Q$dates[is.na(df.Q$val)], collapse=", ")))

# gap-fill with linear interpolation
df.Q$val <- na.approx(df.Q$val, maxgap=14)

## unit conversions
# cfs to m3/d
df.Q$discharge_m3d <- df.Q$val*cfs.to.m3d

# calculate baseflow using Nathan & McMahon digital filter
bf <- BaseflowSeparation(df.Q$discharge_m3d)
df.Q$baseflow_m3d <- bf[,1]
df.Q$quickflow_m3d <- bf[,2]

# year and water year
df.Q$date <- ymd(df.Q$date)
df.Q$year <- year(df.Q$date)
df.Q$water.year <- year(df.Q$date+days(sum(days_in_month(c(10,11,12)))))
df.Q$month <- month(df.Q$date)
df.Q$DOY <- yday(df.Q$date)

# summarize by year and month
df.Q.yr.mo <-
  summarize(group_by(df.Q, water.year, month),
            discharge_m3d.mean = mean(discharge_m3d),
            discharge_m3d.cum = sum(discharge_m3d),
            baseflow_m3d.mean = mean(baseflow_m3d),
            baseflow_m3d.cum = sum(baseflow_m3d))

# summarize by month
df.Q.mo <-
  summarize(group_by(df.Q.yr.mo, month),
            baseflow_m3d = mean(baseflow_m3d.mean),
            baseflow_m3d.std = sd(baseflow_m3d.mean),
            baseflow_m3d.min = min(baseflow_m3d.mean)) %>% 
  transform(baseflow_m3d.ribbon.min = baseflow_m3d - baseflow_m3d.std,
            baseflow_m3d.ribbon.max = baseflow_m3d + baseflow_m3d.std)

# summarize by DOY
df.Q.DOY <- summarize(group_by(df.Q, DOY),
                      discharge_m3d.mean = mean(discharge_m3d, na.rm=T),
                      baseflow_m3d.mean = mean(baseflow_m3d, na.rm=T))

## plots
# mean annual hydrograph
p.Q.DOY <-
  ggplot() +
  geom_line(data=df.Q, aes(x=DOY, y=discharge_m3d/86400, group=year), color=col.gray, alpha=0.15) +
  geom_line(data=df.Q.DOY, aes(x=DOY, y=discharge_m3d.mean/86400), color="black", size=2) +
  geom_line(data=df.Q.DOY, aes(x=DOY, y=baseflow_m3d.mean/86400), color=col.cat.blu, size=2) +
  annotate("text", x=60, y=0.9, label="Individual Years", color=col.gray) +
  annotate("text", x=140, y=35, label="Mean Annual\nHydrograph", color="black") +
  annotate("text", x=140, y=0.2, label="Mean Annual\nBaseflow Hydrograph", color=col.cat.blu) +
  scale_y_continuous(name="Discharge [m3/s]", trans="log10", breaks=(10^seq(-2,3)), labels=c("0.01", "0.1", "1", "10", "100", "1000")) +
  scale_x_continuous(name="Day of Year", limits=c(0,366), breaks=seq(0,360,90), expand=c(0,0))

## mean monthly trends
# mean month trends table
df.Q.mo.trend <- data.frame(month=seq(1,12),
                            R2 = NaN,
                            p = NaN,
                            slope = NaN)
for (mo in seq(1,12)){
  df.Q.mo.trend$R2[mo] <- summary(lm(baseflow_m3d.mean ~ water.year, subset(df.Q.yr.mo, month==mo)))$r.squared
  df.Q.mo.trend$p[mo] <- lmp(lm(baseflow_m3d.mean ~ water.year, subset(df.Q.yr.mo, month==mo)))
  df.Q.mo.trend$slope[mo] <- coef(lm(baseflow_m3d.mean ~ water.year, subset(df.Q.yr.mo, month==mo)))[2]
}

# significance label
df.Q.mo.trend$sig.label <- ""
df.Q.mo.trend$sig.label[df.Q.mo.trend$p<0.05] <- "*"
df.Q.mo.trend$sig.label[df.Q.mo.trend$p<0.01] <- "**"
df.Q.mo.trend$sig.label[df.Q.mo.trend$p<0.001] <- "***"
labs.mo.sig <- labs.mo
for (mo in seq(1,12)){
  labs.mo.sig[mo] <- paste0(labs.mo[mo], df.Q.mo.trend$sig.label[mo])
}

p.baseflow.mo.trend <-
  ggplot(df.Q.yr.mo, aes(x=water.year, y=baseflow_m3d.mean*days_in_month(month)/86400, group=month)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_point(color=col.cat.blu) +
  stat_smooth(method="lm", color="black") +
  facet_wrap(~month, scales="free_y", labeller=as_labeller(labs.mo.sig)) +
  scale_x_continuous(name="Water Year", breaks=seq(1950, 2010, 20)) +
  scale_y_continuous(name="Monthly Baseflow [m3/s]")

## save plots
ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_StreamflowData.png"),
       grid.arrange(p.Q.DOY+labs(title="(a)") + theme(plot.title=element_text(hjust=0.01, vjust=-7, face="plain"),
                                                      plot.margin=unit(c(-5,1,0,1), "mm")), 
                    p.baseflow.mo.trend+labs(title="(b)") + theme(plot.title=element_text(hjust=0, vjust=-7, face="plain"),
                                                                  plot.margin=unit(c(-7,1,0,1), "mm")), 
                    ncol=1, heights=c(0.75,1)),
       width=190, height=160, units="mm")