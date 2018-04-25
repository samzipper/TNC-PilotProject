## Figure_MetData.R
#' This figure plots meteorological data for the Navarro River Watershed.
#' Very similar to script Navarro_MetData.R

source(file.path("src", "paths+packages.R"))

## load raw data downloaded from NOAA CDO
# input units are metric
df.raw <- read.csv("data/MetData_GHCND+CoCoRaHS.csv", stringsAsFactors=F)

## Extract useful columns and summarize by date
df <- 
  data.frame(station = df.raw$STATION,
             date = ymd(df.raw$DATE),
             prec = df.raw$PRCP,
             Tmax = df.raw$TMAX,
             Tmin = df.raw$TMIN) %>% 
  group_by(date) %>% 
  summarize(prec = mean(prec, na.rm=T),
            Tmax = mean(Tmax, na.rm=T),
            Tmin = mean(Tmin, na.rm=T))
df$year <- year(df$date)
df$month <- month(df$date)

## gap-filling < 7 days via linear interpolation
max.gap.days <- 7
df$prec <- as.numeric(na.approx(df$prec, maxgap=max.gap.days, na.rm=F))
df$Tmax <- as.numeric(na.approx(df$Tmax, maxgap=max.gap.days, na.rm=F))
df$Tmin <- as.numeric(na.approx(df$Tmin, maxgap=max.gap.days, na.rm=F))

## summarize to monthly
df.mo <- 
  df %>% 
  group_by(year, month) %>% 
  summarize(mo.days = mean(days_in_month(month)),
            date.mid = mean(date),
            prec.days = sum(is.finite(prec)),
            prec.sum = sum(prec, na.rm=T),
            Tmax.days = sum(is.finite(Tmax)),
            Tmax.mean = mean(Tmax, na.rm=T),
            Tmin.days = sum(is.finite(Tmin)),
            Tmin.mean = mean(Tmin, na.rm=T))

# make month a factor
df.mo$month <- factor(df.mo$month, levels=seq(1,12), labels=labs.mo)

## calculate monthly PET using Hargreaves formulation (temperature-based)
df.mo$PET <- hargreaves(df.mo$Tmin.mean, df.mo$Tmax.mean, lat=station.outlet.lat, na.rm=T)
df.mo$defc <- df.mo$PET - df.mo$prec.sum

## select months to keep for prec and temp
missing.day.threshold <- 3
df.mo$prec.keep <- (df.mo$mo.days - df.mo$prec.days) < missing.day.threshold
df.mo$Tmax.keep <- (df.mo$mo.days - df.mo$Tmax.days) < missing.day.threshold
df.mo$Tmin.keep <- (df.mo$mo.days - df.mo$Tmin.days) < missing.day.threshold

## summarize everything to monthly means
df.mo.mean <- 
  df.mo %>% 
  group_by(month) %>% 
  summarize(prec = mean(prec.sum, na.rm=T),
            Tmax = mean(Tmax.mean, na.rm=T),
            PET = mean(PET, na.rm=T),
            defc = PET - prec) %>% 
  melt(id=c("month"))

## make some summary plots
p.mo.mean <-
  ggplot(df.mo.mean, aes(x=month, y=value)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity") +
  facet_wrap(~variable, ncol=2, scales="free_y", 
             labeller=as_labeller(c("prec"="Precipitation [mm/mo]",
                                    "Tmax"="Max Daily Temperature [C]",
                                    "PET"="Potential ET [mm/mo]",
                                    "defc"="Precip. Deficit (PET - Precip) [mm/mo]"))) +
  scale_y_continuous(name="Long-Term Monthly Mean") +
  scale_x_discrete(name="Month")
ggsave(file.path("figures+tables", "Figure_MetData.png"), p.mo.mean,
       width=190, height=95, units="mm")