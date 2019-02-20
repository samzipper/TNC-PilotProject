## Navarro_Residential_04_SummarizeDepletion.R
#' This script is intended to look at output from Navarro_Residential_03_DepletionApportionment.R

source(file.path("src", "paths+packages.R"))

#### load residential data
# depletion by segment associated with each well - output from Navarro_Residential_03_DepletionBySegment.R
df.res <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Residential_03_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(time_yrs = time_days/365,
            depletion_m3s = depletion_m3d*86400)

# calculate year and month
df.res$year.dec <- df.res$time_days/365 + 1
df.res$year <- floor(df.res$year.dec)
df.res$month <- round((df.res$year.dec - df.res$year)*12)+1

## summarize by HouseNum and time
df.res.house.sum <-
  df.res %>% 
  group_by(HouseNum, time_days, year.dec) %>% 
  summarize(depletion_m3d_sum = sum(depletion_m3d))

## summarize by time
df.res.sum <-
  df.res %>% 
  group_by(time_days, year.dec) %>% 
  summarize(depletion_m3d_sum = sum(depletion_m3d))


## plots
ggplot(df.res.house.sum, aes(x=year.dec, y=depletion_m3d_sum, group=factor(HouseNum))) +
  geom_line(alpha=0.25)

## plots
ggplot(df.res.sum, aes(x=year.dec, y=depletion_m3d_sum)) +
  geom_line()
