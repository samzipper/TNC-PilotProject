## Navarro_Residential_03_DepletionBySegment.R
#' This script is intended to calculate depletion for all stream segments.
#' 
#' It requires output from Navarro_Residential_02_DepletionApportionment.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

## load depletion apportionment output
df.out <- read.csv(file.path(dir.TNC, "DerivedData", "Navarro_Residential_02_DepletionApportionment.csv"), stringsAsFactors=F)

## define pumping rates
# load in data from state water board
df.pump.northcoast <- 
  file.path("data", "WaterUse", "uw_supplier_data040219.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(date = mdy(Reporting_Month),
            Year = year(mdy(Reporting_Month)),
            MonthNum = month(mdy(Reporting_Month))) %>% 
  subset(R_GPCD_Reported < 5000)  # get rid of some outliers

ggplot(df.pump.northcoast, aes(x=R_GPCD_Calculated, y=R_GPCD_Reported)) + geom_point() + geom_abline(intercept=0, slope=1)
ggplot(df.pump.northcoast, aes(x=MonthNum, y=R_GPCD_Reported)) + geom_point()
ggplot(df.pump.northcoast, aes(x=MonthNum, y=R_GPCD_Calculated)) + geom_point()

# summarize by month
house.size <- 2.65  # assumed residents per house; this middle of range from Mendocino County Water Agency, 2010 report (2.5-2.8 people)

df.pump.region <-
  df.pump.northcoast %>% 
  dplyr::group_by(MonthNum) %>% 
  dplyr::summarize(WaterUseMean_GPCD = mean(R_GPCD_Reported),
                   WaterUseStd_GPCD = sd(R_GPCD_Reported),
                   WaterUseMean_GalHouseDay = mean(R_GPCD_Reported*house.size),
                   WaterUseStd_GalHouseDay = sd(R_GPCD_Reported*house.size)) %>% 
  transform(domain = "Region")

df.pump.closest <- 
  df.pump.northcoast %>% 
  subset(Supplier_Name == "Healdsburg  City of") %>% 
  dplyr::group_by(MonthNum) %>% 
  dplyr::summarize(WaterUseMean_GPCD = mean(R_GPCD_Reported),
                   WaterUseStd_GPCD = sd(R_GPCD_Reported),
                   WaterUseMean_GalHouseDay = mean(R_GPCD_Reported*house.size),
                   WaterUseStd_GalHouseDay = sd(R_GPCD_Reported*house.size)) %>% 
  transform(domain = "Healdsburg")

df.pump <- rbind(df.pump.region, df.pump.closest)

write.csv(df.pump, file.path("results", "Residential_WaterUse.csv"), row.names=F)

ggplot(df.pump, aes(x=MonthNum, y=WaterUseMean_GalHouseDay, color=domain)) + geom_point() +geom_line()

# subset to only healdsburg
df.pump <- subset(df.pump, domain=="Healdsburg")

df.pump$MonthLengthDays <- lubridate::days_in_month(df.pump$MonthNum)
df.pump$m3HouseDay <- df.pump$WaterUseMean_GalHouseDay*gal.to.m3

# set up long data frames for intermittent_pumping script
t.max.yrs <- ceiling(max(df.out$time_days/365))
df.pump.long <- 
  df.pump %>% 
  dplyr::select(MonthNum, MonthLengthDays, m3HouseDay) %>% 
  replicate(t.max.yrs, ., simplify = FALSE) %>% 
  dplyr::bind_rows() %>% 
  transform(EndOfMonthDays = cumsum(MonthLengthDays))
df.pump.long$StartOfMonthDays <- c(1, df.pump.long$EndOfMonthDays[1:(t.max.yrs*12)-1]+1)

# start of each pumping period based on difference
i.starts <- which(c(1, diff(df.pump.long$m3HouseDay)) != 0)
i.ends <- c((i.starts-1)[2:length(i.starts)], i.starts[length(i.starts)])

df.pump.compressed <-
  data.frame(StartOfMonthDays = df.pump.long$StartOfMonthDays[i.starts],
             EndOfMonthDays = df.pump.long$EndOfMonthDays[i.ends],
             m3HouseDay = df.pump.long$m3HouseDay[i.ends])

## unique well-seg combos
df.combos <- 
  df.out %>% 
  dplyr::select(SegNum, HouseNum, dist_wellToStream_m, S_bulk, Tr_bulk_m2d, lmda_m2d) %>% 
  unique()

## calculate depletion for all combos
start.flag.Qs <- T
for (i in 1:dim(df.combos)[1]){
  # identify well-seg combo
  seg <- df.combos$SegNum[i]
  house <- df.combos$HouseNum[i]
  
  # get times
  output_t_days <- df.out$time_days[df.out$SegNum==seg & df.out$HouseNum==house]
  output_frac <- df.out$frac_depletion[df.out$SegNum==seg & df.out$HouseNum==house]
  
  # different calculation depending on outdoor vs indoor plants
  Qs <- intermittent_pumping(t = output_t_days,
                             starts = df.pump.compressed$StartOfMonthDays,
                             stops  = df.pump.compressed$EndOfMonthDays,
                             rates  = df.pump.compressed$m3HouseDay,
                             method = "hunt",
                             d = df.combos$dist_wellToStream_m[i],
                             S = df.combos$S_bulk[i],
                             Tr = df.combos$Tr_bulk_m2d[i],
                             lmda = df.combos$lmda_m2d[i])
  
  # compile output
  df.depletion <- data.frame(SegNum = seg,
                             HouseNum = house,
                             time_days = output_t_days,
                             Qs = Qs) %>% 
    subset(Qs > 1e-6)
  
  if (start.flag.Qs){
    df.Qs <- df.depletion
    start.flag.Qs <- F
  } else {
    df.Qs <- rbind(df.Qs, df.depletion)
  }
  
  print(paste0("Depletion ", i, " of ", dim(df.combos)[1], " complete, ", Sys.time()))
  
}

# combine and save
df.Qs %>% 
  left_join(df.out, by=c("SegNum", "HouseNum", "time_days")) %>% 
  transform(depletion_m3d = Qs*frac_depletion) %>% 
  dplyr::select(SegNum, HouseNum, time_days, Qs, frac_depletion, depletion_m3d) %>% 
  subset(depletion_m3d >= 1e-6) %>% 
  dfDigits(x=., digits=7) %>% 
  write.csv(file.path(dir.TNC, "DerivedData", "Navarro_Residential_03_DepletionBySegment.csv"),
            row.names=F)
