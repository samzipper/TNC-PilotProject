## Navarro_Residential_03_DepletionBySegment.R
#' This script is intended to calculate depletion for all stream segments.
#' 
#' It requires output from Navarro_Residential_02_DepletionApportionment.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

## load depletion apportionment output
df.out <- read.csv(file.path(dir.TNC, "DerivedData", "Navarro_Residential_02_DepletionApportionment.csv"), stringsAsFactors=F)

## define pumping rates: 250 gpm in winter, 500 gpm in summer
df.pump <- data.frame(
  Month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  MeanWaterUse_GalHouseDay = c(250, 250, 250, 250, 500, 500, 500, 500, 500, 500, 250, 250)
)
df.pump$Month <- factor(df.pump$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump$MonthNum <- match(df.pump$Month, month.abb)
df.pump$MonthLengthDays <- lubridate::days_in_month(df.pump$MonthNum)
df.pump$m3HouseDay <- df.pump$MeanWaterUse_GalHouseDay*gal.to.m3

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
                             Qs = Qs)
  
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
