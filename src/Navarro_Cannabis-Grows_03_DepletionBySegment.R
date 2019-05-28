## Navarro_Cannabis-Grows_03_DepletionBySegment.R
#' This script is intended to calculate depletion for all stream segments.
#' 
#' It requires output from Navarro_Cannabis-Grows_02_DepletionApportionment.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

## load depletion apportionment output
df.out <- read.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_02_DepletionApportionment.csv"), stringsAsFactors=F)

### what dates do you want output for? (units: number of days since start of pumping)
DOYs.all <- unique(df.out$time_days)

### load and pre-process data
# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>%
  sf::st_transform(crs.MODFLOW)

# load stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

# grow locations from TNC (filtered and transformed in Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R)
sf.grows <-
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

## for each well: generate a 50-year pumping schedule
start.flag.g <- T
for (g in unique(sf.grows$GrowNum)){
  MonthName <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

  # extract data in gal/month
  df.pump.g <- 
    # extract data from Chris
    data.frame(
      GrowNum = g,
      MonthNum = seq(1,12),
      MonthLengthDays = lubridate::days_in_month(seq(1,12)),
      WaterUseMean_GalMonth = as.numeric(sf.grows[sf.grows$GrowNum==g, paste0(MonthName, "WU_Estimate")])[1:12],  # 1:12 to remove geometry column
      WaterUseUpperCI_GalMonth = as.numeric(sf.grows[sf.grows$GrowNum==g, paste0(MonthName, "WU_Upper95")])[1:12]  # 1:12 to remove geometry column
    ) %>% 
    # convert to m3/day
    transform(
      WaterUseMean_m3d = WaterUseMean_GalMonth*gal.to.m3/MonthLengthDays,
      WaterUseUpperCI_m3d = WaterUseUpperCI_GalMonth*gal.to.m3/MonthLengthDays
    ) %>% 
    # calculate standard deviation
    transform(WaterUseStd_m3d = (WaterUseUpperCI_m3d-WaterUseMean_m3d)/2)
  
  # add to overall output data frame
  if (start.flag.g){
    df.pump <- df.pump.g[,c("GrowNum", "MonthNum", "WaterUseMean_m3d", "WaterUseStd_m3d")]
    start.flag.g <- F
  } else {
    df.pump <- rbind(df.pump, df.pump.g[,c("GrowNum", "MonthNum", "WaterUseMean_m3d", "WaterUseStd_m3d")])
  }
}

# save pumping output
df.pump %>% 
  write.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_03_PumpingRate.csv"),
            row.names=F)

## unique well-seg combos
df.combos <- 
  df.out %>% 
  dplyr::select(SegNum, GrowNum, dist_wellToStream_m, S_bulk, Tr_bulk_m2d, lmda_m2d) %>% 
  unique()

## calculate depletion for all combos
start.flag.Qs <- T
counter <- 0
for (grow in unique(df.combos$GrowNum)){
  # subset combos for only this cultivation site
  df.combos.g <- subset(df.combos, GrowNum==grow)
  
  # set up long data frame for intermittent_pumping script
  t.max.yrs <- max(yrs.model)
  df.pump.long <- 
    subset(df.pump, GrowNum==grow) %>% 
    transform(MonthLengthDays = lubridate::days_in_month(MonthNum)) %>% 
    dplyr::select(MonthNum, MonthLengthDays, WaterUseMean_m3d, WaterUseStd_m3d) %>% 
    replicate(t.max.yrs, ., simplify = FALSE) %>% 
    dplyr::bind_rows() %>% 
    transform(EndOfMonthDays = cumsum(MonthLengthDays))
  df.pump.long$StartOfMonthDays <- c(1, df.pump.long$EndOfMonthDays[1:(t.max.yrs*12)-1]+1)
  
  # see if there are any back-to-back months with same pumping rate which can be compressed - unlikely but couldn't hurt...
  i.starts <- which(c(1, diff(df.pump.long$WaterUseMean_m3d)) != 0)
  i.ends <- c((i.starts-1)[2:length(i.starts)], i.starts[length(i.starts)])
  df.pump.grow <-
    data.frame(StartOfMonthDays = df.pump.long$StartOfMonthDays[i.starts],
               EndOfMonthDays = df.pump.long$EndOfMonthDays[i.ends],
               WaterUseMean_m3d = df.pump.long$WaterUseMean_m3d[i.ends])
  
  for (i in 1:dim(df.combos.g)[1]){
    # identify well-seg combo
    seg <- df.combos.g$SegNum[i]
    
    # get times
    output_t_days <- df.out$time_days[df.out$SegNum==seg & df.out$GrowNum==grow]
    output_frac <- df.out$frac_depletion[df.out$SegNum==seg & df.out$GrowNum==grow]
    
    # pump
    Qs <- intermittent_pumping(t = output_t_days,
                               starts = df.pump.grow$StartOfMonthDays,
                               stops  = df.pump.grow$EndOfMonthDays,
                               rates  = df.pump.grow$WaterUseMean_m3d,
                               method = "hunt",
                               d = df.combos.g$dist_wellToStream_m[i],
                               S = df.combos.g$S_bulk[i],
                               Tr = df.combos.g$Tr_bulk_m2d[i],
                               lmda = df.combos.g$lmda_m2d[i])
    
    # compile output
    df.depletion <- data.frame(SegNum = seg,
                               GrowNum = grow,
                               time_days = output_t_days,
                               Qs = Qs) %>% 
      subset(Qs > 1e-6)
    
    if (start.flag.Qs){
      df.Qs <- df.depletion
      start.flag.Qs <- F
    } else {
      df.Qs <- rbind(df.Qs, df.depletion)
    }
    
    counter <- counter + 1
    print(paste0("Depletion ", counter, " of ", dim(df.combos)[1], " complete, ", Sys.time()))
    
  }
}

# combine and save
df.Qs %>% 
  left_join(df.out, by=c("SegNum", "GrowNum", "time_days")) %>% 
  transform(depletion_m3d = Qs*frac_depletion) %>% 
  dplyr::select(SegNum, GrowNum, time_days, Qs, frac_depletion, depletion_m3d) %>% 
  subset(depletion_m3d >= 1e-6) %>% 
  dfDigits(x=., digits=7) %>% 
  write.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_03_DepletionBySegment.csv"),
            row.names=F)
