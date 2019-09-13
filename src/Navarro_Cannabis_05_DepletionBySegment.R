## Navarro_Cannabis_05_DepletionBySegment.R
#' This script is intended to calculate depletion for stream segments.
#' 
#' It requires output from Navarro_Cannabis_04_DepletionApportionment.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

# pumping factors for sensitivity analysis (to multiply mean monthly pumping rate by)
pump_factors <- c(0.5, 1, 1.5)

## load depletion apportionment output
df.out <- read.csv(file.path("results", "Navarro_Cannabis_04_DepletionApportionment.csv"), stringsAsFactors=F)

### what dates do you want output for? (units: number of days since start of pumping)
# convert years and dates to DOY since pumping started
yrs.model <- c(1, 10, 50) 
DOYs.model <- yday(c("2017-09-15"))
DOYs.all <- rep(DOYs.model, times=length(yrs.model)) + rep(365*(yrs.model-1), each=length(DOYs.model))
t.max.yrs <- max(yrs.model)

### load and pre-process data
# domain boundary shapefile
sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

# load stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F)

# synthetic pumping wells
sf.wel <- 
  file.path("modflow", "input", "iwel.txt") %>% 
  read.table(header=T) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(sf.streams))

## generate a 50-year pumping schedule
# grow locations from TNC (filtered and transformed in Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R)
#  these are used to determine mean monthly water requirements
sf.grows <-
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

MonthName <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# extract data in gal/month
df.pump <- 
  # extract data from Chris
  data.frame(
    MonthNum = seq(1,12),
    MonthLengthDays = lubridate::days_in_month(seq(1,12)),
    WaterUseMean_GalMonth = c(
      mean(as.numeric(sf.grows$JanWU_Estimate)),
      mean(as.numeric(sf.grows$FebWU_Estimate)),
      mean(as.numeric(sf.grows$MarWU_Estimate)),
      mean(as.numeric(sf.grows$AprWU_Estimate)),
      mean(as.numeric(sf.grows$MayWU_Estimate)),
      mean(as.numeric(sf.grows$JunWU_Estimate)),
      mean(as.numeric(sf.grows$JulWU_Estimate)),
      mean(as.numeric(sf.grows$AugWU_Estimate)),
      mean(as.numeric(sf.grows$SepWU_Estimate)),
      mean(as.numeric(sf.grows$OctWU_Estimate)),
      mean(as.numeric(sf.grows$NovWU_Estimate)),
      mean(as.numeric(sf.grows$DecWU_Estimate)))
    ) %>% 
  # convert to m3/day
  transform(WaterUseMean_m3d = WaterUseMean_GalMonth*gal.to.m3/MonthLengthDays)

# conver to 50 yr long
df.pump.average <- 
  df.pump %>% 
  dplyr::select(MonthNum, MonthLengthDays, WaterUseMean_m3d) %>% 
  replicate(t.max.yrs, ., simplify = FALSE) %>% 
  dplyr::bind_rows() %>% 
  transform(EndOfMonthDays = cumsum(MonthLengthDays))
df.pump.average$StartOfMonthDays <- c(1, df.pump.average$EndOfMonthDays[1:(t.max.yrs*12)-1]+1)

## unique well-seg combos
df.combos <- 
  df.out %>% 
  dplyr::select(SegNum, WellNum, dist_wellToStream_m, S_bulk, Tr_bulk_m2d, lmda_m2d) %>% 
  unique()

## calculate depletion for all combos
start.flag.Qs <- T
for (i in 1:dim(df.combos)[1]){
  # identify well-seg combo
  seg <- df.combos$SegNum[i]
  w <- df.combos$WellNum[i]
  
  # get times
  output_t_days <- df.out$time_days[df.out$SegNum==seg & df.out$WellNum==w]
  output_frac <- df.out$frac_depletion[df.out$SegNum==seg & df.out$WellNum==w]
  
  for (p in pump_factors){
    
    # use 'average' grow for depletion calculation
    Qs <- intermittent_pumping(t = output_t_days,
                               starts = df.pump.average$StartOfMonthDays,
                               stops  = df.pump.average$EndOfMonthDays,
                               rates  = df.pump.average$WaterUseMean_m3d*p,
                               method = "hunt",
                               d = df.combos$dist_wellToStream_m[i],
                               S = df.combos$S_bulk[i],
                               Tr = df.combos$Tr_bulk_m2d[i],
                               lmda = df.combos$lmda_m2d[i])
    
    # compile output
    df.depletion <- data.frame(SegNum = seg,
                               WellNum = w,
                               time_days = output_t_days,
                               Qs = Qs,
                               pump_factor = p)
    
    if (start.flag.Qs){
      df.Qs <- df.depletion
      start.flag.Qs <- F
    } else {
      df.Qs <- rbind(df.Qs, df.depletion)
    }
  }
  
  print(paste0("Depletion ", i, " of ", dim(df.combos)[1], " complete, ", Sys.time()))
  
}

# combine and save
df.Qs %>% 
  left_join(df.out, by=c("SegNum", "WellNum", "time_days")) %>% 
  transform(depletion_m3d = Qs*frac_depletion) %>% 
  dplyr::select(SegNum, WellNum, time_days, pump_factor, Qs, frac_depletion, depletion_m3d) %>% 
  dfDigits(x=., digits=7) %>% 
  write.csv(file.path("results", "Navarro_Cannabis_05_DepletionBySegment.csv"),
            row.names=F)
