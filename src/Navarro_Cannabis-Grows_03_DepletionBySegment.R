## Navarro_Cannabis-Grows_03_DepletionBySegment.R
#' This script is intended to calculate depletion for all stream segments.
#' 
#' It requires output from Navarro_Cannabis-Grows_02_DepletionApportionment.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

## load depletion apportionment output
df.out <- read.csv(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_02_DepletionApportionment.csv"), stringsAsFactors=F)

### what dates do you want output for? (units: number of days since start of pumping)
# convert years and dates to DOY since pumping started
yrs.model <- c(1, 2, 3, 4, 5, 10, 15, 20, 30, 50) 
DOYs.model <- yday(c("2017-01-15", "2017-02-14", "2017-03-15", "2017-04-15", "2017-05-15", "2017-06-15",
                     "2017-07-15", "2017-08-15", "2017-09-15", "2017-10-15", "2017-11-15", "2017-12-15"))
DOYs.all <- rep(DOYs.model, times=length(yrs.model)) + rep(365*(yrs.model-1), each=length(DOYs.model))

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
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.shp"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW) %>% 
  rename(greenhouse=greenhs,
         growsize=growsiz,
         GrowNum=GrowNum,
         elev_m=elev_m,
         dtb_m=dtb_m,
         wte_m=wte_m,
         lon=lon,
         lat=lat,
         screenTopDepth_m=scrnTD_,
         screenBotDepth_m=scrnBD_)

## load pumping rates - this is proprietary from TNC, cannot be shared (until Wilson et al paper published)
df.pump <- read.csv(file.path(dir.TNC, "CannabisMonthlyWaterUse_WilsonEtAl.csv"), stringsAsFactors=F)
df.pump$Month <- factor(df.pump$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
df.pump$Setting <- factor(df.pump$Setting, levels=c("Outdoor", "Greenhouse"))
df.pump$MonthNum <- match(df.pump$Month, month.abb)
df.pump$MonthLengthDays <- lubridate::days_in_month(df.pump$MonthNum)
df.pump$m3PlantDay <- df.pump$MeanWaterUse_GalPlantDay*gal.to.m3

# set up long data frames for intermittent_pumping script; separate for outdoor and greenhouse
t.max.yrs <- max(yrs.model)
df.pump.outdoor <- 
  subset(df.pump, Setting=="Outdoor") %>% 
  dplyr::select(MonthNum, MonthLengthDays, m3PlantDay) %>% 
  replicate(t.max.yrs, ., simplify = FALSE) %>% 
  dplyr::bind_rows() %>% 
  transform(EndOfMonthDays = cumsum(MonthLengthDays))
df.pump.outdoor$StartOfMonthDays <- c(1, df.pump.outdoor$EndOfMonthDays[1:(t.max.yrs*12)-1]+1)

df.pump.greenhouse <- 
  subset(df.pump, Setting=="Greenhouse") %>% 
  dplyr::select(MonthNum, MonthLengthDays, m3PlantDay) %>% 
  replicate(t.max.yrs, ., simplify = FALSE) %>% 
  dplyr::bind_rows() %>% 
  transform(EndOfMonthDays = cumsum(MonthLengthDays))
df.pump.greenhouse$StartOfMonthDays <- c(1, df.pump.greenhouse$EndOfMonthDays[1:(t.max.yrs*12)-1]+1)

## subdivide sf.streams into points for web squared calculation
# define point spacing and figure out how many points to make
pt.spacing <- 25  # [m]
shp.streams <- as(sf.streams, "Spatial")
n.pts <- round(gLength(shp.streams)/pt.spacing)
set.seed(1)
shp.streams.pts <- spsample(shp.streams, n=n.pts, type="regular")
df.streams.pts <- as.data.frame(shp.streams.pts)
colnames(df.streams.pts) <- c("lon", "lat")
shp.streams.buffer <- buffer(shp.streams, 0.1, dissolve=F) # figure out what SegNum each point corresponds to
int <- intersect(shp.streams.pts, shp.streams.buffer)
df.streams.pts <- cbind(df.streams.pts, int@data)
sf.streams.pts <- sf::st_as_sf(df.streams.pts, coords=c("lon", "lat"), crs=st_crs(sf.streams))

## unique well-seg combos
df.combos <- 
  df.out %>% 
  left_join(sf.grows[,c("GrowNum", "greenhouse", "plants", "growsize")], by="GrowNum") %>% 
  dplyr::select(SegNum, GrowNum, dist_wellToStream_m, S_bulk, Tr_bulk_m2d, lmda_m2d, greenhouse, growsize, plants) %>% 
  unique()

## calculate depletion for all combos
start.flag.Qs <- T
for (i in 1:dim(df.combos)[1]){
  # identify well-seg combo
  seg <- df.combos$SegNum[i]
  grow <- df.combos$GrowNum[i]
  
  # get times
  output_t_days <- df.out$time_days[df.out$SegNum==seg & df.out$GrowNum==grow]
  output_frac <- df.out$frac_depletion[df.out$SegNum==seg & df.out$GrowNum==grow]
  
  # different calculation depending on outdoor vs indoor plants
  if (df.combos$greenhouse[i] == 1){
    Qs <- intermittent_pumping(t = output_t_days,
                               starts = df.pump.greenhouse$StartOfMonthDays,
                               stops  = df.pump.greenhouse$EndOfMonthDays,
                               rates  = df.pump.greenhouse$m3PlantDay*df.combos$plants[i],
                               method = "hunt",
                               d = df.combos$dist_wellToStream_m[i],
                               S = df.combos$S_bulk[i],
                               Tr = df.combos$Tr_bulk_m2d[i],
                               lmda = df.combos$lmda_m2d[i])
  } else {
    Qs <- intermittent_pumping(t = output_t_days,
                               starts = df.pump.outdoor$StartOfMonthDays,
                               stops  = df.pump.outdoor$EndOfMonthDays,
                               rates  = df.pump.outdoor$m3PlantDay*df.combos$plants[i],
                               method = "hunt",
                               d = df.combos$dist_wellToStream_m[i],
                               S = df.combos$S_bulk[i],
                               Tr = df.combos$Tr_bulk_m2d[i],
                               lmda = df.combos$lmda_m2d[i])
  }
  
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
  
  print(paste0("Depletion ", i, " of ", dim(df.combos)[1], " complete, ", Sys.time()))
  
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
