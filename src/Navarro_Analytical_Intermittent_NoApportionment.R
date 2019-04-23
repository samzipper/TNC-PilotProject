## Navarro_Analytical_Intermittent_NoApportionment.R
#' This script is intended to calculate streamflow depletion from a
#' variety of analytical solutions for different pumping wells 
#' based only on the stream reach closest to each well (e.g. 
#' depletion apportionment equations are ignored).
#' 
#' The well locations are created using the script MODFLOW_Navarro_InputPrepData.R
#' The depletion apportionment calculations are created using the script Navarro_DepletionApportionment.R
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

# Define some parameters --------------------------------------------------------

## choose stream boundary condition and modflow version - 
## this is just used to get the well numbers that were run,
## so as long as SFR and RIV have the same wells pumped it
## doesn't matter which you use
stream_BC <- "RIV"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"

## load MODFLOW depletion output
df.MODFLOW <- 
  file.path("modflow", "HTC", "Navarro", "Intermittent", stream_BC, modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d)

## what timesteps are you interested in?
ts.all <- unique(df.MODFLOW$Time)

## when is pumping occurring? need starts and stops
# since we are using 365 day years (ignoring leap years)
# we just have to figure out June 1 and Nov 1 for first
# year and then add 365
ts.pump.starts <- seq(from=(sum(days_in_month(seq(1,5))) + 1),
                      by=365, 
                      length.out=10)
ts.pump.stops <- seq(from=(sum(days_in_month(seq(1,10))) + 1),
                     by=365, 
                     length.out=10)

## pumping rate - should be same as MODFLOW
Qw <- -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

## depletion threshold to not save results?
f.thres <- 0.0001  # =0.01%

## various model parameters
# units: [m] and [d]
# flow parameters
hk <- 1e-12*1e7*86400  # horizontal K [m/d], convert k [m-2] to K [m/s] to K [m/d]
ss <- 1e-5             # specific storage
sy <- 0.10             # specific yield (using 50% of domain mean porosity for now)
vka <- 10              # anisotropy
vk <- hk/vka           # calculate vertical K [m/d] based on horizontal K and anisotropy

## streambed parameters
depth <- 5  # river depth?
riverbed_K <- hk/10
riverbed_thickness <- 1

# Prep input data ---------------------------------------------------------

## load depletion apportionment output for AdjacentOnly - 
## this is used to figure out which reach is closest to each well
df.apportion <- 
  "Navarro_DepletionApportionment_AdjacentOnly_AllMethods+Wells+Reaches.csv" %>% 
  file.path("results", .) %>% 
  read.csv(stringsAsFactors=F) %>% 
  dplyr::select(WellNum, SegNum, distToWell.min.m) %>% 
  subset(WellNum %in% df.MODFLOW$WellNum) %>% 
  # get the closest stream to each well
  dplyr::group_by(WellNum) %>% 
  dplyr::filter(distToWell.min.m == min(distToWell.min.m))

## load well input data
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wte.csv")), header=F)

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$wte_m <- m.wte[as.matrix(df.wel[,c("row", "col")])+1]

# add elevation to df.apportion
df.apportion <- left_join(df.apportion, df.wel[,c("WellNum", "wte_m", "ztop_m")], by="WellNum")

## stream elevation needed for Reeves approximation of Hunt lambda
df.apportion <- 
  read.table(file.path("modflow", "input", "isfr_ReachData.txt"), stringsAsFactors=F, header=T) %>% 
  group_by(SegNum, SFR_NSEG) %>% 
  summarize(totalStreamLength_m = sum(length_m),
            streambed_elev_m = median(elev_m_min)-depth) %>% 
  left_join(df.apportion, ., by=c("SegNum")) 

## stream width needed for conductance
# load stream shapefile
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# calculate stream width
shp.streams@data$width_m <- WidthFromDA(DA=shp.streams@data$TotDASqKM, w.min=1, w.max=100)
df.apportion <- left_join(df.apportion, shp.streams@data[,c("SegNum", "width_m")], by="SegNum")

# calculate aquifer vertical thickness for transmissivity
screen_length <- 50  # [m] - should be same as script MODFLOW_Navarro-SteadyState.py
df.apportion$thickness_m <- abs(df.apportion$wte_m-df.apportion$streambed_elev_m)    # reeves et al- uses vertical distance between top of well screen and streambed 
df.apportion$thickness_m[df.apportion$thickness_m < screen_length] <- screen_length  # if vertical distance is < screen length, use screen length

# analytical calculations: intermittent pumping -----------------------------------------------------

start.flag <- T
for (i in 1:dim(df.apportion)[1]){
  for (analytical in c("glover", "hunt")){
    
    # calculate depletion fraction for each individual segment
    if (analytical=="glover"){
      Qs <- intermittent_pumping(t=ts.all, 
                                 starts=ts.pump.starts, 
                                 stops=ts.pump.stops, 
                                 rates=rep(Qw, length(ts.pump.starts)),
                                 method=analytical, 
                                 d=df.apportion$distToWell.min.m[i], 
                                 S=sy,
                                 Tr=(hk*df.apportion$thickness_m[i]))
    } else if (analytical=="hunt"){
      Qs <- 
        streambed_conductance(w=df.apportion$width_m[i], Kriv=vk, briv=riverbed_thickness) %>% 
        intermittent_pumping(t=ts.all, 
                             starts=ts.pump.starts, 
                             stops=ts.pump.stops, 
                             rates=rep(Qw, length(ts.pump.starts)),
                             method=analytical, 
                             d=df.apportion$distToWell.min.m[i], 
                             S=sy,
                             Tr=(hk*df.apportion$thickness_m[i]),
                             lmda=.)
    }
    
    # make data frame for this well-seg combo
    df.i <- data.frame(SegNum = df.apportion$SegNum[i],
                       WellNum = df.apportion$WellNum[i],
                       distToWell.min.m = df.apportion$distToWell.min.m[i],
                       Time = ts.all,
                       analytical = analytical,
                       Qf.NoApport = Qs/Qw)
    
    # add to overall data frame
    if (start.flag){
      df.out <- subset(df.i, abs(Qf.NoApport) > f.thres)
      start.flag <- F
    } else {
      df.out <- rbind(df.out, 
                      subset(df.i, abs(Qf.NoApport) > f.thres))
    }
    
    # status update
    print(paste0(analytical, " ", i, " complete"))
    
  } # end of analytical loop
} # end of i loop

## save output csv
df.out %>% 
  subset(Qf.NoApport > f.thres) %>% 
  dplyr::select(WellNum, SegNum, Time, analytical, 
                distToWell.min.m, Qf.NoApport) %>% 
  write.csv(., file.path("results", "Depletion_Analytical_Intermittent_NoApportionment_AllMethods+Wells+Reaches.csv"), row.names=F)

ggplot(df.out, aes(x=Time, y=Qf.NoApport, group=WellNum)) +
  geom_line() +
  facet_grid(. ~ analytical)
