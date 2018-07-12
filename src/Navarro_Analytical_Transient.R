## Navarro_Analytical_Transient.R
#' This script is intended to calculate streamflow depletion from a
#' variety of analytical solutions for differnet pumping wells and 
#' distribute it via depletion apportionment for a bunch of stream reaches.
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

## what depletion apportionment output do you want?
apportionment_name <- "_LocalArea"      # output from Navarro_DepletionApportionment_LocalArea.R
#apportionment_name <- "_AdjacentOnly"   # output from Navarro_DepletionApportionment_AdjacentOnly.R
#apportionment_name <- "_MaskDryStreams" # output from Navarro_DepletionApportionment_MaskDryStreams.R

## Make sure these are the same as your MODFLOW script!

## choose stream boundary condition and modflow version - 
## this is just used to get the well numbers that were run,
## so as long as SFR and RIV have the same wells pumped it
## doesn't matter which you use
stream_BC <- "RIV"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"
timeType  <- "Transient" # "SteadyState" or "Transient"

## load MODFLOW depletion output
df.MODFLOW <- 
  file.path("modflow", "HTC", "Navarro", timeType, stream_BC, modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d)

## what timesteps are you interested in?
ts.all <- unique(df.MODFLOW$Time)
ts.pump.start <- sum(days_in_month(seq(1,4))) + 1

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

## load depletion apportionment output
df.apportion <- 
  paste0("Navarro_DepletionApportionment", apportionment_name, "_AllMethods+Wells+Reaches.csv") %>% 
  file.path("results", .) %>% 
  read.csv(stringsAsFactors=F)

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

# figure out which segments are part of Navarro
segs.navarro <- shp.streams@data$SegNum[shp.streams@data$TerminalPa==outlet.TerminalPa]
shp.streams@data$width_m <- WidthFromDA(DA=shp.streams@data$TotDASqKM, w.min=1, w.max=100)
df.apportion <- left_join(df.apportion, shp.streams@data[,c("SegNum", "width_m")], by="SegNum")

# analytical calculations: continuous pumping -----------------------------------------------------

# get rid of any segments with < 0.0001 depletion apportionment, or
# WellNum that you don't have MODFLOW results for
f.thres <- 0.0001  # =0.01%
df.apportion <- 
  subset(df.apportion, 
         (f.InvDist > f.thres | f.InvDistSq > f.thres |
            f.Web > f.thres | f.WebSq > f.thres | f.TPoly > f.thres) &
           WellNum %in% unique(df.MODFLOW$WellNum) & 
           SegNum %in% segs.navarro)

# loop through timesteps
start.flag <- T
for (ts in ts.all){
  for (analytical in c("glover", "hunt")){
    # add column for timestep and method
    df.apportion$Time <- ts
    df.apportion$analytical <- analytical
    
    # calculate aquifer vertical thickness for transmissivity
    screen_length <- 50  # [m] - should be same as script MODFLOW_Navarro-SteadyState.py
    df.apportion$thickness_m <- abs(df.apportion$wte_m-df.apportion$streambed_elev_m)    # reeves et al- uses vertical distance between top of well screen and streambed 
    df.apportion$thickness_m[df.apportion$thickness_m < screen_length] <- screen_length  # if vertical distance is < screen length, use screen length
    
    # calculate depletion fraction for each individual segment
    if (analytical=="glover"){
      df.apportion$Qf <- glover(t=(ts-ts.pump.start), d=df.apportion$distToWell.min.m, S=sy, 
                                Tr=(hk*df.apportion$thickness_m))
    } else if (analytical=="hunt"){
      df.apportion$Qf <- 
        streambed_conductance(w=df.apportion$width_m, Kriv=vk, briv=riverbed_thickness) %>% 
        hunt(t=(ts-ts.pump.start), d=df.apportion$distToWell.min.m, S=sy, Tr=(hk*df.apportion$thickness_m), lmda=.)
    }
    
    # weight segments using depletion apportionment fractions
    df.apportion$Qf.InvDist <- df.apportion$Qf*df.apportion$f.InvDist
    df.apportion$Qf.InvDistSq <- df.apportion$Qf*df.apportion$f.InvDistSq
    df.apportion$Qf.Web <- df.apportion$Qf*df.apportion$f.Web
    df.apportion$Qf.WebSq <- df.apportion$Qf*df.apportion$f.WebSq
    df.apportion$Qf.TPoly <- df.apportion$Qf*df.apportion$f.TPoly
    
    # add to overall data frame
    if (start.flag){
      df.out <- subset(df.apportion, Qf.InvDist > f.thres | Qf.InvDistSq > f.thres |
                         Qf.Web > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres)
      start.flag <- F
    } else {
      df.out <- rbind(df.out, 
                      subset(df.apportion, Qf.InvDist > f.thres | Qf.InvDistSq > f.thres |
                               Qf.Web > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres))
    }
    
    # status update
    print(paste0(analytical, " ", ts, " complete"))
    
  } # end of analytical loop
} # end of ts loop

## save output csv
df.out %>% 
  dplyr::select(SegNum, WellNum, Time, analytical, 
                Qf.InvDist, Qf.InvDistSq, Qf.Web, Qf.WebSq, Qf.TPoly) %>% 
  write.csv(., file.path("results", paste0("Depletion_Analytical_", timeType, apportionment_name, "_AllMethods+Wells+Reaches.csv")), row.names=F)

## plot cumulative depletion for each well at each timestep
df.sum <- 
  df.out %>% 
  group_by(Time, WellNum, analytical) %>% 
  summarize(Qf.InvDist = sum(Qf.InvDist),
            Qf.InvDistSq = sum(Qf.InvDistSq),
            Qf.Web = sum(Qf.Web),
            Qf.WebSq = sum(Qf.WebSq),
            Qf.TPoly = sum(Qf.TPoly)) %>% 
  melt(id=c("Time", "WellNum", "analytical"))

ggplot(df.sum, aes(x=Time, y=value, group=WellNum)) +
  geom_line() +
  facet_grid(variable ~ analytical)
