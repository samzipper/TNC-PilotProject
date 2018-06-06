## Navarro_Analytical_Transient_DetermineSimulationLength.R
#' This script is intended to calculate streamflow depletion from a
#' variety of analytical solutions for differnet pumping wells and 
#' distribute it via depletion apportionment for a bunch of stream reaches.
#' 
#' The goal is to calculate how long to run transient MODFLOW simulations.
#' I will use the cutoff as the time at which 80% of the analytical solutions
#' have asymptoted (rounded up). Asymptote is defined as <1 pp increase in 
#' total capture fraction from one year to the next.
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

# Define some parameters --------------------------------------------------------

## what times to test [days]
ts.all <- c(1, 60, 180, seq(from=365, by=365*5, length.out=11))

## Make sure these are the same as your MODFLOW script!
## various model parameters
# units: [m] and [d]
# flow parameters
hk <- 1e-12*1e7*86400  # horizontal K [m/d], convert k [m-2] to K [m/s] to K [m/d]
ss <- 1e-5             # specific storage
sy <- 0.10             # specific yield (using 50% of domain mean porosity for now)
vka <- 10              # anisotropy
vk <- hk/vka           # calculate vertical K [m/d] based on horizontal K and anisotropy
screen_length <- 50    # [m] - should be same as script MODFLOW_Navarro-SteadyState.py

## streambed parameters - should be same as MODFLOW model
depth <- 5  # river depth?
riverbed_K <- hk/10
riverbed_thickness <- 1

# Prep input data ---------------------------------------------------------

## load depletion apportionment output
df.apportion <- 
  read.csv(file.path("results","Navarro_DepletionApportionment_AllMethods+Wells+Reaches.csv"),
           stringsAsFactors=F)

## load well input data
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario
# using SFR because it will likely be slower to asymptote (deeper WTD = less transmissivity)
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", "SFR", "mfnwt", "wte.csv")), header=F)

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

## load analytical depletion apportionment equations
require(streamDepletr)

# analytical calculations: continuous pumping -----------------------------------------------------

# get rid of any segments with < 0.0001 depletion apportionment, or
# WellNum that you don't have MODFLOW results for
f.thres <- 0.0001  # =0.01%
df.apportion <- 
  subset(df.apportion, 
         (f.InvDist > f.thres | f.InvDistSq > f.thres |
            f.Web > f.thres | f.WebSq > f.thres | f.TPoly > f.thres))

# loop through timesteps
start.flag <- T
for (ts in ts.all){
  
    # add column for timestep and method
    df.apportion$Time <- ts

    # calculate aquifer vertical thickness for transmissivity
    df.apportion$thickness_m <- abs(df.apportion$wte_m-df.apportion$streambed_elev_m)    # reeves et al- uses vertical distance between top of well screen and streambed 
    df.apportion$thickness_m[df.apportion$thickness_m < screen_length] <- screen_length  # if vertical distance is < screen length, use screen length
    
    # calculate depletion fraction for each individual segment
    # loop through analytical models - use Hunt only because it is slower depletion
    #df.apportion$Qf <- 
    #  streambed_conductance(w=df.apportion$width_m, Kriv=vk, briv=riverbed_thickness) %>% 
    #  hunt(t=ts, d=df.apportion$distToWell.min.m, S=sy, Tr=(hk*df.apportion$thickness_m), lmda=.)
    df.apportion$Qf <- 
      glover(t=ts, d=df.apportion$distToWell.min.m, S=sy, Tr=(hk*df.apportion$thickness_m))
    
    # weight segments using depletion apportionment fractions
    df.apportion$Qf.InvDist <- df.apportion$Qf*df.apportion$f.InvDist
    df.apportion$Qf.InvDistSq <- df.apportion$Qf*df.apportion$f.InvDistSq
    df.apportion$Qf.Web <- df.apportion$Qf*df.apportion$f.Web
    df.apportion$Qf.WebSq <- df.apportion$Qf*df.apportion$f.WebSq
    df.apportion$Qf.TPoly <- df.apportion$Qf*df.apportion$f.TPoly
    
    # add to overall data frame
    if (start.flag){
      #df.out <- subset(df.apportion, Qf.InvDist > f.thres | Qf.InvDistSq > f.thres |
      #                   Qf.Web > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres)
      df.out <- df.apportion
      start.flag <- F
    } else {
      df.out <- rbind(df.out, 
                      #subset(df.apportion, Qf.InvDist > f.thres | Qf.InvDistSq > f.thres |
                      #         Qf.Web > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres))
                      df.apportion)
    }
    
    # status update
    print(paste0(ts, " complete: ", Sys.time()))
    
} # end of ts loop

## plot cumulative depletion for each well at each timestep
df.sum <- 
  df.out %>% 
  group_by(Time, WellNum) %>% 
  # calculate cumulative capture fraction (within Navarro only)
  summarize(Qf.InvDist = sum(Qf.InvDist),
            Qf.InvDistSq = sum(Qf.InvDistSq),
            Qf.Web = sum(Qf.Web),
            Qf.WebSq = sum(Qf.WebSq),
            Qf.TPoly = sum(Qf.TPoly)) %>% 
  melt(id=c("Time", "WellNum"))

ggplot(df.sum, aes(x=Time/365, y=value, group=WellNum)) +
  geom_line(alpha=0.1) +
  facet_wrap(~variable, ncol=2) +
  scale_x_continuous(name="Years", expand=c(0,0)) +
  scale_y_continuous(name="Capture Fraction") +
  ggsave(file.path("results", "Navarro_Analytical_Transient_DetermineSimulationLength_CaptureFractionTs.png"),
         width=6, height=6, units="in")

df.sum %>% 
  subset(Time==max(ts.all)) %>% 
  ggplot(aes(x=value)) +
  geom_histogram(binwidth=0.05) +
  facet_wrap(~variable, ncol=2) +
  scale_x_continuous(name=paste0("Capture Fraction after ", max(ts.all)/365, " years"), expand=c(0,0)) +
  scale_y_continuous(name="Number of Wells") +
  ggsave(file.path("results", "Navarro_Analytical_Transient_DetermineSimulationLength_CaptureFractionHist.png"),
         width=6, height=6, units="in")

# calculate change from previous timestep
start.flag <- T
change.thres <- 0.01
for (ts in ts.all[-1]){
  df.sum.ts <- 
    df.sum %>% 
    subset(Time==ts) %>% 
    set_colnames(c("Time", "WellNum", "method", "Qf.current"))
  
  df.sum.ts_1 <- 
    df.sum %>% 
    subset(Time==(ts-365)) %>% 
    dplyr::select(WellNum, variable, value) %>% 
    set_colnames(c("WellNum", "method", "Qf.previous"))
  
  df.ts <-
    full_join(df.sum.ts, df.sum.ts_1, by=c("WellNum", "method")) %>% 
    transform(Qf.change = Qf.current - Qf.previous) %>% 
    transform(change.lt.thres = Qf.change < change.thres)
  df.ts$change.lt.thres[is.na(df.ts$change.lt.thres) & 
                          df.ts$Qf.current < change.thres] <- T
  df.ts$change.lt.thres[is.na(df.ts$change.lt.thres)] <- F
  
  if (start.flag){
    df.ts.out <- df.ts
    start.flag <- F
  } else {
    df.ts.out <- rbind(df.ts.out, df.ts)
  }
}

ggplot(df.ts.out, aes(x=Time/365, y=Qf.change, group=WellNum)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept=change.thres, color=col.gray) +
  facet_wrap(~method, ncol=2)

# calculate percent under change threshold for each timestep and method
df.ts.summary <- 
  df.ts.out %>% 
  group_by(Time, method) %>% 
  summarize(n.lt.thres = sum(change.lt.thres),
            n.total = sum(is.finite(Qf.current)),
            prc.lt.thres = n.lt.thres/n.total)

ggplot(df.ts.summary, aes(x=Time/365, y=prc.lt.thres, color=method)) +
  geom_hline(yintercept=0.8, color=col.gray) +
  geom_line()
