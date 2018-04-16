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

source("src/paths+packages.R")

# Define some parameters --------------------------------------------------------
## Make sure these are the same as your MODFLOW script!

## choose stream boundary condition and modflow version
stream_BC <- "RIV"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"

## what timesteps are you interested in?
ts.all <- c(5,10,100)

## various model parameters
# units: [m] and [d]
# flow parameters
hk <- 1e-11*1e7*86400  # horizontal K [m/d], convert k [m-2] to K [m/s] to K [m/d]
ss <- 1e-5             # specific storage
sy <- 0.10             # specific yield (using 50% of domain mean porosity for now)
vka <- 1.              # anisotropy
vk <- hk/vka           # calculate vertical K [m/d] based on horizontal K and anisotropy

# define bottom elevation
zbot <- -100

# streambed parameters
depth <- 4  # river depth?
riverbed_K <- hk/10
river_width <- 10
riverbed_thickness <- 1

# Prep input data ---------------------------------------------------------

## load depletion apportionment output
df.apportion <- 
  read.csv(file.path("results","Navarro_DepletionApportionment_AllMethods+Wells+Reaches.csv"),
           stringsAsFactors=F)

## load well input data
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario
m.head <- as.matrix(read.table(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "head.txt")))

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$head_m <- m.head[as.matrix(df.wel[,c("row", "col")])+1]

# add elevation to df.apportion
df.apportion <- left_join(df.apportion, df.wel[,c("WellNum", "head_m", "ztop_m")], by="WellNum")

## stream elevation needed for Reeves approximation of Hunt lambda
df.apportion <- 
  read.table(file.path("modflow", "input", "isfr_ReachData.txt"), stringsAsFactors=F, header=T) %>% 
  group_by(SegNum, SFR_NSEG) %>% 
  summarize(totalStreamLength_m = sum(length_m),
            streambed_elev_m = median(elev_m_min)-depth) %>% 
  left_join(df.apportion, ., by=c("SegNum"))


## load analytical depletion apportionment equations
# this is a different repository (StreamflowDepletionModels) so source them from GitHub
Hunt1999 <- 
  source_github("https://raw.githubusercontent.com/szipper/StreamflowDepletionModels/master/src/Hunt1999.R")
Hunt1999_lmda <- 
  source_github("https://raw.githubusercontent.com/szipper/StreamflowDepletionModels/master/src/Hunt1999_lmda.R")
GloverBalmer1954 <- 
  source_github("https://raw.githubusercontent.com/szipper/StreamflowDepletionModels/master/src/GloverBalmer1954.R")

# analytical calculations: continuous pumping -----------------------------------------------------

start.flag <- T
for (ts in ts.all){
  for (analytical in c("Glover", "Hunt")){
    # add column for timestep and method
    df.apportion$ts <- ts
    df.apportion$analytical <- analytical
    
    # calculate depletion fraction for each individual segment
    if (analytical=="Glover"){
      df.apportion$Qf <- GloverBalmer1954(d=df.apportion$distToWell.min.m, S=sy, Tr=(hk*(df.apportion$head_m-zbot)), t=ts)
    } else if (analytical=="Hunt"){
      df.apportion$Qf <- 
        Hunt1999_lmda(Kv=vk, w=river_width, b=abs(df.apportion$ztop_m - df.apportion$streambed_elev_m)) %>% 
        Hunt1999(d=df.apportion$distToWell.min.m, S=sy, Tr=(hk*(df.apportion$head_m-zbot)), t=ts, lmda=.)
    }
    
    # weight segments using depletion apportionment fractions
    df.apportion$Qf.InvDist <- df.apportion$Qf*df.apportion$f.InvDist
    df.apportion$Qf.InvDistSq <- df.apportion$Qf*df.apportion$f.InvDistSq
    df.apportion$Qf.Web <- df.apportion$Qf*df.apportion$f.Web
    df.apportion$Qf.WebSq <- df.apportion$Qf*df.apportion$f.WebSq
    df.apportion$Qf.TPoly <- df.apportion$Qf*df.apportion$f.TPoly
    
    # add to overall data frame
    if (start.flag){
      df.out <- df.apportion
      start.flag <- F
    } else {
      df.out <- rbind(df.out, df.apportion)
    }
    
  } # end of analytical loop
} # end of ts loop