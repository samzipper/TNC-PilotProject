## Navarro_DepletionApportionment+Analytical_Transient.R
#' This script is intended to calculate depletion and apportion it
#' among different stream reaches for a bunch of wells. All of this is
#' in one script, instead of several, because the steps are interlinked:
#'  (1) Load MODFLOW results and figure out timesteps and wells for comparison 
#'      (MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
#'  (2) For each well and timestep, figure out the max distance at which depletion occurs
#'  (3) Calculate depletion apportionment using only stream reaches within that distance
#'  (4) Calculate analytical depletion for each reach and distribute using 
#'      apportionment output
#' 
#' The well locations are created using the script MODFLOW_Navarro_InputPrepData.R
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

#### (0) Define some parameters
## Make sure these are the same as your MODFLOW script!

## choose stream boundary condition and modflow version
stream_BC <- "RIV"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"

## pumping start time in MODFLOW model?
ts.pump.start <- sum(days_in_month(seq(1,4))) + 1

## depletion threshold to retain stream reaches?
Qf.thres <- 0.01   # =0.1%

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
screen_length <- 50    # well screen length [m]

## streambed parameters
depth <- 5  # river depth?
riverbed_K <- hk/10
riverbed_thickness <- 1

#### (0) Prep spatial data

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

# make a spatial points data frame
xy <- df.wel[,c("lon", "lat")]
spdf.wel <- SpatialPointsDataFrame(coords = xy, data = df.wel,
                                   proj4string = CRS(crs.MODFLOW))


## load stream data - created in MODFLOW_Navarro_InputPrepData.R
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# combine into 1 line per stream feature (which is defined by TerminalPa)
shp.streams.dissolve <- gLineMerge(shp.streams)

## convert stream lines to points
# define point spacing and figure out how many points to make
pt.spacing <- 10  # [m]
n.pts <- round(gLength(shp.streams)/pt.spacing)

# sample points
shp.streams.pts <- spsample(shp.streams, n=n.pts, type="regular")
df.streams.pts <- as.data.frame(shp.streams.pts)
colnames(df.streams.pts) <- c("lon", "lat")

# figure out what SegNum each point corresponds to
shp.streams.buffer <- buffer(shp.streams, 0.1, dissolve=F)
int <- intersect(shp.streams.pts, shp.streams.buffer)
df.streams.pts <- cbind(df.streams.pts, int@data)

## stream elevation and width needed for Reeves approximation of Hunt lambda
shp.streams@data$width_m <- WidthFromDA(DA=shp.streams@data$TotDASqKM, w.min=1, w.max=100)  # empirical relationship for Navarro
df.stream.elev <- 
  read.table(file.path("modflow", "input", "isfr_ReachData.txt"), stringsAsFactors=F, header=T) %>% 
  group_by(SegNum, SFR_NSEG) %>% 
  summarize(totalStreamLength_m = sum(length_m),
            streambed_elev_m = median(elev_m_min)-depth) %>% 
  left_join(., shp.streams@data[,c("SegNum", "width_m")], by="SegNum")

#### (1) Load MODFLOW results and figure out timesteps and wells for comparison 
####     (MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)

## load MODFLOW results
df.MODFLOW.RIV <- 
  file.path("modflow","HTC", "Navarro", "Transient", "RIV", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,
            stream_BC = "RIV")

df.MODFLOW.SFR <- 
  file.path("modflow","HTC", "Navarro", "Transient", "SFR", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,
            stream_BC = "SFR")

# combine into single data frame
df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)

# convert Time to time since pumping started
df.MODFLOW$Time <- df.MODFLOW$Time - ts.pump.start
df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)

## load output from steady-state, no pumping scenario - this is used to define the screen interval
## so that screen interval is consistent between MODFLOW and analytical
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wte.csv")), header=F)

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$wte_m <- m.wte[as.matrix(df.wel[,c("row", "col")])+1]

#### (2) For each well and timestep, figure out the max distance at which depletion occurs
# figure out wells and timesteps
wells.all <- sort(unique(df.MODFLOW$WellNum))
times.all <- sort(unique(df.MODFLOW$Time))

start.flag <- T
for (wel in wells.all){
  ## for a given well, calculate the distance to each point
  i.wel <- which(df.wel$WellNum==wel)
  
  # get distance to all stream points
  df.wel.dist <- data.frame(SegNum = df.streams.pts$SegNum,
                            distToWell.m = round(sqrt((df.streams.pts$lon-df.wel$lon[i.wel])^2 + (df.streams.pts$lat-df.wel$lat[i.wel])^2), 2))
  
  # grab the lat/lon for these points
  df.wel.dist$lon <- df.streams.pts$lon
  df.wel.dist$lat <- df.streams.pts$lat
  
  # join stream elevation
  df.wel.dist <- left_join(df.wel.dist, df.stream.elev[,c("SegNum", "streambed_elev_m", "totalStreamLength_m", "width_m")], by="SegNum")
  
  # estimate streambed conductance following Reeves et al (2009)
  df.wel.dist$lmda <- streambed_conductance(w = df.wel.dist$width_m,
                                            Kriv = hk*0.1,
                                            briv = riverbed_thickness)
  
  # calculate aquifer vertical thickness for transmissivity
  df.wel.dist$thickness_m <- abs(df.wel$wte_m[i.wel]-df.wel.dist$streambed_elev_m)    # reeves et al- uses vertical distance between top of well screen and streambed 
  df.wel.dist$thickness_m[df.wel.dist$thickness_m < screen_length] <- screen_length   # if vertical distance is < screen length, use screen length
  
  ## Thiessen Polygon apportionment method is not distance-based so it just has one result
  df.tpoly <- apportion.tpoly(reach=df.wel.dist$SegNum, 
                              dist=df.wel.dist$distToWell.m, 
                              lon=df.wel.dist$lon, 
                              lat=df.wel.dist$lat, 
                              wel.lon=df.wel$lon[i.wel],
                              wel.lat=df.wel$lat[i.wel],
                              wel.num=wel,
                              local.area.m=10000,
                              coord.ref=CRS(crs.MODFLOW),
                              col.names=c("SegNum", "f.TPoly"))
  
  for (t in times.all){
    for (analytical in c("glover", "hunt")){
      ## figure out the maximum distance at which depletion > Qf.thres would occur
      # for transmissivity calculation, use the thickness of the aquifer at the closest stream
      # (this is bounded to not be less than screen_length so typically that's what it is)
      closest.thickness <- df.wel.dist$thickness_m[min(which.min(df.wel.dist$distToWell.m))]
      max.dist <- depletion_max_distance(Qf_thres = Qf.thres,
                                         d_max    = 10000,
                                         method   = analytical,
                                         t  = t,
                                         S  = sy,
                                         Tr = hk*closest.thickness,
                                         lmda = min(df.wel.dist$lmda))  # lmda will be ignored if analytical==glover
      
      #### (3) Calculate depletion apportionment using only stream reaches within that distance
      # subset to stream reaches within that distance
      df.wel.dist.t <- subset(df.wel.dist, distToWell.m <= max.dist)
      
      # if there are any streams within that distance, apportion using distance-based methods
      if (dim(df.wel.dist.t)[1] > 0){
        df.id <- 
          data.frame(reach = df.wel.dist.t$SegNum, 
                     dist = df.wel.dist.t$distToWell.m) %>% 
          apportion_inverse(reach_dist = ., w=1) %>% 
          set_colnames(c("SegNum", "f.InvDist"))
        
        df.idsq <- 
          data.frame(reach = df.wel.dist.t$SegNum, 
                     dist = df.wel.dist.t$distToWell.m) %>% 
          apportion_inverse(reach_dist = ., w=2) %>% 
          set_colnames(c("SegNum", "f.InvDistSq"))
        
        df.web <-
          data.frame(reach = df.wel.dist.t$SegNum, 
                     dist = df.wel.dist.t$distToWell.m) %>% 
          apportion_web(reach_dist = ., w=1) %>% 
          set_colnames(c("SegNum", "f.Web"))
        
        df.websq <-
          data.frame(reach = df.wel.dist.t$SegNum, 
                     dist = df.wel.dist.t$distToWell.m) %>% 
          apportion_web(reach_dist = ., w=2) %>% 
          set_colnames(c("SegNum", "f.WebSq"))
        
        # combine into single data frame
        df.apportion.wel.t <- 
          full_join(df.id, df.idsq, by="SegNum") %>% 
          full_join(x=., y=df.web, by="SegNum") %>% 
          full_join(x=., y=df.websq, by="SegNum") %>% 
          full_join(x=., y=df.tpoly, by="SegNum") %>% 
          transform(WellNum = wel,
                    Time = t)
        
        # any NAs means that no depletion is calculated; set to 0
        df.apportion.wel.t[is.na(df.apportion.wel.t)] <- 0
        
        # add column for minimum distance to well from anywhere on this reach
        df.apportion.wel.t <- 
          group_by(df.wel.dist, SegNum) %>% 
          summarize(distToWell.min.m = min(distToWell.m),
                    thickness.max.m = max(thickness_m),
                    lmda.min = min(lmda)) %>% 
          left_join(x=df.apportion.wel.t, y=., by=c("SegNum"))
        
        #### (4) Calculate analytical depletion and distribute using apportionment output
        if (analytical=="glover"){
          df.apportion.wel.t$Qf <- glover(t  = t, 
                                          d  = df.apportion.wel.t$distToWell.min.m, 
                                          S  = sy, 
                                          Tr = (hk*df.apportion.wel.t$thickness.max.m))
          df.apportion.wel.t$analytical <- analytical
          
        } else if (analytical=="hunt"){
          df.apportion.wel.t$Qf <- hunt(t  = t, 
                                        d  = df.apportion.wel.t$distToWell.min.m, 
                                        S  = sy, 
                                        Tr = (hk*df.apportion.wel.t$thickness.max.m),
                                        lmda = df.apportion.wel.t$lmda.min)
          df.apportion.wel.t$analytical <- analytical
        }
        
        # weight segments using depletion apportionment fractions
        df.apportion.wel.t$Qf.InvDist   <- df.apportion.wel.t$Qf*df.apportion.wel.t$f.InvDist
        df.apportion.wel.t$Qf.InvDistSq <- df.apportion.wel.t$Qf*df.apportion.wel.t$f.InvDistSq
        df.apportion.wel.t$Qf.Web       <- df.apportion.wel.t$Qf*df.apportion.wel.t$f.Web
        df.apportion.wel.t$Qf.WebSq     <- df.apportion.wel.t$Qf*df.apportion.wel.t$f.WebSq
        df.apportion.wel.t$Qf.TPoly     <- df.apportion.wel.t$Qf*df.apportion.wel.t$f.TPoly
        
        # add to overall data frame
        if (start.flag){
          df.out <- subset(df.apportion.wel.t, Qf.InvDist > f.thres | Qf.InvDistSq > f.thres |
                             Qf.Web > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres)
          start.flag <- F
        } else {
          df.out <- rbind(df.out, 
                          subset(df.apportion.wel.t, Qf.InvDist > f.thres | Qf.InvDistSq > f.thres |
                                   Qf.Web > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres))
        }  # end of start.flag check
      }  # end of dist check
    }  # end of analytical loop
    
    # status update
    print(paste0("well ", wel, ", time ", t, " complete, ", Sys.time()))
    
  }  # end of t loop
}  # end of wel loop

## save output
# round to 5 digits to reduce file size
df.out$Time <- round(df.out$Time, 1)
df.out$distToWell.min.m <- round(df.out$distToWell.min.m, 2)
df.out$f.InvDist <- round(df.out$f.InvDist, 5)
df.out$f.InvDistSq <- round(df.out$f.InvDistSq, 5)
df.out$f.Web <- round(df.out$f.Web, 5)
df.out$f.WebSq <- round(df.out$f.WebSq, 5)
df.out$f.TPoly <- round(df.out$f.TPoly, 5)
df.out$Qf.InvDist <- round(df.out$Qf.InvDist, 5)
df.out$Qf.InvDistSq <- round(df.out$Qf.InvDistSq, 5)
df.out$Qf.Web <- round(df.out$Qf.Web, 5)
df.out$Qf.WebSq <- round(df.out$Qf.WebSq, 5)
df.out$Qf.TPoly <- round(df.out$Qf.TPoly, 5)

# convert Time back to model time
df.out$Time <- df.out$Time + ts.pump.start

df.out %>% 
  dplyr::select(SegNum, WellNum, Time, analytical, 
                f.InvDist, f.InvDistSq, f.Web, f.WebSq, f.TPoly, 
                Qf.InvDist, Qf.InvDistSq, Qf.Web, Qf.WebSq, Qf.TPoly) %>% 
  write.csv(., file.path("results", "Depletion_Analytical_Dynamic_AllMethods+Wells+Reaches.csv"), row.names=F)
