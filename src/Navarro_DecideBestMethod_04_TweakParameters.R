## Navarro_DecideBestMethod_04_TweakParameters.R
#' This script is intended to tweak the parameters of the best depletion apportionment (web) +
#' analytical (hunt) + search radius (adjacent + dynamic). The parameters we will tweak are:
#'   -Exponent for web weighting (scripts _01-_03 use 1 or 2)
#'   -Percent depletion threshold for dynamic stream selection (scripts _01-_03 use 1%)
#' The saved results will be analytical model output, which will be compared to MODFLOW in a separate script.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

#### (0) Define some parameters
## define bounds of sensitivity analysis
web.exp.all  <- seq(1, 3, 0.25)
Qf.thres.all <- c(0.0001, 0.001, 0.01)

## choose stream boundary condition and modflow version used for SS head values to define screen interval
stream_BC <- "RIV"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"

## when is pumping occurring? 
# for continuous: just start date
ts.pump.start.transient <- sum(days_in_month(seq(1,4))) + 1

# for intermittent: need starts and stops
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
pt.spacing <- 5  # [m]
n.pts <- round(gLength(shp.streams)/pt.spacing)

# sample points
set.seed(1)
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

## load MODFLOW results to get timesteps and wells
for (timeType in c("Intermittent", "Transient")) {
  df.MODFLOW.RIV <- 
    file.path("modflow","HTC", "Navarro", timeType, "RIV", modflow_v, "Depletion_MODFLOW.csv") %>% 
    read.csv(stringsAsFactors=F) %>% 
    transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,
              stream_BC = "RIV",
              pump = timeType,
              stringsAsFactors=F)
  
  df.MODFLOW.SFR <- 
    file.path("modflow","HTC", "Navarro", "Intermittent", "SFR", modflow_v, "Depletion_MODFLOW.csv") %>% 
    read.csv(stringsAsFactors=F) %>% 
    transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,
              stream_BC = "SFR",
              pump = timeType,
              stringsAsFactors=F)
  
  # combine into single data frame
  if (timeType == "Intermittent") {
    df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)
  } else {
    df.MODFLOW <- rbind(df.MODFLOW, df.MODFLOW.RIV, df.MODFLOW.SFR)
  }
}

# convert Time to time since pumping started
df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)

## load output from steady-state, no pumping scenario - this is used to define the screen interval
## so that screen interval is consistent between MODFLOW and analytical
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wte.csv")), header=F)

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$wte_m <- m.wte[as.matrix(df.wel[,c("row", "col")])+1]

#### (2) For each well and timestep, figure out the max distance at which depletion occurs
# figure out wells and timesteps
wells.all <- sort(unique(df.MODFLOW$WellNum))
times.all.intermittent <- 
  df.MODFLOW$Time[df.MODFLOW$pump=="Intermittent"] %>% 
  unique() %>% 
  sort()
times.all.transient <- 
  df.MODFLOW$Time[df.MODFLOW$pump=="Transient"] %>% 
  unique() %>% 
  sort()

# determine time since pumping started
df.times <- 
  rbind(data.frame(time.since.pump.start = (times.all.transient - ts.pump.start.transient),
                   pump = "Transient"),
        data.frame(time.since.pump.start = (times.all.intermittent - min(ts.pump.starts)),
                   pump = "Intermittent"))
times.since.pump.all <- 
  df.times$time.since.pump.start %>% 
  unique() %>% 
  sort()

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
  rdll <- data.frame(reach = df.wel.dist$SegNum, 
                     dist = df.wel.dist$distToWell.m,
                     lon = df.wel.dist$lon,
                     lat = df.wel.dist$lat)
  
  df.tpoly <- apportion_polygon(reach_dist_lon_lat = rdll,
                                wel_lon  = df.wel$lon[i.wel],
                                wel_lat  = df.wel$lat[i.wel],
                                crs = CRS(crs.MODFLOW)) %>% 
    set_colnames(c("SegNum", "f.TPoly"))
  
  # loop through Qf.thres
  for (Qf.thres in Qf.thres.all) {
    # loop through times
    max.dist.prev <- 0  # initialize
    for (t in times.since.pump.all){
      ## figure out the maximum distance at which depletion > Qf.thres would occur
      # for transmissivity calculation, use the thickness of the aquifer at the closest stream
      # (this is bounded to not be less than screen_length so typically that's what it is)
      closest.thickness <- df.wel.dist$thickness_m[min(which.min(df.wel.dist$distToWell.m))]
      max.dist <- depletion_max_distance(Qf_thres = Qf.thres,
                                         d_min    = max.dist.prev,
                                         d_max    = Inf,
                                         method   = "hunt",
                                         t  = t,
                                         S  = sy,
                                         Tr = hk*closest.thickness,
                                         lmda = min(df.wel.dist$lmda))  # lmda will be ignored if method="glover"
      max.dist.prev <- max.dist  # update to make search quicker next timestep
      
      if (Qf.thres==Qf.thres.all[1] & t==times.since.pump.all[1]){
        df.max.dist <- data.frame(time.since.pump = t,
                                  max.dist = max.dist,
                                  Qf.thres = Qf.thres)
      } else {
        df.max.dist <- rbind(df.max.dist, 
                             data.frame(time.since.pump = t,
                                        max.dist = max.dist,
                                        Qf.thres = Qf.thres))
      }
      
    }  # end of time loop
  }  # end of Qf.thres loop
  
  #### (3) Calculate depletion apportionment using only stream reaches within that distance
  ####     and adjacent catchments
  # get unique distance values
  max.dist.all <- sort(unique(df.max.dist$max.dist))
  
  for (dist in max.dist.all) {
    # get stream reaches within that distance and/or adjacent
    df.wel.dist.t <- subset(df.wel.dist, (distToWell.m <= dist) | SegNum %in% df.tpoly$SegNum)
    
    # apportion with different exponents
    for (web.exp in web.exp.all) {
      df.web.exp <-
        data.frame(reach = df.wel.dist.t$SegNum, 
                   dist  = df.wel.dist.t$distToWell.m) %>% 
        apportion_web(reach_dist = ., w=web.exp) %>% 
        set_colnames(c("SegNum", "f.Web")) %>% 
        transform(web.exp = web.exp,
                  max.dist = dist)
      
      if (web.exp==web.exp.all[1] & dist==max.dist.all[1]) {
        df.web <- df.web.exp
      } else {
        df.web <- rbind(df.web, df.web.exp)
      }
    }  # end of web.exp loop
  }  # end of dist loop
  
  # add column for minimum distance to well from anywhere on this reach
  df.web <- 
    group_by(df.wel.dist, SegNum) %>% 
    summarize(distToWell.min.m = min(distToWell.m),
              thickness.max.m = max(thickness_m),
              lmda.min = min(lmda)) %>% 
    left_join(x=df.web, y=., by=c("SegNum"))
  
  # combine data frames
  df.apportionment <- 
    full_join(df.max.dist, df.web, by=c("max.dist")) 
  
  ## check to make sure sum to 1
  #df.sum <- 
  #  df.apportionment %>% 
  #  group_by(time.since.pump, Qf.thres, web.exp) %>% 
  #  summarize(f.sum = sum(f.Web))
  #unique(df.sum$f.sum)
  
  ## now: for each unique combination of distToWell and lmda, run analytical model
  dist.thick.lmda.all <- unique(df.apportionment[,c("SegNum", "distToWell.min.m", "thickness.max.m", "lmda.min")])
  for (i in 1:dim(dist.thick.lmda.all)[1]){
    
    # get apportionment for this SegNum
    df.apportionment.seg <- subset(df.apportionment, 
                                   SegNum==dist.thick.lmda.all$SegNum[i] &
                                     distToWell.min.m==dist.thick.lmda.all$distToWell.min.m[i] & 
                                     thickness.max.m==dist.thick.lmda.all$thickness.max.m[i] &
                                     lmda.min==dist.thick.lmda.all$lmda.min[i])
    
    times.seg.intermittent <- unique(df.apportionment.seg$time.since.pump) + min(ts.pump.starts)
    times.seg.intermittent <- times.seg.intermittent[times.seg.intermittent %in% times.all.intermittent]
    
    times.seg.transient <- unique(df.apportionment.seg$time.since.pump) + ts.pump.start.transient
    times.seg.transient <- times.seg.transient[times.seg.transient %in% times.all.transient]
    
    ## calculate intermittent depletion
    inter <- F
    if (length(times.seg.intermittent) > 0) {
      inter <- T
      
      df.intermittent <- 
        data.frame(time.since.pump = (times.seg.intermittent - min(ts.pump.starts)),
                   SegNum = dist.thick.lmda.all$SegNum[i],
                   distToWell.min.m = dist.thick.lmda.all$distToWell.min.m[i],
                   thickness.max.m = dist.thick.lmda.all$thickness.max.m[i],
                   lmda.min = dist.thick.lmda.all$lmda.min[i],
                   Qf.analytical = (intermittent_pumping(t = times.seg.intermittent, 
                                                         starts = ts.pump.starts, 
                                                         stops = ts.pump.stops, 
                                                         rates = rep(Qw, length(times.seg.intermittent)),
                                                         method = "hunt", 
                                                         d = dist.thick.lmda.all$distToWell.min.m[i], 
                                                         S = sy,
                                                         Tr = (hk*dist.thick.lmda.all$thickness.max.m[i]),
                                                         lmda = dist.thick.lmda.all$lmda.min[i])/Qw),
                   pump = "Intermittent",
                   stringsAsFactors=F) %>% 
        right_join(subset(df.apportionment.seg, time.since.pump %in% (times.seg.intermittent - min(ts.pump.starts))),
                   by=c("time.since.pump", "SegNum", "distToWell.min.m", "thickness.max.m", "lmda.min")) %>% 
        unique()
      df.intermittent$Time <- df.intermittent$time.since.pump + min(ts.pump.starts)
    }
    
    ## calculate transient depletion
    trans <- F
    if (length(times.seg.transient) > 0) {
      trans <- T
      
      df.transient <- 
        data.frame(time.since.pump = (times.seg.transient - ts.pump.start.transient),
                   SegNum = dist.thick.lmda.all$SegNum[i],
                   distToWell.min.m = dist.thick.lmda.all$distToWell.min.m[i],
                   thickness.max.m = dist.thick.lmda.all$thickness.max.m[i],
                   lmda.min = dist.thick.lmda.all$lmda.min[i],
                   Qf.analytical = hunt(t  = times.seg.transient, 
                                        d = dist.thick.lmda.all$distToWell.min.m[i], 
                                        S  = sy, 
                                        Tr = (hk*dist.thick.lmda.all$thickness.max.m[i]),
                                        lmda = dist.thick.lmda.all$lmda.min[i]),
                   pump = "Transient",
                   stringsAsFactors=F) %>% 
        right_join(subset(df.apportionment.seg, time.since.pump %in% (times.seg.transient - ts.pump.start.transient)),
                   by=c("time.since.pump", "SegNum", "distToWell.min.m", "thickness.max.m", "lmda.min")) %>% 
        unique()
      df.transient$Time <- df.transient$time.since.pump + ts.pump.start.transient
    }
    
    if (inter & trans) {
      df.analytical.dist <- rbind(df.intermittent[,c("Time", "SegNum", "f.Web", "Qf.analytical", "pump", "web.exp", "Qf.thres")], 
                                  df.transient[,c("Time", "SegNum", "f.Web", "Qf.analytical", "pump", "web.exp", "Qf.thres")])
    } else if (inter) {
      df.analytical.dist <- df.intermittent[,c("Time", "SegNum", "f.Web", "Qf.analytical", "pump", "web.exp", "Qf.thres")]
    } else if (trans) {
      df.analytical.dist <- df.transient[,c("Time", "SegNum", "f.Web", "Qf.analytical", "pump", "web.exp", "Qf.thres")]
    } else {
      stop("No analytical results")
    }
    
    if (i==1){
      df.analytical <- df.analytical.dist
    } else {
      df.analytical <- rbind(df.analytical, 
                             df.analytical.dist)
    }
  }  # end of dist.lmda.all loop
  
  ## check to make sure everything sums to 1
  #df.sum <- 
  #  df.analytical %>% 
  #  group_by(Time, Qf.thres, web.exp, pump) %>% 
  #  summarize(f.sum = sum(f.Web))
  #unique(df.sum$f.sum)
  
  ## add WellNum
  df.analytical$WellNum <- wel
  
  ## calculate depletion fraction
  df.analytical$Qf <- df.analytical$Qf.analytical*df.analytical$f.Web
  
  if (wel == wells.all[1]) {
    df.all <- subset(df.analytical, Qf >= f.thres)
  } else {
    df.all <- rbind(df.all, subset(df.analytical, Qf >= f.thres))
  }
  
  # status update
  print(paste0("Well ", wel, " complete, ", Sys.time()))
  
}  # end of wel loop

## save output
# round to 5 digits to reduce file size
df.all$Time <- round(df.all$Time, 1)
df.all$WellNum <- round(df.all$WellNum, 0)
df.all$SegNum <- round(df.all$SegNum, 0)
df.all$web.exp <- round(df.all$web.exp, 2)
df.all$Qf.analytical <- round(df.all$Qf.analytical, 5)
df.all$Qf.thres <- signif(df.all$Qf.thres, 2)
df.all$f.Web <- round(df.all$f.Web, 5)
df.all$Qf <- round(df.all$Qf, 5)

# save
df.all %>% 
  subset(pump=="Intermittent") %>% 
  dplyr::select(Time, WellNum, SegNum, Qf.thres, web.exp, f.Web, Qf.analytical, Qf) %>% 
  write.csv(file.path("results", "Navarro_DecideBestMethod_04_TweakParameters-Intermittent.csv"), row.names=F)

df.all %>% 
  subset(pump=="Transient") %>% 
  dplyr::select(Time, WellNum, SegNum, Qf.thres, web.exp, f.Web, Qf.analytical, Qf) %>% 
  write.csv(file.path("results", "Navarro_DecideBestMethod_04_TweakParameters-Transient.csv"), row.names=F)
