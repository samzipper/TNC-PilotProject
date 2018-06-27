## Navarro_Analytical_Transient_OneWell.R
#' This script is intended to calculate streamflow depletion from a
#' variety of analytical solutions for one pumping well and 
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

# Define some parameters --------------------------------------------------------
## Make sure these are the same as your MODFLOW script!

## choose stream boundary condition and modflow version
stream_BC <- "RIV"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"

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

# what wells do you want to plot?
wells.plot <- c(151)

# what stream segments are adjacent?
segs.plot <- c(71, 72, 73, seq(116, 123), 180, 181, 272, 274, 276, 277)

# MODFLOW data processing -------------------------------------------------

df.MODFLOW.RIV <- 
  file.path("modflow","HTC", "Navarro", "Transient", "RIV", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,
            stream_BC = "RIV") %>% 
  subset(WellNum %in% wells.plot)

ts.pump.start <- sum(days_in_month(seq(1,4))) + 1
df.MODFLOW.RIV$Time <- df.MODFLOW.RIV$Time-ts.pump.start

# calculate total capture fraction
df.MODFLOW.sum <- 
  df.MODFLOW.RIV %>% 
  group_by(Time, WellNum) %>% 
  summarize(depletion.prc.max = max(depletion.prc.modflow),
            depletion.prc.sum = sum(depletion.prc.modflow))

# Calculate depletion apportionment ---------------------------------------

## load well locations
df.wel <- 
  read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T) %>% 
  subset(WellNum %in% wells.plot)

# make a spatial points data frame
xy <- df.wel[,c("lon", "lat")]
spdf.wel <- SpatialPointsDataFrame(coords = xy, data = df.wel,
                                   proj4string = CRS(crs.MODFLOW))


## load stream data - created in MODFLOW_Navarro_InputPrepData.R
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# subset to only adjacent streams
shp.streams <- subset(shp.streams, SegNum %in% segs.plot)

# data frame for ggplots
df.streams <- tidy(shp.streams, id=SegNum)

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

# get distance to all stream points
df.wel.dist <- data.frame(SegNum = df.streams.pts$SegNum,
                          distToWell.m = round(sqrt((df.streams.pts$lon-df.wel$lon)^2 + (df.streams.pts$lat-df.wel$lat)^2), 2))

# grab the lat/lon for these points
df.wel.dist$lon <- df.streams.pts$lon
df.wel.dist$lat <- df.streams.pts$lat

local.area.m <- 4616.931  # for 100 m resolution
local.area.m <- local.area.m*5

# calculate depletion apportionment fractions for different methods
df.id <- apportion.inv.dist(reach=df.wel.dist$SegNum, 
                            dist=df.wel.dist$distToWell.m, 
                            w=1, col.names=c("SegNum", "f.InvDist"))

df.idsq <- apportion.inv.dist(reach=df.wel.dist$SegNum, 
                              dist=df.wel.dist$distToWell.m, 
                              w=2, col.names=c("SegNum", "f.InvDistSq"))

df.web <- apportion.web.dist(reach=df.wel.dist$SegNum, 
                             dist=df.wel.dist$distToWell.m, 
                             w=1, col.names=c("SegNum", "f.Web"))

df.websq <- apportion.web.dist(reach=df.wel.dist$SegNum, 
                               dist=df.wel.dist$distToWell.m, 
                               w=2, col.names=c("SegNum", "f.WebSq"))

df.tpoly <- apportion.tpoly(reach=df.wel.dist$SegNum, 
                            dist=df.wel.dist$distToWell.m, 
                            lon=df.wel.dist$lon, 
                            lat=df.wel.dist$lat, 
                            wel.lon=df.wel$lon,
                            wel.lat=df.wel$lat,
                            wel.num=wells.plot,
                            local.area.m=local.area.m,
                            coord.ref=CRS(crs.MODFLOW),
                            col.names=c("SegNum", "f.TPoly"))

# combine into single data frame
df.apportion <- 
  full_join(df.id, df.idsq, by="SegNum") %>% 
  full_join(x=., y=df.web, by="SegNum") %>% 
  full_join(x=., y=df.websq, by="SegNum") %>% 
  full_join(x=., y=df.tpoly, by="SegNum")
df.apportion$WellNum <- wells.plot

# add column for minimum distance to well from anywhere on this reach
df.apportion <- 
  group_by(df.wel.dist, SegNum) %>% 
  summarize(distToWell.min.m = min(distToWell.m)) %>% 
  left_join(x=df.apportion, y=., by=c("SegNum"))

# Prep input data ---------------------------------------------------------

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

# figure out which segments are part of Navarro
shp.streams@data$width_m <- WidthFromDA(DA=shp.streams@data$TotDASqKM, w.min=1, w.max=100)
df.apportion <- left_join(df.apportion, shp.streams@data[,c("SegNum", "width_m")], by="SegNum")

## load analytical depletion apportionment equations
# this is a different repository (StreamflowDepletionModels) so source them from GitHub
#devtools::install_github("szipper/streamDepletr")
require(streamDepletr)

# analytical calculations: continuous pumping -----------------------------------------------------

# what timesteps do you want to count?
ts.all <- unique(c(1:9 %o% 10^(0:8), unique(df.MODFLOW.sum$Time)))

# extract WellNum you want to plot
df.apportion.wel <- 
  subset(df.apportion, WellNum %in% wells.plot)

# loop through timesteps
start.flag <- T
for (ts in ts.all){
  # add column for timestep and method
  df.apportion.wel$Time <- ts
  
  # calculate aquifer vertical thickness for transmissivity
  screen_length <- 50  # [m] - should be same as script MODFLOW_Navarro-SteadyState.py
  df.apportion.wel$thickness_m <- abs(df.apportion.wel$wte_m-df.apportion.wel$streambed_elev_m)    # reeves et al- uses vertical distance between top of well screen and streambed 
  df.apportion.wel$thickness_m[df.apportion.wel$thickness_m < screen_length] <- screen_length  # if vertical distance is < screen length, use screen length
  
  # riverbed thickness - same as script MODFLOW_Navarro-SteadyState.py
  riverbed_thickness <- 1
  
  # calculate depletion fraction for each individual segment
  df.apportion.wel$Qf <- glover(t=ts, d=df.apportion.wel$distToWell.min.m, S=sy, 
                                Tr=(hk*df.apportion.wel$thickness_m))
  
  # weight segments using depletion apportionment fractions
  df.apportion.wel$Qf.InvDist <- df.apportion.wel$Qf*df.apportion.wel$f.InvDist
  df.apportion.wel$Qf.InvDistSq <- df.apportion.wel$Qf*df.apportion.wel$f.InvDistSq
  df.apportion.wel$Qf.Web <- df.apportion.wel$Qf*df.apportion.wel$f.Web
  df.apportion.wel$Qf.WebSq <- df.apportion.wel$Qf*df.apportion.wel$f.WebSq
  df.apportion.wel$Qf.TPoly <- df.apportion.wel$Qf*df.apportion.wel$f.TPoly
  
  # add to overall data frame
  if (start.flag){
    df.out <- df.apportion.wel
    start.flag <- F
  } else {
    df.out <- rbind(df.out, 
                    df.apportion.wel)
  }
  
  # status update
  print(paste0(ts, " complete"))
  
} # end of ts loop

## calculate cumulative depletion for each well at each timestep
df.sum <- 
  df.out %>% 
  group_by(Time, WellNum) %>% 
  summarize(Qf.max = max(Qf, na.rm=T),
            Qf.InvDist = sum(Qf.InvDist, na.rm=T),
            Qf.InvDistSq = sum(Qf.InvDistSq, na.rm=T),
            Qf.Web = sum(Qf.Web, na.rm=T),
            Qf.WebSq = sum(Qf.WebSq, na.rm=T),
            Qf.TPoly = sum(Qf.TPoly, na.rm=T)) %>% 
  melt(id=c("Time", "WellNum"))

# Compare MODFLOW and analytical at segment level -------------------------

df.out.join <- 
  full_join(df.out, df.MODFLOW.RIV, by=c("SegNum", "Time", "WellNum")) %>% 
  subset(Time %in% unique(df.MODFLOW.RIV$Time)) %>% 
  dplyr::select(SegNum, distToWell.min.m, Time, Qf.InvDist, Qf.InvDistSq, Qf.Web, Qf.WebSq, Qf.TPoly, depletion.prc.modflow) %>% 
  melt(id=c("SegNum", "distToWell.min.m", "Time", "depletion.prc.modflow"),
       value.name="depletion.prc.analytical", variable.name="method")

df.out.join$depletion.prc.analytical[is.na(df.out.join$depletion.prc.analytical)] <- 0
df.out.join$depletion.prc.modflow[is.na(df.out.join$depletion.prc.modflow)] <- 0

# calculate fit
df.fit <- 
  df.out.join %>% 
  group_by(Time, method) %>% 
  summarize(KGE.overall = KGE(depletion.prc.analytical, depletion.prc.modflow),
            RMSE.overall = rmse(depletion.prc.analytical, depletion.prc.modflow))

# make plots --------------------------------------------------------------

# depletion fraction
p.TimeToCapture <- 
  ggplot() +
  geom_hline(yintercept=c(0, 1), color=col.gray) +
  geom_line(data=df.sum, aes(x=Time, y=value, color=variable)) +
  geom_point(data=df.MODFLOW.sum, aes(x=Time, y=depletion.prc.sum)) +
  geom_point(data=df.MODFLOW.sum, aes(x=Time, y=depletion.prc.max), color="blue", shape=21) +
  scale_x_log10(name="Time [days]", limits=c(1, max(ts.all)), expand=c(0,0)) +
  scale_y_continuous(name="Capture Fraction", expand=c(0,0)) +
  scale_color_manual(name="Approach", values=c("Qf.max"="black", pal.method.Qf),
                     labels=c("Qf.max"="Glover Depletion\nfor Closest Segment\n", labs.method)) +
  labs(title=paste0("Well ", wells.plot)) +
  theme(panel.grid=element_line(color=col.gray),
        legend.position=c(0.01,0.99),
        legend.justification=c(0,1),
        legend.background=element_blank()) +
  ggsave(file.path("results", "Navarro_Analytical_Transient_OneWell_AdjacentOnly_TimeToCapture.png"))

# fit
p.KGE <- 
  ggplot() +
  geom_hline(yintercept=c(1), color=col.gray) +
  geom_line(data=df.fit, aes(x=Time, y=KGE.overall, color=method)) +
  scale_x_continuous(name="Time [days]", limits=c(1, max(df.fit$Time)), expand=c(0,0)) +
  scale_y_continuous(name="KGE") +
  scale_color_manual(name="Approach", values=c(pal.method.Qf), labels=labs.method, guide=F) +
  labs(title=paste0("Well ", wells.plot)) +
  ggsave(file.path("results", "Navarro_Analytical_Transient_OneWell_AdjacentOnly_KGE.png"))

p.RMSE <- 
  ggplot() +
  geom_line(data=df.fit, aes(x=Time, y=RMSE.overall, color=method)) +
  scale_x_continuous(name="Time [days]", limits=c(1, max(df.fit$Time)), expand=c(0,0)) +
  scale_y_continuous(name="RMSE [pp]") +
  scale_color_manual(name="Approach", values=c(pal.method.Qf), labels=labs.method, guide=F) +
  labs(title=paste0("Well ", wells.plot)) +
  ggsave(file.path("results", "Navarro_Analytical_Transient_OneWell_AdjacentOnly_RMSE.png"))

ggsave(file.path("results", "Navarro_Analytical_Transient_OneWell_AdjacentOnly.png"),
       grid.arrange(p.TimeToCapture, p.RMSE, ncol=2),
       width=8, height=6)
