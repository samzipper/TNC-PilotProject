## Navarro_CompareMODFLOWtoDepletionApportionment_IndividualWells.R
#' This script is intended to compare estimated transient streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
#' and geometric depletion apportionment methods (from script Navarro_Analytical_Transient.R)
#' for a couple selected wells.

source(file.path("src", "paths+packages.R"))

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## choose modflow version
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"
masked <- F           # use masked depletion apportionment? should only be T for SFR
stream_BC <- "RIV"    # stream boundary condition to use for setting steady-state head (screen interval)
timeType  <- "Intermittent" # "Transient" or "Intermittent"

# which methods to analyze?
methods.plot <- c("Qf.InvDistSq", "Qf.WebSq", "Qf.TPoly")

## what depletion apportionment output do you want?
#apportionment_name <- "_LocalArea"      # output from Navarro_DepletionApportionment_LocalArea.R run through Navarro_Analytical_Transient.R
#apportionment_name <- "_AdjacentOnly"   # output from Navarro_DepletionApportionment_AdjacentOnly.R run through Navarro_Analytical_Transient.R
apportionment_name <- "_WholeDomain"    # output from Navarro_DepletionApportionment_WholeDomain.R run through Navarro_Analytical_Transient.R
#apportionment_name <- "_Dynamic"        # output from Navarro_DepletionApportionment+Analytical_Transient.R

#### (0) Prep spatial data

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario - this is used to define the screen interval
## so that screen interval is consistent between MODFLOW and analytical
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wte.csv")), header=F)

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$wte_m <- m.wte[as.matrix(df.wel[,c("row", "col")])+1]

# make a spatial points data frame
xy <- df.wel[,c("lon", "lat")]
spdf.wel <- SpatialPointsDataFrame(coords = xy, data = df.wel,
                                   proj4string = CRS(crs.MODFLOW))

## load stream data - created in MODFLOW_Navarro_InputPrepData.R
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# domain boundary shapefile
shp <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro")
shp.UTM <- spTransform(shp, crs.MODFLOW)

shp.adj <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

## prep polygon boundaries for plots
df.basin <- tidy(shp.UTM)
df.basin.adj <- tidy(shp.adj.UTM)
df.riv <- tidy(shp.streams)

#### (1) Load depletion apportionment results results and combine with MODFLOW
####     (Navarro_DepletionApportionment+Analytical_Transient.R)

# remember: this has been trimmed to only stream segments in Navarro,
# and any segments with depletion <= 0.0001 (0.01%) have been removed
df.analytical <- 
  paste0("Depletion_Analytical_", timeType, apportionment_name, "_AllMethods+Wells+Reaches.csv") %>% 
  file.path("results", .) %>% 
  read.csv(stringsAsFactors=F) %>% 
  dplyr::select(c("SegNum", "WellNum", "analytical", "Time", methods.plot))

#### (2) Choose which wells to investigate - goal is to have wells spanning a range of 
####     depletion values (high to low) within the 10 year timespan.
df.analytical.sum <- 
  df.analytical %>% 
  subset(analytical=="glover") %>% 
  melt(id=c("SegNum", "WellNum", "analytical", "Time"),
       value.name = "depletion.prc", variable.name="method") %>% 
  group_by(WellNum, method, Time) %>% 
  summarize(Qf.sum = sum(depletion.prc))

# map of wells with depletion at final timestep
df.wel <-
  df.wel %>% 
  subset(WellNum %in% df.analytical.sum$WellNum) %>% 
  left_join(subset(df.analytical.sum, Time==max(df.analytical.sum$Time) & method=="Qf.InvDistSq")[,c("WellNum", "Qf.sum")],
            by="WellNum") %>%
  replace_na(list("Qf.sum"=0))

wels.plot <- c(365, 393, 421)

ggplot() +
  geom_point(data=df.wel, aes(x=lon, y=lat, size=Qf.sum)) +
  geom_point(data=subset(df.wel, WellNum %in% wels.plot), aes(x=lon, y=lat, size=Qf.sum), color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue")

ggplot(df.analytical.sum, aes(x=Time, y=Qf.sum, group=WellNum)) +
  geom_line() +
  # choose a few wells to put on with red lines
  geom_line(data=subset(df.analytical.sum, WellNum %in% wels.plot), aes(color=factor(WellNum)), size=2) +
  facet_wrap(~method)

#### (3) Load MODFLOW results and figure out timesteps and wells for comparison 
####     (MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
df.MODFLOW.RIV <- 
  file.path("modflow","HTC", "Navarro", timeType, "RIV", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
  transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
            stream_BC = "RIV")

df.MODFLOW.SFR <- 
  file.path("modflow","HTC", "Navarro", timeType, "SFR", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
  transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
            stream_BC = "SFR")

df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)

# calculate sum
df.MODFLOW.sum <-
  df.MODFLOW %>% 
  subset(stream_BC=="RIV") %>% 
  group_by(WellNum, Time) %>% 
  summarize(Qf.sum = sum(depletion.prc.modflow))

## combine MODFLOW with analytical
## Time has long decimals; round before merging to ensure time match
df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
df.analytical$Time <- round(df.analytical$Time, 1)

# remove results with depletion <= f.thres (should be same as Navarro_Analytical_Transient.R)
f.thres <- 0.0001  # =0.01%
df <- 
  full_join(subset(df.analytical, WellNum %in% wels.plot), 
            subset(df.MODFLOW, WellNum %in% wels.plot)[,c("stream_BC", "SegNum", "WellNum", "Time", "depletion.prc.modflow")], 
            by=c("SegNum", "WellNum", "Time")) %>% 
  melt(id=c("stream_BC", "SegNum", "WellNum", "Time", "analytical", "depletion.prc.modflow"),
       value.name="depletion.prc", variable.name="method") %>%  
  subset(depletion.prc.modflow > f.thres | depletion.prc > f.thres) %>% 
  # for now: just look at glover and RIV
  subset(analytical=="glover" & stream_BC=="RIV")

# missing values should be 0 (0s were filtered out in previous scripts)
df$depletion.prc.modflow[is.na(df$depletion.prc.modflow)] <- 0
df$depletion.prc[is.na(df$depletion.prc)] <- 0

# for each well: extract only the three most affected reaches at the end of the time
for (wel in wels.plot){
  df.w <- 
    df %>% 
    subset(WellNum==wel & Time==max(df$Time)) %>% 
    subset(method=="Qf.TPoly") %>% # because you are using MODFLOW for selection, choose an arbitrary apportionment method
    arrange(desc(depletion.prc.modflow))
  
  if (wel==wels.plot[1]){
    segs.keep <- df.w$SegNum[1:3]
    wel.keep <- rep(wel, 3)
  } else {
    segs.keep <- c(segs.keep, df.w$SegNum[1:3])
    wel.keep <- c(wel.keep, rep(wel, 3))
  }
}

df.plot <- 
  df %>% 
  subset((WellNum==wels.plot[1] &
           SegNum %in% segs.keep[wel.keep==wels.plot[1]]) |
           (WellNum==wels.plot[2] &
              SegNum %in% segs.keep[wel.keep==wels.plot[2]]) |
           (WellNum==wels.plot[3] &
              SegNum %in% segs.keep[wel.keep==wels.plot[3]]))

#### (4) Plot comparison through time
ggplot() +
  geom_line(data=df.plot, aes(x=Time, y=depletion.prc.modflow, group=SegNum, color=factor(SegNum))) +
  geom_line(data=df.plot, aes(x=Time, y=depletion.prc, group=SegNum, color=factor(SegNum)), linetype="dashed") +
  facet_grid(WellNum ~ method, scales="free",
             labeller=as_labeller(c(labels.method.Qf, "365"="Well 365", "393"="Well 393", "421"="Well 421"))) +
  scale_x_continuous(name="Time [days]", expand=c(0,0)) +
  scale_y_continuous(name="Capture Fraction") +
  scale_color_discrete(name="Stream Reach")

#### (5) save plot
p.map <- 
  ggplot() +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_point(data=df.wel, aes(x=lon, y=lat, size=Qf.sum), color=col.gray) +
  geom_point(data=subset(df.wel, WellNum %in% wels.plot), aes(x=lon, y=lat, size=Qf.sum, color=factor(WellNum))) +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  scale_size_continuous(name="Capture Fraction\nat t=3650 days") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_color_discrete(name="Well") +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5)) +
  NULL

p.sum <-
  ggplot() +
    geom_line(data=subset(df.MODFLOW.sum, WellNum %in% wels.plot), 
              aes(x=Time, y=Qf.sum, group=WellNum, color=factor(WellNum)), linetype="dashed") +
  geom_line(data=subset(df.analytical.sum, WellNum %in% wels.plot), 
            aes(x=Time, y=Qf.sum, group=WellNum, color=factor(WellNum))) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Total Capture Fraction") +
  scale_color_discrete(name="Well", guide=F) +
  facet_wrap(~method, labeller=as_labeller(labels.method.Qf)) +
  NULL

p.depletion <-
  ggplot() +
  geom_line(data=df.plot, aes(x=Time, y=depletion.prc.modflow, group=SegNum, color=factor(SegNum)), linetype="dashed") +
  geom_line(data=df.plot, aes(x=Time, y=depletion.prc, group=SegNum, color=factor(SegNum))) +
  facet_grid(WellNum ~ method, scales="free",
             labeller=as_labeller(c(labels.method.Qf, "365"="Well 365", "393"="Well 393", "421"="Well 421"))) +
  scale_x_continuous(name="Time [days]", expand=c(0,0)) +
  scale_y_continuous(name="Capture Fraction") +
  scale_color_brewer(name="Stream Reach", type="qual",
                     palette=6) +
  theme(legend.position="bottom") +
  NULL

ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_IndividualWells_", timeType, ".png")),
       grid.arrange(p.map, 
                    p.sum, 
                    p.depletion, 
                    ncol=1, heights=c(1, 1, 1.5)),
       width=190, height=240, units="mm")
