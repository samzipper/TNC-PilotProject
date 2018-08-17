## Navarro_DecideBestMethod_06_ErrorAnalysisBestMethod.R
#' This script is intended to conduct an error analysis of the 
#' best method identified in the previous scripts.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

#### (0) Prep various inputs, parameters, etc
## what is the best method?
web.exp.best <- 1.75   # winner from Navarro_DecideBestMethod_05_CompareModels.R
Qf.thres.best <- 0.01  # winner from Navarro_DecideBestMethod_05_CompareModels.R
BC.plot <- "RIV"

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

# which MODFLOW version? 
modflow_v <- "mfnwt"  # (no other option)

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

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

## Prep spatial data

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario - this is used to define the screen interval
## so that screen interval is consistent between MODFLOW and analytical
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", "RIV", modflow_v, "wte.csv")), header=F)

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

# (1) Load results -------------------------------------------------

start.flag <- T
for (timeType in c("Transient", "Intermittent")) {
  
  ## open MODFLOW results
  df.MODFLOW.RIV <-
    file.path("modflow","HTC", "Navarro", timeType, "RIV", modflow_v, "Depletion_MODFLOW.csv") %>%
    read.csv(stringsAsFactors=F) %>%
    #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
    transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
              stream_BC = "RIV",
              stringsAsFactors=F)
  
  df.MODFLOW.SFR <-
    file.path("modflow","HTC", "Navarro", timeType, "SFR", modflow_v, "Depletion_MODFLOW.csv") %>%
    read.csv(stringsAsFactors=F) %>%
    #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
    transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
              stream_BC = "SFR",
              stringsAsFactors=F)
  
  df.MODFLOW <-
    rbind(df.MODFLOW.RIV, df.MODFLOW.SFR) %>%
    subset(depletion.prc.modflow > f.thres)
  df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
  
  df.MODFLOW.max <-
    df.MODFLOW %>% 
    group_by(stream_BC, WellNum, Time) %>%
    filter(depletion.prc.modflow==max(depletion.prc.modflow))
  
  ## Load output from Navarro_DecideBestMethod_04_TweakParameters.R
  df.analytical <- 
    paste0("Navarro_DecideBestMethod_04_TweakParameters-", timeType, ".csv") %>% 
    file.path("results", .) %>% 
    read.csv(stringsAsFactors=F) %>% 
    subset(Qf.thres==Qf.thres.best & web.exp==web.exp.best) %>% 
    dplyr::select(c("SegNum", "WellNum", "Time", "Qf")) %>%
    set_colnames(c("SegNum", "WellNum", "Time", "depletion.prc")) %>% 
    subset(depletion.prc > f.thres)
  df.analytical$Time <- round(df.analytical$Time, 1)
  
  df.analytical.max <-
    df.analytical %>%
    group_by(WellNum, Time) %>%
    filter(depletion.prc==max(depletion.prc))
  
  # combine
  df <-
    df.MODFLOW %>%
    subset(stream_BC == BC.plot) %>%
    left_join(df.analytical, by=c("WellNum", "Time", "SegNum")) %>%
    replace_na(list("stream_BC"=BC.plot, 
                    "depletion.prc"=0, 
                    "depletion.prc.modflow" = 0)) %>%
    transform(pump = timeType,
              stringsAsFactors=F) %>% 
    dplyr::select(WellNum, Time, SegNum, depletion.prc.modflow, depletion.prc, pump)
  
  df.max <-
    left_join(df.MODFLOW.max[,c("SegNum", "WellNum", "Time", "stream_BC", "depletion.prc.modflow")], 
              df.analytical.max[,c("SegNum", "WellNum", "Time")],
              by=c("WellNum", "Time"), suffix=c(".modflow", ".analytical")) %>% 
    subset(stream_BC == BC.plot) %>% 
    # add depletion in the most affected modflow segment
    left_join(df.analytical, by=c("WellNum", "Time", "SegNum.modflow"="SegNum")) %>% 
    replace_na(list("SegNum.analytical"=9999, "depletion.prc"=0)) %>%
    transform(pump = timeType,
              stringsAsFactors=F)
  
  if (start.flag) {
    df.all <- df
    df.max.all <- df.max
    start.flag <- F
  } else {
    df.all <- rbind(df.all, df)
    df.max.all <- rbind(df.max.all, df.max)
  }
  
  # status update
  print(paste(timeType, "complete"))
  
}  # end of timeType loop

## calculate match percent
df.fit.match <-
  df.max.all %>% 
  group_by(Time, stream_BC, pump) %>% 
  summarize(n.reach = sum(is.finite(SegNum.modflow)),
            n.match = sum(SegNum.modflow==SegNum.analytical),
            n.noMatch.noAnalytical = sum(SegNum.analytical==9999),
            MSE.match = MSE(depletion.prc, depletion.prc.modflow),
            KGE.match = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
  transform(prc.match = n.match/n.reach,
            prc.noAnalytical = n.noMatch.noAnalytical/n.reach)

df.fit.match.gt5 <-
  df.max.all %>% 
  subset(depletion.prc.modflow > 0.05) %>% 
  group_by(Time, stream_BC, pump) %>% 
  summarize(n.reach = sum(is.finite(SegNum.modflow)),
            n.match = sum(SegNum.modflow==SegNum.analytical),
            n.noMatch.noAnalytical = sum(SegNum.analytical==9999),
            MSE.match = MSE(depletion.prc, depletion.prc.modflow),
            KGE.match = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
  transform(prc.match = n.match/n.reach,
            prc.noAnalytical = n.noMatch.noAnalytical/n.reach)

## calculate overall fit through time
df.fit.all <-
  df.all %>% 
  group_by(pump, Time) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))

## calculate overall fit by well
df.fit.well <-
  df.all %>% 
  group_by(pump, WellNum) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>% 
  left_join(df.wel, by="WellNum")

## calculate match percent by well
df.match.well <-
  df.max.all %>% 
  group_by(WellNum, stream_BC, pump) %>% 
  summarize(n.reach = sum(is.finite(SegNum.modflow)),
            n.noMatch.noAnalytical = sum(SegNum.analytical==9999),
            n.match = sum(SegNum.modflow==SegNum.analytical)) %>%
  transform(prc.match = n.match/n.reach,
            prc.noAnalytical = n.noMatch.noAnalytical/n.reach)
  
## fit by distance to well
# load horizontal distance
df.dist <-
  file.path("results", "Navarro_DepletionApportionment_WholeDomain_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  left_join(df.wel[,c("WellNum", "ztop_m")], by=c("WellNum"))

# load streambed info to calculate vertical distance
df.stream.elev <- 
  read.table(file.path("modflow", "input", "isfr_ReachData.txt"), stringsAsFactors=F, header=T) %>% 
  group_by(SegNum) %>% 
  summarize(totalStreamLength_m = sum(length_m),
            streambed_elev_m = median(elev_m_min)-depth)

df.dist <- 
  left_join(df.dist, df.stream.elev, by=c("SegNum")) %>% 
  transform(dist.vert.m = (ztop_m - streambed_elev_m))

# join with all info
df.all.dist <-
  df.all %>% 
  subset(pump=="Transient") %>% 
  left_join(df.dist[,c("SegNum", "WellNum", "distToWell.min.m", "dist.vert.m")], by=c("SegNum", "WellNum"))

dist.breaks <- quantile(df.all.dist$distToWell.min.m, seq(0,1,0.1))
df.all.dist$dist.cut <- cut(df.all.dist$distToWell.min.m, dist.breaks, include.lowest=T)
df.fit.dist <- 
  df.all.dist %>% 
  group_by(pump, dist.cut) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"),
            dist.m.median = median(distToWell.min.m))

dist.vert.breaks <- quantile(df.all.dist$dist.vert.m, seq(0,1,0.1))
df.all.dist$dist.vert.cut <- cut(df.all.dist$dist.vert.m, dist.vert.breaks, include.lowest=T)
df.fit.dist.vert <- 
  df.all.dist %>% 
  group_by(pump, dist.vert.cut) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"),
            dist.m.median = median(dist.vert.m))

## calculate overall fit by reach length
df.fit.length <-
  df.stream.elev %>% 
  transform(length.cut = cut(totalStreamLength_m, 
                             quantile(totalStreamLength_m, seq(0,1,0.1)), 
                             include.lowest=T)) %>%  # cut into 10 percentile chunks
  #transform(length.cut = cut(totalStreamLength_m, c(0, 0.1, 0.5, 1, 5, max(totalStreamLength_m))*1000)) %>%  # cut at same length intervals as Nanaimo paper
  left_join(df.all, ., by=c("SegNum")) %>% 
  group_by(pump, length.cut) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"),
            length.m.median = median(totalStreamLength_m))

## calculate fit by time of year
df.fit.mo <- 
  df.all %>% 
  transform(year = floor(Time/365)+1, 
            month = month(ymd("2012-12-31") + days(floor(Time) %% 365))) %>% 
  group_by(pump, month) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))

df.fit.yr.mo <- 
  df.all %>% 
  transform(year = floor(Time/365)+1, 
            month = month(ymd("2012-12-31") + days(floor(Time) %% 365))) %>% 
  group_by(pump, year, month) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))

# Make some plots ---------------------------------------------------------

## correct identification of most affected segment through time
ggplot() +
  geom_line(data=df.fit.match, aes(x=Time, y=prc.match, color=pump)) +
  geom_line(data=df.fit.match.gt5, aes(x=Time, y=prc.match, color=pump), linetype="longdash") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="% of Wells where Most-Affected\nReach is Correctly Identified", 
                     limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.25)) +
  scale_color_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Match-Prc.png"),
         width=4, height=4, units="in")

## comparison of depletion at max affected segment
ggplot(df.max.all, aes(x=depletion.prc, y=depletion.prc.modflow, color=Time)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_wrap(~pump) +
  scale_x_continuous(name="Analytical Depletion", 
                     limits=c(0,1), expand=c(0,0), labels = scales::percent) +
  scale_y_continuous(name="MODFLOW Depletion in Most Affected Reach", 
                     limits=c(0,1), expand=c(0,0), labels = scales::percent) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Scatter-AnalyticalVsMODFLOW.png"),
         width=6, height=4, units="in")

## fit through time
ggplot(df.fit.all, aes(x=Time, y=KGE.overall, color=pump)) +
  geom_line() +
  scale_x_continuous() +
  scale_y_continuous(name="Overall KGE") +
  scale_color_discrete(name=NULL) +
  theme(legend.position=c(0.99, 0.01),
        legend.justification=c(1,0))

## error distribution through time
# ggtern currently broken with ggplot 3.0.0... try using this Ternary package instead
png(file.path("results", "Navarro_DecideBestMethod_06_Ternary-Time.png"),
    height=6, width=6, units="in", res=150)
TernaryPlot(alab='% MSE due to Bias', blab='% MSE due to Variability', clab='% MSE due to Correlation')
TernaryPoints(df.fit.all[df.fit.all$pump=="Transient", c("MSE.bias.norm", "MSE.var.norm", "MSE.cor.norm")])
dev.off()

df.fit.all[df.fit.all$pump=="Transient", c("MSE.bias.norm", "MSE.var.norm", "MSE.cor.norm")] %>% head()

## error by well
# fit for different wells
df.fit.well %>% 
  ggplot() + 
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_point(aes(x=lon, y=lat, size=KGE.overall)) +
  facet_wrap(~pump, ncol=2) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5)) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-Map-KGE.png"),
         width=8, height=8, units="in")

df.fit.well %>% 
  ggplot() + 
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_point(aes(x=lon, y=lat, size=MSE.overall)) +
  facet_wrap(~pump, ncol=2) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5)) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-Map-MSE.png"),
         width=8, height=8, units="in")

df.match.well %>% 
  left_join(df.wel, by="WellNum") %>% 
  subset(pump=="Transient") %>% 
  ggplot() + 
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_point(aes(x=lon, y=lat, size=1-prc.match)) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_size_continuous(name="% of Time Most Affected\nReach is Identified WRONG") +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5)) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-Map-Match.png"),
         width=6, height=4, units="in")

# error scatterplots
df.fit.well %>% 
  ggplot(aes(x=ztop_m, y=KGE.overall, color=pump)) +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Well Land Surface Elevation [m]") +
  scale_y_continuous(name="Well KGE") +
  scale_color_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-KGEvElev.png"),
         width=4, height=4, units="in")

df.fit.well %>% 
  ggplot(aes(x=(ztop_m-wte_m), y=KGE.overall, color=pump)) +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Well Water Table Depth [m]") +
  scale_y_continuous(name="Well KGE") +
  scale_color_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-KGEvWTD.png"),
         width=4, height=4, units="in")

df.match.well %>% 
  left_join(df.wel, by="WellNum") %>% 
  ggplot(aes(x=(ztop_m-wte_m), y=1-prc.match, color=pump)) +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Well Water Table Depth [m]") +
  scale_y_continuous(name="% of Time Most Affected\nReach is Identified WRONG") +
  scale_color_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-MatchvWTD.png"),
         width=4, height=4, units="in")

df.match.well %>% 
  left_join(df.wel, by="WellNum") %>% 
  ggplot(aes(x=ztop_m, y=1-prc.match, color=pump)) +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Well Elevation [m]") +
  scale_y_continuous(name="% of Time Most Affected\nReach is Identified WRONG") +
  scale_color_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-MatchvElev.png"),
         width=4, height=4, units="in")

df.dist %>% 
  group_by(WellNum) %>% 
  filter(distToWell.min.m==min(distToWell.min.m)) %>% 
  left_join(df.match.well, ., by="WellNum") %>%
  ggplot(aes(x=distToWell.min.m, y=1-prc.match, color=pump)) +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Distance to Closest Reach") +
  scale_y_continuous(name="% of Time Most Affected\nReach is Identified WRONG") +
  scale_color_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Well-MatchvDist.png"),
         width=4, height=4, units="in")

## error by horizontal distance
df.fit.dist %>% 
  ggplot(aes(x=factor(round(dist.m.median, 1)), y=KGE.overall)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity") +
  scale_x_discrete(name="Quantile Median Horizontal Distance from Well to Stream") +
  scale_y_continuous(name="KGE") +
  coord_cartesian(ylim=c(-1,1)) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Dist-KGEbyQuantile.png"),
         width=6, height=4, units="in")

df.fit.dist %>% 
  ggplot(aes(x=factor(round(dist.m.median, 1)), y=MSE.overall)) +
  geom_bar(stat="identity") +
  scale_x_discrete(name="Quantile Median Horizontal Distance from Well to Stream") +
  scale_y_continuous(name="MSE") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Dist-MSEbyQuantile.png"),
         width=6, height=4, units="in")

## error by vertical distance
df.fit.dist.vert %>% 
  ggplot(aes(x=factor(round(dist.m.median, 1)), y=KGE.overall)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity") +
  scale_x_discrete(name="Quantile Median Vertical Distance from Well to Stream") +
  scale_y_continuous(name="KGE") +
  coord_cartesian(ylim=c(-1,1)) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_DistVert-KGEbyQuantile.png"),
         width=6, height=4, units="in")

df.fit.dist.vert %>% 
  ggplot(aes(x=factor(round(dist.m.median, 1)), y=MSE.overall)) +
  geom_bar(stat="identity") +
  scale_x_discrete(name="Quantile Median Vertical Distance from Well to Stream") +
  scale_y_continuous(name="MSE") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_DistVert-MSEbyQuantile.png"),
         width=6, height=4, units="in")

## error by reach length
df.fit.length %>% 
  ggplot(aes(x=length.cut, y=KGE.overall, fill=pump)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="Stream Reach Length [m]") +
  scale_y_continuous(name="KGE") +
  scale_fill_discrete(name=NULL) +
  coord_cartesian(ylim=c(-1,1)) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Length-KGE.png"),
         width=6, height=4, units="in")

df.fit.length %>% 
  ggplot(aes(x=length.cut, y=MSE.overall, fill=pump)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="Stream Reach Length [m]") +
  scale_y_continuous(name="MSE") +
  scale_fill_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Length-MSE.png"),
         width=6, height=4, units="in")

## error by time of year
df.fit.mo %>% 
  ggplot(aes(x=month, y=KGE.overall, color=pump)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() + geom_point() +
  scale_x_continuous(name="Month", breaks=seq(1,12)) +
  scale_y_continuous(name="KGE") +
  scale_fill_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Time-Mo-KGE.png"),
         width=6, height=4, units="in")

df.fit.mo %>% 
  ggplot(aes(x=month, y=MSE.overall, color=pump)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() + geom_point() +
  scale_x_continuous(name="Month", breaks=seq(1,12)) +
  scale_y_continuous(name="MSE") +
  scale_fill_discrete(name=NULL) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", "Navarro_DecideBestMethod_06_Time-Mo-MSE.png"),
         width=6, height=4, units="in")
