## Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW.R
#' This script is intended to compare depletion calculated using the script
#' Navarro_DepletionApportionment+Analytical_Transient.R to MODFLOW results.
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

# round Time to 1 digit
df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)

#### (2) Load depletion apportionment results results and combine with MODFLOW
####     (Navarro_DepletionApportionment+Analytical_Transient.R)

## load analytical results
df.out <- 
  file.path("results", "Depletion_Apportionment+Analytical_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors=F)

# combine and melt
df <- 
  df.out %>% 
  dplyr::select(SegNum, WellNum, Time, analytical, 
                Qf.InvDist, Qf.InvDistSq, Qf.Web, Qf.WebSq, Qf.TPoly) %>% 
  full_join(., df.MODFLOW[,c("stream_BC", "SegNum", "WellNum", "Time", "depletion.prc.modflow")], by=c("Time", "WellNum", "SegNum")) %>% 
  melt(id=c("stream_BC", "SegNum", "WellNum", "Time", "analytical", "depletion.prc.modflow"),
       value.name="depletion.prc", variable.name="method")

# missing values should be 0 (0s were filtered out in previous scripts)
df$depletion.prc.modflow[is.na(df$depletion.prc.modflow)] <- 0
df$depletion.prc[is.na(df$depletion.prc)] <- 0

tmp <- subset(df, is.na(analytical))
tmp$analytical <- "glover"
df <- rbind(df, tmp)
df$analytical[is.na(df$analytical)] <- "hunt"

tmp <- subset(df, is.na(stream_BC))
tmp$stream_BC <- "RIV"
df <- rbind(df, tmp)
df$stream_BC[is.na(df$stream_BC)] <- "SFR"

#### (3) Calculate derived statistics e.g. fit
# all points with any change
df.fit <- 
  df %>% 
  group_by(stream_BC, analytical, method, Time) %>% 
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

# only points with >5% depletion in either MODFLOW or analytical
df.nreach.big <- 
  df %>% 
  subset(depletion.prc > 0.05 | depletion.prc.modflow > 0.05) %>% 
  group_by(stream_BC, analytical, method, Time) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc))) %>% 
  subset(n.reach >= 3) %>% 
  transform(groups = paste(stream_BC, analytical, method, Time, sep="_"))

df.fit.big <-
  df %>% 
  transform(groups = paste(stream_BC, analytical, method, Time, sep="_")) %>% 
  subset(groups %in% df.nreach.big$groups) %>% 
  group_by(stream_BC, analytical, method, Time) %>% 
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

## for each timestep: what % of wells does analytical correctly identify the most affected reach?
df.most.affected.analytical <- 
  df.out %>% 
  group_by(WellNum, Time, analytical) %>% 
  summarize(Qf.InvDist.max = max(Qf.InvDist), 
            Qf.InvDist.max.seg = SegNum[which.max(Qf.InvDist)],
            Qf.InvDistSq.max = max(Qf.InvDistSq), 
            Qf.InvDistSq.max.seg = SegNum[which.max(Qf.InvDistSq)],
            Qf.Web.max = max(Qf.Web), 
            Qf.Web.max.seg = SegNum[which.max(Qf.Web)],
            Qf.WebSq.max = max(Qf.WebSq), 
            Qf.WebSq.max.seg = SegNum[which.max(Qf.WebSq)],
            Qf.TPoly.max = max(Qf.TPoly), 
            Qf.TPoly.max.seg = SegNum[which.max(Qf.TPoly)])

df.most.affected.MODFLOW <- 
  df.MODFLOW %>% 
  group_by(WellNum, Time, stream_BC) %>% 
  summarize(Qf.MODFLOW.max = max(depletion.prc.modflow), 
            Qf.MODFLOW.max.seg = SegNum[which.max(depletion.prc.modflow)])

# melt and combine
df.most.affected.prc <- 
  df.most.affected.analytical %>% 
  dplyr::select(WellNum, Time, analytical, Qf.InvDist.max, Qf.InvDistSq.max, 
                Qf.Web.max, Qf.WebSq.max, Qf.TPoly.max) %>% 
  melt(id=c("WellNum", "Time", "analytical"),
       variable.name="method", value.name="depletion.prc") %>% 
  full_join(., df.most.affected.MODFLOW, by=c("WellNum", "Time")) %>% 
  subset(!is.na(analytical) & !is.na(stream_BC))

df.most.affected.seg <- 
  df.most.affected.analytical %>% 
  dplyr::select(WellNum, Time, analytical, Qf.InvDist.max.seg, Qf.InvDistSq.max.seg, 
                Qf.Web.max.seg, Qf.WebSq.max.seg, Qf.TPoly.max.seg) %>% 
  melt(id=c("WellNum", "Time", "analytical"),
       variable.name="method", value.name="analytical.max.seg") %>% 
  full_join(., df.most.affected.MODFLOW, by=c("WellNum", "Time"))

# for each timestep, analytical model, stream_BC, and method, calculate the % of WellNum where max.seg is the same
df.most.affected.seg.match <-
  df.most.affected.seg %>% 
  subset(!is.na(analytical) & !is.na(stream_BC)) %>% 
  group_by(Time, analytical, stream_BC, method) %>% 
  summarize(n.well = sum(is.finite(WellNum)),
            n.match = sum(analytical.max.seg==Qf.MODFLOW.max.seg)) %>% 
  transform(match.prc = n.match/n.well)

## fit by well
df.fit.well <- 
  df %>% 
  group_by(stream_BC, analytical, method, WellNum) %>% 
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

## calculate sum of depletion.prc for each well/timestep/method
df.sum <-
  df %>% 
  group_by(stream_BC, Time, WellNum, analytical, method) %>% 
  summarize(Qf.MODFLOW = sum(depletion.prc.modflow),
            Qf.analytical = sum(depletion.prc))

df.fit.sum <- 
  df.sum %>% 
  group_by(stream_BC, analytical, method) %>% 
  summarize(n.reach = sum(is.finite(Qf.analytical)),
            cor = cor(Qf.analytical, Qf.MODFLOW, method="pearson"),
            bias = pbias(Qf.analytical, Qf.MODFLOW),
            R2 = R2(Qf.analytical, Qf.MODFLOW),
            MSE.bias = MSE.bias(Qf.analytical, Qf.MODFLOW),
            MSE.var = MSE.var(Qf.analytical, Qf.MODFLOW),
            MSE.cor = MSE.cor(Qf.analytical, Qf.MODFLOW),
            MSE.bias.norm = MSE.bias.norm(Qf.analytical, Qf.MODFLOW),
            MSE.var.norm = MSE.var.norm(Qf.analytical, Qf.MODFLOW),
            MSE.cor.norm = MSE.cor.norm(Qf.analytical, Qf.MODFLOW),
            MSE.overall = MSE(Qf.analytical, Qf.MODFLOW),
            KGE.overall = KGE(Qf.analytical, Qf.MODFLOW, method="2012"))

df.fit.sum.WellNum <- 
  df.sum %>% 
  group_by(stream_BC, analytical, method, WellNum) %>% 
  summarize(n.reach = sum(is.finite(Qf.analytical)),
            cor = cor(Qf.analytical, Qf.MODFLOW, method="pearson"),
            bias = pbias(Qf.analytical, Qf.MODFLOW),
            R2 = R2(Qf.analytical, Qf.MODFLOW),
            MSE.bias = MSE.bias(Qf.analytical, Qf.MODFLOW),
            MSE.var = MSE.var(Qf.analytical, Qf.MODFLOW),
            MSE.cor = MSE.cor(Qf.analytical, Qf.MODFLOW),
            MSE.bias.norm = MSE.bias.norm(Qf.analytical, Qf.MODFLOW),
            MSE.var.norm = MSE.var.norm(Qf.analytical, Qf.MODFLOW),
            MSE.cor.norm = MSE.cor.norm(Qf.analytical, Qf.MODFLOW),
            MSE.overall = MSE(Qf.analytical, Qf.MODFLOW),
            KGE.overall = KGE(Qf.analytical, Qf.MODFLOW, method="2012"))

#### (4) make plots
ts.all <- unique(df$Time)
ts.plot <- ts.all[c(6, 73, 365)]

# which methods to analyze?
methods.plot <- c("Qf.InvDistSq", "Qf.WebSq", "Qf.TPoly")

p.ts.KGE <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=method)) +
  geom_hline(yintercept=0, color=col.gray) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Kling-Gupta Efficiency") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.KGE.png")),
         width=8, height=8, units="in")

p.ts.match.prc <-
  df.most.affected.seg.match %>% 
  transform(method = gsub(".max.seg", "", method)) %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=match.prc, color=method)) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="% of Wells with Most-Depleted Segment Correct", 
                     labels = scales::percent, limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.match.prc.png")),
         width=8, height=8, units="in")

p.ts.MSE <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=MSE.overall, color=method)) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Mean Squared Error", labels=scales::percent) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") + 
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.MSE.png")),
         width=8, height=8, units="in")

p.ts.bias <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=bias, color=method)) +
  geom_hline(yintercept=0, color=col.gray) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Percent Bias") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.bias.png")),
         width=8, height=8, units="in")

p.ts.cor <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=cor, color=method)) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Correlation") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.cor.png")),
       width=8, height=8, units="in")

## scatterplots at different timesteps
p.depletion.scatter <-
  df %>% 
  subset(Time %in% ts.plot & 
           method %in% methods.plot) %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(Time ~ stream_BC+analytical) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") + 
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.depletion.scatter.png")),
         width=8, height=8, units="in")

p.depletion.scatter.limits <-
  df %>% 
  subset(Time %in% ts.plot & 
           method %in% methods.plot) %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(Time ~ stream_BC+analytical) +
  scale_x_continuous(name="Analytical Depletion Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", limits=c(0,1)) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.depletion.scatter.limits.png")),
         width=8, height=8, units="in")

p.match.scatter <- 
  df.most.affected.prc %>% 
  transform(method = gsub(".max", "", method)) %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=depletion.prc, y=Qf.MODFLOW.max, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point(shape=21, alpha=0.25) +
  facet_grid(stream_BC ~ analytical) +
  scale_x_continuous(name="Analytical Depletion Fraction, most affected segment only", limits=c(0,1), expand=c(0,0)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction, most affected segment only", limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") + 
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.match.scatter.png")),
         width=8, height=8, units="in")

# capture fraction is independent of depletion apportionment equation so only plot one
p.sum.scatter <-
  df.sum %>% 
  subset(method %in% methods.plot[1]) %>% 
  ggplot(aes(x=Qf.analytical, y=Qf.MODFLOW, color=Time-ts.pump.start)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(stream_BC ~ analytical, scales="free_y", 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Analytical Capture Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Capture Fraction") +
  scale_color_continuous(name="Time Since Start\nOf Pumping [days]") +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.sum.scatter.png")),
         width=6, height=6, units="in")

p.sum.scatter.limits <-
  df.sum %>% 
  subset(method %in% methods.plot[1]) %>% 
  ggplot(aes(x=Qf.analytical, y=Qf.MODFLOW, color=Time-ts.pump.start)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(stream_BC ~ analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Analytical Capture Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Capture Fraction", limits=c(0,1)) +
  scale_color_continuous(name="Time Since Start\nOf Pumping [days]") +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.sum.scatter.limits.png")),
         width=6, height=6, units="in")

# only points with depletion > 5% (big)
p.ts.big.KGE <-
  df.fit.big %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=method)) +
  geom_hline(yintercept=0, color=col.gray) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Kling-Gupta Efficiency") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.big.KGE.png")),
         width=8, height=8, units="in")

p.ts.big.MSE <-
  df.fit.big %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=MSE.overall, color=method)) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Mean Squared Error", labels=scales::percent) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.big.MSE.png")),
         width=8, height=8, units="in")

p.ts.big.bias <-
  df.fit.big %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=bias, color=method)) +
  geom_hline(yintercept=0, color=col.gray) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Percent Bias") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.big.bias.png")),
         width=8, height=8, units="in")

p.ts.big.cor <-
  df.fit.big %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=cor, color=method)) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Correlation") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.ts.big.cor.png")),
         width=8, height=8, units="in")

## fit for different wells
p.well.map.KGE <-
  df.fit.well %>% 
  subset(stream_BC=="RIV" & analytical=="glover" & method %in% methods.plot) %>% 
  ggplot() + 
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_point(aes(x=lon, y=lat, size=KGE.overall)) +
  facet_wrap(~method, ncol=2) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5)) +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.well.map.KGE.png")),
         width=8, height=8, units="in")

p.well.map.MSE <-
  df.fit.well %>% 
  subset(stream_BC=="RIV" & analytical=="glover" & method %in% methods.plot) %>% 
  ggplot() + 
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_point(aes(x=lon, y=lat, size=MSE.overall)) +
  facet_wrap(~method, ncol=2) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5)) +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.well.map.MSE.png")),
         width=8, height=8, units="in")

# error scatterplots
p.well.scatter.elev <-
  df.fit.well %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=ztop_m, y=KGE.overall, color=method)) +
  geom_point() +
  facet_grid(stream_BC ~ analytical, scales="free_y", 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Well Land Surface Elevation [m]") +
  scale_y_continuous(name="Well KGE") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.well.scatter.elev.png")),
         width=8, height=8, units="in")

p.well.scatter.WTD <-
  df.fit.well %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=(ztop_m-wte_m), y=KGE.overall, color=method)) +
  geom_point() +
  facet_grid(stream_BC ~ analytical, scales="free_y", 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Well Water Table Depth [m]") +
  scale_y_continuous(name="Well KGE") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.well.scatter.WTD.png")),
         width=8, height=8, units="in")

### some RIV only plots for presentation
p.depletion.scatter.limits <-
  df %>% 
  subset(Time %in% ts.plot & 
           method %in% methods.plot & 
           stream_BC=="RIV") %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point(shape=21) +
  facet_grid(analytical ~ Time) +
  scale_x_continuous(name="Analytical Depletion Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", limits=c(0,1)) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.depletion.scatter.limits_RIV.png")),
         width=8, height=6, units="in")

p.sum.scatter.limits <-
  df.sum %>% 
  subset(method %in% methods.plot[1] & 
           stream_BC=="RIV") %>% 
  ggplot(aes(x=Qf.analytical, y=Qf.MODFLOW, color=Time-ts.pump.start)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(. ~ analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Analytical Capture Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Capture Fraction", limits=c(0,1)) +
  scale_color_continuous(name="Time Since Start\nOf Pumping [days]") +
  theme(legend.position="bottom") +
  ggsave(file.path("results", paste0("Navarro_DepletionApportionment+Analytical_Transient_CompareToMODFLOW_p.sum.scatter.limits_RIV.png")),
         width=6, height=4, units="in")