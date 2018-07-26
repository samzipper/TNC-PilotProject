## Navarro_CompareMODFLOWtoDepletionApportionment_SteadyStateMODFLOWtoTransientAnalytical.R
#' This script is intended to compare estimated steady-state streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_RIV-SummarizeLeakage.py)
#' to transient analytical results calculated over a very long time period to see if they
#' ever match up.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

## what depletion apportionment output do you want?
apportionment_name <- "_LocalArea"      # output from Navarro_DepletionApportionment_LocalArea.R
#apportionment_name <- "_AdjacentOnly"   # output from Navarro_DepletionApportionment_AdjacentOnly.R
#apportionment_name <- "_WholeDomain"   # output from Navarro_DepletionApportionment_WholeDomain.R

#### First: process steady-state MODFLOW data
## MODFLOW - steady state results
stream_BC <- "RIV"
modflow_v <- "mfnwt"

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## threshold for analysis
f.thres <- 0.001

## load stream shapefile
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# figure out which segments are part of Navarro
segs.navarro <- shp.streams@data$SegNum[shp.streams@data$TerminalPa==outlet.TerminalPa]

# leakage
df_MODFLOW <- 
  file.path("modflow","HTC", "Navarro", "SteadyState", stream_BC, modflow_v, paste0(stream_BC, "-SummarizeLeakage.csv")) %>% 
  read.csv(stringsAsFactors=F)

## if it's SFR: need to get SegNum
if (stream_BC=="SFR"){
  df_MODFLOW <- 
    read.table(file.path("modflow", "input", "isfr_SegmentData.txt"), 
               stringsAsFactors=F, header=T) %>% 
    left_join(df_MODFLOW, ., by="SFR_NSEG")
  
}

## calculate MODFLOW depletion
# make separate column for no pumping depletion
df_MODFLOW<- 
  df_MODFLOW %>% 
  subset(WellNum==0, select=c("SegNum", "leakage", "MNW_net")) %>% 
  set_colnames(c("SegNum", "leakage_NoPump", "MNW_NoPump")) %>% 
  left_join(subset(df_MODFLOW, WellNum != 0), ., by=c("SegNum")) %>% 
  arrange(WellNum, SegNum) %>% 
  subset(SegNum %in% segs.navarro)

# calculate change in leakage and MNW
# convention: negative value means flow out of aquifer (gaining stream, pumping well)
# so if (leakage_NoPump - leakage) < 0, stream is depleted due to pumping
df_MODFLOW$Qw_m3.d <- df_MODFLOW$MNW_net - df_MODFLOW$MNW_NoPump
df_MODFLOW$depletion_m3.d <- df_MODFLOW$leakage_NoPump - df_MODFLOW$leakage

## calculate MODFLOW depletion
#df_MODFLOW$depletion.prc.modflow <- df_MODFLOW$depletion_m3.d/df_MODFLOW$Qw_m3.d
df_MODFLOW$depletion.prc.modflow <- df_MODFLOW$depletion_m3.d/Qw

df_MODFLOW <- subset(df_MODFLOW, depletion.prc.modflow > f.thres)

#### load depletion apportionment output
# choose your apportionment method - LocalArea and WholeDomain have best performance for SS comparison
df.apportionment <- 
  file.path("results", paste0("Navarro_DepletionApportionment", apportionment_name, "_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  dplyr::select(WellNum, SegNum, distToWell.min.m, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "depletion.prc.modflow")],
            by=c("SegNum", "WellNum"))
df.apportionment[is.na(df.apportionment)] <- 0

## combine into single melted data frame
df <- 
  melt(df.apportionment, id=c("WellNum", "SegNum", "distToWell.min.m", "depletion.prc.modflow"),
       value.name="depletion.prc", variable.name="method")

# check for NAs
sum(is.na(df))

# get rid of reaches with no depletion
df <- subset(df, depletion.prc.modflow > f.thres | depletion.prc > f.thres)

# calculate fit statistics for steady-state
df.fit.SS <- 
  df %>% 
  group_by(method) %>% 
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

#### run transient analytical for long time period
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

## load RIV and WEL input to figure out thickness for transmissivity
## load well input data
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", "RIV", "mfnwt", "wte.csv")), header=F)

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$wte_m <- m.wte[as.matrix(df.wel[,c("row", "col")])+1]

# add elevation to df.apportion
df.apportionment <- left_join(df.apportionment, df.wel[,c("WellNum", "wte_m", "ztop_m")], by="WellNum")

## stream elevation needed for Reeves approximation of Hunt lambda
df.apportionment <- 
  read.table(file.path("modflow", "input", "isfr_ReachData.txt"), stringsAsFactors=F, header=T) %>% 
  group_by(SegNum, SFR_NSEG) %>% 
  summarize(totalStreamLength_m = sum(length_m),
            streambed_elev_m = median(elev_m_min)-depth) %>% 
  left_join(df.apportionment, ., by=c("SegNum")) 

# calculate aquifer vertical thickness for transmissivity
screen_length <- 50  # [m] - should be same as script MODFLOW_Navarro-SteadyState.py
df.apportionment$thickness_m <- abs(df.apportionment$wte_m-df.apportionment$streambed_elev_m)    # reeves et al- uses vertical distance between top of well screen and streambed 
df.apportionment$thickness_m[df.apportionment$thickness_m < screen_length] <- screen_length  # if vertical distance is < screen length, use screen length

## timesteps to calculate
ts.all <- unique(c(1:10 %o% 10^(2:5)))
start.flag <- T
for (ts in ts.all){
  Qf <- glover(ts, d = df.apportionment$distToWell.min.m, S = sy, Tr = hk*df.apportionment$thickness_m)
  
  if (start.flag){
    df.transient <- data.frame(time = ts,
                               WellNum = df.apportionment$WellNum,
                               SegNum = df.apportionment$SegNum,
                               f.InvDistSq = df.apportionment$f.InvDistSq,
                               f.WebSq = df.apportionment$f.WebSq,
                               f.TPoly = df.apportionment$f.TPoly,
                               Qf = Qf)
    start.flag <- F
  } else {
    df.transient <- rbind(df.transient, 
                          data.frame(time = ts,
                                     WellNum = df.apportionment$WellNum,
                                     SegNum = df.apportionment$SegNum,
                                     f.InvDistSq = df.apportionment$f.InvDistSq,
                                     f.WebSq = df.apportionment$f.WebSq,
                                     f.TPoly = df.apportionment$f.TPoly,
                                     Qf = Qf))
  }
  
  print(paste0(ts, " complete"))
}

# combine transient with MODFLOW
df.transient.modflow <- 
  full_join(df.transient, df_MODFLOW[,c("WellNum", "SegNum", "depletion.prc.modflow")], by=c("WellNum", "SegNum")) %>% 
  transform(Qf.InvDistSq = Qf*f.InvDistSq,
            Qf.WebSq = Qf*f.WebSq,
            Qf.TPoly = Qf*f.TPoly) %>% 
  dplyr::select(time, WellNum, SegNum, Qf.InvDistSq, Qf.WebSq, Qf.TPoly, depletion.prc.modflow) %>% 
  melt(id=c("time", "WellNum", "SegNum", "depletion.prc.modflow"),
       value.name = "depletion.prc", variable.name="method") %>% 
  replace_na(list(depletion.prc = 0, depletion.prc.modflow = 0)) %>% 
  subset(depletion.prc.modflow > f.thres | depletion.prc > f.thres)

# calculate fit statistics for transient
df.fit.transient <- 
  df.transient.modflow %>% 
  group_by(method, time) %>% 
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

## plot
ggplot(df.fit.transient, aes(x=time, y=KGE.overall, color=method)) +
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=min(df.fit.transient$time), xmax=3650, 
           fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=df.fit.SS$KGE.overall[df.fit.SS$method=="f.InvDistSq"], 
             linetype="dashed", color=pal.method["f.InvDistSq"]) +
  geom_hline(yintercept=df.fit.SS$KGE.overall[df.fit.SS$method=="f.WebSq"], 
             linetype="dashed", color=pal.method["f.WebSq"]) +
  geom_hline(yintercept=df.fit.SS$KGE.overall[df.fit.SS$method=="f.TPoly"], 
             linetype="dashed", color=pal.method["f.TPoly"]) +
  geom_line() +
  scale_color_manual(name="Depletion\nApportionment", 
                     values=pal.method.Qf, labels=labels.method.Qf) +
  scale_x_log10(name="Analytical Timestep [days]", expand=c(0,0)) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(name="KGE") +
  #labs(title="Comparing Steady-State MODFLOW to Transient Analytical (Glover)", 
  #     subtitle="Dashed Line = KGE for Depletion Apportionment Only, no analytical") +
  theme(legend.position=c(0.99,0.01),
        legend.justification=c(1,0),
        plot.margin=unit(c(1,5,1,1), "mm")) +
  ggsave(file.path("results", "Navarro_CompareMODFLOWtoDepletionApportionment_SteadyStateMODFLOWtoTransientAnalytical.png"),
         width=190, height=120, units="mm")
